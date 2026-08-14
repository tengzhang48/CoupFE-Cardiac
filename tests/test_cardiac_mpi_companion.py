from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest

import post
from distributed_local_pressure import RankLocalInvalidDeformationError
from distributed_solver import (
    DistributedPetscSnesSolver,
    DistributedSnesSolveError,
    FGMRES_ASM_ILU1_PROFILE,
)


ROOT = Path(__file__).resolve().parents[1]
BENCHMARK = ROOT / "examples" / "cardiac_benchmark"


def _require_parallel_runtime():
    try:
        from petsc4py import PETSc
        from mpi4py import MPI
    except ImportError as error:
        pytest.fail(f"MPI gate requested without its dependencies: {error}")
    return PETSc, MPI


def _mpi_launcher():
    sibling = Path(sys.executable).with_name("mpiexec")
    if sibling.is_file():
        return str(sibling)
    return shutil.which("mpiexec") or shutil.which("mpirun")


def _subprocess_environment():
    """Expose application modules to child interpreters launched by MPI."""
    environment = os.environ.copy()
    inherited = environment.get("PYTHONPATH")
    entries = [str(BENCHMARK)]
    if inherited:
        entries.append(inherited)
    environment["PYTHONPATH"] = os.pathsep.join(entries)
    environment["OMP_NUM_THREADS"] = "1"
    return environment


class _ScalarBatch:
    def __init__(self, target=1.0):
        self.target = float(target)
        self.commits = 0
        self.clears = 0
        self.residual_evaluations = 0
        self.tangent_evaluations = 0

    def element_residual_batch(self, coordinates, displacement, increment):
        self.residual_evaluations += 1
        return np.asarray(displacement, dtype=float) - self.target

    def element_tangent_batch(self, coordinates, displacement, increment):
        self.tangent_evaluations += 1
        return np.ones((len(displacement), 1, 1), dtype=float)

    def clear_cache(self):
        self.clears += 1

    def commit(self):
        self.commits += 1


class _PositiveLogBatch(_ScalarBatch):
    def element_residual_batch(self, coordinates, displacement, increment):
        value = float(displacement[0, 0])
        if value <= 0.0:
            raise RankLocalInvalidDeformationError(9, value)
        return np.array([[np.log(value) + 2.0]])

    def element_tangent_batch(self, coordinates, displacement, increment):
        value = float(displacement[0, 0])
        if value <= 0.0:
            raise RankLocalInvalidDeformationError(9, value)
        return np.array([[[1.0 / value]]])


class _IndependentCheckFailureBatch(_ScalarBatch):
    def element_residual_batch(self, coordinates, displacement, increment):
        residual = super().element_residual_batch(
            coordinates, displacement, increment
        )
        if self.clears >= 2:
            residual = residual + 1.0e-4
        return residual


class _TwoDofBatch:
    def __init__(self):
        self.matrix = np.diag([1.0, 10.0])
        self.target = np.ones(2)
        self.commits = 0

    def clear_cache(self):
        pass

    def element_residual_batch(self, coordinates, displacement, increment):
        values = np.asarray(displacement[0], dtype=float)
        return (self.matrix @ values - self.target)[None, :]

    def element_tangent_batch(self, coordinates, displacement, increment):
        return self.matrix[None, :, :]

    def commit(self):
        self.commits += 1


@pytest.mark.mpi
def test_distributed_snes_backtracks_collectively_and_commits_once():
    PETSc, _MPI = _require_parallel_runtime()
    if PETSc.COMM_WORLD.getSize() != 1:
        pytest.skip("this focused callback test owns the one-rank COMM_WORLD case")
    batch = _PositiveLogBatch()
    application_solver = DistributedPetscSnesSolver(
        1,
        np.array([[0]], dtype=int),
        np.zeros((1, 1)),
        batch,
        np.zeros(1),
        dof_per_node=1,
        u0=np.ones(1),
    )
    try:
        displacement, diagnostics = application_solver.solve_step(t=0.1, dt=0.1)
        np.testing.assert_allclose(
            displacement, [np.exp(-2.0)], rtol=1.0e-9, atol=1.0e-12
        )
        assert diagnostics.function_domain_rejections >= 1
        assert diagnostics.last_function_domain_error == (
            "distributed cardiac trial rejected: global_element=9"
        )
        assert diagnostics.snes_converged_reason > 0
        assert diagnostics.ksp_converged_reason > 0
        assert diagnostics.final_residual_norm < (
            diagnostics.residual_acceptance_threshold
        )
        assert batch.commits == 1
        assert application_solver.accepted_steps == 1
    finally:
        application_solver.close()


@pytest.mark.mpi
def test_independent_residual_failure_does_not_commit_any_history():
    PETSc, _MPI = _require_parallel_runtime()
    if PETSc.COMM_WORLD.getSize() != 1:
        pytest.skip("this focused commit test owns the one-rank COMM_WORLD case")
    batch = _IndependentCheckFailureBatch()
    application_solver = DistributedPetscSnesSolver(
        1,
        np.array([[0]], dtype=int),
        np.zeros((1, 1)),
        batch,
        np.zeros(1),
        dof_per_node=1,
    )
    try:
        with pytest.raises(
            DistributedSnesSolveError,
            match="independent residual rule before state commit",
        ):
            application_solver.solve_step(t=0.1, dt=0.1)
        assert batch.commits == 0
        assert application_solver.accepted_steps == 0
        np.testing.assert_array_equal(application_solver._u_prev, [0.0])
        np.testing.assert_array_equal(application_solver._v_prev, [0.0])
    finally:
        application_solver.close()


@pytest.mark.mpi
def test_iterative_limit_reports_failed_ksp_without_committing_state():
    PETSc, _MPI = _require_parallel_runtime()
    if PETSc.COMM_WORLD.getSize() != 1:
        pytest.skip("this focused KSP failure test owns one-rank COMM_WORLD")
    batch = _TwoDofBatch()
    application_solver = DistributedPetscSnesSolver(
        2,
        np.array([[0, 1]], dtype=int),
        np.zeros((1, 2)),
        batch,
        np.zeros(2),
        dof_per_node=1,
        linear_solver_profile=FGMRES_ASM_ILU1_PROFILE,
    )
    ksp = application_solver._snes.getKSP()
    ksp.getPC().setType("none")
    ksp.setTolerances(rtol=1.0e-12, atol=1.0e-14, max_it=1)
    ksp.setErrorIfNotConverged(False)
    try:
        with pytest.raises(DistributedSnesSolveError) as caught:
            application_solver.solve_step(t=0.125, dt=0.025)
        diagnostics = caught.value.diagnostics
        assert diagnostics is not None
        assert diagnostics.time == 0.125
        assert diagnostics.dt == 0.025
        assert diagnostics.snes_converged_reason < 0
        assert diagnostics.ksp_converged_reason < 0
        assert diagnostics.ksp_iterations_per_solve[-1] == 1
        assert diagnostics.ksp_reasons_per_solve[-1] < 0
        assert len(diagnostics.ksp_residual_histories) == len(
            diagnostics.ksp_reasons_per_solve
        )
        assert diagnostics.ksp_residual_histories[-1]
        assert batch.commits == 0
        assert application_solver.accepted_steps == 0
        np.testing.assert_array_equal(application_solver._u_prev, np.zeros(2))
        np.testing.assert_array_equal(application_solver._v_prev, np.zeros(2))
    finally:
        application_solver.close()


@pytest.mark.mpi
def test_distributed_snes_context_persists_and_commits_once_per_step():
    PETSc, _MPI = _require_parallel_runtime()
    if PETSc.COMM_WORLD.getSize() != 1:
        pytest.skip("this focused reuse test owns the one-rank COMM_WORLD case")
    batch = _ScalarBatch(target=1.0)
    application_solver = DistributedPetscSnesSolver(
        1,
        np.array([[0]], dtype=int),
        np.zeros((1, 1)),
        batch,
        np.zeros(1),
        dof_per_node=1,
    )
    try:
        snes_identity = id(application_solver._snes)
        first, _ = application_solver.solve_step(t=0.1, dt=0.1)
        batch.target = 2.0
        second, _ = application_solver.solve_step(t=0.2, dt=0.1)
        np.testing.assert_allclose(first, [1.0])
        np.testing.assert_allclose(second, [2.0])
        assert id(application_solver._snes) == snes_identity
        assert batch.commits == 2
        assert application_solver.accepted_steps == 2
    finally:
        application_solver.close()


@pytest.mark.mpi
def test_exact_equilibrium_needs_no_tangent_or_linear_solve():
    PETSc, _MPI = _require_parallel_runtime()
    if PETSc.COMM_WORLD.getSize() != 1:
        pytest.skip("this focused zero-iteration test owns one-rank COMM_WORLD")
    batch = _ScalarBatch(target=0.0)
    application_solver = DistributedPetscSnesSolver(
        1,
        np.array([[0]], dtype=int),
        np.zeros((1, 1)),
        batch,
        np.zeros(1),
        dof_per_node=1,
    )
    try:
        displacement, diagnostics = application_solver.solve_step(t=0.1, dt=0.1)
        np.testing.assert_array_equal(displacement, [0.0])
        assert diagnostics.snes_converged_reason > 0
        assert diagnostics.ksp_converged_reason == 0
        assert diagnostics.nonlinear_iterations == 0
        assert diagnostics.linear_iterations == 0
        assert batch.tangent_evaluations == 0
        assert batch.commits == 1
    finally:
        application_solver.close()


@pytest.mark.mpi
@pytest.mark.parametrize("ranks", [2, 4])
def test_collective_invalid_trial_uses_lowest_global_element(ranks):
    _require_parallel_runtime()
    launcher = _mpi_launcher()
    if launcher is None:
        pytest.fail("MPI gate requested but no mpiexec/mpirun is available")
    code = r'''
import numpy as np
from mpi4py import MPI
from distributed_local_pressure import RankLocalInvalidDeformationError
from distributed_solver import DistributedPetscSnesSolver

comm = MPI.COMM_WORLD
rank, size = comm.rank, comm.size
ndof = 8
global_ids = np.arange(rank, ndof, size, dtype=int)

class Batch:
    def __init__(self):
        self.commits = 0
    def clear_cache(self):
        pass
    def commit(self):
        self.commits += 1
    def element_residual_batch(self, coordinates, displacement, increment):
        values = displacement[:, 0]
        bad = np.flatnonzero(values <= 0.0)
        if len(bad):
            local = int(bad[np.argmin(global_ids[bad])])
            raise RankLocalInvalidDeformationError(global_ids[local], values[local])
        return (np.log(values) + 2.0)[:, None]
    def element_tangent_batch(self, coordinates, displacement, increment):
        values = displacement[:, 0]
        bad = np.flatnonzero(values <= 0.0)
        if len(bad):
            local = int(bad[np.argmin(global_ids[bad])])
            raise RankLocalInvalidDeformationError(global_ids[local], values[local])
        return (1.0 / values)[:, None, None]

batch = Batch()
solver = DistributedPetscSnesSolver(
    ndof, global_ids[:, None], np.zeros((len(global_ids), 1)), batch,
    np.zeros(ndof), dof_per_node=1, u0=np.ones(ndof),
)
try:
    displacement, diagnostics = solver.solve_step(t=0.1, dt=0.1)
    assert np.max(np.abs(displacement - np.exp(-2.0))) < 1.0e-9
    assert diagnostics.function_domain_rejections >= 1
    assert diagnostics.last_function_domain_error.endswith("global_element=0")
    assert batch.commits == 1
    if rank == 0:
        print("COLLECTIVE_INVALID_OK", flush=True)
finally:
    solver.close()
'''
    environment = _subprocess_environment()
    completed = subprocess.run(
        [launcher, "-n", str(ranks), sys.executable, "-c", code],
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=180,
    )
    assert completed.returncode == 0, completed.stdout
    assert "COLLECTIVE_INVALID_OK" in completed.stdout


@pytest.mark.mpi
def test_tiny_case_b_mpi_ranks_match_retained_serial_companion(tmp_path):
    _require_parallel_runtime()
    launcher = _mpi_launcher()
    if launcher is None:
        pytest.fail("MPI gate requested but no mpiexec/mpirun is available")
    environment = _subprocess_environment()

    common = [
        "--case",
        "B",
        "--mesh-topology",
        "polar-ring",
        "--nt",
        "1",
        "--nmu",
        "2",
        "--ntheta",
        "4",
        "--apex-offset",
        "0.2",
        "--dt",
        "0.002",
        "--tend",
        "0.004",
    ]
    serial_output = tmp_path / "serial.npz"
    serial_command = [
        sys.executable,
        str(BENCHMARK / "run.py"),
        *common,
        "--formulation",
        "local-pressure",
        "--mass",
        "lumped",
        "--nonlinear-solver",
        "petsc-snes",
        "--build-dir",
        str(tmp_path / "serial-build"),
        "--out",
        str(serial_output),
    ]
    serial = subprocess.run(
        serial_command,
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=900,
    )
    assert serial.returncode == 0, serial.stdout
    assert serial_output.is_file()

    with np.load(serial_output, allow_pickle=False) as archive:
        assert str(archive["mass_representation"]) == "lumped_row_sum"
        reference = {
            name: archive[name].copy()
            for name in (
                "U_peak",
                "u0",
                "u1",
                "det_f_gauss_peak",
                "element_pressure_peak_pa",
            )
        }

    for ranks in (1, 2, 4):
        output = tmp_path / f"mpi-{ranks}.npz"
        build_root = tmp_path / f"mpi-{ranks}-build"
        command = [
            launcher,
            "-n",
            str(ranks),
            sys.executable,
            str(BENCHMARK / "run_mpi.py"),
            *common,
            "--build-dir",
            str(build_root),
            "--out",
            str(output),
        ]
        completed = subprocess.run(
            command,
            cwd=ROOT,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=900,
        )
        assert completed.returncode == 0, completed.stdout
        assert output.is_file()
        assert sorted(path.name for path in build_root.iterdir()) == [
            f"rank-{rank:05d}" for rank in range(ranks)
        ]
        with np.load(output, allow_pickle=False) as archive:
            validated = post._load_validated_result(
                archive, output, requested_case="step_0B"
            )
            assert validated["solver_name"] == "petsc-snes-mpi"
            assert validated["mass_representation"] == "lumped_row_sum"
            for name, expected in reference.items():
                np.testing.assert_allclose(
                    archive[name], expected, rtol=2.0e-11, atol=2.0e-13
                )
            assert int(archive["mpi_world_size"]) == ranks
            assert int(archive["mpi_ranks"]) == ranks
            assert int(archive["completed_steps"]) == 2
            assert str(archive["driver"]) == (
                "examples/cardiac_benchmark/run_mpi.py"
            )
            assert str(archive["nonlinear_solver"]) == "petsc-snes-mpi"
            assert str(archive["mass_representation"]) == "lumped_row_sum"
            assert str(archive["fiber_sampling"]) == "cg1_gram_schmidt"
            assert str(archive["fiber_sampling_option"]) == "cg1"
            assert str(archive["tbar_definition"]) == "analytic_parametric"
            assert str(archive["tbar_source_filename"]) == ""
            assert str(archive["tbar_source_sha256"]) == ""
            assert str(archive["tbar_metadata_filename"]) == ""
            assert str(archive["tbar_metadata_sha256"]) == ""
            assert str(archive["tbar_metadata_schema"]) == ""
            assert not bool(archive["isotropic"])
            assert float(archive["material_eta_pa_s"]) == 100.0
            assert bool(archive["viscous_term_active"])
            assert str(archive["parameter_variant"]) == "benchmark_eta"
            assert str(archive["mesh_topology"]) == "polar_ring"
            assert int(np.sum(archive["mpi_local_element_counts"])) == 8
            configuration = json.loads(str(archive["solver_configuration_json"]))
            assert configuration["communicator"] == "PETSc.COMM_WORLD"
            assert configuration["ranks"] == ranks
            assert configuration["function_domain_rejection_api"] == (
                "nonfinite residual for PETSc BT"
            )
            assert configuration["snes_type"] == "newtonls"
            assert configuration["line_search_type"] == "bt"
            assert configuration["ksp_type"] == "preonly"
            assert configuration["pc_type"] == "lu"
            assert configuration["configured_factor_solver_type"] == (
                "superlu_dist"
            )
            assert configuration["element_evaluation_mode"] == "joint"
            assert isinstance(
                configuration["compiled_material_residual_only_available"], bool
            )
            assert str(archive["element_evaluation_mode"]) == "joint"
            assert bool(archive["compiled_material_residual_only_available"]) == (
                configuration["compiled_material_residual_only_available"]
            )


@pytest.mark.mpi
def test_active_case_a_mpi_ranks_match_serial_companion(tmp_path):
    """Exercise nonzero active tension, not only the zero-load Case A prefix."""
    _require_parallel_runtime()
    launcher = _mpi_launcher()
    if launcher is None:
        pytest.fail("MPI gate requested but no mpiexec/mpirun is available")
    environment = _subprocess_environment()
    common = [
        "--case",
        "A",
        "--mesh-topology",
        "polar-ring",
        "--nt",
        "1",
        "--nmu",
        "2",
        "--ntheta",
        "4",
        "--apex-offset",
        "0.2",
        "--dt",
        "0.01",
        "--tend",
        "0.2",
    ]
    serial_output = tmp_path / "serial-case-a.npz"
    serial = subprocess.run(
        [
            sys.executable,
            str(BENCHMARK / "run.py"),
            *common,
            "--formulation",
            "local-pressure",
            "--mass",
            "lumped",
            "--nonlinear-solver",
            "petsc-snes",
            "--element-evaluation",
            "split",
            "--build-dir",
            str(tmp_path / "serial-case-a-build"),
            "--out",
            str(serial_output),
        ],
        cwd=ROOT,
        env=environment,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=900,
    )
    assert serial.returncode == 0, serial.stdout
    with np.load(serial_output, allow_pickle=False) as archive:
        assert str(archive["mass_representation"]) == "lumped_row_sum"
        reference = {
            name: archive[name].copy()
            for name in ("U_peak", "u0", "u1", "det_f_gauss_peak")
        }
        assert float(np.max(archive["tau"])) > 0.0
        np.testing.assert_array_equal(archive["pres"], 0.0)

    for ranks in (1, 2, 4):
        output = tmp_path / f"mpi-case-a-{ranks}.npz"
        completed = subprocess.run(
            [
                launcher,
                "-n",
                str(ranks),
                sys.executable,
                str(BENCHMARK / "run_mpi.py"),
                *common,
                "--element-evaluation",
                "split",
                "--build-dir",
                str(tmp_path / f"mpi-case-a-{ranks}-build"),
                "--out",
                str(output),
            ],
            cwd=ROOT,
            env=environment,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=900,
        )
        assert completed.returncode == 0, completed.stdout
        with np.load(output, allow_pickle=False) as archive:
            validated = post._load_validated_result(
                archive, output, requested_case="step_0A"
            )
            assert validated["solver_name"] == "petsc-snes-mpi"
            assert validated["mass_representation"] == "lumped_row_sum"
            for name, expected in reference.items():
                np.testing.assert_allclose(
                    archive[name], expected, rtol=3.0e-11, atol=3.0e-13
                )
            assert str(archive["case"]) == "A"
            assert str(archive["mass_representation"]) == "lumped_row_sum"
            assert str(archive["fiber_sampling_option"]) == "cg1"
            assert str(archive["tbar_definition"]) == "analytic_parametric"
            assert float(archive["material_eta_pa_s"]) == 100.0
            assert str(archive["parameter_variant"]) == "benchmark_eta"
            assert str(archive["mesh_topology"]) == "polar_ring"
            assert int(archive["mpi_world_size"]) == ranks
            assert str(archive["element_evaluation_mode"]) == "split"
            assert float(np.max(archive["tau"])) > 0.0
            np.testing.assert_array_equal(archive["pres"], 0.0)


@pytest.mark.mpi
def test_mpi_driver_input_failure_writes_no_archive(tmp_path):
    _require_parallel_runtime()
    output = tmp_path / "must-not-exist.npz"
    completed = subprocess.run(
        [
            sys.executable,
            str(BENCHMARK / "run_mpi.py"),
            "--apex-offset",
            "0",
            "--out",
            str(output),
        ],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=30,
    )
    assert completed.returncode != 0
    assert not output.exists()
