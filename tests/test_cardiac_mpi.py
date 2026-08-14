from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import solver
from coupfe.operators.base import Residual, Tangent
from local_pressure import InvalidDeformationError


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = sorted((ROOT / "examples" / "mpi_smoke").glob("distributed_cardiac_*.py"))


class _ScalarOperator:
    def __init__(self):
        self.commits = 0

    def residual(self, U, state, t, dt):
        return Residual(np.array([0]), np.array([float(U[0]) - 1.0]))

    def tangent(self, U, state, t, dt):
        return Tangent(np.array([0]), np.array([0]), np.array([1.0]))

    def commit(self, U, state, t, dt):
        self.commits += 1
        return state


class _PositiveLogOperator(_ScalarOperator):
    """Scalar domain analogue whose first full Newton trial is invalid."""

    def residual(self, U, state, t, dt):
        value = float(U[0])
        if value <= 0.0:
            raise InvalidDeformationError("logarithm trial is outside x > 0")
        return Residual(np.array([0]), np.array([np.log(value) + 2.0]))

    def tangent(self, U, state, t, dt):
        value = float(U[0])
        if value <= 0.0:
            raise InvalidDeformationError("logarithm tangent is outside x > 0")
        return Tangent(np.array([0]), np.array([0]), np.array([1.0 / value]))


def _mpi_launcher():
    # Conda PETSc builds are commonly linked against the environment's MPICH;
    # a system OpenMPI launcher can be present first on PATH but is not ABI
    # compatible. Prefer the launcher beside the selected Python interpreter.
    sibling = Path(sys.executable).with_name("mpiexec")
    if sibling.is_file():
        return str(sibling)
    return shutil.which("mpiexec") or shutil.which("mpirun")


@pytest.mark.mpi
def test_application_petsc_snes_callbacks_converge_and_reuse_context():
    try:
        __import__("petsc4py")
    except ImportError as error:
        pytest.fail(f"PETSc gate requested without petsc4py: {error}")

    operator = _ScalarOperator()
    application_solver = solver.PetscSnesSolver()
    first, _, first_diagnostics = application_solver.solve(
        [operator], np.zeros(1), None, 1, {}, t=0.1, dt=0.1
    )
    snes_identity = id(application_solver._snes)
    second, _, second_diagnostics = application_solver.solve(
        [operator], np.zeros(1), None, 1, {}, t=0.2, dt=0.1
    )

    np.testing.assert_allclose(first, [1.0])
    np.testing.assert_allclose(second, [1.0])
    assert id(application_solver._snes) == snes_identity
    assert operator.commits == 2
    for diagnostics in (first_diagnostics, second_diagnostics):
        assert diagnostics.snes_converged_reason > 0
        assert diagnostics.ksp_converged_reason > 0
        assert diagnostics.final_residual_norm <= (
            diagnostics.residual_acceptance_threshold
        )
    configuration = application_solver.configuration()
    assert configuration["line_search_type"] == "bt"
    assert configuration["line_search_configuration_api"] in {
        "SNES.getLineSearch",
        "namespaced PETSc option",
    }
    assert configuration["factor_solver_type"] != "unknown"
    application_solver.close()


@pytest.mark.mpi
def test_application_petsc_snes_backtracks_nonfinite_residual_domain_rejection():
    try:
        __import__("petsc4py")
    except ImportError as error:
        pytest.fail(f"PETSc gate requested without petsc4py: {error}")

    operator = _PositiveLogOperator()
    application_solver = solver.PetscSnesSolver()
    displacement, _, diagnostics = application_solver.solve(
        [operator], np.ones(1), None, 1, {}, t=0.1, dt=0.1
    )

    np.testing.assert_allclose(displacement, [np.exp(-2.0)], rtol=1.0e-9)
    assert diagnostics.function_domain_rejections >= 1
    assert diagnostics.last_function_domain_error == (
        "logarithm trial is outside x > 0"
    )
    assert diagnostics.snes_converged_reason > 0
    assert diagnostics.final_residual_norm <= diagnostics.residual_acceptance_threshold
    assert np.all(np.isfinite(diagnostics.residual_history))
    assert operator.commits == 1
    assert application_solver.configuration()["function_domain_rejection_api"] == (
        "nonfinite residual for PETSc BT"
    )
    assert json.loads(json.dumps(diagnostics.as_dict(), allow_nan=False))[
        "last_function_domain_error"
    ] == "logarithm trial is outside x > 0"
    application_solver.close()


@pytest.mark.mpi
@pytest.mark.parametrize("script", SCRIPTS, ids=lambda path: path.stem)
def test_distributed_cardiac_scripts_match_their_declared_reference(script, tmp_path):
    try:
        __import__("petsc4py")
        __import__("mpi4py")
    except ImportError as error:
        pytest.fail(f"MPI gate requested without its dependencies: {error}")
    launcher = _mpi_launcher()
    if launcher is None:
        pytest.fail("MPI gate requested but no mpiexec/mpirun is available")

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "1"
    env["CARDIAC_SCALING_DIR"] = str(tmp_path)
    for ranks in (1, 2, 4):
        command = [launcher, "-n", str(ranks), sys.executable, str(script)]
        if script.name == "distributed_cardiac_scaling.py":
            command.extend(["1", "8", "12"])
        completed = subprocess.run(
            command,
            cwd=tmp_path,
            env=env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=900,
        )
        assert completed.returncode == 0, completed.stdout
        if script.name != "distributed_cardiac_scaling.py":
            assert "FAIL" not in completed.stdout
            assert "OK" in completed.stdout
