"""Focused gates for the closed, like-for-like Cardiac MPI configuration."""
from __future__ import annotations

import copy

import numpy as np
import pytest
import scipy.sparse as sp

import geometry
import run as serial_driver
import run_mpi
import distributed_solver
from consistent_mass import consistent_mass_coo
from distributed_mass import (
    OwnedRowMassMatrix,
    assemble_owned_consistent_mass,
    balanced_row_range,
)
from distributed_material import RankLocalStandardMaterialBatch
from material import CardiacHex8


def _small_closed_mesh():
    return geometry.build_closed_mesh(
        n_t=1,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
        flip_helix=True,
    )


def _serial_mass_matrix(mesh, density=1000.0):
    ndof = 3 * mesh.n_node
    rows, columns, values = consistent_mass_coo(
        mesh.nodes, mesh.elems, density, dof_per_node=3
    )
    return sp.coo_matrix(
        (values, (rows, columns)), shape=(ndof, ndof)
    ).tocsr()


def test_owned_row_consistent_mass_equals_serial_matrix_and_action_1_2_4():
    mesh = _small_closed_mesh()
    ndof = 3 * mesh.n_node
    serial = _serial_mass_matrix(mesh)
    vector = np.linspace(-0.23, 0.31, ndof)
    element_dofs = (
        mesh.elems[:, :, None] * 3 + np.arange(3)[None, None, :]
    ).reshape(mesh.n_elem, -1)

    for size in (1, 2, 4):
        rank_rows = []
        records = []
        for rank in range(size):
            start, stop = balanced_row_range(ndof, rank, size)
            owned = assemble_owned_consistent_mass(
                mesh.nodes,
                mesh.elems,
                1000.0,
                ndof,
                start,
                stop,
                dof_per_node=3,
            )
            expected_touching = np.flatnonzero(
                np.any(
                    (element_dofs >= start) & (element_dofs < stop), axis=1
                )
            )
            np.testing.assert_array_equal(
                owned.touching_element_ids, expected_touching
            )
            np.testing.assert_array_equal(
                owned.action(vector), (serial @ vector)[start:stop]
            )
            rank_rows.append(owned.matrix)
            records.append(owned.metadata())

        reconstructed = sp.vstack(rank_rows, format="csr")
        difference = reconstructed - serial
        assert difference.nnz == 0
        np.testing.assert_array_equal(reconstructed @ vector, serial @ vector)
        summary = run_mpi._validate_mass_partition_records(records, ndof)
        assert summary["partition"] == "owned-row-csr-all-touching-elements"
        assert sum(stop - start for start, stop in summary["row_ranges"]) == ndof


def test_mass_partition_metadata_fails_closed_on_a_missing_row():
    records = [
        {
            "representation": "consistent_q1_hex8",
            "partition": "owned-row-csr-all-touching-elements",
            "row_start": 0,
            "row_end": 3,
            "owned_rows": 3,
            "global_columns": 7,
            "local_nnz": 8,
            "touching_elements": 2,
        },
        {
            "representation": "consistent_q1_hex8",
            "partition": "owned-row-csr-all-touching-elements",
            "row_start": 4,
            "row_end": 7,
            "owned_rows": 3,
            "global_columns": 7,
            "local_nnz": 8,
            "touching_elements": 2,
        },
    ]
    with pytest.raises(ValueError, match="contiguous global partition"):
        run_mpi._validate_mass_partition_records(records, 7)


class _StateOnlyElement:
    def __init__(self, count, state_size):
        self.svars = np.zeros((count, state_size), dtype=float)
        self.svars_trial = self.svars.copy()


@pytest.mark.parametrize("sampling", ["cg1", "gp-direct"])
@pytest.mark.parametrize("size", [1, 2, 4])
def test_rank_fiber_state_initialization_equals_serial(sampling, size):
    mesh = _small_closed_mesh()
    problem = CardiacHex8()
    schema = problem._mat.state_schema
    per_gp = sum(record["size"] for record in schema.values())
    # Nonlinear but bounded tbar data makes this test exercise the supplied
    # field path rather than silently falling back to mesh.param[:, 0].
    tbar = np.asarray(mesh.param[:, 0], dtype=float) ** 2

    serial = _StateOnlyElement(mesh.n_elem, 8 * per_gp)
    serial_driver._init_fiber_svars(
        serial,
        mesh,
        schema,
        per_gp,
        sampling=sampling,
        tbar=tbar,
        asign=-1.0,
    )

    reconstructed = np.empty_like(serial.svars)
    for identifiers in np.array_split(np.arange(mesh.n_elem), size):
        local = _StateOnlyElement(len(identifiers), 8 * per_gp)
        run_mpi._init_fiber_svars_subset(
            local,
            mesh,
            schema,
            per_gp,
            identifiers,
            sampling=sampling,
            tbar=tbar,
            asign=-1.0,
        )
        reconstructed[identifiers] = local.svars

    np.testing.assert_array_equal(reconstructed, serial.svars)


class _DeterministicMaterialElement:
    has_element_r_batch = True

    def __init__(self):
        self.props = np.array([1.75, -0.25])
        self.commits = 0

    def element_r_batch(self, coordinates, displacement, increment):
        return displacement + 0.5 * increment + self.props[0]

    def element_rk_batch(self, coordinates, displacement, increment):
        residual = self.element_r_batch(coordinates, displacement, increment)
        tangent = np.broadcast_to(np.eye(24), (len(displacement), 24, 24)).copy()
        return residual, tangent

    def commit(self):
        self.commits += 1


@pytest.mark.parametrize("evaluation_mode", ["joint", "split"])
@pytest.mark.parametrize("size", [1, 2, 4])
def test_standard_material_partition_blocks_equal_monolithic_operator(
    evaluation_mode, size
):
    mesh = _small_closed_mesh()
    coordinates = mesh.nodes[mesh.elems]
    identifiers = np.arange(mesh.n_elem, dtype=np.int64)
    rng = np.random.default_rng(20260803)
    displacement = rng.normal(0.0, 1.0e-5, (mesh.n_elem, 24))
    increment = rng.normal(0.0, 2.0e-6, (mesh.n_elem, 24))

    monolithic = RankLocalStandardMaterialBatch(
        _DeterministicMaterialElement(),
        coordinates,
        global_element_ids=identifiers,
        evaluation_mode=evaluation_mode,
    )
    expected_residual = monolithic.element_residual_batch(
        coordinates, displacement, increment
    )
    expected_tangent = monolithic.element_tangent_batch(
        coordinates, displacement, increment
    )
    expected_det_f = monolithic.deformation_jacobians(displacement)

    residual = np.empty_like(expected_residual)
    tangent = np.empty_like(expected_tangent)
    determinant = np.empty_like(expected_det_f)
    for local_ids in np.array_split(identifiers, size):
        local = RankLocalStandardMaterialBatch(
            _DeterministicMaterialElement(),
            coordinates[local_ids],
            global_element_ids=local_ids,
            evaluation_mode=evaluation_mode,
        )
        residual[local_ids] = local.element_residual_batch(
            coordinates[local_ids],
            displacement[local_ids],
            increment[local_ids],
        )
        tangent[local_ids] = local.element_tangent_batch(
            coordinates[local_ids],
            displacement[local_ids],
            increment[local_ids],
        )
        determinant[local_ids] = local.deformation_jacobians(
            displacement[local_ids]
        )

    np.testing.assert_array_equal(residual, expected_residual)
    np.testing.assert_array_equal(tangent, expected_tangent)
    np.testing.assert_array_equal(determinant, expected_det_f)


@pytest.mark.parametrize("case", ["A", "B"])
@pytest.mark.parametrize(
    ("formulation", "implementation"),
    [
        (
            "std-kappa",
            run_mpi.DistributedPetscSnesSolver.CLOSED_STD_KAPPA_IMPLEMENTATION,
        ),
        (
            "local-pressure",
            run_mpi.DistributedPetscSnesSolver.CLOSED_LOCAL_PRESSURE_IMPLEMENTATION,
        ),
        (
            "local-pressure-paper",
            run_mpi.DistributedPetscSnesSolver.
            CLOSED_LOCAL_PRESSURE_PAPER_IMPLEMENTATION,
        ),
    ],
)
def test_closed_step0_parser_retains_exact_fixed_configuration(
    case, formulation, implementation
):
    parser = run_mpi._parser()
    arguments = parser.parse_args(
        [
            "--case",
            case,
            "--mesh-topology",
            "closed-multiblock",
            "--formulation",
            formulation,
            "--mass",
            "consistent",
            "--material-eta",
            "100",
            "--fiber-sampling",
            "gp-direct",
            "--nt",
            "2",
            "--ncore",
            "20",
            "--nradial",
            "17",
            "--core-half-width",
            "0.36",
            "--dt",
            "0.001",
            "--tend",
            "0.32",
            "--load-horizon",
            "1.0",
            "--tbar-laplace",
            "field.npy",
        ]
    )
    run_mpi._validate_arguments(parser, arguments)
    assert arguments.case == case
    assert arguments.mesh_topology == "closed-multiblock"
    assert arguments.formulation == formulation
    assert arguments.mass == "consistent"
    assert arguments.material_eta == 100.0
    assert arguments.fiber_sampling == "gp-direct"
    assert arguments.load_horizon == 1.0
    assert arguments.tbar_laplace == "field.npy"
    assert arguments.mpi_implementation == implementation


@pytest.mark.parametrize(
    ("profile", "node_aligned"),
    [
        (distributed_solver.DIRECT_SUPERLU_DIST_PROFILE, False),
        (distributed_solver.FGMRES_ASM_LU_PROFILE, False),
        (distributed_solver.FGMRES_ASM_ILU1_PROFILE, False),
        (distributed_solver.FGMRES_GAMG_RIGID_PROFILE, True),
        (distributed_solver.FGMRES_GAMG_RIGID_REBUILD_PROFILE, True),
    ],
)
def test_linear_solver_profiles_are_explicit_without_changing_physics_contract(
    profile, node_aligned
):
    parser = run_mpi._parser()
    arguments = parser.parse_args(
        [
            "--case",
            "B",
            "--mesh-topology",
            "closed-multiblock",
            "--formulation",
            "std-kappa",
            "--mass",
            "consistent",
            "--fiber-sampling",
            "gp-direct",
            "--tbar-laplace",
            "field.npy",
            "--linear-solver-profile",
            profile,
        ]
    )
    run_mpi._validate_arguments(parser, arguments)
    assert arguments.linear_solver_profile == profile
    assert arguments.mpi_implementation == (
        run_mpi.DistributedPetscSnesSolver.CLOSED_STD_KAPPA_IMPLEMENTATION
    )
    assert (
        distributed_solver.linear_solver_requires_node_aligned_ownership(profile)
        is node_aligned
    )


@pytest.mark.parametrize("case", ["A", "B"])
def test_fixed_load_horizon_gives_exact_serial_mpi_rank_gate_prefix(case):
    short = serial_driver._benchmark_load_histories(
        case, 0.001, 0.32, load_horizon=1.0
    )
    production = serial_driver._benchmark_load_histories(
        case, 0.001, 1.0, load_horizon=1.0
    )
    mpi_short = run_mpi.serial_driver._benchmark_load_histories(
        case, 0.001, 0.32, load_horizon=1.0
    )

    assert short[3] == production[3] == mpi_short[3] == 1.0
    for short_values, mpi_values, full_values in zip(
        short[:3], mpi_short[:3], production[:3]
    ):
        np.testing.assert_array_equal(short_values, mpi_values)
        np.testing.assert_array_equal(short_values, full_values[: len(short_values)])


@pytest.mark.parametrize(
    ("tend", "horizon"),
    [(0.32, 0.31), (0.32, 0.9995), (0.32, float("nan"))],
)
def test_mpi_parser_rejects_invalid_load_horizon(tend, horizon):
    parser = run_mpi._parser()
    arguments = parser.parse_args(
        ["--tend", str(tend), "--load-horizon", str(horizon)]
    )
    with pytest.raises(SystemExit):
        run_mpi._validate_arguments(parser, arguments)


@pytest.mark.parametrize("eta", [0.0, 99.0, 101.0, float("nan")])
def test_mpi_parser_rejects_physical_viscosity_changes(eta):
    parser = run_mpi._parser()
    arguments = parser.parse_args(["--material-eta", str(eta)])
    with pytest.raises(SystemExit):
        run_mpi._validate_arguments(parser, arguments)


def test_owned_mass_matrix_rejects_nonfinite_metadata_or_action():
    matrix = sp.eye(2, format="csr")
    owned = OwnedRowMassMatrix(matrix, 0, 2, 2, np.array([0]))
    with pytest.raises(ValueError, match="finite vector"):
        owned.action([0.0, np.nan])

    broken = copy.deepcopy(matrix)
    broken.data[0] = np.nan
    with pytest.raises(ValueError, match="finite and nonnegative"):
        OwnedRowMassMatrix(broken, 0, 2, 2, np.array([0]))


@pytest.mark.mpi
def test_distributed_solver_uses_off_diagonal_owned_mass_rows():
    from petsc4py import PETSc

    from distributed_solver import DistributedPetscSnesSolver

    if PETSc.COMM_WORLD.getSize() != 1:
        pytest.skip("focused owned-row solver test uses one-rank COMM_WORLD")

    class Batch:
        evaluation_mode = "split"
        material_residual_only_available = True

        def __init__(self):
            self.commits = 0

        def clear_cache(self):
            pass

        def element_residual_batch(self, coordinates, displacement, increment):
            target = np.array([1.0, 2.0])
            return displacement - target[:, None]

        def element_tangent_batch(self, coordinates, displacement, increment):
            return np.ones((2, 1, 1))

        def commit(self):
            self.commits += 1

    matrix = sp.csr_matrix([[2.0, 0.5], [0.5, 3.0]])
    mass = OwnedRowMassMatrix(matrix, 0, 2, 2, np.array([0, 1]))
    batch = Batch()
    solver = DistributedPetscSnesSolver(
        2,
        np.array([[0], [1]], dtype=int),
        np.zeros((2, 1)),
        batch,
        mass,
        dof_per_node=1,
    )
    try:
        dt = 0.5
        displacement, diagnostics = solver.solve_step(t=dt, dt=dt)
        expected = np.linalg.solve(
            np.eye(2) + matrix.toarray() / dt**2,
            np.array([1.0, 2.0]),
        )
        np.testing.assert_allclose(
            displacement, expected, rtol=2.0e-13, atol=2.0e-15
        )
        assert diagnostics.snes_converged_reason > 0
        assert diagnostics.final_residual_norm < (
            diagnostics.residual_acceptance_threshold
        )
        assert batch.commits == 1
        configuration = solver.configuration()
        assert configuration["mass_representation"] == "consistent_q1_hex8"
        assert configuration["mass_partition"] == (
            "owned-row-csr-all-touching-elements"
        )
        assert configuration["mass_local_nnz"] == 4
    finally:
        solver.close()


@pytest.mark.mpi
@pytest.mark.parametrize(
    "profile",
    [
        distributed_solver.FGMRES_ASM_LU_PROFILE,
        distributed_solver.FGMRES_ASM_ILU1_PROFILE,
        distributed_solver.FGMRES_GAMG_RIGID_PROFILE,
        distributed_solver.FGMRES_GAMG_RIGID_REBUILD_PROFILE,
    ],
)
def test_iterative_profiles_solve_and_report_tiny_vector_problem(profile):
    from petsc4py import PETSc

    if PETSc.COMM_WORLD.getSize() != 1:
        pytest.skip("focused iterative-profile test uses one-rank COMM_WORLD")

    class Batch:
        evaluation_mode = "split"
        material_residual_only_available = True

        def __init__(self, target):
            self.target = target
            self.commits = 0

        def clear_cache(self):
            pass

        def element_residual_batch(self, coordinates, displacement, increment):
            return displacement - self.target[None, :]

        def element_tangent_batch(self, coordinates, displacement, increment):
            return np.eye(24)[None, :, :]

        def commit(self):
            self.commits += 1

    nodes = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
    )
    target = np.linspace(-1.0, 1.0, 24) * 1.0e-3
    batch = Batch(target)
    solver = distributed_solver.DistributedPetscSnesSolver(
        24,
        np.arange(24).reshape(1, 24),
        nodes[None, :, :],
        batch,
        np.ones(24),
        linear_solver_profile=profile,
        node_coordinates=nodes,
    )
    try:
        displacement, diagnostics = solver.solve_step(t=0.1, dt=0.1)
        np.testing.assert_allclose(
            displacement, target / 101.0, rtol=2.0e-12, atol=2.0e-15
        )
        assert diagnostics.snes_converged_reason > 0
        assert diagnostics.linear_iterations > 0
        assert diagnostics.ksp_reasons_per_solve
        assert all(reason > 0 for reason in diagnostics.ksp_reasons_per_solve)
        assert batch.commits == 1
        configuration = solver.configuration()
        assert configuration["linear_solver_profile"] == profile
        assert configuration["ksp_type"] == "fgmres"
        if profile in {
            distributed_solver.FGMRES_GAMG_RIGID_PROFILE,
            distributed_solver.FGMRES_GAMG_RIGID_REBUILD_PROFILE,
        }:
            assert configuration["pc_type"] == "gamg"
            assert configuration["matrix_block_size"] == 3
            assert configuration["near_nullspace_mode_count"] == 6
            if profile == distributed_solver.FGMRES_GAMG_RIGID_REBUILD_PROFILE:
                assert configuration["ksp_max_it"] == 400
                assert configuration["gmres_restart"] == 100
                assert configuration["ksp_error_if_not_converged"] is False
                assert any(
                    key.endswith("pc_gamg_reuse_interpolation")
                    and value == "false"
                    for key, value in configuration["petsc_options"].items()
                )
        else:
            assert configuration["pc_type"] == "asm"
            assert configuration["near_nullspace_mode_count"] == 0
    finally:
        solver.close()
