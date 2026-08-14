"""Focused gates for the source-matched generalized-alpha MPI increment."""
from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

import geometry
import material
import run as serial_driver
import run_mpi
from distributed_local_pressure import FusedQ1P0Batch, RankLocalQ1P0PressureBatch
from distributed_mass import OwnedRowMassMatrix
from distributed_solver import DistributedPetscSnesSolver
from generalized_alpha import (
    GeneralizedAlphaParameters,
    SOURCE_MATCHED_GENERALIZED_ALPHA,
)


def _closed_arguments(case="A", formulation="local-pressure"):
    return [
        "--case",
        case,
        "--mesh-topology",
        "closed-multiblock",
        "--formulation",
        formulation,
        "--mass",
        "consistent",
        "--fiber-sampling",
        "gp-direct",
        "--tbar-laplace",
        "field.npy",
        "--integrator",
        "generalized-alpha",
    ]


def test_source_parameters_and_all_stage_chain_rules_are_exact():
    parameters = SOURCE_MATCHED_GENERALIZED_ALPHA
    assert parameters == GeneralizedAlphaParameters(
        alpha_m=0.2, alpha_f=0.4, gamma=0.7, beta=0.36
    )
    dt = 0.013
    u_n = np.array([-0.04, 0.02, 0.01])
    v_n = np.array([0.3, -0.2, 0.1])
    a_n = np.array([-0.5, 0.4, 0.2])
    u_np1 = np.array([-0.031, 0.015, 0.014])

    a_np1 = parameters.acceleration(u_np1, u_n, v_n, a_n, dt)
    v_np1 = parameters.velocity(a_np1, v_n, a_n, dt)
    np.testing.assert_allclose(
        parameters.force_stage(u_n, u_np1), 0.4 * u_n + 0.6 * u_np1
    )
    np.testing.assert_allclose(
        parameters.force_stage(v_n, v_np1), 0.4 * v_n + 0.6 * v_np1
    )
    np.testing.assert_allclose(
        parameters.acceleration_stage(a_n, a_np1), 0.2 * a_n + 0.8 * a_np1
    )
    assert parameters.load_time(0.2, dt) == pytest.approx(0.2 - 0.4 * dt)

    direction = np.array([0.7, -0.4, 0.2])
    epsilon = 1.0e-7

    def stages(displacement):
        acceleration = parameters.acceleration(
            displacement, u_n, v_n, a_n, dt
        )
        velocity = parameters.velocity(acceleration, v_n, a_n, dt)
        return (
            parameters.acceleration_stage(a_n, acceleration),
            parameters.force_stage(v_n, velocity),
        )

    plus = stages(u_np1 + epsilon * direction)
    minus = stages(u_np1 - epsilon * direction)
    acceleration_fd = (plus[0] - minus[0]) / (2.0 * epsilon)
    velocity_fd = (plus[1] - minus[1]) / (2.0 * epsilon)
    np.testing.assert_allclose(
        acceleration_fd,
        parameters.inertia_tangent(dt) * direction,
        rtol=2.0e-10,
    )
    np.testing.assert_allclose(
        velocity_fd,
        parameters.force_velocity_tangent(dt) * direction,
        rtol=2.0e-10,
    )


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"gamma": 0.6}, "gamma"),
        ({"beta": 0.35}, "beta"),
        (
            {"alpha_m": 0.1, "gamma": 0.8, "beta": 0.4225},
            "source-matched",
        ),
        ({"alpha_f": 1.0, "gamma": 1.3, "beta": 0.81}, "alpha"),
    ],
)
def test_non_source_generalized_alpha_parameters_fail_closed(kwargs, message):
    with pytest.raises(ValueError, match=message):
        GeneralizedAlphaParameters(**kwargs)


def test_mpi_load_grid_matches_simula_shifted_ode_not_endpoint_interpolation():
    dt = 0.001
    times, tau, pressure, horizon, evaluation_times = run_mpi._mpi_load_histories(
        "A", dt, 0.32, 1.0, integrator="generalized-alpha"
    )
    _solve, schedule, expected_horizon = serial_driver._load_schedule_grid(
        dt, 0.32, 1.0
    )
    shifted = schedule[1:] - 0.4 * dt
    expected_tau = serial_driver.tau_of_t(
        shifted, t_span=(0.0, expected_horizon)
    )[: len(times) - 1]

    assert horizon == expected_horizon == 1.0
    np.testing.assert_array_equal(pressure, np.zeros_like(pressure))
    np.testing.assert_array_equal(evaluation_times[1:], shifted[: len(times) - 1])
    np.testing.assert_array_equal(tau[1:], expected_tau)
    assert tau[0] == 0.0
    # Radau's adaptive path changes if it is incorrectly started at the first
    # shifted sample instead of the source's fixed t=0 initial condition.
    wrong_span = serial_driver.tau_of_t(shifted)[: len(times) - 1]
    assert np.max(np.abs(tau[1:] - wrong_span)) > 1.0

    be = run_mpi._mpi_load_histories("A", dt, 0.32, 1.0, integrator="be")
    retained = serial_driver._benchmark_load_histories("A", dt, 0.32, 1.0)
    for observed, expected in zip(be[:4], retained):
        np.testing.assert_array_equal(observed, expected)
    np.testing.assert_array_equal(be[4], be[0])


@pytest.mark.parametrize(
    "formulation", ["std-kappa", "local-pressure", "local-pressure-paper"]
)
@pytest.mark.parametrize("case", ["A", "B"])
def test_parser_labels_closed_step0_cases_as_source_matched(formulation, case):
    parser = run_mpi._parser()
    arguments = parser.parse_args(
        _closed_arguments(formulation=formulation, case=case)
    )
    run_mpi._validate_arguments(parser, arguments)
    expected = {
        "std-kappa": (
            run_mpi.DistributedPetscSnesSolver.
            CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION
        ),
        "local-pressure": (
            run_mpi.DistributedPetscSnesSolver.
            CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION
        ),
        "local-pressure-paper": (
            run_mpi.DistributedPetscSnesSolver.
            CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION
        ),
    }[formulation]
    assert arguments.mpi_implementation == expected


def test_parser_pins_generalized_alpha_to_source_one_second_load_horizon():
    parser = run_mpi._parser()
    arguments = parser.parse_args(
        _closed_arguments() + ["--tend", "0.2"]
    )
    run_mpi._validate_arguments(parser, arguments)
    assert arguments.load_horizon == 1.0

    arguments = parser.parse_args(
        _closed_arguments()
        + ["--tend", "0.2", "--load-horizon", "0.2"]
    )
    with pytest.raises(SystemExit):
        run_mpi._validate_arguments(parser, arguments)


def test_solver_rejects_a_generalized_alpha_label_on_a_be_material_batch():
    class Batch:
        time_integrator = "be"

        def element_residual_batch(self, coordinates, displacement, increment):
            return np.zeros_like(displacement)

        def element_tangent_batch(self, coordinates, displacement, increment):
            return np.zeros((1, 1, 1))

        def commit(self):
            pass

    with pytest.raises(ValueError, match="material batch time integrators disagree"):
        DistributedPetscSnesSolver(
            1,
            np.array([[0]]),
            np.zeros((1, 1)),
            Batch(),
            np.ones(1),
            dof_per_node=1,
            integrator="generalized-alpha",
            implementation=(
                DistributedPetscSnesSolver.
                CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION
            ),
        )


def test_solver_rejects_local_pressure_law_and_implementation_mismatch():
    class Batch:
        time_integrator = "be"
        pressure_law = "paper"

    with pytest.raises(ValueError, match="local-pressure law disagree"):
        DistributedPetscSnesSolver(
            1,
            np.array([[0]]),
            np.zeros((1, 1)),
            Batch(),
            np.ones(1),
            dof_per_node=1,
            implementation=(
                DistributedPetscSnesSolver.CLOSED_LOCAL_PRESSURE_IMPLEMENTATION
            ),
        )


class _ZeroMaterial:
    def __init__(self, count):
        self.props = np.zeros(1)
        self.count = count

    def element_r_batch(self, coordinates, displacement, increment):
        return np.zeros((self.count, 24))

    def element_rk_batch(self, coordinates, displacement, increment):
        return (
            np.zeros((self.count, 24)),
            np.zeros((self.count, 24, 24)),
        )

    def commit(self):
        pass


@pytest.mark.parametrize("pressure_law", ["log", "paper"])
@pytest.mark.parametrize("evaluation_mode", ["joint", "split"])
def test_condensed_local_pressure_uses_alpha_f_stage_and_endpoint_chain_rule(
    evaluation_mode, pressure_law,
):
    mesh = geometry.build_mesh(
        n_t=1, n_mu=2, n_theta=4, flip_helix=True, apex_offset=0.2
    )
    identifiers = np.array([1, 5], dtype=np.int64)
    coordinates = mesh.nodes[mesh.elems[identifiers]]
    rng = np.random.default_rng(20260803)
    old = rng.normal(scale=5.0e-6, size=(2, 24))
    endpoint = old + rng.normal(scale=2.0e-6, size=(2, 24))
    increment = endpoint - old
    alpha_f = SOURCE_MATCHED_GENERALIZED_ALPHA.alpha_f
    stage = alpha_f * old + (1.0 - alpha_f) * endpoint

    fused = FusedQ1P0Batch(
        _ZeroMaterial(len(identifiers)),
        coordinates,
        bulk_modulus=1.0e6,
        global_element_ids=identifiers,
        evaluation_mode=evaluation_mode,
        pressure_stage_alpha_f=alpha_f,
        pressure_law=pressure_law,
    )
    direct = RankLocalQ1P0PressureBatch(
        coordinates,
        bulk_modulus=1.0e6,
        global_element_ids=identifiers,
        pressure_law=pressure_law,
    )
    residual = fused.element_residual_batch(coordinates, endpoint, increment)
    tangent = fused.element_tangent_batch(coordinates, endpoint, increment)
    np.testing.assert_allclose(
        residual,
        direct.element_residual_batch(stage),
        rtol=2.0e-14,
        atol=2.0e-10,
    )
    np.testing.assert_allclose(
        tangent,
        (1.0 - alpha_f) * direct.element_tangent_batch(stage),
        rtol=2.0e-14,
        atol=2.0e-8,
    )


def test_material_uses_velocity_consistent_green_lagrange_rate_and_updates_state():
    parameters = SOURCE_MATCHED_GENERALIZED_ALPHA
    dt = 0.001
    fiber = np.array([0.8, 0.6, 0.0])
    sheet = np.array([-0.6, 0.8, 0.0])
    ff, ss, fssym = material.struct_tensors(fiber, sheet)
    grad_u_old = np.array(
        [[0.012, -0.003, 0.001], [0.002, -0.007, 0.003], [0.0, 0.001, 0.004]]
    )
    grad_v_old = np.array(
        [[0.08, -0.02, 0.01], [0.01, -0.04, 0.02], [0.0, 0.01, 0.03]]
    )
    grad_a_old = np.array(
        [[0.5, -0.1, 0.02], [0.04, -0.2, 0.03], [0.0, 0.02, 0.1]]
    )
    grad_u_new = np.array(
        [[0.014, -0.002, 0.0], [0.003, -0.006, 0.004], [0.001, 0.0, 0.005]]
    )
    F_new = np.eye(3) + grad_u_new

    with_viscosity = material.HolzapfelOgdenActiveGeneralizedAlpha()
    without_viscosity = material.HolzapfelOgdenActiveGeneralizedAlpha()
    without_viscosity.eta = 0.0
    arguments = (
        F_new.astype(complex),
        ff.astype(complex),
        ss.astype(complex),
        fssym.astype(complex),
        grad_u_old.astype(complex),
        grad_v_old.astype(complex),
        grad_a_old.astype(complex),
        dt,
    )
    stress, state = with_viscosity.stress_PK1(*arguments)
    elastic_stress, _ = without_viscosity.stress_PK1(*arguments)

    grad_a_new = (
        grad_u_new
        - grad_u_old
        - dt * grad_v_old
        - (0.5 - parameters.beta) * dt**2 * grad_a_old
    ) / (parameters.beta * dt**2)
    grad_v_new = grad_v_old + dt * (
        (1.0 - parameters.gamma) * grad_a_old
        + parameters.gamma * grad_a_new
    )
    grad_u_stage = parameters.alpha_f * grad_u_old + (
        1.0 - parameters.alpha_f
    ) * grad_u_new
    grad_v_stage = parameters.alpha_f * grad_v_old + (
        1.0 - parameters.alpha_f
    ) * grad_v_new
    F_stage = np.eye(3) + grad_u_stage
    strain_rate = 0.5 * (
        F_stage.T @ grad_v_stage + grad_v_stage.T @ F_stage
    )
    expected_viscous_pk1 = F_stage @ (with_viscosity.eta * strain_rate)

    np.testing.assert_allclose(
        np.asarray(stress - elastic_stress).real,
        expected_viscous_pk1,
        rtol=3.0e-13,
        atol=2.0e-7,
    )
    np.testing.assert_allclose(np.asarray(state["grad_u_prev"]).real, grad_u_new)
    np.testing.assert_allclose(np.asarray(state["grad_v_prev"]).real, grad_v_new)
    np.testing.assert_allclose(np.asarray(state["grad_a_prev"]).real, grad_a_new)
    assert material.verify_generalized_alpha_material_tangent(verbose=False)


@pytest.mark.slow
def test_compiled_generalized_alpha_hex8_directional_tangent_and_commit(tmp_path):
    from coupfe.runtime.compiled_element import CompiledElement

    coordinates = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
        dtype=float,
    )[None, :, :]
    problem = material.CardiacHex8GeneralizedAlpha()
    schema = problem._mat.state_schema
    per_gp = sum(entry["size"] for entry in schema.values())
    properties = np.asarray(problem._mat.props_array, dtype=float)
    _source, module = material.build_generalized_alpha_kernel(
        tmpdir=str(tmp_path), module_name="cardiac_ga_focused_test"
    )
    element = CompiledElement(
        module,
        props=properties,
        dof_per_node=3,
        n_svars=8 * per_gp,
        mcrd=3,
        n_elem=1,
        dt=0.001,
        backend="native",
        state_schema=schema,
    )

    gradient = np.array(
        [[2.0e-3, -4.0e-4, 2.0e-4], [3.0e-4, -1.0e-3, 5.0e-4], [0.0, 2.0e-4, 8.0e-4]]
    )
    endpoint = np.einsum("eai,ji->eaj", coordinates, gradient).reshape(1, 24)
    rng = np.random.default_rng(20260803)
    direction = rng.normal(size=(1, 24))
    direction /= np.linalg.norm(direction)
    residual, tangent = element.element_rk_batch(
        coordinates, endpoint, endpoint
    )
    committed_before = element.svars.copy()

    epsilon = 2.0e-8
    plus = endpoint + epsilon * direction
    minus = endpoint - epsilon * direction
    residual_plus = element.element_rk_batch(coordinates, plus, plus)[0]
    residual_minus = element.element_rk_batch(coordinates, minus, minus)[0]
    finite_difference = (residual_plus - residual_minus) / (2.0 * epsilon)
    np.testing.assert_allclose(
        finite_difference,
        np.einsum("eij,ej->ei", tangent, direction),
        rtol=2.0e-5,
        atol=2.0e-3,
    )

    # Refresh trial state at the accepted endpoint before the transactional copy.
    accepted_residual = element.element_rk_batch(
        coordinates, endpoint, endpoint
    )[0]
    trial = element.svars_trial.copy()
    np.testing.assert_allclose(accepted_residual, residual, rtol=2.0e-13)
    np.testing.assert_array_equal(element.svars, committed_before)
    assert np.linalg.norm(
        trial[:, schema["grad_u_prev"]["offset"] :]
    ) > 0.0
    element.commit()
    np.testing.assert_array_equal(element.svars, trial)


@pytest.mark.slow
@pytest.mark.parametrize(
    ("formulation", "pressure_law"),
    [("local-pressure", "log"), ("local-pressure-paper", "paper")],
)
def test_mpi_builder_fuses_source_staged_material_and_local_pressure(
    tmp_path, formulation, pressure_law
):
    mesh = geometry.build_closed_mesh(
        n_t=1,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
        flip_helix=True,
    )
    identifiers = np.array([0, mesh.n_elem - 1], dtype=np.int64)
    coordinates = mesh.nodes[mesh.elems[identifiers]]
    batch, element, tension_index, bulk_modulus = (
        run_mpi._build_rank_local_batch(
            mesh,
            identifiers,
            coordinates,
            0.001,
            tmp_path,
            0,
            formulation=formulation,
            element_evaluation="joint",
            fiber_sampling="gp-direct",
            fiber_tbar=np.asarray(mesh.param[:, 0]),
            fiber_asign=-1.0,
            integrator="generalized-alpha",
        )
    )
    assert batch.time_integrator == "generalized-alpha-source-matched"
    assert batch.pressure_stage_alpha_f == 0.4
    assert batch.pressure_law == pressure_law
    assert bulk_modulus == 1.0e6
    assert element.props[tension_index] == 0.0
    schema = material.CardiacHex8GeneralizedAlpha()._mat.state_schema
    assert set(schema) == {
        "ff",
        "ss",
        "fssym",
        "grad_u_prev",
        "grad_v_prev",
        "grad_a_prev",
    }

    displacement = np.zeros((len(identifiers), 24))
    residual, tangent = batch.element_rk_batch(
        coordinates, displacement, displacement
    )
    np.testing.assert_allclose(residual, 0.0, atol=2.0e-10)
    assert tangent.shape == (len(identifiers), 24, 24)
    assert np.all(np.isfinite(tangent))
    np.testing.assert_allclose(
        batch.deformation_jacobians(displacement), 1.0, atol=2.0e-15
    )
    batch.commit()


@pytest.mark.mpi
def test_distributed_solver_stages_consistent_inertia_and_robin_then_commits():
    from petsc4py import PETSc

    if PETSc.COMM_WORLD.getSize() != 1:
        pytest.skip("focused generalized-alpha solver test uses one-rank COMM_WORLD")

    class LoadBatch:
        evaluation_mode = "split"
        material_residual_only_available = True
        time_integrator = "generalized-alpha-source-matched"
        pressure_stage_alpha_f = 0.4
        pressure_law = "log"

        def __init__(self):
            self.force = 0.0
            self.commits = 0

        def clear_cache(self):
            pass

        def element_residual_batch(self, coordinates, displacement, increment):
            return np.array([[-self.force]])

        def element_tangent_batch(self, coordinates, displacement, increment):
            return np.zeros((1, 1, 1))

        def commit(self):
            self.commits += 1

    class ScalarRobin:
        def __init__(self, stiffness, damping):
            self.Kmat = sp.csr_matrix([[stiffness]])
            self.Cmat = sp.csr_matrix([[damping]])
            self.dofs = np.array([0])
            self.u_prev = np.zeros(1)

        def commit(self, displacement, state, t, dt):
            self.u_prev = np.asarray(displacement, dtype=float).copy()
            return state

    mass_value = 2.3
    stiffness = 7.0
    damping = 0.45
    dt = 0.08
    mass = OwnedRowMassMatrix(
        sp.csr_matrix([[mass_value]]), 0, 1, 1, np.array([0])
    )
    batch = LoadBatch()
    batch.material_dt = dt
    robin = ScalarRobin(stiffness, damping)
    solver = DistributedPetscSnesSolver(
        1,
        np.array([[0]], dtype=int),
        np.zeros((1, 1)),
        batch,
        mass,
        dof_per_node=1,
        robin=robin,
        integrator="generalized-alpha",
        implementation=(
            DistributedPetscSnesSolver.
            CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION
        ),
    )
    parameters = SOURCE_MATCHED_GENERALIZED_ALPHA

    def expected_endpoint(u_n, v_n, a_n, force):
        def residual(endpoint):
            endpoint = np.array([endpoint])
            acceleration = parameters.acceleration(
                endpoint, u_n, v_n, a_n, dt
            )
            velocity = parameters.velocity(acceleration, v_n, a_n, dt)
            return float(
                mass_value
                * parameters.acceleration_stage(a_n, acceleration)[0]
                + stiffness * parameters.force_stage(u_n, endpoint)[0]
                + damping * parameters.force_stage(v_n, velocity)[0]
                - force
            )

        zero = residual(0.0)
        slope = residual(1.0) - zero
        return -zero / slope

    try:
        expected_u = np.zeros(1)
        expected_v = np.zeros(1)
        expected_a = np.zeros(1)
        for step, force in enumerate((3.2, -0.7), start=1):
            batch.force = force
            endpoint = expected_endpoint(
                expected_u, expected_v, expected_a, force
            )
            new_u = np.array([endpoint])
            new_a = parameters.acceleration(
                new_u, expected_u, expected_v, expected_a, dt
            )
            new_v = parameters.velocity(new_a, expected_v, expected_a, dt)

            displacement, diagnostics = solver.solve_step(
                t=step * dt, dt=dt
            )
            np.testing.assert_allclose(displacement, new_u, rtol=3.0e-13)
            np.testing.assert_allclose(solver._v_prev, new_v, rtol=3.0e-13)
            np.testing.assert_allclose(solver._a_prev, new_a, rtol=3.0e-13)
            assert diagnostics.snes_converged_reason > 0
            expected_u, expected_v, expected_a = new_u, new_v, new_a

        assert batch.commits == 2
        configuration = solver.configuration()
        assert configuration["time_integrator"] == "generalized-alpha"
        assert configuration["generalized_alpha"]["alpha_m"] == 0.2
        assert configuration["generalized_alpha"]["alpha_f"] == 0.4
        assert configuration["acceleration_stage"] == "1-alpha_m"
        assert configuration["force_stage"] == "1-alpha_f"
        assert configuration["local_pressure_stage_alpha_f"] == 0.4
        assert configuration["compiled_material_dt"] == dt
        with pytest.raises(
            ValueError, match="differs from compiled material dt"
        ):
            solver.solve_step(t=3.0 * dt, dt=2.0 * dt)
        assert batch.commits == 2
    finally:
        solver.close()
