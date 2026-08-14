from __future__ import annotations

import hashlib
import json
import pickle
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp

import activation
import fiber_crosscheck
import geometry
import local_pressure
import material
import post
import pressure
import result_io
import run as cardiac_run
import solver
from newmark import NewmarkInertia
from robin import RobinOperator, _assemble
from coupfe.operators.base import Residual, Tangent


ROOT = Path(__file__).resolve().parents[1]
RETAINED_CASE_A = (
    ROOT
    / "examples"
    / "cardiac_benchmark"
    / "results"
    / "archive"
    / "truncated_polar"
    / "case_a"
    / "case_a_reduced.json"
)


def _assembled_residual(operator, displacement, dt=1.0):
    contribution = operator.residual(displacement, None, 0.0, dt)
    result = np.zeros(operator.ndof)
    np.add.at(result, contribution.gdofs, contribution.values)
    return result


class _ScalarOperator:
    def __init__(self, *, stubborn=False):
        self.stubborn = stubborn
        self.commits = 0

    def residual(self, U, state, t, dt):
        value = 1.0 if self.stubborn else float(U[0] - 1.0)
        return Residual(np.array([0]), np.array([value]))

    def tangent(self, U, state, t, dt):
        return Tangent(np.array([0]), np.array([0]), np.array([1.0]))

    def commit(self, U, state, t, dt):
        self.commits += 1
        return state


class _TrialTrackingScalarOperator(_ScalarOperator):
    def __init__(self):
        super().__init__()
        self.trial = None

    def residual(self, U, state, t, dt):
        self.trial = float(U[0])
        return super().residual(U, state, t, dt)

    def commit(self, U, state, t, dt):
        assert self.trial == pytest.approx(float(U[0]))
        return super().commit(U, state, t, dt)


class _NonfiniteAtRootOperator(_ScalarOperator):
    def residual(self, U, state, t, dt):
        if np.isclose(float(U[0]), 1.0):
            return Residual(np.array([0]), np.array([float("nan")]))
        return super().residual(U, state, t, dt)


class _FakeVec:
    def __init__(self, values):
        self.values = np.asarray(values, dtype=float)

    def getArray(self, readonly=False):
        return self.values


class _FakeKsp:
    def __init__(self, reason):
        self.reason = reason

    def getConvergedReason(self):
        return self.reason


class _FakeSnes:
    def __init__(
        self,
        result,
        *,
        snes_reason=2,
        ksp_reason=4,
        nonlinear_iterations=1,
        linear_iterations=1,
        function_norm=0.0,
    ):
        self.result = np.asarray(result, dtype=float)
        self.snes_reason = snes_reason
        self.ksp = _FakeKsp(ksp_reason)
        self.nonlinear_iterations = nonlinear_iterations
        self.linear_iterations = linear_iterations
        self.function_norm = function_norm

    def solve(self, _rhs, vector):
        vector.values = self.result.copy()

    def getConvergedReason(self):
        return self.snes_reason

    def getIterationNumber(self):
        return self.nonlinear_iterations

    def getLinearSolveIterations(self):
        return self.linear_iterations

    def getFunctionNorm(self):
        return self.function_norm

    def getKSP(self):
        return self.ksp


def _fake_petsc_solver(monkeypatch, result, **snes_options):
    application_solver = solver.PetscSnesSolver()
    vector = _FakeVec(np.zeros(1))
    snes = _FakeSnes(result, **snes_options)

    def ensure_context(_operators, _U0, _state, ndof, *, t, dt):
        assert ndof == 1
        application_solver._ndof = ndof
        application_solver._x = vector
        application_solver._snes = snes

    monkeypatch.setattr(application_solver, "_ensure_context", ensure_context)
    return application_solver


def _write_qualified_result(path, *, case="A", **overrides):
    times = np.array([0.0, 0.5, 1.0])
    history = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, -1.0], [2.0, 0.0, -2.0]])
    payload = {
        "result_schema": "coupfe-cardiac-result-v1",
        "case": case,
        "integrator": "newmark",
        "formulation": "hex8_fbar",
        "fiber_sampling": "cg1_gram_schmidt",
        "point_sampling": "global_delaunay_tetra",
        "viscous_rate": "backward_difference",
        "dt": 0.5,
        "t_end": 1.0,
        "apex_offset": 0.2,
        "n_t": 1,
        "n_mu": 2,
        "n_theta": 4,
        "flip_helix": True,
        "converged": True,
        "completed_steps": 2,
        "expected_steps": 2,
        "app_revision": "a" * 40,
        "app_tree_state": "clean",
        "app_source_kind": "git-checkout",
        "core_revision": "454f73ce2de284262b214a2b37bd676c6aca3c0a",
        "core_tree_state": "clean",
        "core_source_kind": "git-checkout",
        "core_source_url": "https://github.com/tengzhang48/CoupFE.git",
        "times": times,
        "u0": history,
        "u1": history,
    }
    payload.update(overrides)
    np.savez(path, **payload)


def test_activation_and_pressure_histories_match_benchmark_scales():
    times = np.linspace(0.0, 1.0, 1001)
    tension = activation.tau_of_t(times)
    cavity_pressure = activation.p_of_t(times)

    assert tension.shape == times.shape
    assert cavity_pressure.shape == times.shape
    assert np.all(np.isfinite(tension))
    assert np.all(np.isfinite(cavity_pressure))
    assert tension[0] == 0.0
    assert cavity_pressure[0] == 0.0
    assert 1.15e5 < tension.max() < 1.22e5
    assert 1.55e4 < cavity_pressure.max() < 1.68e4


def test_material_kernel_rejects_unknown_formulation_before_codegen():
    with pytest.raises(ValueError, match="fbar_mechanics.*standard"):
        material.build_kernel(formulation="not-a-formulation")


def test_retained_case_a_result_is_traceable_and_semantically_consistent():
    evidence_bytes = RETAINED_CASE_A.read_bytes()
    assert hashlib.sha256(evidence_bytes).hexdigest() == (
        "3b97d6cede4c98a2d79e373f676029d8a4286f9a96e986942c9991f48f2c58bf"
    )
    evidence = json.loads(evidence_bytes)

    assert evidence["schema"] == "coupfe-cardiac-retained-result-v1"
    assert evidence["case"] == "A"
    assert evidence["evidence_label"] == "checked reduced execution demonstration"
    assert evidence["license"] == "CC-BY-4.0"

    source = evidence["source"]
    assert source["app_revision"] == "44cbfed9e09d4150203faae3087f2e4617d1fc47"
    assert source["core_revision"] == "454f73ce2de284262b214a2b37bd676c6aca3c0a"
    assert source["app_tree_state"] == source["core_tree_state"] == "clean"
    assert source["core_url"] == "https://github.com/tengzhang48/CoupFE.git"
    source_hashes = source["source_files_sha256"]
    assert "examples/cardiac_benchmark/run.py" in source_hashes
    assert "examples/cardiac_benchmark/result_io.py" in source_hashes
    for relative_name, expected_hash in source_hashes.items():
        relative_path = Path(relative_name)
        assert relative_path.parts[:2] == ("examples", "cardiac_benchmark")
        assert len(expected_hash) == 64
        int(expected_hash, 16)

    artifacts = evidence["artifacts"]
    stdout_name = artifacts["normalized_stdout"]
    assert Path(stdout_name).name == stdout_name
    stdout_path = RETAINED_CASE_A.parent / stdout_name
    stdout_bytes = stdout_path.read_bytes()
    assert hashlib.sha256(stdout_bytes).hexdigest() == artifacts[
        "normalized_stdout_sha256"
    ]
    transcript = stdout_bytes.decode("utf-8")
    assert "mesh: 24 nodes, 8 hexes, ndof=72" in transcript
    assert "step  500 t=1.000s" in transcript
    assert "finished 500/500 steps" in transcript
    assert "elapsed " not in transcript
    assert artifacts["generated_npz_distributed"] is False
    assert len(artifacts["generated_npz_sha256"]) == 64
    int(artifacts["generated_npz_sha256"], 16)

    configuration = evidence["configuration"]
    time_config = configuration["time"]
    history = evidence["retained_history"]
    samples = history["samples"]
    stride = int(history["source_step_stride"])
    assert time_config["integrator"] == "be"
    assert time_config["completed_steps"] == time_config["expected_steps"] == 500
    assert stride == 10
    assert len(samples) == time_config["expected_steps"] // stride + 1
    assert history["sample_interval_s"] == pytest.approx(
        stride * time_config["dt_s"]
    )

    times = np.array([sample["time_s"] for sample in samples])
    tension = np.array([sample["active_tension_pa"] for sample in samples])
    cavity_pressure = np.array([sample["pressure_pa"] for sample in samples])
    u0 = np.asarray([sample["u0_m"] for sample in samples], dtype=float)
    u1 = np.asarray([sample["u1_m"] for sample in samples], dtype=float)
    np.testing.assert_allclose(
        times,
        np.linspace(0.0, time_config["t_end_s"], len(samples)),
        rtol=0.0,
        atol=1.0e-15,
    )
    tolerances = evidence["regression_tolerances"]
    assert all(float(value) > 0.0 for value in tolerances.values())
    np.testing.assert_allclose(
        tension,
        activation.tau_of_t(times),
        rtol=tolerances["load_relative"],
        atol=tolerances["load_absolute_pa"],
    )
    np.testing.assert_array_equal(cavity_pressure, np.zeros_like(cavity_pressure))
    assert u0.shape == u1.shape == (len(samples), 3)
    assert np.all(np.isfinite(u0)) and np.all(np.isfinite(u1))
    np.testing.assert_array_equal(u0[0], np.zeros(3))
    np.testing.assert_array_equal(u1[0], np.zeros(3))

    result = evidence["result"]
    peak_index = int(np.argmax(tension))
    assert result["converged"] is True
    assert times[peak_index] == pytest.approx(result["peak_time_s"])
    assert stride * peak_index == result["peak_step"]
    assert tension[peak_index] == pytest.approx(result["peak_active_tension_pa"])
    np.testing.assert_allclose(
        u0[peak_index], result["u0_at_peak_m"], rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(
        u1[peak_index], result["u1_at_peak_m"], rtol=0.0, atol=0.0
    )
    np.testing.assert_allclose(u0[-1], result["u0_final_m"], rtol=0.0, atol=0.0)
    np.testing.assert_allclose(u1[-1], result["u1_final_m"], rtol=0.0, atol=0.0)
    assert np.linalg.norm(u0[peak_index]) > np.linalg.norm(u0[-1])
    assert np.linalg.norm(u1[peak_index]) > np.linalg.norm(u1[-1])


def test_fiber_formula_crosscheck_is_an_executable_gate():
    unflipped, flipped = fiber_crosscheck.compare_fiber_formulas()
    assert len(flipped) == 480
    assert np.mean(flipped > 0.999) > 0.99
    assert np.mean(unflipped > 0.999) < 0.5


def test_noncollapsed_mesh_has_positive_volume_and_orthonormal_fibers():
    mesh = geometry.build_mesh(n_t=1, n_mu=4, n_theta=8, apex_offset=0.12)
    volumes = np.array([geometry._hex_volume(mesh.nodes[element]) for element in mesh.elems])
    assert np.all(volumes > 0.0)
    assert np.all(np.isfinite(mesh.nodes))

    for frame in zip(mesh.fiber, mesh.sheet, mesh.normal):
        basis = np.stack(frame)
        np.testing.assert_allclose(basis @ basis.T, np.eye(3), atol=2e-14)


def test_structural_tensors_reorthogonalize_interpolated_directions():
    ff, ss, fssym = material.struct_tensors(
        np.array([1.0, 0.2, 0.0]), np.array([0.4, 1.0, 0.3])
    )
    assert np.trace(ff @ ss) == pytest.approx(0.0, abs=1.0e-14)
    assert np.trace(ff) == pytest.approx(1.0)
    assert np.trace(ss) == pytest.approx(1.0)
    assert np.all(np.isfinite(fssym))

    model = material.HolzapfelOgdenActive()
    stress, _ = model.stress_PK1(
        np.eye(3) + 0j,
        ff + 0j,
        ss + 0j,
        fssym + 0j,
        np.zeros((3, 3), dtype=complex),
        1.0e-3,
    )
    np.testing.assert_allclose(np.asarray(stress).real, 0.0, atol=1.0e-10)


def test_material_stress_is_derivative_of_declared_benchmark_energy():
    """Check PK1 against the cited smooth-switch energy, not against itself."""
    model = material.HolzapfelOgdenActive()
    fiber = np.array([0.8, 0.6, 0.0])
    sheet = np.array([-0.6, 0.8, 0.0])
    ff, ss, fssym = material.struct_tensors(fiber, sheet)
    deformation = np.array(
        [[1.015, 0.012, -0.004], [0.003, 0.992, 0.008], [0.0, -0.002, 1.006]]
    )
    identity = np.eye(3)
    strain = 0.5 * (deformation.T @ deformation - identity)

    def energy(F):
        C = F.T @ F
        J = np.linalg.det(F)
        I1_bar = J ** (-2.0 / 3.0) * np.trace(C)
        I4f = np.trace(C @ ff)
        I4s = np.trace(C @ ss)
        I8fs = np.trace(C @ fssym)

        def directional(I4, a, b):
            switch = 1.0 / (1.0 + np.exp(-model.k_sw * (I4 - 1.0)))
            return (
                a
                / (2.0 * b)
                * switch
                * (np.exp(b * (I4 - 1.0) ** 2) - 1.0)
            )

        return (
            model.a / (2.0 * model.b) * (np.exp(model.b * (I1_bar - 3.0)) - 1.0)
            + directional(I4f, model.a_f, model.b_f)
            + directional(I4s, model.a_s, model.b_s)
            + model.a_fs
            / (2.0 * model.b_fs)
            * (np.exp(model.b_fs * I8fs**2) - 1.0)
            + 0.25 * model.kappa * (J * J - 1.0 - 2.0 * np.log(J))
        )

    step = 1.0e-20
    energy_gradient = np.empty((3, 3))
    for row in range(3):
        for column in range(3):
            perturbed = deformation.astype(complex)
            perturbed[row, column] += 1j * step
            energy_gradient[row, column] = np.imag(energy(perturbed)) / step

    stress, _ = model.stress_PK1(
        deformation.astype(complex),
        ff.astype(complex),
        ss.astype(complex),
        fssym.astype(complex),
        strain.astype(complex),
        1.0e-3,
    )
    np.testing.assert_allclose(
        np.asarray(stress).real, energy_gradient, rtol=2.0e-11, atol=2.0e-8
    )


def test_robin_matrix_is_symmetric_and_has_analytic_uniform_reaction():
    nodes = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    facets = np.array([[0, 1, 2, 3]])
    alpha = 3.0
    operator = RobinOperator(nodes, 12, [(facets, alpha, 0.0, "full")])
    displacement = np.zeros(12)
    displacement[0::3] = 2.0e-3

    reaction = _assembled_residual(operator, displacement)
    np.testing.assert_allclose(reaction.reshape(-1, 3).sum(axis=0), [alpha * 2e-3, 0, 0])
    np.testing.assert_allclose(operator.Kmat.toarray(), operator.Kmat.toarray().T, atol=1e-14)


def test_robin_smoothed_mode_is_bit_identical_on_a_flat_facet():
    # on a flat facet every node receives the same facet normal, so the
    # area-weighted smoothed normal field equals the per-Gauss facet normal
    # and the assembled matrices are bit-identical to the default mode.
    nodes = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    facets = np.array([[0, 1, 2, 3]])
    k_facet, c_facet, _, _, _, _ = _assemble(nodes, [(facets, 3.0, 5.0, "normal")], 12)
    k_smooth, c_smooth, _, _, _, _ = _assemble(
        nodes, [(facets, 3.0, 5.0, "normal-smoothed")], 12
    )
    np.testing.assert_array_equal(k_facet.toarray(), k_smooth.toarray())
    np.testing.assert_array_equal(c_facet.toarray(), c_smooth.toarray())


def test_robin_smoothed_k_is_symmetric_and_psd_on_the_epicardium():
    mesh = geometry.build_mesh(n_t=2, n_mu=10, n_theta=16, apex_offset=0.2)
    ndof = 3 * mesh.n_node
    operator = RobinOperator(
        mesh.nodes,
        ndof,
        [
            (mesh.facets_base, 1.0e5, 5.0e3, "full"),
            (mesh.facets_epi, 1.0e8, 5.0e3, "normal-smoothed"),
        ],
    )
    k = operator.Kmat.toarray()
    np.testing.assert_allclose(k, k.T, atol=1e-9)
    eig = np.linalg.eigvalsh(k)
    assert eig.min() >= -1.0e-6 * max(1.0, abs(eig.max()))


def test_robin_smoothed_mode_collapses_the_long_axis_twist_restraint():
    # the candidate selector for the snap-window family is the spurious
    # epicardial restraint of long-axis (rx) rotation.  the smoothed normal
    # field must drive that quadratic form toward zero, leaving it as a
    # mechanism diagnostic rather than a benchmark-faithful operator.
    mesh = geometry.build_mesh(n_t=2, n_mu=20, n_theta=17, apex_offset=0.2)
    ndof = 3 * mesh.n_node

    def rx_quadratic_form(mode: str) -> float:
        operator = RobinOperator(
            mesh.nodes,
            ndof,
            [
                (mesh.facets_base, 1.0e5, 5.0e3, "full"),
                (mesh.facets_epi, 1.0e8, 5.0e3, mode),
            ],
        )
        column = np.zeros(ndof)
        for node in range(mesh.n_node):
            column[3 * node : 3 * node + 3] = np.cross(
                np.array([1.0, 0.0, 0.0]), mesh.nodes[node]
            )
        column /= np.linalg.norm(column)
        return float(column @ operator.Kmat @ column)

    facet_rx = rx_quadratic_form("normal")
    smoothed_rx = rx_quadratic_form("normal-smoothed")
    assert smoothed_rx < 0.05 * facet_rx


def test_robin_rejects_unknown_projection_mode():
    nodes = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]])
    facets = np.array([[0, 1, 2, 3]])
    with pytest.raises(ValueError):
        RobinOperator(nodes, 12, [(facets, 1.0, 1.0, "not-a-mode")])


def test_cardiac_run_exposes_epicardial_normal_mode_option():
    parser = cardiac_run._parser()
    parsed = parser.parse_args(
        ["--case", "A", "--epicardial-normal-mode", "nodal-smoothed", "--out", "x.npz"]
    )
    assert parsed.epicardial_normal_mode == "nodal-smoothed"
    assert parser.parse_args(["--case", "A"]).epicardial_normal_mode == "facet"
    with pytest.raises(SystemExit):
        parser.parse_args(["--case", "A", "--epicardial-normal-mode", "bogus"])


def test_follower_pressure_orientation_batch_and_tangent_broken_control():
    nodes = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    facets = np.array([[0, 1, 2, 3]])
    interior = np.array([[0.5, 0.5, 1.0]])
    operator = pressure.FollowerPressureOperator(
        nodes, 12, facets, p=2000.0, interior=interior
    )

    zero = np.zeros(12)
    result = _assembled_residual(operator, zero).reshape(-1, 3).sum(axis=0)
    np.testing.assert_allclose(result, [0.0, 0.0, -2000.0], atol=2e-12)

    scalar = pressure._facet_residual(nodes, zero, operator.p, operator._orient[0])
    batch = pressure._facet_residual_batch(
        nodes[None, :, :], zero.reshape(1, 4, 3), operator.p, operator._orient
    )[0]
    np.testing.assert_allclose(batch, scalar, atol=0.0)

    # Deliberately reverse the normal: the entire load must reverse.
    broken = pressure._facet_residual_batch(
        nodes[None, :, :], zero.reshape(1, 4, 3), operator.p, -operator._orient
    )[0]
    np.testing.assert_allclose(broken, -batch, atol=0.0)

    rng = np.random.default_rng(7)
    displacement = rng.normal(scale=2e-3, size=12)
    tangent = operator.tangent(displacement, None, 0.0, 1.0)
    matrix = sp.coo_matrix(
        (tangent.values, (tangent.rows, tangent.cols)), shape=(12, 12)
    ).toarray()
    eps = 1e-7
    finite_difference = np.column_stack(
        [
            (
                _assembled_residual(operator, displacement + eps * np.eye(12)[column])
                - _assembled_residual(operator, displacement - eps * np.eye(12)[column])
            )
            / (2 * eps)
            for column in range(12)
        ]
    )
    np.testing.assert_allclose(matrix, finite_difference, rtol=2e-8, atol=2e-7)


def test_follower_pressure_requires_parent_element_centroids():
    nodes = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    with pytest.raises(ValueError, match="parent-element centroids are required"):
        pressure.FollowerPressureOperator(
            nodes, 12, np.array([[0, 1, 2, 3]]), p=1.0
        )


def test_newmark_inertia_contract():
    operator = NewmarkInertia(np.array([2.0, 0.0, 3.0]), 3)
    displacement = np.array([1e-3, 9.0, -2e-3])
    dt = 0.1
    residual = operator.residual(displacement, None, 0.0, dt)
    np.testing.assert_array_equal(residual.gdofs, [0, 2])
    np.testing.assert_allclose(residual.values, [0.8, -2.4])
    tangent = operator.tangent(displacement, None, 0.0, dt)
    np.testing.assert_allclose(tangent.values, [800.0, 1200.0])
    np.testing.assert_allclose(operator.velocity(displacement, dt), 20.0 * displacement)
    assert operator.velocity_tangent(dt) == pytest.approx(20.0)


def test_robin_dashpot_uses_newmark_velocity_and_consistent_tangent():
    nodes = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
    )
    facets = np.array([[0, 1, 2, 3]])
    kinematics = NewmarkInertia(np.zeros(12), 12)
    operator = RobinOperator(
        nodes,
        12,
        [(facets, 3.0, 2.0, "full")],
        kinematics=kinematics,
    )
    displacement = np.linspace(-1.0e-3, 1.0e-3, 12)
    dt = 0.1
    tangent = operator.tangent(displacement, None, 0.0, dt)
    matrix = sp.coo_matrix(
        (tangent.values, (tangent.rows, tangent.cols)), shape=(12, 12)
    ).toarray()
    eps = 1.0e-7
    finite_difference = np.column_stack(
        [
            (
                _assembled_residual(operator, displacement + eps * np.eye(12)[column], dt)
                - _assembled_residual(operator, displacement - eps * np.eye(12)[column], dt)
            )
            / (2.0 * eps)
            for column in range(12)
        ]
    )
    np.testing.assert_allclose(matrix, finite_difference, rtol=2.0e-10, atol=2.0e-10)


def test_checked_newton_rejects_exhaustion_before_physical_commit():
    operator = _ScalarOperator(stubborn=True)
    with pytest.raises(RuntimeError, match="did not converge before state commit"):
        solver.checked_newton_solve(
            [operator], np.zeros(1), None, 1, {}, rtol=1.0e-10, maxit=2
        )
    assert operator.commits == 0


def test_checked_newton_accepts_a_residual_converged_iterate():
    operator = _ScalarOperator()
    displacement, _, _ = solver.checked_newton_solve(
        [operator], np.zeros(1), None, 1, {}, rtol=1.0e-10, maxit=3
    )
    np.testing.assert_allclose(displacement, [1.0])
    assert operator.commits == 1


def test_checked_newton_rejects_a_noncontractive_relative_tolerance():
    with pytest.raises(ValueError, match="0 < rtol < 1"):
        solver.checked_newton_solve(
            [_ScalarOperator(stubborn=True)],
            np.zeros(1),
            None,
            1,
            {},
            rtol=2.0,
            maxit=1,
        )


@pytest.mark.parametrize("reason", [0, -3])
def test_petsc_snes_rejects_nonpositive_reason_before_commit(monkeypatch, reason):
    operator = _ScalarOperator()
    application_solver = _fake_petsc_solver(
        monkeypatch, [1.0], snes_reason=reason
    )

    with pytest.raises(solver.SnesSolveError, match="did not report convergence"):
        application_solver.solve([operator], np.zeros(1), None, 1, {})

    assert operator.commits == 0
    assert application_solver.last_diagnostics.snes_converged_reason == reason


@pytest.mark.parametrize("reason", [-4, 0])
def test_petsc_snes_rejects_unconverged_linear_solve_before_commit(
    monkeypatch, reason
):
    operator = _ScalarOperator()
    application_solver = _fake_petsc_solver(
        monkeypatch, [1.0], snes_reason=2, ksp_reason=reason,
        linear_iterations=1
    )

    with pytest.raises(solver.SnesSolveError, match="linear solve did not report"):
        application_solver.solve([operator], np.zeros(1), None, 1, {})

    assert operator.commits == 0


def test_petsc_snes_allows_zero_ksp_reason_at_iteration_zero(monkeypatch):
    operator = _ScalarOperator()
    application_solver = _fake_petsc_solver(
        monkeypatch,
        [1.0],
        snes_reason=2,
        ksp_reason=0,
        nonlinear_iterations=0,
        linear_iterations=0,
    )

    _displacement, _committed, diagnostics = application_solver.solve(
        [operator], np.ones(1), None, 1, {}
    )

    assert operator.commits == 1
    assert diagnostics.ksp_converged_reason == 0


@pytest.mark.parametrize("result", [[float("nan")], [float("inf")], [1.0, 2.0]])
def test_petsc_snes_rejects_invalid_iterate_before_commit(
    monkeypatch, result
):
    operator = _ScalarOperator()
    application_solver = _fake_petsc_solver(monkeypatch, result)

    with pytest.raises(solver.SnesSolveError, match="malformed displacement"):
        application_solver.solve([operator], np.zeros(1), None, 1, {})

    assert operator.commits == 0


def test_petsc_snes_rejects_nonfinite_initial_residual_before_petsc(monkeypatch):
    operator = _NonfiniteAtRootOperator()
    application_solver = solver.PetscSnesSolver()
    monkeypatch.setattr(
        application_solver,
        "_ensure_context",
        lambda *_args, **_kwargs: pytest.fail("PETSc must not be initialized"),
    )

    with pytest.raises(solver.SnesSolveError, match="non-finite initial residual"):
        application_solver.solve([operator], np.ones(1), None, 1, {})

    assert operator.commits == 0


def test_petsc_snes_rejects_false_step_tolerance_convergence(monkeypatch):
    operator = _ScalarOperator(stubborn=True)
    application_solver = _fake_petsc_solver(monkeypatch, [2.0], snes_reason=3)

    with pytest.raises(solver.SnesSolveError, match="recovered residual rule"):
        application_solver.solve([operator], np.zeros(1), None, 1, {})

    assert operator.commits == 0
    diagnostics = application_solver.last_diagnostics
    assert diagnostics.final_residual_norm > diagnostics.residual_acceptance_threshold


def test_petsc_snes_rejects_nonfinite_final_residual(monkeypatch):
    operator = _NonfiniteAtRootOperator()
    application_solver = _fake_petsc_solver(monkeypatch, [1.0])

    with pytest.raises(solver.SnesSolveError, match="recovered residual rule"):
        application_solver.solve([operator], np.zeros(1), None, 1, {})

    assert operator.commits == 0
    assert np.isinf(application_solver.last_diagnostics.final_residual_norm)


def test_petsc_snes_rejects_nonfinite_petsc_diagnostics(monkeypatch):
    operator = _ScalarOperator()
    application_solver = _fake_petsc_solver(
        monkeypatch, [1.0], function_norm=float("nan")
    )

    with pytest.raises(solver.SnesSolveError, match="non-finite solver diagnostics"):
        application_solver.solve([operator], np.zeros(1), None, 1, {})

    assert operator.commits == 0


@pytest.mark.parametrize(("time_value", "dt_value"), [(float("nan"), 1.0), (1.0, 0.0)])
def test_petsc_snes_rejects_invalid_time_before_petsc(
    monkeypatch, time_value, dt_value
):
    application_solver = solver.PetscSnesSolver()
    monkeypatch.setattr(
        application_solver,
        "_ensure_context",
        lambda *_args, **_kwargs: pytest.fail("PETSc must not be initialized"),
    )
    with pytest.raises(ValueError, match="time"):
        application_solver.solve(
            [_ScalarOperator()],
            np.zeros(1),
            None,
            1,
            {},
            t=time_value,
            dt=dt_value,
        )


def test_petsc_snes_refreshes_trial_then_commits_once(monkeypatch):
    operator = _TrialTrackingScalarOperator()
    application_solver = _fake_petsc_solver(monkeypatch, [1.0])

    displacement, committed, diagnostics = application_solver.solve(
        [operator], np.zeros(1), None, 1, {}, t=0.2, dt=0.01
    )

    np.testing.assert_allclose(displacement, [1.0])
    assert committed == [None]
    assert operator.commits == 1
    assert diagnostics.time == pytest.approx(0.2)
    assert diagnostics.dt == pytest.approx(0.01)
    assert diagnostics.initial_residual_norm == pytest.approx(1.0)
    assert diagnostics.final_residual_norm == pytest.approx(0.0)
    assert diagnostics.residual_acceptance_threshold == pytest.approx(1.0e-9)
    assert diagnostics.snes_converged_reason == 2
    assert diagnostics.ksp_converged_reason == 4
    assert diagnostics.function_domain_rejections == 0
    assert diagnostics.last_function_domain_error is None
    assert json.loads(json.dumps(diagnostics.as_dict(), allow_nan=False))[
        "last_function_domain_error"
    ] is None


def test_petsc_function_callback_turns_only_deformation_domain_error_into_infinity():
    class DomainFailureOperator(_ScalarOperator):
        def residual(self, U, state, t, dt):
            raise local_pressure.InvalidDeformationError("invalid trial det(F)")

    application_solver = solver.PetscSnesSolver()
    fake_snes = _FakeSnes([0.0])
    application_solver._active = {
        "operators": (DomainFailureOperator(),),
        "state": None,
        "t": 0.1,
        "dt": 0.1,
        "ndof": 1,
    }
    target = _FakeVec([123.0])
    application_solver._form_function(fake_snes, _FakeVec([0.0]), target)

    assert np.isposinf(target.values[0])
    assert application_solver._function_domain_rejections == 1
    assert application_solver._last_function_domain_error == "invalid trial det(F)"
    assert application_solver._function_domain_rejection_api == (
        "nonfinite residual for PETSc BT"
    )


@pytest.mark.parametrize("jacobian", [False, True])
def test_petsc_callbacks_do_not_mask_unrelated_operator_errors(jacobian):
    class UnexpectedFailureOperator(_ScalarOperator):
        def residual(self, U, state, t, dt):
            if jacobian:
                return super().residual(U, state, t, dt)
            raise ValueError("unrelated residual failure")

        def tangent(self, U, state, t, dt):
            if not jacobian:
                return super().tangent(U, state, t, dt)
            raise ValueError("unrelated tangent failure")

    application_solver = solver.PetscSnesSolver()
    fake_snes = _FakeSnes([0.0])
    application_solver._active = {
        "operators": (UnexpectedFailureOperator(),),
        "state": None,
        "t": 0.1,
        "dt": 0.1,
        "ndof": 1,
    }
    with pytest.raises(ValueError, match="unrelated"):
        if jacobian:
            application_solver._form_jacobian(
                fake_snes, _FakeVec([0.0]), object(), object()
            )
        else:
            application_solver._form_function(
                fake_snes, _FakeVec([0.0]), _FakeVec([0.0])
            )
    assert application_solver._function_domain_rejections == 0


def test_petsc_snes_rejects_dirichlet_data_before_petsc(monkeypatch):
    application_solver = solver.PetscSnesSolver()
    monkeypatch.setattr(
        application_solver,
        "_ensure_context",
        lambda *_args, **_kwargs: pytest.fail("PETSc must not be initialized"),
    )

    with pytest.raises(NotImplementedError, match="empty Dirichlet"):
        application_solver.solve(
            [_ScalarOperator()], np.zeros(1), None, 1, {0: 0.0}
        )


def test_recovered_petsc_snes_settings_are_explicit_and_stable():
    configuration = solver.PetscSnesSolver().configuration()

    assert configuration["settings_source"] == (
        "recovered 2026-06-27 Case B development adapter"
    )
    assert configuration["snes_type"] == "newtonls"
    assert configuration["line_search_type"] == "bt"
    assert configuration["ksp_type"] == "preonly"
    assert configuration["pc_type"] == "lu"
    assert configuration["function_domain_rejection_api"] == (
        "nonfinite residual for PETSc BT"
    )
    assert configuration["rtol"] == pytest.approx(1.0e-9)
    assert configuration["atol"] == pytest.approx(1.0e-10)
    assert configuration["stol"] == pytest.approx(1.0e-12)
    assert configuration["max_it"] == 60


def test_distributed_result_check_rejects_false_success_metadata():
    with pytest.raises(RuntimeError, match="final residual tolerance"):
        solver.require_distributed_success(
            np.zeros(2),
            {"rnorm": 1.0e-3, "ksp_diverged": False},
            tol=1.0e-9,
            context="test solve",
        )
    with pytest.raises(RuntimeError, match="diverged PETSc"):
        solver.require_distributed_success(
            np.zeros(2),
            {"rnorm": 0.0, "ksp_diverged": True},
            tol=1.0e-9,
            context="test solve",
        )
    with pytest.raises(ValueError, match="tol must be finite and positive"):
        solver.require_distributed_success(
            np.zeros(2),
            {"rnorm": 0.0, "ksp_diverged": False},
            tol=0.0,
            context="test solve",
        )


def test_time_grid_requires_the_requested_end_point():
    np.testing.assert_allclose(cardiac_run._time_grid(0.2, 1.0), np.linspace(0, 1, 6))
    with pytest.raises(ValueError, match="integer multiple"):
        cardiac_run._time_grid(0.3, 1.0)


def test_serial_cli_keeps_joint_default_and_exposes_explicit_split():
    parser = cardiac_run._parser()
    assert parser.parse_args([]).element_evaluation == "joint"
    assert parser.parse_args(
        ["--element-evaluation", "split"]
    ).element_evaluation == "split"
    with pytest.raises(SystemExit):
        parser.parse_args(["--element-evaluation", "auto"])


def test_local_pressure_rejects_collapsed_apex_before_build(capsys):
    with pytest.raises(SystemExit) as error:
        cardiac_run.main(
            ["--formulation", "local-pressure", "--apex-offset", "0",
             "--mesh-topology", "polar-ring"]
        )
    assert error.value.code == 2
    assert "requires a nondegenerate open apex" in capsys.readouterr().err


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        (["--material-eta", "50"], "must equal the paper value 100 Pa s"),
        (["--isotropic"], "historical forensic control"),
        (
            ["--viscous-evidence-out", "eta-split.npz"],
            "eta-split diagnostic is a historical forensic control",
        ),
        (
            ["--viscous-evidence-start", "0.20"],
            "parameter sweeps are not run",
        ),
    ],
)
def test_benchmark_driver_rejects_physical_parameter_sweeps(
    arguments, message, capsys
):
    with pytest.raises(SystemExit) as error:
        cardiac_run.main(arguments)
    assert error.value.code == 2
    assert message in capsys.readouterr().err


def test_reference_directory_resolution_has_no_machine_local_fallback(tmp_path, monkeypatch):
    data = tmp_path / "benchmark_article_data" / "results_time_curves" / "data"
    data.mkdir(parents=True)
    (data / "x_nonblinded_step_0A_group_demo.pickle").touch()
    assert post.resolve_reference_dir(tmp_path) == data

    monkeypatch.delenv("CARDIAC_BENCHMARK_DATA_DIR", raising=False)
    try:
        post.resolve_reference_dir()
    except FileNotFoundError as error:
        assert "CARDIAC_BENCHMARK_DATA_DIR" in str(error)
    else:
        raise AssertionError("missing reference configuration must fail closed")


def test_reference_comparison_reads_completed_u0_u1_histories(tmp_path, capsys):
    data = tmp_path / "results_time_curves" / "data"
    data.mkdir(parents=True)
    times = np.array([0.0, 0.5, 1.0])
    displacement = {
        point: {
            "ux": np.array([0.0, 1.0, 2.0]),
            "uy": np.array([0.0, 0.0, 0.0]),
            "uz": np.array([0.0, -1.0, -2.0]),
        }
        for point in ("p0", "p1")
    }
    payload = pickle.dumps({"time": times, "displacement": displacement})
    for suffix in post.REFERENCE_MANIFEST_SUFFIXES:
        (data / post._reference_filename("step_0A", suffix)).write_bytes(payload)
    (data / post._reference_filename("step_0A", "simvascular")).write_bytes(
        payload
    )

    result = tmp_path / "result.npz"
    history = np.stack(
        [displacement["p0"][name] for name in ("ux", "uy", "uz")], axis=1
    )
    _write_qualified_result(result, times=times, u0=history, u1=history)
    post.main([str(result), "--case", "step_0A", "--reference-dir", str(tmp_path)])
    assert "ours=0.000" in capsys.readouterr().out


def test_reference_comparison_rejects_case_mismatch_before_loading_pickle(
    tmp_path, monkeypatch
):
    result = tmp_path / "result.npz"
    _write_qualified_result(result, case="B")
    monkeypatch.setattr(
        post,
        "load_reference",
        lambda *_args, **_kwargs: pytest.fail("reference must not load first"),
    )

    with pytest.raises(RuntimeError, match="wrong benchmark"):
        post.main([str(result), "--case", "step_0A", "--reference-dir", str(tmp_path)])


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"times": np.array([0.1, 0.5, 1.0])}, "time grid is inconsistent"),
        ({"times": np.array([0.0, 0.4, 1.0])}, "time grid is inconsistent"),
        ({"t_end": 1.5}, "time grid is inconsistent"),
        ({"n_mu": 2.5}, "invalid mesh field"),
    ],
)
def test_reference_comparison_rejects_inconsistent_result_coordinates(
    tmp_path, monkeypatch, overrides, message
):
    result = tmp_path / "result.npz"
    _write_qualified_result(result, **overrides)
    monkeypatch.setattr(
        post,
        "load_reference",
        lambda *_args, **_kwargs: pytest.fail("reference must not load first"),
    )

    with pytest.raises(RuntimeError, match=message):
        post.main([str(result), "--reference-dir", str(tmp_path)])


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("app_tree_state", "dirty", "requires a clean source tree"),
        ("core_revision", "unknown", "full 40-hex Git revision"),
    ],
)
def test_reference_comparison_rejects_unqualified_source_state(
    tmp_path, monkeypatch, field, value, message
):
    result = tmp_path / "result.npz"
    _write_qualified_result(result, **{field: value})
    monkeypatch.setattr(
        post,
        "load_reference",
        lambda *_args, **_kwargs: pytest.fail("reference must not load first"),
    )

    with pytest.raises(RuntimeError, match=message):
        post.main([str(result), "--reference-dir", str(tmp_path)])


def test_incomplete_solve_cannot_create_a_result_archive(tmp_path):
    output = tmp_path / "partial.npz"
    try:
        result_io.save_completed(
            output,
            completed_steps=2,
            expected_steps=3,
            times=np.arange(4),
        )
    except RuntimeError as error:
        assert "completed 2/3" in str(error)
    else:
        raise AssertionError("partial solve must not be writable as a completed result")
    assert not output.exists()


def test_driver_solver_failure_cannot_reach_result_writer(
    tmp_path, monkeypatch
):
    class Element:
        props = np.zeros(1)

    class FailingSnesSolver:
        def solve(self, *_args, **_kwargs):
            raise solver.SnesSolveError("forced divergence")

        def close(self):
            pass

        def configuration(self):
            return {"name": "petsc-snes"}

    monkeypatch.setattr(
        cardiac_run,
        "build_group",
        lambda *_args, **_kwargs: (object(), Element(), 0, None),
    )
    monkeypatch.setattr(cardiac_run, "RobinOperator", lambda *_args, **_kwargs: object())
    monkeypatch.setattr(
        cardiac_run,
        "audit_robin",
        lambda *_args, **_kwargs: {"passed": True},
    )
    monkeypatch.setattr(cardiac_run, "PetscSnesSolver", FailingSnesSolver)
    monkeypatch.setattr(
        cardiac_run,
        "locate_hex8_point",
        lambda *_args, **_kwargs: type(
            "Location",
            (),
            {
                "element_index": 0,
                "reconstruction_error": 0.0,
                "natural_coordinates": (0.0, 0.0, 0.0),
                "weights": (0.125,) * 8,
            },
        )(),
    )
    monkeypatch.setattr(
        cardiac_run,
        "save_completed",
        lambda *_args, **_kwargs: pytest.fail("failed solve must not write output"),
    )
    output = tmp_path / "diverged.npz"

    with pytest.raises(solver.SnesSolveError, match="forced divergence"):
        cardiac_run.main(
            [
                "--case", "A",
                "--nt", "1",
                "--nmu", "2",
                "--ntheta", "4",
                "--mesh-topology", "polar-ring",
                "--dt", "0.5",
                "--tend", "0.5",
                "--apex-offset", "0.2",
                "--nonlinear-solver", "petsc-snes",
                "--build-dir", str(tmp_path / "build"),
                "--out", str(output),
            ]
        )

    assert not output.exists()


def test_result_write_failure_cannot_replace_or_leave_a_partial_archive(
    tmp_path, monkeypatch
):
    output = tmp_path / "result.npz"
    output.write_bytes(b"previous complete result")

    def fail_after_partial_write(path, **_payload):
        Path(path).write_bytes(b"partial")
        raise OSError("simulated write failure")

    monkeypatch.setattr(result_io.np, "savez", fail_after_partial_write)
    with pytest.raises(OSError, match="simulated write failure"):
        result_io.save_completed(
            output,
            completed_steps=2,
            expected_steps=2,
            times=np.arange(3),
        )

    assert output.read_bytes() == b"previous complete result"
    assert list(tmp_path.glob(".result.npz.*.tmp.npz")) == []


def test_runtime_metadata_records_revisions_without_source_paths(monkeypatch):
    app_revision = "a" * 40
    core_revision = "b" * 40
    monkeypatch.setenv("COUPFE_CARDIAC_APP_REVISION", app_revision)
    monkeypatch.setenv("COUPFE_CARDIAC_TREE_STATE", "clean")
    monkeypatch.setenv("COUPFE_CARDIAC_SOURCE_KIND", "asserted")
    monkeypatch.setenv("COUPFE_CORE_REVISION", core_revision)
    monkeypatch.setenv("COUPFE_CORE_TREE_STATE", "clean")
    monkeypatch.setenv("COUPFE_CORE_SOURCE_KIND", "asserted")

    metadata = cardiac_run._runtime_metadata()

    assert metadata["result_schema"] == "coupfe-cardiac-result-v1"
    assert metadata["driver"] == "examples/cardiac_benchmark/run.py"
    assert metadata["app_revision"] == app_revision
    assert metadata["app_tree_state"] == "clean"
    assert metadata["app_source_kind"] == "asserted"
    assert metadata["core_revision"] == core_revision
    assert metadata["core_tree_state"] == "clean"
    assert metadata["core_source_kind"] == "asserted"
    assert metadata["core_source_url"] == "https://github.com/tengzhang48/CoupFE.git"
    assert metadata["python_version"]
    assert metadata["numpy_version"]
    assert metadata["scipy_version"]
    assert "source_path" not in metadata
    assert "command_argv" not in metadata


def test_pep610_core_identity_qualifies_a_standard_vcs_install(tmp_path, monkeypatch):
    installed_core = tmp_path / "site-packages" / "coupfe" / "__init__.py"
    installed_core.parent.mkdir(parents=True)
    installed_core.write_text("", encoding="utf-8")

    class Distribution:
        def read_text(self, name):
            assert name == "direct_url.json"
            return (
                '{"url":"https://github.com/tengzhang48/CoupFE.git",'
                '"vcs_info":{"vcs":"git","commit_id":"' + "4" * 40 + '"}}'
            )

        def locate_file(self, name):
            assert name == "coupfe/__init__.py"
            return installed_core

    monkeypatch.setattr(cardiac_run.importlib.metadata, "distribution", lambda _name: Distribution())
    identity = cardiac_run._installed_vcs_identity(module_file=installed_core)
    assert identity == {
        "revision": "4" * 40,
        "tree_state": "installed",
        "source_kind": "pep610-vcs",
        "source_url": "https://github.com/tengzhang48/CoupFE.git",
    }


def test_pep610_identity_rejects_a_shadowing_module(tmp_path, monkeypatch):
    distributed = tmp_path / "site-packages" / "coupfe" / "__init__.py"
    shadow = tmp_path / "checkout" / "coupfe" / "__init__.py"
    distributed.parent.mkdir(parents=True)
    shadow.parent.mkdir(parents=True)
    distributed.write_text("", encoding="utf-8")
    shadow.write_text("", encoding="utf-8")

    class Distribution:
        def read_text(self, _name):
            return (
                '{"url":"https://github.com/tengzhang48/CoupFE.git",'
                '"vcs_info":{"vcs":"git","commit_id":"' + "4" * 40 + '"}}'
            )

        def locate_file(self, _name):
            return distributed

    monkeypatch.setattr(cardiac_run.importlib.metadata, "distribution", lambda _name: Distribution())
    assert cardiac_run._installed_vcs_identity(module_file=shadow) is None


def test_core_checkout_lookup_does_not_mistake_nested_app_venv_for_core(tmp_path):
    app = tmp_path / "cardiac-app"
    (app / ".git").mkdir(parents=True)
    installed_core = app / ".venv" / "lib" / "python" / "site-packages" / "coupfe"
    installed_core.mkdir(parents=True)
    (installed_core / "__init__.py").write_text("", encoding="utf-8")

    assert cardiac_run._checkout_root(installed_core) == app
    assert cardiac_run._checkout_root(installed_core, require_core_package=True) is None


def test_core_checkout_lookup_does_not_bind_a_nested_core_venv(tmp_path):
    checkout = tmp_path / "core-checkout"
    (checkout / ".git").mkdir(parents=True)
    (checkout / "coupfe").mkdir()
    (checkout / "coupfe" / "__init__.py").write_text("", encoding="utf-8")
    installed_core = (
        checkout / ".venv" / "lib" / "python" / "site-packages" / "coupfe"
    )
    installed_core.mkdir(parents=True)
    (installed_core / "__init__.py").write_text("", encoding="utf-8")

    assert cardiac_run._checkout_root(
        installed_core, require_core_package=True
    ) is None
    assert cardiac_run._checkout_root(
        checkout / "coupfe", require_core_package=True
    ) == checkout
