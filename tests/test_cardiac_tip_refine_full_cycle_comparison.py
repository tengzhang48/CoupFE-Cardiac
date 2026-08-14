from __future__ import annotations

import hashlib
import json
import pickle
from pathlib import Path

import numpy as np
import pytest

import compare_tip_refine_full_cycle as comparison


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _histories(times: np.ndarray, shift_m: float = 0.0):
    p0 = np.column_stack((0.010 * times, 0.002 * times, -0.003 * times))
    p1 = np.column_stack((-0.005 * times, 0.001 * times, -0.002 * times))
    if shift_m:
        p0 = p0 + np.column_stack(
            (shift_m * times, -0.5 * shift_m * times, 0.25 * shift_m * times)
        )
        p1 = p1 + np.column_stack(
            (-0.25 * shift_m * times, shift_m * times, 0.5 * shift_m * times)
        )
    return p0, p1


def _write_run(
    path: Path,
    *,
    n_t: int,
    t_end: float,
    tip_refine: float = 6.0,
    landmark_p0: tuple[float, float, float] = (0.025, 0.030, 0.0),
    history_shift_m: float = 0.0,
    runtime_source_sha256: str = "c" * 64,
    semantic_overrides: dict | None = None,
) -> None:
    completed_steps = int(round(t_end / 1.0e-3))
    times = np.arange(completed_steps + 1, dtype=float) * 1.0e-3
    p0, p1 = _histories(times, history_shift_m)
    semantic = {
        "benchmark_activation_parameters_json": '{"Ta":0.0}',
        "benchmark_identity_scope": "benchmark-step-and-configuration-id",
        "benchmark_load_contract": "pressure-only",
        "benchmark_material_parameters_json": '{"eta":100.0,"kappa":1000000.0}',
        "benchmark_peak_load_definition": "argmax(abs(active_tension_pa+pressure_pa))",
        "benchmark_pressure_parameters_json": '{"sigma_mid":16000.0}',
        "benchmark_reproduction_profile": "not-applicable",
        "fiber_direction_reconstruction": "toolkit-physical-coordinate-u-v-v1",
        "generalized_alpha_stage_contract": "simula-source-matched-v1",
        "local_pressure_volume_law": "linear-reference-volume-mean-log-j-v1",
        "material_kernel_formulation": "standard",
        "material_model_id": (
            "holzapfel-ogden-smooth-switch-complete-energy-derivative-v1"
        ),
        "parameter_variant": "benchmark_eta",
        "point_sampling": "hex8_reference_isoparametric",
        "tbar_definition": "laplace_presolved",
        "viscous_rate": "velocity_consistent_green_lagrange_at_alpha_f_stage",
    }
    semantic.update(semantic_overrides or {})
    np.savez(
        path,
        result_schema="coupfe-cardiac-result-v1",
        completed_steps=completed_steps,
        expected_steps=completed_steps,
        converged=True,
        case="B",
        benchmark_step=0,
        benchmark_configuration_id="benchmark-1-step-0-case-B-pressure-only",
        benchmark_active_stress_enabled=False,
        benchmark_pressure_enabled=True,
        mpi_enabled=True,
        mpi_ranks=8,
        mpi_linear_solver_profile="fgmres-gamg-rigid-rebuild",
        integrator="generalized-alpha",
        generalized_alpha_alpha_m=0.2,
        generalized_alpha_alpha_f=0.4,
        generalized_alpha_gamma=0.7,
        generalized_alpha_beta=0.36,
        formulation="hex8_local_pressure_p0_condensed_logj",
        material_kappa_pa=0.0,
        local_pressure_bulk_modulus_pa=1.0e6,
        mass_representation="consistent_q1_hex8",
        material_eta_pa_s=100.0,
        density=1000.0,
        mesh_topology="closed_multiblock_disk",
        a_top=1.0e5,
        b_top=5.0e3,
        a_epi=1.0e8,
        b_epi=5.0e3,
        core_half_width=0.36,
        apex_offset=0.0,
        perturb=0.0,
        flip_helix=True,
        isotropic=False,
        n_t=n_t,
        n_core=20,
        n_radial=17,
        tip_refine=tip_refine,
        dt=1.0e-3,
        t_end=t_end,
        load_horizon=1.0,
        fiber_sampling_option="gp-direct",
        fiber_sampling="gp_direct_rule",
        viscous_term_active=True,
        element_evaluation_mode="joint",
        app_tree_state="clean",
        core_tree_state="clean",
        app_revision=("a" if n_t == 2 else "d") * 40,
        core_revision="b" * 40,
        benchmark_runtime_source_sha256=runtime_source_sha256,
        p0=np.array(landmark_p0),
        p1=np.array([0.000, 0.030, 0.0]),
        nodes=np.array([[float(n_t), 0.0, 0.0]]),
        elems=np.full((1, 8), n_t, dtype=np.int64),
        times=times,
        u0=p0,
        u1=p1,
        **comparison.EXPECTED_TBAR_IDENTITY_BY_NT[n_t],
        **semantic,
    )


def _write_fenics(path: Path, monkeypatch) -> None:
    path.mkdir()
    times = comparison.EXPECTED_FENICS_TIMES
    p0, p1 = _histories(times)
    p0[:, 0] += 1.0e-3
    p1[:, 2] += 2.0e-3
    (path / "parameters.json").write_text("{}\n", encoding="utf-8")
    np.save(path / "time_stamps.npy", times)
    np.save(path / "componentwise_displacement_up0.npy", p0)
    np.save(path / "componentwise_displacement_up1.npy", p1)
    monkeypatch.setattr(
        comparison,
        "EXPECTED_FENICS_SHA256",
        {name: _sha256(path / name) for name in comparison.EXPECTED_FENICS_SHA256},
    )


def _write_teams(path: Path, monkeypatch) -> None:
    path.mkdir()
    times = np.linspace(0.0, 1.0, 101)
    p0, p1 = _histories(times)
    hashes = {}
    team_names = tuple(comparison.TEAM_FILE_SHA256)
    for index, team in enumerate(team_names):
        offset = -1.0e-3 + 2.0e-3 * index / (len(team_names) - 1)
        record = {
            "time": times,
            "displacement": {
                "p0": {
                    component: p0[:, ci] + offset
                    for ci, component in enumerate(("ux", "uy", "uz"))
                },
                "p1": {
                    component: p1[:, ci] + offset
                    for ci, component in enumerate(("ux", "uy", "uz"))
                },
            },
        }
        name = f"{comparison.TEAM_FILENAME_PREFIX}{team}.pickle"
        target = path / name
        target.write_bytes(pickle.dumps(record, protocol=4))
        hashes[team] = _sha256(target)
    alias = path / f"{comparison.TEAM_FILENAME_PREFIX}{comparison.TEAM_ALIAS}.pickle"
    selected = path / (
        f"{comparison.TEAM_FILENAME_PREFIX}{comparison.TEAM_ALIAS_TARGET}.pickle"
    )
    alias.write_bytes(selected.read_bytes())
    monkeypatch.setattr(comparison, "TEAM_FILE_SHA256", hashes)


def _inputs(tmp_path: Path, monkeypatch):
    two = tmp_path / "two_layer.npz"
    four = tmp_path / "four_layer.npz"
    prefix = tmp_path / "four_layer_prefix.npz"
    fenics = tmp_path / "fenics"
    teams = tmp_path / "teams"
    _write_run(two, n_t=2, t_end=1.0)
    _write_run(four, n_t=4, t_end=1.0)
    _write_run(prefix, n_t=4, t_end=0.32)
    monkeypatch.setattr(comparison, "EXPECTED_COUPFE_SHA256", _sha256(two))
    monkeypatch.setattr(
        comparison, "EXPECTED_COUPFE_FOUR_LAYER_SHA256", _sha256(four)
    )
    monkeypatch.setattr(
        comparison,
        "EXPECTED_COUPFE_FOUR_LAYER_PREFIX_SHA256",
        _sha256(prefix),
    )
    _write_fenics(fenics, monkeypatch)
    _write_teams(teams, monkeypatch)

    def validated_result(archive, result_path, requested_case=None):
        assert requested_case == "step_0B"
        times = np.asarray(archive["times"], dtype=float)
        n_t = int(np.asarray(archive["n_t"]).item())
        return {
            "times": times,
            "histories": {
                "u0": np.asarray(archive["u0"], dtype=float),
                "u1": np.asarray(archive["u1"], dtype=float),
            },
            "det_f_gauss_peak": np.ones((1, 8), dtype=float),
            "solver_diagnostics": [
                {
                    "function_domain_rejections": 0,
                    "nonlinear_iterations": 2,
                }
                for _ in range(len(times) - 1)
            ],
            "n_peak": len(times) - 1,
            "mesh": {
                "nodes": np.array([[float(n_t), 0.0, 0.0]]),
                "elements": np.array([[n_t]], dtype=np.int64),
            },
            "sampling_metadata": {"p0": {"element": n_t}, "p1": {"element": n_t}},
            "runtime_versions": {
                "python_version": "3.10.8",
                "numpy_version": "1.26.4",
                "scipy_version": "1.15.2",
                "coupfe_version": "0.0.1",
            },
            "solver_configuration": {
                "petsc_version": "3.18.4",
                "petsc4py_version": "3.18.4",
            },
        }

    monkeypatch.setattr(comparison.post, "_load_validated_result", validated_result)
    return two, four, prefix, fenics, teams


def _metric_args(inputs):
    two, four, prefix, fenics, teams = inputs
    return two, four, prefix, fenics, teams


def _cli_args(inputs, report: Path, figure: Path):
    two, four, prefix, fenics, teams = inputs
    return [
        "--coupfe-npz",
        str(two),
        "--coupfe-four-layer",
        str(four),
        "--coupfe-four-layer-prefix",
        str(prefix),
        "--fenics-dir",
        str(fenics),
        "--teams-dir",
        str(teams),
        "--out-metrics",
        str(report),
        "--out-figure",
        str(figure),
    ]


def test_metrics_are_v2_source_bound_and_use_exact_shared_times(tmp_path, monkeypatch):
    inputs = _inputs(tmp_path, monkeypatch)

    metrics = comparison.compute_metrics(*_metric_args(inputs))

    assert metrics["schema"].endswith("-v2")
    assert set(metrics["runs"]) == {"two_layer", "four_layer"}
    for role in metrics["runs"]:
        assert metrics["runs"][role]["fenics"]["p0"][
            "vector_rmse_mm"
        ] == pytest.approx(1.0)
        assert metrics["runs"][role]["fenics"]["p1"][
            "vector_rmse_mm"
        ] == pytest.approx(2.0)
        assert metrics["runs"][role]["fenics"]["p0"][
            "last_shared_time_s"
        ] == 0.999
        assert metrics["runs"][role]["ten_team_envelope"]["p0"][
            "all_components_inside_fraction"
        ] == pytest.approx(1.0)
    assert metrics["cross_mesh"]["windows"]["snap_window"]["time_window_s"] == [
        0.2,
        0.32,
    ]
    assert metrics["cross_mesh"]["windows"]["full_cycle"]["points"]["p0"][
        "vector_rmse_mm"
    ] == 0.0
    decay = metrics["cross_mesh"]["cross_mesh_decay_diagnostic"]["points"]["p0"]
    assert decay["transient_max_difference_mm"] == 0.0
    assert decay["ratios"]["endpoint_over_transient_max"] == 0.0
    assert "classification" not in json.dumps(
        metrics["cross_mesh"]["cross_mesh_decay_diagnostic"]
    )
    assert metrics["cross_mesh"]["four_layer_prefix_continuity"]["passed"]
    assert metrics["inputs"]["coupfe"]["four_layer"]["validation_evidence"][
        "function_domain_rejections"
    ] == 0
    assert metrics["inputs"]["coupfe"]["four_layer"]["tbar_identity"][
        "tbar_source_sha256"
    ] == comparison.EXPECTED_TBAR_IDENTITY_BY_NT[4]["tbar_source_sha256"]
    assert len(metrics["inputs"]["ten_team_dataset"]["files"]) == 10
    assert "parameters.json" in metrics["inputs"]["fenics_files"]
    assert str(tmp_path) not in json.dumps(metrics)


def test_renderer_is_deterministic_accessible_and_distinguishes_runs(
    tmp_path, monkeypatch
):
    pytest.importorskip("matplotlib")
    inputs = _inputs(tmp_path, monkeypatch)

    first = comparison.render_figure(*_metric_args(inputs))
    second = comparison.render_figure(*_metric_args(inputs))
    decoded = first.decode("utf-8")

    assert first == second
    assert all(line == line.rstrip() for line in decoded.splitlines())
    assert 'role="img"' in decoded
    assert 'aria-labelledby="figure-title figure-description"' in decoded
    assert comparison.FIGURE_SUBTITLE in decoded
    assert "CoupFE 2x20x17 tip_refine=6.0" in decoded
    assert "CoupFE 4x20x17 tip_refine=6.0" in decoded
    assert "local FEniCS" in decoded
    assert "#d97706" in decoded
    assert "#2563eb" in decoded
    assert str(tmp_path) not in decoded


def test_uz_minus_5mm_events_require_downward_then_later_upward_crossing():
    times = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
    displacement = np.zeros((len(times), 3))
    displacement[:, 2] = np.array([0.0, -0.004, -0.006, -0.004, 0.0])

    events = comparison._uz_minus_5mm_events(times, displacement)

    assert events == pytest.approx({"downward_s": 0.15, "upward_s": 0.25})

    starts_below = displacement.copy()
    starts_below[:, 2] = np.array([-0.006, -0.004, -0.003, -0.002, -0.001])
    assert comparison._uz_minus_5mm_events(times, starts_below) == {
        "downward_s": None,
        "upward_s": None,
    }


def test_expected_tbar_identities_are_literal_source_bindings():
    assert comparison.EXPECTED_TBAR_IDENTITY_BY_NT[2]["tbar_source_sha256"] == (
        "d848b6cafa0e74c6e8cf56ddd825e8b3d8c91fd8490237ed8b204539f7d3cbeb"
    )
    assert comparison.EXPECTED_TBAR_IDENTITY_BY_NT[4]["tbar_source_sha256"] == (
        "1578362593495b6fe48d6a2fd2e1332150121be4d6b361915d04f3980d78da8f"
    )


def test_comparison_rejects_wrong_run_configuration(tmp_path, monkeypatch):
    inputs = _inputs(tmp_path, monkeypatch)
    two, four, prefix, fenics, teams = inputs
    _write_run(two, n_t=2, t_end=1.0, tip_refine=2.5)
    monkeypatch.setattr(comparison, "EXPECTED_COUPFE_SHA256", _sha256(two))

    with pytest.raises(ValueError, match="tip_refine"):
        comparison.compute_metrics(two, four, prefix, fenics, teams)


@pytest.mark.parametrize(
    ("field", "changed"),
    [
        ("a_top", 2.0e5),
        ("flip_helix", False),
        ("isotropic", True),
        ("fiber_sampling", "changed_rule"),
    ],
)
def test_comparison_rejects_changed_boundary_or_fiber_contract(
    tmp_path, monkeypatch, field, changed
):
    two, four, prefix, fenics, teams = _inputs(tmp_path, monkeypatch)
    with np.load(four, allow_pickle=False) as archive:
        record = {key: np.asarray(archive[key]) for key in archive.files}
    record[field] = np.asarray(changed)
    np.savez(four, **record)
    monkeypatch.setattr(
        comparison, "EXPECTED_COUPFE_FOUR_LAYER_SHA256", _sha256(four)
    )

    with pytest.raises(ValueError, match=field):
        comparison.compute_metrics(two, four, prefix, fenics, teams)


def test_comparison_rejects_wrong_laplace_field_identity(tmp_path, monkeypatch):
    two, four, prefix, fenics, teams = _inputs(tmp_path, monkeypatch)
    with np.load(four, allow_pickle=False) as archive:
        record = {key: np.asarray(archive[key]) for key in archive.files}
    record["tbar_source_sha256"] = np.asarray("f" * 64)
    np.savez(four, **record)
    monkeypatch.setattr(
        comparison, "EXPECTED_COUPFE_FOUR_LAYER_SHA256", _sha256(four)
    )

    with pytest.raises(ValueError, match="tbar_source_sha256"):
        comparison.compute_metrics(two, four, prefix, fenics, teams)


def test_comparison_rejects_invalid_rejection_counts_even_if_sum_is_zero(
    tmp_path, monkeypatch
):
    inputs = _inputs(tmp_path, monkeypatch)

    def invalid_diagnostics(archive, result_path, requested_case=None):
        times = np.asarray(archive["times"], dtype=float)
        diagnostics = [
            {"function_domain_rejections": 0, "nonlinear_iterations": 2}
            for _ in range(len(times) - 1)
        ]
        diagnostics[0]["function_domain_rejections"] = -1
        diagnostics[1]["function_domain_rejections"] = 1
        return {
            "times": times,
            "histories": {
                "u0": np.asarray(archive["u0"], dtype=float),
                "u1": np.asarray(archive["u1"], dtype=float),
            },
            "det_f_gauss_peak": np.ones((1, 8), dtype=float),
            "solver_diagnostics": diagnostics,
        }

    monkeypatch.setattr(
        comparison.post, "_load_validated_result", invalid_diagnostics
    )
    with pytest.raises(ValueError, match="invalid domain-rejection"):
        comparison.compute_metrics(*_metric_args(inputs))


def test_comparison_rejects_swapped_mesh_roles_even_when_hashes_are_rebound(
    tmp_path, monkeypatch
):
    two, four, prefix, fenics, teams = _inputs(tmp_path, monkeypatch)
    monkeypatch.setattr(comparison, "EXPECTED_COUPFE_SHA256", _sha256(four))
    monkeypatch.setattr(
        comparison, "EXPECTED_COUPFE_FOUR_LAYER_SHA256", _sha256(two)
    )

    with pytest.raises(ValueError, match=r"two_layer.*n_t"):
        comparison.compute_metrics(four, two, prefix, fenics, teams)


def test_comparison_rejects_semantic_contract_mismatch(tmp_path, monkeypatch):
    two, four, prefix, fenics, teams = _inputs(tmp_path, monkeypatch)
    _write_run(
        four,
        n_t=4,
        t_end=1.0,
        semantic_overrides={"benchmark_material_parameters_json": '{"eta":0.0}'},
    )
    monkeypatch.setattr(
        comparison, "EXPECTED_COUPFE_FOUR_LAYER_SHA256", _sha256(four)
    )

    with pytest.raises(ValueError, match="semantic contract mismatch"):
        comparison.compute_metrics(two, four, prefix, fenics, teams)


def test_comparison_rejects_runtime_source_identity_mismatch(tmp_path, monkeypatch):
    two, four, prefix, fenics, teams = _inputs(tmp_path, monkeypatch)
    _write_run(four, n_t=4, t_end=1.0, runtime_source_sha256="e" * 64)
    monkeypatch.setattr(
        comparison, "EXPECTED_COUPFE_FOUR_LAYER_SHA256", _sha256(four)
    )

    with pytest.raises(ValueError, match="runtime source identity mismatch"):
        comparison.compute_metrics(two, four, prefix, fenics, teams)


def test_comparison_rejects_four_layer_prefix_continuity_failure(
    tmp_path, monkeypatch
):
    two, four, prefix, fenics, teams = _inputs(tmp_path, monkeypatch)
    with np.load(four, allow_pickle=False) as archive:
        record = {key: np.asarray(archive[key]) for key in archive.files}
    changed = np.asarray(record["u0"]).copy()
    changed[100, 1] += 2.0e-12  # 2e-9 mm, above the declared 1e-9 mm tolerance.
    record["u0"] = changed
    np.savez(four, **record)
    monkeypatch.setattr(
        comparison, "EXPECTED_COUPFE_FOUR_LAYER_SHA256", _sha256(four)
    )

    with pytest.raises(ValueError, match="retained-prefix continuity"):
        comparison.compute_metrics(two, four, prefix, fenics, teams)


def test_comparison_rejects_four_layer_prefix_mesh_mismatch(
    tmp_path, monkeypatch
):
    two, four, prefix, fenics, teams = _inputs(tmp_path, monkeypatch)
    with np.load(prefix, allow_pickle=False) as archive:
        record = {key: np.asarray(archive[key]) for key in archive.files}
    changed_nodes = np.asarray(record["nodes"]).copy()
    changed_nodes[0, 0] += 1.0e-12
    record["nodes"] = changed_nodes
    np.savez(prefix, **record)
    monkeypatch.setattr(
        comparison,
        "EXPECTED_COUPFE_FOUR_LAYER_PREFIX_SHA256",
        _sha256(prefix),
    )

    with pytest.raises(ValueError, match="mesh or point sampling"):
        comparison.compute_metrics(two, four, prefix, fenics, teams)


def test_comparison_rejects_runtime_environment_mismatch(
    tmp_path, monkeypatch
):
    inputs = _inputs(tmp_path, monkeypatch)
    original = comparison._load_coupfe

    def changed_environment(*args, **kwargs):
        times, histories, identity = original(*args, **kwargs)
        if kwargs["role"] == "four_layer":
            identity["environment_identity"] = {
                **identity["environment_identity"],
                "petsc_version": "changed",
            }
        return times, histories, identity

    monkeypatch.setattr(comparison, "_load_coupfe", changed_environment)
    with pytest.raises(ValueError, match="runtime environment mismatch"):
        comparison.compute_metrics(*_metric_args(inputs))


def test_comparison_rejects_changed_reference_input(tmp_path, monkeypatch):
    inputs = _inputs(tmp_path, monkeypatch)
    two, four, prefix, fenics, teams = inputs
    values = np.load(fenics / "componentwise_displacement_up0.npy")
    values[0, 0] += 1.0e-6
    np.save(fenics / "componentwise_displacement_up0.npy", values)

    with pytest.raises(ValueError, match="SHA-256"):
        comparison.compute_metrics(two, four, prefix, fenics, teams)


def test_comparison_rejects_changed_retained_prefix(tmp_path, monkeypatch):
    two, four, prefix, fenics, teams = _inputs(tmp_path, monkeypatch)
    prefix.write_bytes(prefix.read_bytes() + b"changed")

    with pytest.raises(ValueError, match=r"four_layer_prefix.*SHA-256"):
        comparison.compute_metrics(two, four, prefix, fenics, teams)


def test_comparison_rejects_wrong_physical_landmark(tmp_path, monkeypatch):
    two, four, prefix, fenics, teams = _inputs(tmp_path, monkeypatch)
    _write_run(
        two,
        n_t=2,
        t_end=1.0,
        landmark_p0=(0.024, 0.030, 0.0),
    )
    monkeypatch.setattr(comparison, "EXPECTED_COUPFE_SHA256", _sha256(two))

    with pytest.raises(ValueError, match="physical landmark p0"):
        comparison.compute_metrics(two, four, prefix, fenics, teams)


def test_comparison_rejects_mutated_team_alias(tmp_path, monkeypatch):
    two, four, prefix, fenics, teams = _inputs(tmp_path, monkeypatch)
    alias = teams / (
        f"{comparison.TEAM_FILENAME_PREFIX}{comparison.TEAM_ALIAS}.pickle"
    )
    alias.write_bytes(alias.read_bytes() + b"changed")

    with pytest.raises(ValueError, match="alias is not byte-identical"):
        comparison.compute_metrics(two, four, prefix, fenics, teams)


def test_post_validator_failure_propagates(tmp_path, monkeypatch):
    inputs = _inputs(tmp_path, monkeypatch)

    def reject(*args, **kwargs):
        raise RuntimeError("archive audit failed")

    monkeypatch.setattr(comparison.post, "_load_validated_result", reject)
    with pytest.raises(RuntimeError, match="archive audit failed"):
        comparison.compute_metrics(*_metric_args(inputs))


@pytest.mark.parametrize("role_index", [0, 1, 2])
def test_post_validator_is_invoked_for_each_coupfe_role(
    tmp_path, monkeypatch, role_index
):
    inputs = _inputs(tmp_path, monkeypatch)
    target = inputs[role_index]
    delegated_validator = comparison.post._load_validated_result

    def reject_target(archive, result_path, requested_case=None):
        if Path(result_path) == target:
            raise RuntimeError(f"rejected role {role_index}")
        return delegated_validator(archive, result_path, requested_case)

    monkeypatch.setattr(comparison.post, "_load_validated_result", reject_target)
    with pytest.raises(RuntimeError, match=f"rejected role {role_index}"):
        comparison.compute_metrics(*_metric_args(inputs))


def test_render_failure_preserves_existing_output_pair(tmp_path, monkeypatch):
    inputs = _inputs(tmp_path, monkeypatch)
    report = tmp_path / "old.json"
    figure = tmp_path / "old.svg"
    report.write_bytes(b"old report\n")
    figure.write_bytes(b"old figure\n")

    def reject_render(*args, **kwargs):
        raise RuntimeError("render failed")

    monkeypatch.setattr(comparison, "render_figure", reject_render)
    with pytest.raises(RuntimeError, match="render failed"):
        comparison.main(_cli_args(inputs, report, figure))
    assert report.read_bytes() == b"old report\n"
    assert figure.read_bytes() == b"old figure\n"


def test_second_write_failure_restores_existing_output_pair(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    inputs = _inputs(tmp_path, monkeypatch)
    report = tmp_path / "old.json"
    figure = tmp_path / "old.svg"
    report.write_bytes(b"old report\n")
    figure.write_bytes(b"old figure\n")
    atomic_write = comparison._write_bytes_atomic
    calls = 0

    def fail_once_on_second_write(path, payload):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError("second artifact write failed")
        return atomic_write(path, payload)

    monkeypatch.setattr(
        comparison, "_write_bytes_atomic", fail_once_on_second_write
    )
    with pytest.raises(OSError, match="second artifact write failed"):
        comparison.main(_cli_args(inputs, report, figure))
    assert report.read_bytes() == b"old report\n"
    assert figure.read_bytes() == b"old figure\n"


def test_cli_rejects_output_collision_before_loading(tmp_path):
    output = tmp_path / "same.svg"
    inputs = (
        tmp_path / "two.npz",
        tmp_path / "four.npz",
        tmp_path / "prefix.npz",
        tmp_path / "fenics",
        tmp_path / "teams",
    )
    with pytest.raises(ValueError, match="output paths must be distinct"):
        comparison.main(_cli_args(inputs, output, output))


def test_cli_rejects_output_collision_with_any_coupfe_input(tmp_path):
    inputs = (
        tmp_path / "two.npz",
        tmp_path / "four.json",
        tmp_path / "prefix.npz",
        tmp_path / "fenics",
        tmp_path / "teams",
    )
    with pytest.raises(ValueError, match="must not overwrite any CoupFE"):
        comparison.main(_cli_args(inputs, inputs[1], tmp_path / "figure.svg"))


def test_cli_rejects_duplicate_coupfe_input_roles(tmp_path):
    duplicate = tmp_path / "same.npz"
    inputs = (
        duplicate,
        duplicate,
        tmp_path / "prefix.npz",
        tmp_path / "fenics",
        tmp_path / "teams",
    )
    with pytest.raises(ValueError, match="inputs must be distinct"):
        comparison.main(
            _cli_args(inputs, tmp_path / "report.json", tmp_path / "figure.svg")
        )


def test_cli_writes_public_readable_artifacts(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    inputs = _inputs(tmp_path, monkeypatch)
    report = tmp_path / "report.json"
    figure = tmp_path / "figure.svg"

    assert comparison.main(_cli_args(inputs, report, figure)) == 0
    assert report.stat().st_mode & 0o777 == 0o644
    assert figure.stat().st_mode & 0o777 == 0o644
    assert json.loads(report.read_text(encoding="utf-8"))["schema"].endswith("-v2")


def test_cli_expands_user_paths_for_validation_and_loading(tmp_path, monkeypatch):
    pytest.importorskip("matplotlib")
    _inputs(tmp_path, monkeypatch)
    monkeypatch.setenv("HOME", str(tmp_path))
    tilde_inputs = (
        Path("~/two_layer.npz"),
        Path("~/four_layer.npz"),
        Path("~/four_layer_prefix.npz"),
        Path("~/fenics"),
        Path("~/teams"),
    )

    assert comparison.main(
        _cli_args(
            tilde_inputs,
            Path("~/expanded.json"),
            Path("~/expanded.svg"),
        )
    ) == 0
    assert (tmp_path / "expanded.json").is_file()
    assert (tmp_path / "expanded.svg").is_file()


def test_cli_rejects_non_svg_figure_path_before_loading(tmp_path):
    inputs = (
        tmp_path / "two.npz",
        tmp_path / "four.npz",
        tmp_path / "prefix.npz",
        tmp_path / "fenics",
        tmp_path / "teams",
    )
    with pytest.raises(ValueError, match=r"\.svg"):
        comparison.main(
            _cli_args(inputs, tmp_path / "report.json", tmp_path / "figure.png")
        )
