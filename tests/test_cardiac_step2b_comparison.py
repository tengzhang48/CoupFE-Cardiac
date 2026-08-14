from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import pickle

import numpy as np
import pytest

import compare_step2b_case_b as comparison
import plot_step2b_case_b as plotter
import plot_step2b_current_rerun as current_plotter


ROOT = Path(__file__).resolve().parents[1]
OFFICIAL_MANIFEST = (
    ROOT
    / "examples"
    / "cardiac_benchmark"
    / "step2b_case_b_reference_hashes.json"
)
RETAINED_STEP2B_REPORT = (
    ROOT
    / "examples"
    / "cardiac_benchmark"
    / "results"
    / plotter.REPORT_NAME
)


def _base_ux(time_s, point):
    if point == "p0":
        knot_time = [0.00, 0.16, 0.20, 0.235, 0.31, 0.40, 0.48, 0.51, 0.58, 1.00]
        knot_value = [0.00, 0.00, -0.005, -0.015, -0.028, -0.031, -0.031, -0.025, -0.012, 0.00]
    else:
        knot_time = [0.00, 0.16, 0.204, 0.241, 0.30, 0.40, 0.48, 0.51, 0.58, 1.00]
        knot_value = [0.00, 0.00, -0.005, -0.015, -0.0225, -0.0245, -0.0245, -0.0195, -0.009, 0.00]
    return np.interp(time_s, knot_time, knot_value)


def _synthetic_displacement(time_s, team_index):
    scale = 0.982 + 0.004 * team_index
    result = {}
    for point_index, point in enumerate(comparison.POINTS):
        ux = scale * _base_ux(time_s, point)
        normalized = -ux / (0.031 if point == "p0" else 0.0245)
        uy = (0.0026 if point == "p0" else -0.0017) * normalized
        uz = (0.0050 if point == "p0" else 0.0019) * normalized
        if point_index == 1:
            uz *= 1.0 - 0.002 * team_index
        magnitude = np.sqrt(ux * ux + uy * uy + uz * uz)
        result[point] = {"ux": ux, "uy": uy, "uz": uz, "magnitude": magnitude}
    return result


def _write_publisher_data(root: Path):
    root.mkdir(parents=True)
    team_arrays = []
    for team_index, (team_id, _label, filename) in enumerate(
        comparison.PUBLISHER_SELECTION
    ):
        stored_time = comparison.PUBLISHED_TIME_S.copy()
        if team_id == "ambit":
            stored_time[0] = 0.001
        displacement = _synthetic_displacement(
            comparison.PUBLISHED_TIME_S, team_index
        )
        value = {
            "time": stored_time,
            "displacement": displacement,
            "stress": {"synthetic": np.zeros(101)},
            "volume": np.ones(101),
        }
        path = root / filename
        with path.open("wb") as stream:
            pickle.dump(value, stream, protocol=4)
        team_arrays.append(
            np.stack(
                [
                    np.stack(
                        [displacement[point][component] for component in comparison.COMPONENTS]
                    )
                    for point in comparison.POINTS
                ]
            )
        )
    return np.stack(team_arrays)


def _write_manifest(path: Path, data_directory: Path):
    files = []
    for team_id, label, filename in comparison.PUBLISHER_SELECTION:
        payload = (data_directory / filename).read_bytes()
        files.append(
            {
                "team_id": team_id,
                "label": label,
                "filename": filename,
                "sha256": hashlib.sha256(payload).hexdigest(),
                "size_bytes": len(payload),
            }
        )
    path.write_text(
        json.dumps(
            {
                "schema": comparison.HASH_MANIFEST_SCHEMA,
                "source_doi": comparison.OFFICIAL_SOURCE_DOI,
                "selection_source": comparison.OFFICIAL_SELECTION_SOURCE,
                "files": files,
            }
        ),
        encoding="utf-8",
    )


def _write_run(path: Path, common_mean, **overrides):
    native = np.empty((2, 3, 1001))
    for point_index in range(2):
        for component_index in range(3):
            native[point_index, component_index] = np.interp(
                comparison.COUPFE_TIME_S,
                comparison.PUBLISHED_TIME_S,
                common_mean[point_index, component_index],
            )
    tau, pressure = comparison._expected_step2b_load_histories()
    source_metadata = comparison.runtime_source_metadata()
    diagnostics = [
        {
            "time": float(comparison.COUPFE_TIME_S[index]),
            "ranks": 4,
            "final_residual_norm": 0.0,
            "residual_acceptance_threshold": 1.0,
            "snes_converged_reason": 2,
            "active_tension_pa": float(tau[index]),
            "pressure_pa": float(pressure[index]),
        }
        for index in range(1, 1001)
    ]
    payload = {
        "result_schema": comparison.RESULT_SCHEMA,
        "converged": True,
        "completed_steps": 1000,
        "expected_steps": 1000,
        "case": "B",
        "benchmark_step": 2,
        "benchmark_configuration_id": comparison.BENCHMARK_CONFIGURATION_ID,
        "benchmark_load_contract": comparison.BENCHMARK_LOAD_CONTRACT,
        "benchmark_peak_load_definition": comparison.BENCHMARK_PEAK_LOAD_DEFINITION,
        "benchmark_active_stress_enabled": True,
        "benchmark_pressure_enabled": True,
        "dt": 0.001,
        "t_end": 1.0,
        "integrator": "generalized-alpha",
        "generalized_alpha_alpha_m": 0.2,
        "generalized_alpha_alpha_f": 0.4,
        "generalized_alpha_gamma": 0.7,
        "generalized_alpha_beta": 0.36,
        "generalized_alpha_stage_contract": comparison.GENERALIZED_ALPHA_STAGE_CONTRACT,
        "load_evaluation_times_s": np.concatenate(
            ([0.0], comparison.COUPFE_TIME_S[1:] - 0.4 * 0.001)
        ),
        "benchmark_material_parameters_json": json.dumps(
            comparison.EXPECTED_MATERIAL
        ),
        "benchmark_activation_parameters_json": json.dumps(
            comparison.EXPECTED_ACTIVATION
        ),
        "benchmark_pressure_parameters_json": json.dumps(
            comparison.EXPECTED_PRESSURE
        ),
        **source_metadata,
        "driver": "examples/cardiac_benchmark/run_mpi.py",
        "benchmark_reproduction_profile": "paper-source-matched-full-cycle",
        "density": 1000.0,
        "material_kernel_formulation": "standard",
        "material_model_id": comparison.EXPECTED_MATERIAL_MODEL_ID,
        "mass_representation": "consistent_q1_hex8",
        "nonlinear_solver": "petsc-snes-mpi",
        "nonlinear_step_diagnostics_json": json.dumps(diagnostics),
        "mesh_topology": "closed_multiblock_disk",
        "fiber_sampling": "gp_direct_rule",
        "fiber_sampling_option": "gp-direct",
        "tbar_definition": "laplace_presolved",
        "tbar_source_sha256": "1" * 64,
        "tbar_metadata_sha256": "2" * 64,
        "tbar_metadata_schema": "coupfe-cardiac-laplace-tbar-v1",
        "point_sampling": "hex8_reference_isoparametric",
        "viscous_rate": "velocity_consistent_green_lagrange_at_alpha_f_stage",
        "material_eta_pa_s": 100.0,
        "parameter_variant": "benchmark_eta",
        "a_top": 1.0e5,
        "b_top": 5.0e3,
        "a_epi": 1.0e8,
        "b_epi": 5.0e3,
        "perturb": 0.0,
        "apex_offset": 0.0,
        "mpi_enabled": True,
        "viscous_term_active": True,
        "flip_helix": True,
        "isotropic": False,
        "mpi_ranks": 4,
        "mpi_world_size": 4,
        "formulation": "hex8_standard_pointwise_kappa",
        **comparison.FORMULATION_CONTRACTS["hex8_standard_pointwise_kappa"],
        "app_revision": "0" * 40,
        "app_source_kind": "git-checkout",
        "app_tree_state": "dirty",
        "core_revision": comparison.EXPECTED_CORE_REVISION,
        "core_source_url": comparison.EXPECTED_CORE_URL,
        "core_source_kind": "pep610-vcs",
        "core_tree_state": "installed",
        "pre_solve_audit_json": json.dumps(
            {
                name: {"schema": schema, "passed": True}
                for name, schema in {
                    "geometry": "coupfe-cardiac-pre-solve-geometry-v1",
                    "pressure": "coupfe-cardiac-pre-solve-pressure-v1",
                    "robin": "coupfe-cardiac-pre-solve-robin-v1",
                }.items()
            }
        ),
        "times": comparison.COUPFE_TIME_S,
        "u0": native[0].T,
        "u1": native[1].T,
        "tau": tau,
        "pres": pressure,
        "p0": comparison.LANDMARKS_M["p0"],
        "p1": comparison.LANDMARKS_M["p1"],
    }
    payload.update(overrides)
    np.savez(path, **payload)


def _write_current_rerun(path: Path, **overrides):
    payload = {
        **current_plotter.REVIEWED_CORRECTED_RUN_SCALARS,
        "times": comparison.COUPFE_TIME_S,
        "load_evaluation_times_s": np.concatenate(
            ([0.0], comparison.COUPFE_TIME_S[1:] - 0.4 * 0.001)
        ),
        "u0": np.zeros((1001, 3)),
        "u1": np.zeros((1001, 3)),
        "p0": comparison.LANDMARKS_M["p0"],
        "p1": comparison.LANDMARKS_M["p1"],
        "nodes": np.zeros((5403, 3)),
        "elems": np.zeros((3520, 8), dtype=np.int64),
    }
    payload.update(overrides)
    np.savez(path, **payload)


def _trust_current_rerun(monkeypatch, path: Path):
    payload = path.read_bytes()
    monkeypatch.setattr(
        current_plotter, "REVIEWED_CORRECTED_NPZ_NAME", path.name
    )
    monkeypatch.setattr(
        current_plotter,
        "REVIEWED_CORRECTED_NPZ_SHA256",
        hashlib.sha256(payload).hexdigest(),
    )
    monkeypatch.setattr(
        current_plotter, "REVIEWED_CORRECTED_NPZ_SIZE_BYTES", len(payload)
    )


def _trust_synthetic_manifest(monkeypatch, manifest):
    payload = manifest.read_bytes()
    monkeypatch.setattr(
        comparison,
        "OFFICIAL_HASH_MANIFEST_SHA256",
        hashlib.sha256(payload).hexdigest(),
    )
    monkeypatch.setattr(
        comparison, "OFFICIAL_HASH_MANIFEST_SIZE_BYTES", len(payload)
    )


def _synthetic_inputs(tmp_path, monkeypatch):
    data = tmp_path / "publisher"
    teams = _write_publisher_data(data)
    manifest = tmp_path / "hashes.json"
    _write_manifest(manifest, data)
    _trust_synthetic_manifest(monkeypatch, manifest)
    run = tmp_path / "step2b.npz"
    _write_run(run, np.mean(teams, axis=0))
    return data, manifest, run


def test_synthetic_end_to_end_separates_plumbing_and_reproduction(
    tmp_path, monkeypatch
):
    data, manifest, run = _synthetic_inputs(tmp_path, monkeypatch)

    report = comparison.build_report(data, manifest, run)
    json.dumps(report, allow_nan=False)

    assert report["plumbing"]["status"] == "passed"
    assert report["reproduction"]["status"] == "quantified-no-paper-acceptance-threshold"
    errors = report["reproduction"]["trajectory_errors"]
    assert list(errors) == [
        "contraction_event",
        "post_snap_plateau",
        "relaxation",
        "full_history",
    ]
    assert errors["post_snap_plateau"]["all_components"]["rmse_m"] < 1.0e-14
    assert errors["post_snap_plateau"]["points"]["p0"]["vector"]["rmse_vector_m"] < 1.0e-14
    paper_red = report["reproduction"]["paper_relative_discrepancy"]
    assert paper_red["metric"] == "relative discrepancy (benchmark paper Eq. 21)"
    assert paper_red["points"] == pytest.approx({"p0": 0.0, "p1": 0.0})
    assert (
        report["reproduction"]["post_snap_branch"]["points"]["p0"]
        ["shape_rmse_after_0p32_level_alignment_m"]
        < 1.0e-14
    )
    envelopes = (
        report["reproduction"]["post_snap_branch"]["points"]["p0"]
        ["component_envelopes_0p32_to_0p48"]
    )
    assert envelopes["ux"]["samplewise_official_envelope_coverage_fraction"] == 1.0
    assert envelopes["ux"]["plateau_level_sign"][
        "all_coupfe_samples_match_official_unanimous_sign"
    ] is True
    assert report["reproduction"]["post_snap_branch"]["points"]["p0"][
        "plateau_vector_octant"
    ]["matches_official_unanimous_octant"] is True
    events = report["reproduction"]["contraction_events"]
    assert events["fixed_threshold_crossings"]["p0"]["-5"][
        "coupfe_common_grid"
    ]["crossing_count"] == 1
    assert events["fixed_threshold_crossings"]["p0"]["-5"][
        "common_grid_delta_from_simula_s"
    ] is not None
    assert set(events["normalized_drop_crossings"]["p1"]) == {"10", "50", "90"}
    assert report["publisher_mean_std"]["sample_count"] == 101
    assert report["publisher_mean_std"]["team_count"] == 10
    attribution = report["publisher_material_attribution"]
    assert attribution["source_doi"] == comparison.OFFICIAL_SOURCE_DOI
    assert attribution["creators"] == list(comparison.OFFICIAL_SOURCE_CREATORS)
    assert attribution["license"] == "CC-BY-4.0"
    assert "modified/derived data" in attribution["transformation_notice"]
    curves = report["comparison_curves"]
    assert curves["sample_count"] == 101
    assert curves["team_count"] == 10
    assert set(curves["points"]) == {"p0", "p1"}
    assert set(curves["points"]["p0"]) == {"ux", "uy", "uz"}
    p0_x = curves["points"]["p0"]["ux"]
    assert set(p0_x) == {
        "coupfe_m",
        "publisher_mean_m",
        "publisher_min_m",
        "publisher_max_m",
        "publisher_simula_m",
    }
    assert len(p0_x["coupfe_m"]) == 101
    assert p0_x["coupfe_m"] == pytest.approx(p0_x["publisher_mean_m"])
    assert np.all(
        np.asarray(p0_x["publisher_min_m"])
        <= np.asarray(p0_x["publisher_mean_m"])
    )
    assert np.all(
        np.asarray(p0_x["publisher_mean_m"])
        <= np.asarray(p0_x["publisher_max_m"])
    )
    assert "trajectory_errors_vs_named_simula" in report["reproduction"]


def test_official_manifest_has_a_hard_coded_root_of_trust(tmp_path):
    manifest, identity = comparison.load_hash_manifest(OFFICIAL_MANIFEST)
    assert identity["sha256"] == comparison.OFFICIAL_HASH_MANIFEST_SHA256
    assert manifest["source_doi"] == comparison.OFFICIAL_SOURCE_DOI

    changed = tmp_path / "changed-manifest.json"
    changed.write_bytes(OFFICIAL_MANIFEST.read_bytes() + b"\n")
    with pytest.raises(
        comparison.ComparisonInputError, match="reviewed official manifest"
    ):
        comparison.load_hash_manifest(changed)


def test_manifest_and_pickles_are_identified_and_decoded_from_one_read(
    tmp_path, monkeypatch
):
    data, manifest, _run = _synthetic_inputs(tmp_path, monkeypatch)
    tracked = {
        manifest.resolve(),
        *(
            (data / filename).resolve()
            for _team_id, _label, filename in comparison.PUBLISHER_SELECTION
        ),
    }
    open_counts = {path: 0 for path in tracked}
    original_open = Path.open

    def counted_open(path, *args, **kwargs):
        resolved = path.resolve()
        if resolved in open_counts:
            open_counts[resolved] += 1
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", counted_open)
    reference = comparison.load_publisher_reference(data, manifest)

    assert all(count == 1 for count in open_counts.values())
    assert reference["team_times_s"].shape == (10, 101)
    assert reference["team_times_s"][1, 0] == pytest.approx(0.001)
    assert reference["team_times_s"][0, 0] == pytest.approx(0.0)


def test_current_rerun_rejects_nonreviewed_npz_even_with_matching_metadata(
    tmp_path,
):
    run = tmp_path / current_plotter.REVIEWED_CORRECTED_NPZ_NAME
    _write_current_rerun(run)

    with pytest.raises(
        comparison.ComparisonInputError, match="exact reviewed NPZ"
    ):
        current_plotter.load_reviewed_corrected_run(run)


def test_current_rerun_rejects_backward_euler_metadata_after_identity_check(
    tmp_path, monkeypatch
):
    run = tmp_path / current_plotter.REVIEWED_CORRECTED_NPZ_NAME
    _write_current_rerun(run, integrator="backward-euler")
    _trust_current_rerun(monkeypatch, run)

    with pytest.raises(
        comparison.ComparisonInputError, match="field 'integrator'"
    ):
        current_plotter.load_reviewed_corrected_run(run)


def test_current_rerun_accepts_hash_pinned_matching_contract(
    tmp_path, monkeypatch
):
    run = tmp_path / current_plotter.REVIEWED_CORRECTED_NPZ_NAME
    _write_current_rerun(run)
    _trust_current_rerun(monkeypatch, run)

    loaded = current_plotter.load_reviewed_corrected_run(run)

    assert loaded["times"].shape == (1001,)
    assert loaded["ours"]["p0"].shape == (1001, 3)
    assert loaded["identity"]["sha256"] == hashlib.sha256(
        run.read_bytes()
    ).hexdigest()


def test_current_rerun_legacy_report_rejects_schema_change_after_identity_check(
    tmp_path, monkeypatch
):
    source = ROOT / "examples" / "cardiac_benchmark" / "results" / (
        current_plotter.REVIEWED_LEGACY_REPORT_NAME
    )
    changed = json.loads(source.read_text(encoding="utf-8"))
    changed["comparison_curves"]["points"].pop("p1")
    path = tmp_path / current_plotter.REVIEWED_LEGACY_REPORT_NAME
    path.write_text(json.dumps(changed), encoding="utf-8")
    payload = path.read_bytes()
    monkeypatch.setattr(
        current_plotter,
        "REVIEWED_LEGACY_REPORT_SHA256",
        hashlib.sha256(payload).hexdigest(),
    )
    monkeypatch.setattr(
        current_plotter, "REVIEWED_LEGACY_REPORT_SIZE_BYTES", len(payload)
    )

    with pytest.raises(
        comparison.ComparisonInputError, match="exactly p0 and p1"
    ):
        current_plotter.load_reviewed_legacy_report(path)


def test_current_rerun_rejects_nonreviewed_legacy_report(tmp_path):
    source = ROOT / "examples" / "cardiac_benchmark" / "results" / (
        current_plotter.REVIEWED_LEGACY_REPORT_NAME
    )
    path = tmp_path / current_plotter.REVIEWED_LEGACY_REPORT_NAME
    path.write_bytes(source.read_bytes() + b"\n")

    with pytest.raises(
        comparison.ComparisonInputError, match="exact reviewed report"
    ):
        current_plotter.load_reviewed_legacy_report(path)


def test_current_rerun_cli_requires_canonical_report_basename(
    tmp_path, capsys
):
    with pytest.raises(SystemExit, match="2"):
        current_plotter.main(
            [
                "--corrected-npz",
                str(tmp_path / "missing.npz"),
                "--legacy-report",
                str(tmp_path / "missing.json"),
                "--teams-dir",
                str(tmp_path),
                "--out-metrics",
                str(tmp_path / "misleading-name.json"),
                "--out-figure",
                str(tmp_path / "figure.svg"),
            ]
        )
    assert (
        f"--out-metrics basename must be {current_plotter.REPORT_NAME}"
        in capsys.readouterr().err
    )


def test_retained_runtime_source_manifest_is_pinned_and_accepts_old_run(
    tmp_path,
):
    reviewed, identity = comparison.load_reviewed_runtime_source_manifest()
    assert (
        reviewed["benchmark_runtime_source_sha256"]
        == comparison.REVIEWED_RUNTIME_SOURCE_SHA256
    )
    assert (
        identity["sha256"]
        == comparison.REVIEWED_RUNTIME_SOURCE_MANIFEST_SHA256
    )

    current = comparison.runtime_source_metadata()
    reviewed_files = json.loads(
        reviewed["benchmark_runtime_source_manifest_json"]
    )
    current_files = json.loads(current["benchmark_runtime_source_manifest_json"])
    assert {
        name
        for name in reviewed_files
        if reviewed_files[name] != current_files[name]
    } == {
        "examples/cardiac_benchmark/activation.py",
        "examples/cardiac_benchmark/benchmark_parameters.py",
        "examples/cardiac_benchmark/boundary_audit.py",
        "examples/cardiac_benchmark/geometry.py",
        "examples/cardiac_benchmark/robin.py",
        "examples/cardiac_benchmark/run.py",
        "examples/cardiac_benchmark/run_mpi.py",
        "examples/cardiac_benchmark/tbar_laplace.py",
    }

    run = tmp_path / "retained-source.npz"
    _write_run(
        run,
        np.zeros((2, 3, 101)),
        app_revision=comparison.REVIEWED_RESULT_APP_REVISION,
        **reviewed,
    )
    loaded = comparison.load_coupfe_step2b(run)
    assert (
        loaded["run_contract"]["runtime_source_sha256"]
        == comparison.REVIEWED_RUNTIME_SOURCE_SHA256
    )

    for changed_identity in (
        {"app_revision": "0" * 40},
        {"app_tree_state": "clean"},
    ):
        retained_identity = {
            **reviewed,
            "app_revision": comparison.REVIEWED_RESULT_APP_REVISION,
            "app_tree_state": "dirty",
            **changed_identity,
        }
        _write_run(run, np.zeros((2, 3, 101)), **retained_identity)
        with pytest.raises(
            comparison.ComparisonInputError,
            match="reviewed application revision and tree state",
        ):
            comparison.load_coupfe_step2b(run)


def test_retained_runtime_source_manifest_fails_closed_on_byte_change(
    tmp_path,
):
    changed = tmp_path / "changed-runtime-source-hashes.json"
    changed.write_bytes(
        comparison.DEFAULT_RUNTIME_SOURCE_MANIFEST.read_bytes() + b"\n"
    )
    with pytest.raises(
        comparison.ComparisonInputError,
        match="reviewed retained-run manifest",
    ):
        comparison.load_reviewed_runtime_source_manifest(changed)


def test_retained_step2b_plot_report_is_exact_and_exposes_both_error_metrics():
    parsed = plotter._validate_report(RETAINED_STEP2B_REPORT)
    metrics = parsed["metrics"]

    assert metrics["aggregate_relative_l2"] == pytest.approx(
        0.0980377406480768
    )
    assert metrics["p0_relative_l2"] == pytest.approx(0.09055465979197143)
    assert metrics["p1_relative_l2"] == pytest.approx(0.10927093276997242)
    assert metrics["p0_paper_red"] == pytest.approx(0.28300351203509677)
    assert metrics["p1_paper_red"] == pytest.approx(0.35277371701033466)
    assert parsed["points"]["p1"]["uz"]["coupfe_m"][40] < 0.0
    assert parsed["points"]["p1"]["uz"]["publisher_min_m"][40] > 0.0


@pytest.mark.parametrize("failure", ["missing", "unexpected", "changed"])
def test_publisher_selection_and_hashes_fail_closed(tmp_path, monkeypatch, failure):
    data, manifest, _run = _synthetic_inputs(tmp_path, monkeypatch)
    first = data / comparison.PUBLISHER_SELECTION[0][2]
    if failure == "missing":
        first.unlink()
        match = "missing selected"
    elif failure == "unexpected":
        (data / "monoventricular_blinded_B_group_unselected.pkl").write_bytes(b"x")
        match = "unexpected publisher"
    else:
        first.write_bytes(first.read_bytes() + b"changed")
        match = "input changed"

    with pytest.raises(comparison.ComparisonInputError, match=match):
        comparison.load_publisher_reference(data, manifest)


def test_restricted_unpickler_rejects_manifest_approved_global(
    tmp_path, monkeypatch
):
    data, manifest, _run = _synthetic_inputs(tmp_path, monkeypatch)

    class Unsafe:
        def __reduce__(self):
            return eval, ("40 + 2",)

    first = data / comparison.PUBLISHER_SELECTION[0][2]
    with first.open("wb") as stream:
        pickle.dump(Unsafe(), stream, protocol=4)
    _write_manifest(manifest, data)
    _trust_synthetic_manifest(monkeypatch, manifest)

    with pytest.raises(comparison.ComparisonInputError, match="disallowed pickle global"):
        comparison.load_publisher_reference(data, manifest)


@pytest.mark.parametrize(
    ("override", "message"),
    [
        ({"benchmark_step": 0}, "not Benchmark 1 Step 2 Case B"),
        ({"completed_steps": 999}, "complete at 1000/1000"),
        (
            {
                "benchmark_configuration_id": "benchmark-1-step-0-case-B-pressure-only"
            },
            "wrong benchmark configuration identity",
        ),
        ({"integrator": "be"}, "requires generalized-alpha"),
        ({"density": 999.0}, "noncanonical 'density'"),
        (
            {"benchmark_runtime_source_sha256": "0" * 64},
            "reviewed runtime source manifest",
        ),
        (
            {
                "tau": comparison._expected_step2b_load_histories()[0]
                + np.concatenate((np.zeros(300), [1000.0], np.zeros(700)))
            },
            "step diagnostic disagrees",
        ),
    ],
)
def test_coupfe_archive_rejects_incomplete_or_step0_identity(
    tmp_path, override, message
):
    teams = _write_publisher_data(tmp_path / "publisher")
    run = tmp_path / "bad.npz"
    _write_run(run, np.mean(teams, axis=0), **override)

    with pytest.raises(comparison.ComparisonInputError, match=message):
        comparison.load_coupfe_step2b(run)


def _official_data_directory():
    configured = os.environ.get("CARDIAC_BENCHMARK_REFERENCE_DATA")
    if configured:
        return Path(configured)
    return (
        Path.home()
        / "cardiac_benchmark_reference"
        / "benchmark_article_data"
        / "results_time_curves"
        / "data"
    )


def test_local_official_selection_hashes_schema_and_event_oracles():
    data = _official_data_directory()
    if not data.is_dir():
        pytest.skip("local publisher results_time_curves/data directory is absent")

    reference = comparison.load_publisher_reference(data, OFFICIAL_MANIFEST)
    native = np.empty((2, 3, 1001))
    for point_index in range(2):
        for component_index in range(3):
            native[point_index, component_index] = np.interp(
                comparison.COUPFE_TIME_S,
                comparison.PUBLISHED_TIME_S,
                reference["mean_m"][point_index, component_index],
            )
    mirror = {
        "common_displacement_m": reference["mean_m"],
        "displacement_m": native,
        "time_s": comparison.COUPFE_TIME_S,
    }
    events = comparison.contraction_event_metrics(reference, mirror)

    assert reference["teams_m"].shape == (10, 2, 3, 101)
    assert [entry["team_id"] for entry in reference["identities"]] == [
        selected[0] for selected in comparison.PUBLISHER_SELECTION
    ]
    assert reference["excluded_generic_identity"] is not None
    fixed = events["fixed_threshold_crossings"]
    assert fixed["p0"]["-5"]["reference"]["team_mean_s"] == pytest.approx(
        0.19944670252764482
    )
    assert fixed["p0"]["-15"]["reference"]["team_mean_s"] == pytest.approx(
        0.23144924499043493
    )
    assert fixed["p1"]["-5"]["reference"]["team_mean_s"] == pytest.approx(
        0.20316157748621433
    )
    assert fixed["p1"]["-15"]["reference"]["team_mean_s"] == pytest.approx(
        0.24075889669140862
    )
    normalized = events["normalized_drop_crossings"]
    assert normalized["p0"]["90"]["reference"]["team_min_s"] == pytest.approx(
        0.2953122702444225
    )
    assert normalized["p1"]["90"]["reference"]["team_max_s"] == pytest.approx(
        0.3092497603145973
    )
