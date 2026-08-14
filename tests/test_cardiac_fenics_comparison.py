from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

import compare_fenics_case_b as comparison


ROOT = Path(__file__).resolve().parents[1]
KNOWN_FENICS_HASH_MANIFEST = (
    ROOT
    / "examples"
    / "cardiac_benchmark"
    / "fenics_case_b_reference_hashes.example.json"
)


def _fenics_parameters():
    return {
        "activation_parameters": {
            "a_max": 5.0,
            "a_min": -30.0,
            "gamma": 0.005,
            "sigma_0": 150000.0,
            "t_dias": 0.484,
            "t_sys": 0.16,
        },
        "benchmark": 1,
        "case": "b",
        "fiber_parameters": {
            "alpha_endo": -60.0,
            "alpha_epi": 60.0,
            "function_space": "P_2",
        },
        "geometry_path": "/private/campaign/lv_ellipsoid.h5",
        "material_parameters": {
            "a": 59.0,
            "a_f": 18472.0,
            "a_fs": 216.0,
            "a_s": 2481.0,
            "b": 8.023,
            "b_f": 16.026,
            "b_fs": 11.436,
            "b_s": 11.12,
            "eta": 100.0,
            "k": 100.0,
            "kappa": 1000000.0,
        },
        "mesh_parameters": {
            "mesh_size_factor": 1.0,
            "mu_apex_endo": -np.pi,
            "mu_apex_epi": -np.pi,
            "mu_base_endo": -1.2722641256100204,
            "mu_base_epi": -1.318116071652818,
            "psize_ref": 0.005,
            "r_long_endo": 0.09,
            "r_long_epi": 0.097,
            "r_short_endo": 0.025,
            "r_short_epi": 0.035,
        },
        "outdir": "/private/campaign/results",
        "outpath": "/private/campaign/results/result.h5",
        "pressure_parameters": {
            "a_max": 5.0,
            "a_min": -30.0,
            "alpha_mid": 1.0,
            "alpha_pre": 5.0,
            "gamma": 0.005,
            "sigma_mid": 16000.0,
            "sigma_pre": 7000.0,
            "t_dias_pre": 0.484,
            "t_sys_pre": 0.17,
        },
        "problem_parameters": {
            "alpha_epi": 100000000.0,
            "alpha_f": 0.4,
            "alpha_m": 0.2,
            "alpha_top": 100000.0,
            "beta_epi": 5000.0,
            "beta_top": 5000.0,
            "dt": 0.001,
            "function_space": "P_2",
            "p": 0.0,
            "rho": 1000.0,
        },
        "step": 0,
        "timestamp": "2026-06-26T00:41:16.484740",
        "zero_activation": True,
        "zero_pressure": False,
    }


def _histories(times):
    p0 = np.column_stack((0.010 * times, 0.004 * times, -0.025 * times))
    p1 = np.column_stack((0.006 * times, -0.003 * times, -0.020 * times))
    return p0, p1


def _diagnostics(times):
    return [
        {
            "time": float(times[index]),
            "dt": 0.001,
            "initial_residual_norm": 1.0,
            "final_residual_norm": 1.0e-12,
            "residual_acceptance_threshold": 1.0e-10,
            "snes_converged_reason": 2,
            "ksp_converged_reason": 4,
            "nonlinear_iterations": 2,
            "linear_iterations": 2,
        }
        for index in range(1, len(times))
    ]


def _write_inputs(root: Path):
    root.mkdir(parents=True, exist_ok=True)
    times = comparison.COUPFE_TIME.copy()
    fenics_times = comparison.FENICS_TIME.copy()
    fenics_p0, fenics_p1 = _histories(fenics_times)
    run_p0, run_p1 = _histories(times)
    run_p0[:, 0] += 1.0e-3
    run_p1[:, 1] += 2.0e-3

    audit = {
        "geometry": {
            "passed": True,
            "mesh_topology": "closed_multiblock_disk",
            "unclassified_exterior_faces": 0,
            "nonpositive_extended_jacobians": 0,
        },
        "pressure": {"passed": True},
        "robin": {"passed": True},
    }
    solver_configuration = {
        "name": "petsc-snes",
        "snes_type": "newtonls",
        "line_search_type": "bt",
        "ksp_type": "preonly",
        "pc_type": "lu",
        "function_domain_rejection_api": "nonfinite residual for PETSc BT",
        "element_evaluation_mode": "joint",
        "compiled_material_residual_only_available": True,
    }
    run = root / "case_b.npz"
    np.savez(
        run,
        result_schema="coupfe-cardiac-result-v1",
        converged=True,
        completed_steps=1000,
        expected_steps=1000,
        case="B",
        times=times,
        u0=run_p0,
        u1=run_p1,
        pres=comparison.p_of_t(times),
        tau=np.zeros_like(times),
        p0=np.array([0.025, 0.030, 0.0]),
        p1=np.array([0.000, 0.030, 0.0]),
        p0_sampling_reconstruction_error_m=1.0e-16,
        p1_sampling_reconstruction_error_m=2.0e-16,
        app_revision="a" * 40,
        app_tree_state="clean",
        app_source_kind="git-checkout",
        core_revision="b" * 40,
        core_tree_state="clean",
        core_source_kind="git-checkout",
        core_source_url="https://github.com/tengzhang48/CoupFE.git",
        integrator="be",
        nonlinear_solver="petsc-snes",
        formulation="hex8_standard_pointwise_kappa",
        material_kernel_formulation="standard",
        material_model_id=comparison.MATERIAL_MODEL_ID,
        mass_representation="consistent_q1_hex8",
        fiber_sampling="gp_direct_rule",
        fiber_sampling_option="gp-direct",
        mesh_topology="closed_multiblock_disk",
        point_sampling="hex8_reference_isoparametric",
        viscous_rate="backward_difference",
        parameter_variant="benchmark_eta",
        tbar_definition="laplace_presolved",
        tbar_source_filename="tbar.npy",
        tbar_source_sha256="c" * 64,
        tbar_metadata_filename="tbar.meta.json",
        tbar_metadata_sha256="d" * 64,
        tbar_metadata_schema="coupfe-cardiac-laplace-tbar-v1",
        dt=0.001,
        t_end=1.0,
        load_horizon=1.0,
        density=1000.0,
        material_eta_pa_s=100.0,
        material_kappa_pa=1.0e6,
        local_pressure_bulk_modulus_pa=0.0,
        apex_offset=0.0,
        perturb=0.0,
        a_top=1.0e5,
        b_top=5.0e3,
        a_epi=1.0e8,
        b_epi=5.0e3,
        isotropic=False,
        viscous_term_active=True,
        flip_helix=True,
        driver="examples/cardiac_benchmark/run.py",
        element_evaluation_mode="joint",
        compiled_material_residual_only_available=True,
        nodes=np.array(
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
        ),
        elems=np.arange(8, dtype=np.int64).reshape(1, 8),
        det_f_gauss_peak=np.ones((1, 8)),
        pre_solve_audit_json=json.dumps(audit, sort_keys=True),
        solver_configuration_json=json.dumps(solver_configuration, sort_keys=True),
        nonlinear_step_diagnostics_json=json.dumps(_diagnostics(times)),
        n_t=1,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )

    parameters = root / "parameters.json"
    parameters.write_text(json.dumps(_fenics_parameters()), encoding="utf-8")
    fenics_times_path = root / "time_stamps.npy"
    fenics_p0_path = root / "componentwise_displacement_up0.npy"
    fenics_p1_path = root / "componentwise_displacement_up1.npy"
    np.save(fenics_times_path, fenics_times)
    np.save(fenics_p0_path, fenics_p0)
    np.save(fenics_p1_path, fenics_p1)
    return {
        "coupfe-run": run,
        "fenics-parameters": parameters,
        "fenics-times": fenics_times_path,
        "fenics-p0": fenics_p0_path,
        "fenics-p1": fenics_p1_path,
    }


def _identities(paths):
    return {role: comparison._identity(path) for role, path in paths.items()}


def _load(paths):
    identities = _identities(paths)
    coupfe = comparison.load_coupfe_run(paths["coupfe-run"], identities["coupfe-run"])
    fenics = comparison.load_fenics_reference(
        paths["fenics-parameters"],
        paths["fenics-times"],
        paths["fenics-p0"],
        paths["fenics-p1"],
        identities,
    )
    return coupfe, fenics


def _rewrite_npz(path: Path, **updates):
    with np.load(path, allow_pickle=False) as archive:
        payload = {key: archive[key] for key in archive.files}
    payload.update(updates)
    np.savez(path, **payload)


def _remove_npz_field(path: Path, field: str):
    with np.load(path, allow_pickle=False) as archive:
        payload = {key: archive[key] for key in archive.files if key != field}
    np.savez(path, **payload)


def _convert_run_to_mpi(path: Path, ranks: int = 2):
    with np.load(path, allow_pickle=False) as archive:
        payload = {key: archive[key] for key in archive.files}
    diagnostics = json.loads(str(payload["nonlinear_step_diagnostics_json"]))
    for record in diagnostics:
        record["ranks"] = ranks
    configuration = {
        "name": "petsc-snes-mpi",
        "ranks": ranks,
        "line_search_configuration_api": "namespaced PETSc option",
        "petsc4py_version": "3.18.4",
        "petsc_version": "3.18.4",
        "mass_representation": "consistent_q1_hex8",
        "mass_partition": "owned-row-csr-all-touching-elements",
        "mass_owned_row_range": [0, 12],
        "mass_local_nnz": 64,
        "element_evaluation_mode": "joint",
        "compiled_material_residual_only_available": True,
        **comparison.MPI_FIXED_CONFIGURATION,
    }
    payload.update(
        {
            "nonlinear_solver": "petsc-snes-mpi",
            "driver": "examples/cardiac_benchmark/run_mpi.py",
            "solver_configuration_json": json.dumps(configuration, sort_keys=True),
            "nonlinear_step_diagnostics_json": json.dumps(diagnostics),
            "mpi_enabled": True,
            "mpi_ranks": ranks,
            "mpi_world_size": ranks,
            "mpi_local_element_counts": np.array([1, 0], dtype=np.int64),
            "mpi_implementation": comparison.MPI_IMPLEMENTATION,
            "mpi_partition": "coupfe.partition_elements",
            "mpi_build_layout": "isolated-rank-directories",
            "mpi_factor_solver_type": "superlu_dist",
            "mpi_mass_partition": "owned-row-csr-all-touching-elements",
            "mpi_mass_owned_row_ranges": np.array([[0, 12], [12, 24]], dtype=np.int64),
            "mpi_mass_local_nnz": np.array([64, 64], dtype=np.int64),
            "mpi_mass_touching_element_counts": np.array([1, 1], dtype=np.int64),
        }
    )
    np.savez(path, **payload)


def _convert_run_to_generalized_alpha_mpi(path: Path, ranks: int = 2):
    _convert_run_to_mpi(path, ranks=ranks)
    with np.load(path, allow_pickle=False) as archive:
        payload = {key: archive[key] for key in archive.files}
    configuration = json.loads(str(payload["solver_configuration_json"]))
    configuration.update(
        {
            "implementation": comparison.MPI_GENERALIZED_ALPHA_IMPLEMENTATION,
            "time_integrator": "generalized-alpha",
            "material_batch_time_integrator": (
                "generalized-alpha-source-matched"
            ),
            "compiled_material_dt": comparison.EXPECTED_DT_S,
            "acceleration_stage": "1-alpha_m",
            "force_stage": "1-alpha_f",
            "material_viscous_rate": "sym(F_stage^T*grad(v_stage))",
            "nonlinear_initial_guess": "accepted-u_n-like-simula",
            "generalized_alpha": comparison.GENERALIZED_ALPHA_CONFIGURATION,
        }
    )
    load_times = comparison.COUPFE_TIME.copy()
    load_times[0] = 0.0
    load_times[1:] -= 0.4 * comparison.EXPECTED_DT_S
    pressure = np.zeros_like(load_times)
    pressure[1:] = comparison.p_of_t(
        load_times[1:], t_span=(0.0, 1.0)
    )
    selected = comparison.benchmark_configuration(0, "B")
    payload.update(
        comparison.benchmark_metadata(
            selected,
            material_parameters=selected.material_parameters,
            activation_parameters=selected.activation_parameters,
            pressure_parameters=selected.pressure_parameters,
        )
    )
    payload.update(
        {
            "integrator": "generalized-alpha",
            "viscous_rate": (
                "velocity_consistent_green_lagrange_at_alpha_f_stage"
            ),
            "generalized_alpha_alpha_m": 0.2,
            "generalized_alpha_alpha_f": 0.4,
            "generalized_alpha_gamma": 0.7,
            "generalized_alpha_beta": 0.36,
            "generalized_alpha_stage_contract": "simula-source-matched-v1",
            "load_evaluation_times_s": load_times,
            "pres": pressure,
            "tau": np.zeros_like(load_times),
            "solver_configuration_json": json.dumps(
                configuration, sort_keys=True
            ),
            "mpi_implementation": (
                comparison.MPI_GENERALIZED_ALPHA_IMPLEMENTATION
            ),
        }
    )
    np.savez(path, **payload)


def test_report_metrics_common_grid_and_path_sanitization(tmp_path):
    paths = _write_inputs(tmp_path / "inputs")
    coupfe, fenics = _load(paths)
    report, figure = comparison.build_report(coupfe, fenics, None)

    assert report["schema"] == comparison.REPORT_SCHEMA
    common = report["common_grid"]
    assert common["samples"] == 999
    assert common["start_s"] == pytest.approx(0.001)
    assert common["end_s"] == pytest.approx(0.999)
    assert common["coupfe_mapping"] == {
        "p0": {"identity_sampling": True, "method": "identity_index_selection"},
        "p1": {"identity_sampling": True, "method": "identity_index_selection"},
    }

    p0 = report["comparison"]["p0"]
    p1 = report["comparison"]["p1"]
    assert p0["phases"]["full_0p001_to_0p999_s"]["vector_rmse_mm"] == pytest.approx(1.0)
    assert p1["phases"]["full_0p001_to_0p999_s"]["vector_rmse_mm"] == pytest.approx(2.0)
    assert p0["snap_onset"]["candidate_time_s"] == pytest.approx(0.2)
    assert p0["snap_onset"]["reference_time_s"] == pytest.approx(0.2)
    assert p1["snap_onset"]["candidate_time_s"] == pytest.approx(0.25)
    assert p1["snap_onset"]["reference_time_s"] == pytest.approx(0.25)

    expected_relative = np.linalg.norm(
        np.tile([1.0e-3, 0.0, 0.0], (999, 1))
    ) / np.linalg.norm(fenics["histories"]["p0"])
    assert p0["phases"]["full_0p001_to_0p999_s"]["relative_l2"] == pytest.approx(
        expected_relative
    )

    public_json = json.dumps(report)
    assert str(tmp_path) not in public_json
    assert "/private/campaign" not in public_json
    assert report["configuration"]["fenics"]["geometry_filename"] == "lv_ellipsoid.h5"
    assert report["configuration"]["fenics_source_identity"]["producing_revision"] == (
        "not-retained-contemporaneously"
    )
    assert report["evidence_status"]["public_retained_candidate"] is False
    assert figure["title"].startswith("DEVELOPMENT")
    assert figure["retention_note"] == "Development output—not a retained public result."
    assert report["inputs"]["hash_verification"] == {
        "mode": "development",
        "all_input_roles_required": False,
        "all_input_roles_verified": False,
        "verified_inputs": {},
        "unverified_roles": list(comparison.INPUT_ROLES),
        "interpretation": (
            "Development comparison: no caller-expected hashes were supplied; "
            "this report is not a retained public result."
        ),
    }


def test_windows_geometry_path_is_reduced_to_public_basename(tmp_path):
    paths = _write_inputs(tmp_path / "inputs")
    parameters = json.loads(paths["fenics-parameters"].read_text(encoding="utf-8"))
    parameters["geometry_path"] = "C:" + "\\Users\\Alice\\private\\lv_ellipsoid.h5"
    paths["fenics-parameters"].write_text(json.dumps(parameters), encoding="utf-8")

    coupfe, fenics = _load(paths)
    report, _figure = comparison.build_report(coupfe, fenics, None)
    public_json = json.dumps(report)
    assert report["configuration"]["fenics"]["geometry_filename"] == (
        "lv_ellipsoid.h5"
    )
    assert "C:" + "\\Users" not in public_json
    assert "Alice" not in public_json


@pytest.mark.parametrize(
    "field,filename",
    [
        ("tbar_source_filename", "tbar.npy"),
        ("tbar_metadata_filename", "tbar.meta.json"),
    ],
)
def test_windows_tbar_paths_are_rejected_instead_of_recorded(tmp_path, field, filename):
    paths = _write_inputs(tmp_path / field)
    value = "C:" + "\\Users\\Alice\\" + filename
    _rewrite_npz(paths["coupfe-run"], **{field: value})
    with pytest.raises(comparison.ComparisonInputError, match="not a portable basename"):
        _load(paths)


def test_cli_checks_manifest_and_writes_atomic_public_outputs(tmp_path, capsys):
    pytest.importorskip("matplotlib")
    paths = _write_inputs(tmp_path / "inputs")
    manifest = {
        "schema": comparison.HASH_MANIFEST_SCHEMA,
        "files": {
            role: comparison._identity(path)
            for role, path in paths.items()
            if role.startswith("fenics-")
        },
    }
    manifest_path = tmp_path / "known-reference-hashes.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    report_path = tmp_path / "public" / "comparison.json"
    figure_path = tmp_path / "public" / "comparison.png"
    report = comparison.main(
        [
            "--retained",
            "--coupfe-run",
            str(paths["coupfe-run"]),
            "--fenics-parameters",
            str(paths["fenics-parameters"]),
            "--fenics-times",
            str(paths["fenics-times"]),
            "--fenics-p0",
            str(paths["fenics-p0"]),
            "--fenics-p1",
            str(paths["fenics-p1"]),
            "--expected-hashes",
            str(manifest_path),
            "--expect-sha256",
            "coupfe-run=" + comparison._identity(paths["coupfe-run"])["sha256"],
            "--report",
            str(report_path),
            "--figure",
            str(figure_path),
        ]
    )
    assert report_path.is_file()
    assert figure_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert report["inputs"]["expected_hash_manifest"]["filename"] == manifest_path.name
    assert report["inputs"]["coupfe_run"]["expected_sha256_checked"] is True
    assert all(
        record["expected_sha256_checked"] is True
        for record in report["inputs"]["fenics"].values()
    )
    verification = report["inputs"]["hash_verification"]
    assert verification["mode"] == "retained"
    assert verification["all_input_roles_required"] is True
    assert verification["all_input_roles_verified"] is True
    assert verification["unverified_roles"] == []
    assert set(verification["verified_inputs"]) == set(comparison.INPUT_ROLES)
    assert all(
        record["sha256"] == comparison._identity(paths[role])["sha256"]
        for role, record in verification["verified_inputs"].items()
    )
    assert report["evidence_status"]["public_retained_candidate"] is True
    text = report_path.read_text(encoding="utf-8")
    image = figure_path.read_bytes()
    assert str(tmp_path) not in text
    assert str(tmp_path).encode() not in image
    assert b"/private/campaign" not in image
    assert not list((tmp_path / "public").glob(".*.tmp"))
    stdout = capsys.readouterr().out
    assert str(tmp_path) not in stdout
    assert "saved report -> comparison.json" in stdout
    assert "saved figure -> comparison.png" in stdout


def test_retained_mode_rejects_a_missing_expected_role_without_outputs(tmp_path):
    paths = _write_inputs(tmp_path / "inputs")
    manifest = {
        "schema": comparison.HASH_MANIFEST_SCHEMA,
        "files": {
            role: comparison._identity(path)
            for role, path in paths.items()
            if role.startswith("fenics-")
        },
    }
    manifest_path = tmp_path / "partial-hashes.json"
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    report_path = tmp_path / "comparison.json"
    figure_path = tmp_path / "comparison.png"
    with pytest.raises(
        comparison.ComparisonInputError,
        match="retained mode requires expected SHA-256.*coupfe-run",
    ):
        comparison.main(
            [
                "--retained",
                "--coupfe-run", str(paths["coupfe-run"]),
                "--fenics-parameters", str(paths["fenics-parameters"]),
                "--fenics-times", str(paths["fenics-times"]),
                "--fenics-p0", str(paths["fenics-p0"]),
                "--fenics-p1", str(paths["fenics-p1"]),
                "--expected-hashes", str(manifest_path),
                "--report", str(report_path),
                "--figure", str(figure_path),
            ]
        )
    assert not report_path.exists()
    assert not figure_path.exists()


def test_dirty_coupfe_source_is_rejected_without_outputs(tmp_path):
    paths = _write_inputs(tmp_path / "inputs")
    _rewrite_npz(paths["coupfe-run"], app_tree_state="dirty")
    report_path = tmp_path / "report.json"
    figure_path = tmp_path / "figure.png"
    with pytest.raises(comparison.ComparisonInputError, match="requires clean source"):
        comparison.main(
            [
                "--coupfe-run", str(paths["coupfe-run"]),
                "--fenics-parameters", str(paths["fenics-parameters"]),
                "--fenics-times", str(paths["fenics-times"]),
                "--fenics-p0", str(paths["fenics-p0"]),
                "--fenics-p1", str(paths["fenics-p1"]),
                "--report", str(report_path),
                "--figure", str(figure_path),
            ]
        )
    assert not report_path.exists()
    assert not figure_path.exists()


@pytest.mark.parametrize(
    "source_url",
    [
        "https://" + "token@" + "github.com/tengzhang48/CoupFE.git",
        "https://example.invalid/tengzhang48/CoupFE.git",
    ],
)
def test_nonpublic_or_credential_bearing_core_url_is_rejected(tmp_path, source_url):
    paths = _write_inputs(tmp_path / "inputs")
    _rewrite_npz(paths["coupfe-run"], core_source_url=source_url)
    with pytest.raises(comparison.ComparisonInputError, match="unexpected Core source URL"):
        _load(paths)


def test_missing_core_url_is_rejected(tmp_path):
    paths = _write_inputs(tmp_path / "inputs")
    _remove_npz_field(paths["coupfe-run"], "core_source_url")
    with pytest.raises(comparison.ComparisonInputError, match="missing scalar 'core_source_url'"):
        _load(paths)


def test_closed_case_b_mpi_metadata_is_validated_and_reported(tmp_path):
    paths = _write_inputs(tmp_path / "mpi")
    _convert_run_to_mpi(paths["coupfe-run"])
    coupfe, fenics = _load(paths)
    report, _ = comparison.build_report(coupfe, fenics, None)
    mpi = report["configuration"]["coupfe"]["mpi"]
    assert mpi["implementation"] == comparison.MPI_IMPLEMENTATION
    assert mpi["world_size"] == 2
    assert mpi["local_element_counts"] == [1, 0]
    assert mpi["mass_partition"] == {
        "policy": "owned-row-csr-all-touching-elements",
        "owned_row_ranges": [[0, 12], [12, 24]],
        "local_nnz": [64, 64],
        "touching_element_counts": [1, 1],
    }


def test_step0_case_b_generalized_alpha_is_validated_and_labeled(tmp_path):
    paths = _write_inputs(tmp_path / "mpi-ga")
    _convert_run_to_generalized_alpha_mpi(paths["coupfe-run"])
    coupfe, fenics = _load(paths)
    report, figure = comparison.build_report(coupfe, fenics, None)

    configuration = report["configuration"]["coupfe"]
    assert configuration["benchmark_identity"]["step"] == 0
    assert configuration["benchmark_identity"]["load_contract"] == "pressure-only"
    assert configuration["integrator"] == "generalized-alpha"
    assert configuration["time_integrator"] == (
        "source-matched generalized-alpha"
    )
    assert configuration["generalized_alpha"] == {
        "alpha_m": 0.2,
        "alpha_f": 0.4,
        "gamma": 0.7,
        "beta": 0.36,
        "stage_contract": "simula-source-matched-v1",
        "load_time": "t_np1 - alpha_f*dt",
    }
    assert configuration["mpi"]["implementation"] == (
        comparison.MPI_GENERALIZED_ALPHA_IMPLEMENTATION
    )
    assert figure["source_labels"]["coupfe"] == (
        "CoupFE Q1-Hex8, source-matched generalized-alpha"
    )


@pytest.mark.parametrize(
    ("updates", "message"),
    [
        (
            {"generalized_alpha_alpha_f": 0.3},
            "CoupFE generalized-alpha alpha_f",
        ),
        (
            {"mpi_implementation": comparison.MPI_IMPLEMENTATION},
            "time-integration contract",
        ),
        (
            {"benchmark_step": 2},
            "not Benchmark 1 Step 0 Case B",
        ),
    ],
)
def test_step0_case_b_generalized_alpha_broken_controls_are_rejected(
    tmp_path, updates, message
):
    paths = _write_inputs(tmp_path / message.replace(" ", "-"))
    _convert_run_to_generalized_alpha_mpi(paths["coupfe-run"])
    _rewrite_npz(paths["coupfe-run"], **updates)
    with pytest.raises(comparison.ComparisonInputError, match=message):
        _load(paths)


def test_step0_case_b_generalized_alpha_rejects_endpoint_pressure(tmp_path):
    paths = _write_inputs(tmp_path / "ga-endpoint-pressure")
    _convert_run_to_generalized_alpha_mpi(paths["coupfe-run"])
    _rewrite_npz(
        paths["coupfe-run"],
        pres=comparison.p_of_t(comparison.COUPFE_TIME),
    )
    with pytest.raises(
        comparison.ComparisonInputError, match="selected load stage"
    ):
        _load(paths)


@pytest.mark.parametrize("mutation", ["missing-mass", "diagnostic-ranks"])
def test_mpi_metadata_missing_or_inconsistent_is_rejected(tmp_path, mutation):
    paths = _write_inputs(tmp_path / mutation)
    _convert_run_to_mpi(paths["coupfe-run"])
    if mutation == "missing-mass":
        _remove_npz_field(paths["coupfe-run"], "mpi_mass_local_nnz")
        match = "incomplete MPI provenance"
    else:
        with np.load(paths["coupfe-run"], allow_pickle=False) as archive:
            diagnostics = json.loads(str(archive["nonlinear_step_diagnostics_json"]))
        diagnostics[0]["ranks"] = 1
        _rewrite_npz(
            paths["coupfe-run"],
            nonlinear_step_diagnostics_json=json.dumps(diagnostics),
        )
        match = "diagnostic rank count disagrees"
    with pytest.raises(comparison.ComparisonInputError, match=match):
        _load(paths)


def test_incomplete_coupfe_archive_is_rejected(tmp_path):
    paths = _write_inputs(tmp_path / "inputs")
    _rewrite_npz(paths["coupfe-run"], completed_steps=999)
    with pytest.raises(comparison.ComparisonInputError, match="incomplete"):
        _load(paths)


@pytest.mark.parametrize(
    "mutation",
    ["missing", "holzapfel-ogden-smooth-switch-stress-without-switch-derivative-v0"],
)
def test_coupfe_material_model_id_fails_closed(tmp_path, mutation):
    paths = _write_inputs(tmp_path / mutation.replace("/", "_"))
    if mutation == "missing":
        _remove_npz_field(paths["coupfe-run"], "material_model_id")
        match = "missing scalar 'material_model_id'"
    else:
        _rewrite_npz(paths["coupfe-run"], material_model_id=mutation)
        match = "material_model_id"
    with pytest.raises(comparison.ComparisonInputError, match=match):
        _load(paths)


def test_fenics_object_array_is_rejected_without_pickle(tmp_path):
    paths = _write_inputs(tmp_path / "inputs")
    np.save(paths["fenics-p0"], np.array([{"unsafe": True}], dtype=object))
    with pytest.raises(comparison.ComparisonInputError, match="without pickle"):
        _load(paths)


def test_fenics_shape_time_and_nonfinite_data_are_rejected(tmp_path):
    paths = _write_inputs(tmp_path / "shape")
    np.save(paths["fenics-p1"], np.zeros((998, 3)))
    with pytest.raises(comparison.ComparisonInputError, match="shape"):
        _load(paths)

    paths = _write_inputs(tmp_path / "time")
    bad_times = comparison.FENICS_TIME.copy()
    bad_times[500] += 1.0e-5
    np.save(paths["fenics-times"], bad_times)
    with pytest.raises(comparison.ComparisonInputError, match="retained 0.001--0.999"):
        _load(paths)

    paths = _write_inputs(tmp_path / "finite")
    p0, _ = _histories(comparison.FENICS_TIME)
    p0[12, 1] = np.nan
    np.save(paths["fenics-p0"], p0)
    with pytest.raises(comparison.ComparisonInputError, match="non-finite"):
        _load(paths)


def test_fenics_physical_metadata_mismatch_is_rejected(tmp_path):
    paths = _write_inputs(tmp_path / "inputs")
    parameters = json.loads(paths["fenics-parameters"].read_text(encoding="utf-8"))
    parameters["material_parameters"]["eta"] = 99.0
    paths["fenics-parameters"].write_text(json.dumps(parameters), encoding="utf-8")
    with pytest.raises(comparison.ComparisonInputError, match="material_parameters.eta"):
        _load(paths)


@pytest.mark.parametrize(
    "timestamp",
    [
        "/private/campaign/results",
        "2026-06-26",
        "2026-06-26T00:41:16+00:00",
    ],
)
def test_fenics_timestamp_must_match_the_recorded_naive_iso_form(
    tmp_path, timestamp
):
    paths = _write_inputs(tmp_path / "inputs")
    parameters = json.loads(paths["fenics-parameters"].read_text(encoding="utf-8"))
    parameters["timestamp"] = timestamp
    paths["fenics-parameters"].write_text(json.dumps(parameters), encoding="utf-8")
    with pytest.raises(comparison.ComparisonInputError, match="recorded naive ISO-8601"):
        _load(paths)


def test_nonbenchmark_pressure_schedule_is_rejected(tmp_path):
    paths = _write_inputs(tmp_path / "inputs")
    _rewrite_npz(
        paths["coupfe-run"],
        pres=16000.0 * np.sin(np.pi * comparison.COUPFE_TIME) ** 2,
    )
    with pytest.raises(comparison.ComparisonInputError, match="benchmark schedule"):
        _load(paths)


def test_windows_looking_input_basename_is_rejected(tmp_path):
    source = tmp_path / "ordinary.npy"
    np.save(source, np.zeros(1))
    windows_looking = tmp_path / ("C:" + "\\Users\\Alice\\secret.npy")
    source.rename(windows_looking)
    with pytest.raises(comparison.ComparisonInputError, match="portable basename"):
        comparison._identity(windows_looking)


def test_expected_hash_mismatch_and_unknown_role_fail_closed(tmp_path):
    paths = _write_inputs(tmp_path / "inputs")
    identities = _identities(paths)
    with pytest.raises(comparison.ComparisonInputError, match="SHA-256 mismatch"):
        comparison._validate_expected_hashes(
            identities, {"fenics-p0": {"sha256": "0" * 64}}
        )
    with pytest.raises(comparison.ComparisonInputError, match="unknown expected-hash role"):
        comparison.load_expected_hashes(None, ["not-a-role=" + "0" * 64])


def test_known_input_hashes_are_recorded_by_role(tmp_path):
    paths = _write_inputs(tmp_path / "inputs")
    identities = _identities(paths)
    for role, path in paths.items():
        expected = hashlib.sha256(path.read_bytes()).hexdigest()
        assert identities[role] == {
            "filename": path.name,
            "sha256": expected,
            "size_bytes": path.stat().st_size,
        }


def test_checked_known_fenics_manifest_has_exact_safe_file_identities():
    expected, manifest_identity = comparison.load_expected_hashes(
        KNOWN_FENICS_HASH_MANIFEST, []
    )
    assert manifest_identity["filename"] == KNOWN_FENICS_HASH_MANIFEST.name
    assert set(expected) == {
        "fenics-parameters",
        "fenics-times",
        "fenics-p0",
        "fenics-p1",
    }
    assert expected["fenics-parameters"] == {
        "filename": "parameters.json",
        "sha256": "c1cd4c8d2521fd6c28774975843740a8af12568edd1240f5daa133d469e6fb76",
        "size_bytes": 1385,
    }
    assert expected["fenics-times"]["sha256"] == (
        "ddba330b1c8f8c1bb61282e187047f3aa99d0df37b2c4ed2139ea1b0e0ff0f0c"
    )
    assert expected["fenics-p0"]["sha256"] == (
        "4344a4f599a6eabb16159682339a735bff572eaa18eedd1fe2a97ebd3ee7f4a0"
    )
    assert expected["fenics-p1"]["sha256"] == (
        "88a679de2189bc137de5d64186c698f1702e9df20333849377cdfd01aac8bf1e"
    )


def test_onset_uses_fixed_threshold_and_linear_interpolation():
    times = np.array([0.1, 0.2, 0.3])
    history = np.zeros((3, 3))
    history[:, 2] = np.array([-0.004, -0.0045, -0.0055])
    assert comparison.onset_time(times, history) == pytest.approx(0.25)
    history[:, 2] = -0.004
    assert comparison.onset_time(times, history) is None


def test_mapping_reports_identity_or_interpolation_without_extrapolation():
    source_t = np.array([0.0, 0.1, 0.2])
    source_u = np.column_stack((source_t, 2.0 * source_t, -source_t))
    mapped, identity, method = comparison.map_to_common_grid(
        source_t, source_u, np.array([0.1, 0.2])
    )
    assert identity is True
    assert method == "identity_index_selection"
    np.testing.assert_array_equal(mapped, source_u[1:])

    mapped, identity, method = comparison.map_to_common_grid(
        source_t, source_u, np.array([0.05, 0.15])
    )
    assert identity is False
    assert method == "linear_interpolation_without_extrapolation"
    np.testing.assert_allclose(mapped[:, 0], [0.05, 0.15])

    # Decimal 1 ms ticks produced by linspace and arange can differ by one
    # floating-point ulp. They still select the same retained samples and must
    # not be mislabeled as interpolated data.
    source_t = np.linspace(0.0, 1.0, 1001)
    target_t = np.arange(1, 1000, dtype=float) * 0.001
    source_u = np.column_stack((source_t, source_t, source_t))
    mapped, identity, method = comparison.map_to_common_grid(
        source_t, source_u, target_t
    )
    assert identity is True
    assert method == "identity_index_selection"
    np.testing.assert_array_equal(mapped, source_u[1:-1])
