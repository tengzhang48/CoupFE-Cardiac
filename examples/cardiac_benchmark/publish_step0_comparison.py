"""Publish a compact, hash-bound Benchmark 1 Step 0 comparison report.

``post.py`` emits a forensic generic-v2 report containing all accepted-step
and KSP histories.  This tool verifies that report and its external campaign
manifest/stdout, then retains only the configuration, completion evidence,
pre-solve audits, reference identities, 101-point curves, and RED values used
by the public figure.  The NPZ, forensic report, manifest, and stdout remain
external and are bound by filename, size, and SHA-256.

The output is deterministic strict JSON smaller than 100 kB.  It deliberately
omits per-step diagnostics, raw result histories, and generalized-alpha load
evaluation arrays.  It supports explicit recorded Step 0 Cases A/B; retained
release guards still pin the exact selected run.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Dict, Mapping, Tuple

import numpy as np

try:  # package import
    from .activation import p_of_t, tau_of_t
    from .benchmark_parameters import benchmark_configuration, benchmark_metadata
    from . import post
except ImportError:  # direct example/script import
    from activation import p_of_t, tau_of_t
    from benchmark_parameters import benchmark_configuration, benchmark_metadata
    import post


GENERIC_REPORT_SCHEMA = "coupfe-cardiac-reference-comparison-v2"
COMPACT_REPORT_SCHEMA = "coupfe-cardiac-step0-retained-comparison-v1"
CAMPAIGN_MANIFEST_SCHEMA = "coupfe-cardiac-external-run-manifest-v1"
MAX_COMPACT_REPORT_BYTES = 100_000
DT_S = 0.001
T_END_S = 1.0
MPI_RANKS = 8
LINEAR_SOLVER_PROFILE = "fgmres-gamg-rigid-rebuild"
FORMULATION = "hex8_local_pressure_p0_condensed_logj"
GENERALIZED_ALPHA = {
    "alpha_f": 0.4,
    "alpha_m": 0.2,
    "beta": 0.36,
    "gamma": 0.7,
    "stage_contract": "simula-source-matched-v1",
}
MPI_IMPLEMENTATION = (
    "cardiac-owned-distributed-closed-local-pressure-generalized-alpha-step0"
)
MPI_CONTRACT = {
    "A": "closed_case_a_local_pressure_consistent_generalized_alpha",
    "B": "closed_case_b_local_pressure_consistent_generalized_alpha",
}
PRECONDITIONER = (
    "PETSc GAMG aggregation with six rigid-body near-null modes and "
    "interpolation rebuilt for changed matrices"
)
MODEL_PARAMETERS = {
    "base_robin_damping_pa_s_m": 5000.0,
    "base_robin_stiffness_pa_m": 100000.0,
    "density_kg_m3": 1000.0,
    "epicardial_robin_damping_pa_s_m": 5000.0,
    "epicardial_robin_stiffness_pa_m": 100000000.0,
    "local_pressure_bulk_modulus_pa": 1000000.0,
    "local_pressure_volume_law": "linear-reference-volume-mean-log-j-v1",
    "material_kappa_pa": 0.0,
    "material_kernel_formulation": "standard",
    "material_model_id": (
        "holzapfel-ogden-smooth-switch-complete-energy-derivative-v1"
    ),
    "mesh_perturbation_std_m": 0.0,
}
CONFIGURATION_FIELDS = (
    "apex_offset_rad",
    "dt_s",
    "fiber_sampling",
    "fiber_sampling_option",
    "flip_helix",
    "formulation",
    "integrator",
    "isotropic",
    "load_horizon_origin",
    "load_horizon_s",
    "mass_representation",
    "material_eta_pa_s",
    "mesh",
    "method_metadata_origin",
    "model_parameters",
    "nonlinear_solver",
    "parameter_variant",
    "point_sampling",
    "sampling_points",
    "t_end_s",
    "tbar",
    "viscous_rate",
    "viscous_term_active",
)
PHYSICAL_POINTS_M = {"p0": [0.025, 0.030, 0.0], "p1": [0.0, 0.030, 0.0]}


class PublicationError(RuntimeError):
    """Input evidence is not safe to publish under the compact schema."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PublicationError(message)


def _reject_constant(value: str) -> None:
    raise ValueError("non-finite JSON constant {!r}".format(value))


def _read_json(path: Path, label: str) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    try:
        payload = path.read_bytes()
        value = json.loads(
            payload.decode("utf-8"), parse_constant=_reject_constant
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise PublicationError("cannot read strict {}: {}".format(label, error)) from error
    _require(isinstance(value, dict), "{} must contain an object".format(label))
    return value, {
        "filename": path.name,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size_bytes": len(payload),
    }


def _identity(path: Path) -> Tuple[Dict[str, Any], bytes]:
    try:
        payload = path.read_bytes()
    except OSError as error:
        raise PublicationError("cannot read {}: {}".format(path, error)) from error
    return {
        "filename": path.name,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size_bytes": len(payload),
    }, payload


def _finite_array(value: Any, shape: Tuple[int, ...], label: str) -> np.ndarray:
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as error:
        raise PublicationError("{} is not numeric".format(label)) from error
    _require(array.shape == shape, "{} has shape {}, expected {}".format(label, array.shape, shape))
    _require(np.all(np.isfinite(array)), "{} contains non-finite values".format(label))
    return array


def _sha(value: Any) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def _source_identity(value: Any) -> Dict[str, Any]:
    _require(isinstance(value, dict), "source identity must be an object")
    app = value.get("app", {})
    core = value.get("core", {})
    revisions = (app.get("revision"), core.get("revision"))
    _require(
        all(isinstance(item, str) and re.fullmatch(r"[0-9a-f]{40}", item) for item in revisions),
        "source revisions are malformed",
    )
    _require(
        app.get("source_kind") == "git-checkout" and app.get("tree_state") == "clean",
        "application source is not a clean Git checkout",
    )
    _require(
        core.get("source_kind") == "git-checkout"
        and core.get("tree_state") == "clean"
        and core.get("source_url") == post.PUBLIC_CORE_URL,
        "Core source is not the clean public checkout",
    )
    return copy.deepcopy(value)


def _benchmark_identity(value: Any, case: str) -> Tuple[Any, Dict[str, Any]]:
    _require(isinstance(value, dict), "recorded benchmark identity is missing")
    selected = benchmark_configuration(0, case)
    metadata = benchmark_metadata(
        selected,
        material_parameters=value.get("material_parameters", {}),
        activation_parameters=value.get("activation_parameters", {}),
        pressure_parameters=value.get("pressure_parameters", {}),
    )
    expected = {
        "benchmark": 1,
        "step": 0,
        "case": case,
        "configuration_id": selected.identity,
        "identity_scope": metadata["benchmark_identity_scope"],
        "load_contract": selected.load_contract,
        "peak_load_definition": metadata["benchmark_peak_load_definition"],
        "active_stress_enabled": selected.active_stress_enabled,
        "pressure_enabled": selected.pressure_enabled,
        "material_parameters": dict(selected.material_parameters),
        "activation_parameters": dict(selected.activation_parameters),
        "pressure_parameters": dict(selected.pressure_parameters),
        "runtime_source_manifest": json.loads(
            metadata["benchmark_runtime_source_manifest_json"]
        ),
        "runtime_source_sha256": metadata["benchmark_runtime_source_sha256"],
    }
    _require(value == expected, "generic report has the wrong recorded Step 0 identity")
    compact = {
        key: copy.deepcopy(value[key])
        for key in (
            "benchmark",
            "step",
            "configuration_id",
            "identity_scope",
            "load_contract",
            "peak_load_definition",
            "active_stress_enabled",
            "pressure_enabled",
            "material_parameters",
            "activation_parameters",
            "pressure_parameters",
            "runtime_source_sha256",
        )
    }
    compact.update(
        case="step_0" + case,
        status="recorded",
        explicit_archive_identity_fields=True,
        recorded_archive_fields=sorted(post.BENCHMARK_ARCHIVE_FIELDS),
    )
    return selected, compact


def _reference_and_comparison(
    reference: Any, comparison: Any, reference_case: str
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    _require(isinstance(reference, dict) and isinstance(comparison, dict), "comparison evidence is malformed")
    _require(
        reference.get("case") == reference_case
        and reference.get("doi") == post.REFERENCE_DOI
        and reference.get("license") == post.REFERENCE_LICENSE
        and reference.get("published_archive_identity") == post.REFERENCE_ARCHIVE_IDENTITY,
        "external benchmark identity differs",
    )
    grid = _finite_array(reference.get("canonical_grid_s"), (101,), "canonical grid")
    _require(np.allclose(grid, post.CANONICAL_TIME_GRID, rtol=0.0, atol=1.0e-15), "canonical grid differs")
    selection = reference.get("selection", {})
    expected_names = [
        "monoventricular_nonblinded_{}_group_{}.pickle".format(reference_case, suffix)
        for suffix in post.REFERENCE_MANIFEST_SUFFIXES
    ]
    _require(
        selection.get("policy") == post.REFERENCE_SELECTION_POLICY
        and selection.get("selected_count") == 10
        and selection.get("selected_files") == expected_names
        and selection.get("upstream_figures_py_identity") == post.REFERENCE_FIGURES_PY_IDENTITY
        and selection.get("upstream_manifest_variable") == post.REFERENCE_MANIFEST_VARIABLES[reference_case],
        "reference selection differs from upstream figures.py",
    )
    team_files = reference.get("team_files")
    _require(isinstance(team_files, list) and len(team_files) == 10, "ten-team manifest is missing")
    for record, name, team in zip(team_files, expected_names, post.REFERENCE_MANIFEST_SUFFIXES):
        _require(
            isinstance(record, dict)
            and record.get("filename") == name
            and record.get("team") == team
            and _sha(record.get("sha256"))
            and isinstance(record.get("size_bytes"), int)
            and record["size_bytes"] > 0,
            "ten-team manifest identity differs",
        )
    means = reference.get("mean_curves_m", {})
    ours = comparison.get("ours_on_canonical_grid_m", {})
    red = comparison.get("red", {})
    _require(
        comparison.get("metric") == "relative discrepancy (benchmark paper Eq. 21)",
        "comparison metric differs",
    )
    expected_teams = set(post.REFERENCE_MANIFEST_SUFFIXES)
    for point in ("p0", "p1"):
        mean_curve = _finite_array(means.get(point), (101, 3), "mean " + point)
        ours_curve = _finite_array(ours.get(point), (101, 3), "CoupFE " + point)
        record = red.get(point, {})
        observed = record.get("ours")
        _require(
            isinstance(observed, (int, float))
            and not isinstance(observed, bool)
            and math.isclose(observed, post.red(ours_curve, mean_curve), rel_tol=1.0e-12, abs_tol=1.0e-14),
            "{} RED is not reproducible".format(point),
        )
        teams = record.get("teams", {})
        _require(
            set(teams) == expected_teams
            and all(isinstance(item, (int, float)) and not isinstance(item, bool) and math.isfinite(item) and item >= 0.0 for item in teams.values()),
            "{} team RED evidence differs".format(point),
        )
    return copy.deepcopy(reference), copy.deepcopy(comparison)


def _validate_generic(report: Dict[str, Any]) -> Dict[str, Any]:
    _require(
        set(report) == {"bounded_claim", "comparison", "reference", "result", "schema"}
        and report.get("schema") == GENERIC_REPORT_SCHEMA,
        "input is not an uncorrected generic-v2 report",
    )
    result = report.get("result", {})
    case = result.get("case")
    _require(case in {"A", "B"} and result.get("reference_case") == "step_0" + case, "input is not Step 0 Case A/B")
    _require(result.get("result_schema") == post.RESULT_SCHEMA, "result schema differs")
    _require(
        isinstance(result.get("filename"), str)
        and Path(result["filename"]).name == result["filename"]
        and _sha(result.get("sha256"))
        and isinstance(result.get("size_bytes"), int)
        and result["size_bytes"] > 0,
        "external NPZ identity is malformed",
    )
    source = _source_identity(result.get("source_identity"))
    configuration = result.get("configuration", {})
    selected, benchmark = _benchmark_identity(configuration.get("benchmark"), case)

    mesh = configuration.get("mesh", {})
    _require(
        mesh.get("topology") == "closed_multiblock_disk"
        and mesh.get("core_half_width") == 0.36
        and all(isinstance(mesh.get(key), int) and mesh[key] > 0 for key in ("n_t", "n_core", "n_radial", "elements", "nodes"))
        and mesh.get("degrees_of_freedom") == 3 * mesh["nodes"],
        "closed multiblock mesh identity differs",
    )
    exact = {
        "apex_offset_rad": 0.0,
        "dt_s": DT_S,
        "fiber_sampling": "gp_direct_rule",
        "fiber_sampling_option": "gp-direct",
        "flip_helix": True,
        "formulation": FORMULATION,
        "integrator": "generalized-alpha",
        "isotropic": False,
        "load_horizon_origin": "recorded",
        "load_horizon_s": 1.0,
        "mass_representation": "consistent_q1_hex8",
        "material_eta_pa_s": 100.0,
        "method_metadata_origin": "recorded",
        "model_parameters": MODEL_PARAMETERS,
        "nonlinear_solver": "petsc-snes-mpi",
        "parameter_variant": "benchmark_eta",
        "point_sampling": "hex8_reference_isoparametric",
        "t_end_s": T_END_S,
        "viscous_rate": "velocity_consistent_green_lagrange_at_alpha_f_stage",
        "viscous_term_active": True,
    }
    _require(all(configuration.get(key) == value for key, value in exact.items()), "selected local-pressure/GA configuration differs")
    tbar = configuration.get("tbar", {})
    _require(
        tbar.get("definition") == "laplace_presolved"
        and tbar.get("metadata_schema") == post.TBAR_METADATA_SCHEMA
        and _sha(tbar.get("source_sha256"))
        and _sha(tbar.get("metadata_sha256")),
        "Laplace tbar identity differs",
    )

    histories = result.get("retained_histories", {})
    times = _finite_array(histories.get("times_s"), (1001,), "result times")
    _require(np.allclose(times, np.linspace(0.0, 1.0, 1001), rtol=0.0, atol=1.0e-15), "result time grid differs")
    tau = _finite_array(histories.get("active_tension_pa"), (1001,), "active tension")
    pressure = _finite_array(histories.get("cavity_pressure_pa"), (1001,), "pressure")
    for point in ("u0_m", "u1_m"):
        _finite_array(histories.get(point), (1001, 3), point)
    load_times = times.copy()
    load_times[1:] -= GENERALIZED_ALPHA["alpha_f"] * DT_S
    expected_tau = tau_of_t(load_times, t_span=(0.0, 1.0)) if selected.active_stress_enabled else np.zeros_like(times)
    expected_pressure = p_of_t(load_times, t_span=(0.0, 1.0)) if selected.pressure_enabled else np.zeros_like(times)
    _require(
        np.allclose(tau, expected_tau, rtol=0.0, atol=1.0e-10)
        and np.allclose(pressure, expected_pressure, rtol=0.0, atol=1.0e-10),
        "recorded load histories differ from the selected Step 0 mode",
    )
    benchmark["load_history_audit"] = {
        "active_tension_minimum_pa": float(np.min(tau)),
        "active_tension_maximum_pa": float(np.max(tau)),
        "cavity_pressure_minimum_pa": float(np.min(pressure)),
        "cavity_pressure_maximum_pa": float(np.max(pressure)),
    }
    ga = configuration.get("generalized_alpha", {})
    _require(all(ga.get(key) == value for key, value in GENERALIZED_ALPHA.items()), "generalized-alpha parameters differ")
    ga_times = _finite_array(ga.get("load_evaluation_times_s"), (1001,), "GA load times")
    _require(np.allclose(ga_times, load_times, rtol=0.0, atol=1.0e-15), "generalized-alpha load staging differs")

    mpi = result.get("mpi_metadata", {})
    solver = result.get("solver_configuration", {})
    _require(
        mpi.get("world_size") == MPI_RANKS
        and mpi.get("contract") == MPI_CONTRACT[case]
        and mpi.get("implementation") == MPI_IMPLEMENTATION
        and mpi.get("linear_solver_profile") == LINEAR_SOLVER_PROFILE
        and solver.get("ranks") == MPI_RANKS
        and solver.get("linear_solver_profile") == LINEAR_SOLVER_PROFILE
        and solver.get("preconditioner") == PRECONDITIONER
        and solver.get("time_integrator") == "generalized-alpha"
        and solver.get("mass_representation") == "consistent_q1_hex8"
        and solver.get("local_pressure_law") == "log",
        "MPI solver contract differs",
    )
    diagnostics = result.get("nonlinear_step_diagnostics")
    _require(isinstance(diagnostics, list) and len(diagnostics) == 1000, "accepted-step diagnostics are incomplete")
    rejections = 0
    for index, record in enumerate(diagnostics, start=1):
        _require(
            record.get("time") == float(times[index])
            and record.get("dt") == DT_S
            and record.get("ranks") == MPI_RANKS
            and isinstance(record.get("snes_converged_reason"), int)
            and record["snes_converged_reason"] > 0
            and math.isfinite(record.get("final_residual_norm", math.nan))
            and record["final_residual_norm"] <= record.get("residual_acceptance_threshold", -1.0),
            "accepted-step diagnostic {} is invalid".format(index),
        )
        count = record.get("function_domain_rejections", 0)
        _require(isinstance(count, int) and count >= 0, "domain rejection count is invalid")
        rejections += count

    audit = result.get("pre_solve_audit", {})
    geometry = audit.get("geometry", {})
    _require(
        set(audit) == {"geometry", "pressure", "robin"}
        and all(audit[name].get("passed") is True and audit[name].get("failures") == [] for name in audit)
        and geometry.get("mesh_topology") == "closed_multiblock_disk"
        and geometry.get("require_closed") is True
        and geometry.get("nodes") == mesh["nodes"]
        and geometry.get("elements") == mesh["elements"]
        and geometry.get("intentional_traction_free_tip_faces") == 0
        and geometry.get("unexpected_unclassified_exterior_faces") == 0,
        "closed geometry/boundary pre-solve audits differ",
    )
    _require(audit["pressure"].get("pressure_surface_policy") == "closed_analytic_base_projection", "pressure audit policy differs")

    def summary(name: str, shape: list, positive: bool) -> Dict[str, Any]:
        record = result.get(name, {})
        values = [record.get(key) for key in ("minimum", "mean", "maximum")]
        _require(
            record.get("available") is True
            and record.get("shape") == shape
            and record.get("count") == int(np.prod(shape))
            and all(isinstance(item, (int, float)) and math.isfinite(item) for item in values)
            and values[0] <= values[1] <= values[2]
            and (not positive or values[0] > 0.0),
            "{} is invalid".format(name),
        )
        return copy.deepcopy(record)

    peak = result.get("peak", {})
    peak_index = peak.get("index")
    _require(
        peak.get("available") is True
        and isinstance(peak_index, int)
        and 0 <= peak_index < 1001
        and peak.get("time_s") == float(times[peak_index])
        and peak.get("u0_m") == histories["u0_m"][peak_index]
        and peak.get("u1_m") == histories["u1_m"][peak_index],
        "peak evidence differs from retained histories",
    )
    reference, comparison = _reference_and_comparison(
        report.get("reference"), report.get("comparison"), result["reference_case"]
    )
    compact_configuration = {
        key: copy.deepcopy(configuration[key]) for key in CONFIGURATION_FIELDS
    }
    compact_configuration.update(
        generalized_alpha=copy.deepcopy(GENERALIZED_ALPHA),
        sampling_points_m=copy.deepcopy(PHYSICAL_POINTS_M),
        boundary_conditions={
            "base": "full-vector Robin spring and dashpot",
            "cavity_pressure": (
                "Benchmark 1 Step 0B pressure-only load"
                if case == "B"
                else "identically-zero"
            ),
            "epicardium": "normal-only Robin spring and dashpot",
        },
        mpi={
            "contract": mpi["contract"],
            "implementation": mpi["implementation"],
            "linear_solver_profile": LINEAR_SOLVER_PROFILE,
            "preconditioner": PRECONDITIONER,
            "ranks": MPI_RANKS,
            "threads_per_rank": 1,
        },
    )
    return {
        "benchmark_identity": benchmark,
        "bounded_claim": report["bounded_claim"],
        "comparison": comparison,
        "normalized_run_log": copy.deepcopy(result.get("normalized_run_log")),
        "pre_solve_audit": copy.deepcopy(audit),
        "reference": reference,
        "result": {
            "case": case,
            "completion": {
                "completed_steps": 1000,
                "converged": True,
                "det_f_gauss_peak_summary": summary(
                    "det_f_gauss_peak_summary", [mesh["elements"], 8], True
                ),
                "element_pressure_peak_pa_summary": summary(
                    "element_pressure_peak_pa_summary", [mesh["elements"]], False
                ),
                "expected_steps": 1000,
                "function_domain_rejections": rejections,
                "peak": copy.deepcopy(peak),
            },
            "configuration": compact_configuration,
            "filename": result["filename"],
            "reference_case": result["reference_case"],
            "result_schema": result["result_schema"],
            "sha256": result["sha256"],
            "size_bytes": result["size_bytes"],
            "source_identity": source,
        },
    }


def _normalize_stdout(path: Path) -> Tuple[Dict[str, Any], Dict[str, Any], str]:
    raw_identity, payload = _identity(path)
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise PublicationError("raw stdout must be UTF-8") from error
    _require("\x00" not in text, "raw stdout contains a NUL byte")
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    if text and not text.endswith("\n"):
        text += "\n"
    normalized = text.encode("utf-8")
    return raw_identity, {
        "filename": path.name,
        "normalization": "UTF-8; CRLF/CR converted to LF; final LF added when nonempty",
        "sha256": hashlib.sha256(normalized).hexdigest(),
        "size_bytes": len(normalized),
    }, text


def _campaign_records(
    manifest_path: Path,
    attempt_name: str,
    stdout_path: Path,
    generic_identity: Mapping[str, Any],
    validated: Mapping[str, Any],
) -> Tuple[Dict[str, Any], float]:
    manifest, manifest_identity = _read_json(manifest_path, "campaign manifest")
    _require(manifest.get("schema") == CAMPAIGN_MANIFEST_SCHEMA, "campaign manifest schema differs")
    attempt = manifest.get("attempts", {}).get(attempt_name, {})
    raw_log, normalized_log, stdout_text = _normalize_stdout(stdout_path)
    result = validated["result"]
    expected = {
        "completed_steps": 1000,
        "expected_steps": 1000,
        "mpi_ranks": MPI_RANKS,
        "output_sha256": result["sha256"],
        "output_size_bytes": result["size_bytes"],
        "report_sha256": generic_identity["sha256"],
        "report_size_bytes": generic_identity["size_bytes"],
        "log_sha256": raw_log["sha256"],
        "log_size_bytes": raw_log["size_bytes"],
        "t_end_s": T_END_S,
        "status": "completed_and_validated",
    }
    _require(all(attempt.get(key) == value for key, value in expected.items()), "selected manifest attempt identity differs")
    _require(
        Path(attempt.get("output", "")).name == result["filename"]
        and Path(attempt.get("report", "")).name == generic_identity["filename"]
        and Path(attempt.get("log", "")).name == raw_log["filename"],
        "selected manifest attempt paths differ",
    )
    source = result["source_identity"]
    _require(
        manifest.get("application", {}).get("revision") == source["app"]["revision"]
        and manifest.get("application", {}).get("tree_state") == "clean"
        and manifest.get("core", {}).get("revision") == source["core"]["revision"]
        and manifest.get("core", {}).get("tree_state") == "clean",
        "campaign manifest source identity differs",
    )
    mesh = result["configuration"]["mesh"]
    attempt_mesh = attempt.get("mesh", {})
    _require(
        all(
            attempt_mesh.get(key) == mesh[key]
            for key in (
                "n_t",
                "n_core",
                "n_radial",
                "elements",
                "nodes",
                "degrees_of_freedom",
            )
        ),
        "selected manifest attempt mesh differs",
    )
    manifest_configuration = manifest.get("configuration", {})
    _require(
        all(
            manifest_configuration.get(key) == value
            for key, value in {
                "benchmark": 1,
                "case": result["case"],
                "step": 0,
                "dt_s": DT_S,
                "integrator": "generalized-alpha",
                "mass": "consistent",
                "formulation": "local-pressure",
                "mesh_topology": "closed-multiblock",
            }.items()
        ),
        "campaign manifest configuration differs",
    )
    _require(
        manifest.get("runtime", {}).get("selected_mpi_ranks") == MPI_RANKS
        and all(value == "1" for value in manifest.get("runtime", {}).get("thread_caps", {}).values()),
        "campaign manifest rank/thread contract differs",
    )
    _require(validated["normalized_run_log"] == normalized_log, "raw stdout does not match generic report")
    elapsed = attempt.get("elapsed_s")
    _require(isinstance(elapsed, (int, float)) and not isinstance(elapsed, bool) and elapsed > 0.0, "manifest elapsed time is invalid")
    _require(
        re.search(r"step\s+1000\s+t=1\.000s", stdout_text) is not None
        and re.search(r"^saved -> .*{}\s*$".format(re.escape(result["filename"])), stdout_text, re.MULTILINE) is not None
        and re.findall(r"^elapsed\s+([0-9]+(?:\.[0-9]+)?)s\s*$", stdout_text, re.MULTILINE) == [str(elapsed)],
        "raw stdout completion markers differ",
    )
    return {
        "generic_report": {
            **generic_identity,
            "retained_in_repository": False,
            "reason": "The forensic v2 report contains per-step solver diagnostics; only its external identity is retained.",
        },
        "manifest": {
            **manifest_identity,
            "retained_in_repository": False,
            "reason": "The campaign manifest contains machine-local paths; only its external identity is retained.",
        },
        "stdout": {
            **raw_log,
            "normalization": normalized_log["normalization"],
            "normalized_sha256": normalized_log["sha256"],
            "normalized_size_bytes": normalized_log["size_bytes"],
            "retained_in_repository": False,
            "reason": "The raw transcript contains machine-local paths; only its external identity is retained.",
        },
    }, float(elapsed)


def build_compact_report(
    generic_report_path: Path,
    manifest_path: Path,
    stdout_path: Path,
    *,
    manifest_attempt: str = "full",
) -> Dict[str, Any]:
    """Validate external evidence and return the six-key compact report."""
    generic, generic_identity = _read_json(Path(generic_report_path), "generic v2 report")
    validated = _validate_generic(generic)
    records, elapsed = _campaign_records(
        Path(manifest_path),
        manifest_attempt,
        Path(stdout_path),
        generic_identity,
        validated,
    )
    result = validated["result"]
    result["campaign_records"] = records
    result["completion"]["solve_elapsed_s"] = elapsed
    return {
        "benchmark_identity": validated["benchmark_identity"],
        "bounded_claim": validated["bounded_claim"],
        "comparison": validated["comparison"],
        "reference": validated["reference"],
        "result": result,
        "schema": COMPACT_REPORT_SCHEMA,
    }


def encode_compact_report(report: Mapping[str, Any]) -> bytes:
    """Return deterministic bytes and enforce the compact package boundary."""
    try:
        payload = (json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")
    except (TypeError, ValueError) as error:
        raise PublicationError("compact report is not finite JSON") from error
    _require(len(payload) < MAX_COMPACT_REPORT_BYTES, "compact report exceeds 100 kB")
    return payload


def publish_compact_report(
    generic_report_path: Path,
    manifest_path: Path,
    stdout_path: Path,
    output_path: Path,
    *,
    manifest_attempt: str = "full",
) -> Dict[str, Any]:
    report = build_compact_report(
        generic_report_path,
        manifest_path,
        stdout_path,
        manifest_attempt=manifest_attempt,
    )
    payload = encode_compact_report(report)
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name("." + destination.name + ".tmp")
    try:
        temporary.write_bytes(payload)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    return report


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generic-report", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--manifest-attempt", default="full")
    parser.add_argument("--stdout", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    arguments = _arguments()
    report = publish_compact_report(
        arguments.generic_report,
        arguments.manifest,
        arguments.stdout,
        arguments.output,
        manifest_attempt=arguments.manifest_attempt,
    )
    payload = arguments.output.read_bytes()
    print("{} ({} bytes, SHA-256 {})".format(arguments.output, len(payload), hashlib.sha256(payload).hexdigest()))
    print(
        "{} RED: p0={:.7g}, p1={:.7g}".format(
            report["benchmark_identity"]["case"],
            report["comparison"]["red"]["p0"]["ours"],
            report["comparison"]["red"]["p1"]["ours"],
        )
    )


if __name__ == "__main__":
    main()
