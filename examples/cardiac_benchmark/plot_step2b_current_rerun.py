"""Retained Step 2 Case B corrected-setup diagnostic renderer.

Renders one source-bound but provenance-incomplete development archive against
the official ten-team Step 2 Case B envelope and writes the compact metrics
report. The archive records application ``97d4474`` and Core ``454f73c``; it is
not a result from the current release tree. Publisher pickles are accepted only
after the reviewed ten-file hash manifest matches and are decoded by the
restricted NumPy unpickler shared with the release comparator. No machine-local
path is written into an output. The figure is a retained diagnostic, not a
reproduction, validation, or pass claim.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path

import numpy as np

try:  # package import
    from .compare_step2b_case_b import (
        ComparisonInputError,
        DEFAULT_HASH_MANIFEST,
        PUBLISHER_SELECTION,
        load_publisher_reference,
        require,
    )
except ImportError:  # direct example/script import
    from compare_step2b_case_b import (
        ComparisonInputError,
        DEFAULT_HASH_MANIFEST,
        PUBLISHER_SELECTION,
        load_publisher_reference,
        require,
    )

SCHEMA = "coupfe-cardiac-step2b-current-rerun-comparison-v1"
REPORT_NAME = "step2b_current_rerun_comparison.report.json"
FIGURE_TITLE = (
    "Benchmark 1, Step 2 Case B - retained corrected-setup diagnostic"
)
FIGURE_DESCRIPTION = (
    "Six line charts compare a retained CoupFE-Cardiac Step 2 Case B "
    "corrected-setup diagnostic (straight-wall geometry, physical frame, "
    "generalized-alpha Q1/P0 local pressure, tip_refine=6.0) with the official "
    "ten-team Step 2 Case B envelope from zero to one second. Rows show p0 and "
    "p1; columns show x, y, and z displacement in millimetres. The archive "
    "records application 97d4474 and Core 454f73c and is not a current-release "
    "result. The chart is a provenance-incomplete diagnostic, not a "
    "reproduction, validation, or pass claim."
)
VISIBLE_LABELS = (
    "closed t/core/radial 2x20x17 Hex8 Q1/P0 local pressure, tip_refine=6.0",
    "consistent-mass generalized-alpha",
    f"Report: {REPORT_NAME}",
    "benchmark DOI 10.5281/zenodo.14260459",
    "CC-BY-4.0",
)
PLATEAU = (0.32, 0.48)
REVIEWED_CORRECTED_NPZ_NAME = (
    "step2_caseb_local_pressure_ga_nt2_core20_rad17_tip6p0_t1p0.npz"
)
REVIEWED_CORRECTED_NPZ_SHA256 = (
    "63a8de59b7b8b9ab309896ff69989d6ff89f6dfe2532151605486ad67967dd41"
)
REVIEWED_CORRECTED_NPZ_SIZE_BYTES = 15396658
REVIEWED_LEGACY_REPORT_NAME = (
    "step2_case_b_std_kappa_2x20x17_dt0p001.report.json"
)
REVIEWED_LEGACY_REPORT_SHA256 = (
    "098e316daaea369a2a595cf43829d28597e53d2ff5a38cf32388e01c8dfa74aa"
)
REVIEWED_LEGACY_REPORT_SIZE_BYTES = 214309
REVIEWED_CORRECTED_RUN_SCALARS = {
    "result_schema": "coupfe-cardiac-result-v1",
    "driver": "examples/cardiac_benchmark/run_mpi.py",
    "completed_steps": 1000,
    "expected_steps": 1000,
    "converged": True,
    "app_revision": "97d447491498830e19a0b791a35cff4f1d13694e",
    "app_tree_state": "clean",
    "app_source_kind": "git-checkout",
    "core_revision": "454f73ce2de284262b214a2b37bd676c6aca3c0a",
    "core_tree_state": "clean",
    "core_source_kind": "git-checkout",
    "core_source_url": "https://github.com/tengzhang48/CoupFE.git",
    "mpi_enabled": True,
    "mpi_ranks": 8,
    "mpi_world_size": 8,
    "mpi_implementation": (
        "cardiac-owned-distributed-closed-local-pressure-"
        "generalized-alpha-step0"
    ),
    "case": "B",
    "benchmark_step": 2,
    "benchmark_configuration_id": (
        "benchmark-1-step-2-case-B-active-stress-plus-pressure"
    ),
    "benchmark_load_contract": "active-stress-plus-pressure",
    "benchmark_active_stress_enabled": True,
    "benchmark_pressure_enabled": True,
    "dt": 0.001,
    "t_end": 1.0,
    "load_horizon": 1.0,
    "integrator": "generalized-alpha",
    "generalized_alpha_alpha_m": 0.2,
    "generalized_alpha_alpha_f": 0.4,
    "generalized_alpha_gamma": 0.7,
    "generalized_alpha_beta": 0.36,
    "generalized_alpha_stage_contract": "simula-source-matched-v1",
    "formulation": "hex8_local_pressure_p0_condensed_logj",
    "material_kappa_pa": 0.0,
    "local_pressure_bulk_modulus_pa": 1.0e6,
    "local_pressure_volume_law": "linear-reference-volume-mean-log-j-v1",
    "mass_representation": "consistent_q1_hex8",
    "material_kernel_formulation": "standard",
    "material_model_id": (
        "holzapfel-ogden-smooth-switch-complete-energy-derivative-v1"
    ),
    "density": 1000.0,
    "material_eta_pa_s": 100.0,
    "viscous_term_active": True,
    "viscous_rate": "velocity_consistent_green_lagrange_at_alpha_f_stage",
    "mesh_topology": "closed_multiblock_disk",
    "n_t": 2,
    "n_core": 20,
    "n_radial": 17,
    "core_half_width": 0.36,
    "tip_refine": 6.0,
    "apex_offset": 0.0,
    "perturb": 0.0,
    "flip_helix": True,
    "isotropic": False,
    "fiber_sampling": "gp_direct_rule",
    "fiber_sampling_option": "gp-direct",
    "fiber_direction_reconstruction": "toolkit-physical-coordinate-u-v-v1",
    "tbar_definition": "laplace_presolved",
    "point_sampling": "hex8_reference_isoparametric",
}
PUBLISHER_MATERIAL_ATTRIBUTION = {
    "source_title": "A software benchmark for cardiac elastodynamics",
    "creators": ["Arostica Barrera, R.A.", "Bertoglio, Cristobal"],
    "source_doi": "https://doi.org/10.5281/zenodo.14260459",
    "license": "CC-BY-4.0",
    "license_url": "https://creativecommons.org/licenses/by/4.0/",
    "transformation_notice": (
        "CoupFE-Cardiac selected ten Step 2 Case B team curves, linearly "
        "interpolated each curve to the corrected-run time grid, and computed "
        "the team mean, minimum/maximum envelope, containment fractions, and "
        "plateau summaries used in this report and SVG. These are "
        "modified/derived data; raw publisher pickle files are not "
        "redistributed."
    ),
}


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _archive_scalar(archive, key: str, path: Path):
    require(key in archive, f"{path.name} is missing {key!r}")
    try:
        value = np.asarray(archive[key])
    except (TypeError, ValueError) as error:
        raise ComparisonInputError(
            f"cannot read {key!r} from {path.name}"
        ) from error
    require(value.shape == (), f"{path.name} field {key!r} must be scalar")
    return value.item()


def _scalar_matches(actual, expected) -> bool:
    if isinstance(expected, bool):
        return isinstance(actual, bool) and actual is expected
    if isinstance(expected, int):
        return isinstance(actual, int) and not isinstance(actual, bool) and actual == expected
    if isinstance(expected, float):
        return (
            isinstance(actual, (int, float))
            and not isinstance(actual, bool)
            and np.isfinite(float(actual))
            and float(actual) == expected
        )
    return isinstance(actual, str) and actual == expected


def _archive_array(archive, key: str, path: Path, shape) -> np.ndarray:
    require(key in archive, f"{path.name} is missing {key!r}")
    try:
        value = np.asarray(archive[key], dtype=float)
    except (TypeError, ValueError) as error:
        raise ComparisonInputError(
            f"{path.name} field {key!r} must be numeric"
        ) from error
    require(value.shape == shape, f"{path.name} field {key!r} has the wrong shape")
    require(np.all(np.isfinite(value)), f"{path.name} field {key!r} must be finite")
    return value.copy()


def load_reviewed_corrected_run(path: Path) -> dict:
    """Load only the exact retained archive and verify its plotted contract."""
    path = Path(path).expanduser().resolve()
    require(path.is_file(), f"reviewed corrected archive does not exist: {path}")
    require(
        path.name == REVIEWED_CORRECTED_NPZ_NAME,
        "corrected archive does not have the reviewed filename",
    )
    try:
        payload = path.read_bytes()
    except OSError as error:
        raise ComparisonInputError(
            f"cannot read reviewed corrected archive {path.name}: {error}"
        ) from error
    identity = {
        "name": path.name,
        "sha256": _sha256_bytes(payload),
        "size_bytes": len(payload),
    }
    require(
        identity["sha256"] == REVIEWED_CORRECTED_NPZ_SHA256
        and identity["size_bytes"] == REVIEWED_CORRECTED_NPZ_SIZE_BYTES,
        "corrected archive is not the exact reviewed NPZ",
    )
    try:
        context = np.load(io.BytesIO(payload), allow_pickle=False)
    except (OSError, ValueError) as error:
        raise ComparisonInputError(
            f"cannot load reviewed corrected archive {path.name}: {error}"
        ) from error
    try:
        with context as archive:
            observed_scalars = {}
            for key, expected in REVIEWED_CORRECTED_RUN_SCALARS.items():
                actual = _archive_scalar(archive, key, path)
                require(
                    _scalar_matches(actual, expected),
                    f"reviewed corrected archive field {key!r} is {actual!r}, "
                    f"expected {expected!r}",
                )
                observed_scalars[key] = actual

            times = _archive_array(archive, "times", path, (1001,))
            expected_times = np.arange(1001, dtype=float) * 0.001
            require(
                np.allclose(times, expected_times, rtol=0.0, atol=5.0e-13),
                "reviewed corrected archive does not use the exact 1 ms time grid",
            )
            load_times = _archive_array(
                archive, "load_evaluation_times_s", path, (1001,)
            )
            expected_load_times = expected_times.copy()
            expected_load_times[1:] -= 0.4 * 0.001
            require(
                np.allclose(
                    load_times, expected_load_times, rtol=0.0, atol=5.0e-13
                ),
                "reviewed corrected archive has the wrong generalized-alpha "
                "load stage",
            )
            ours = {
                "p0": _archive_array(archive, "u0", path, (1001, 3)),
                "p1": _archive_array(archive, "u1", path, (1001, 3)),
            }
            expected_points = {
                "p0": np.array([0.025, 0.030, 0.0]),
                "p1": np.array([0.000, 0.030, 0.0]),
            }
            for label, expected in expected_points.items():
                point = _archive_array(archive, label, path, (3,))
                require(
                    np.allclose(point, expected, rtol=0.0, atol=1.0e-14),
                    f"reviewed corrected archive has the wrong {label} landmark",
                )
            nodes = _archive_array(archive, "nodes", path, (5403, 3))
            require(
                np.all(np.isfinite(nodes)),
                "reviewed corrected archive has invalid mesh nodes",
            )
            require(
                "elems" in archive
                and np.asarray(archive["elems"]).shape == (3520, 8)
                and np.issubdtype(np.asarray(archive["elems"]).dtype, np.integer),
                "reviewed corrected archive has the wrong Hex8 mesh connectivity",
            )
    except ComparisonInputError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ComparisonInputError(
            f"malformed reviewed corrected archive {path.name}: {error}"
        ) from error

    return {
        "times": times,
        "ours": ours,
        "identity": identity,
        "app_revision": observed_scalars["app_revision"],
        "core_revision": observed_scalars["core_revision"],
    }


def load_reviewed_legacy_report(path: Path) -> dict:
    """Load the exact retained pre-correction report used by this diagnostic."""
    path = Path(path).expanduser().resolve()
    require(path.is_file(), f"reviewed legacy report does not exist: {path}")
    require(
        path.name == REVIEWED_LEGACY_REPORT_NAME,
        "legacy report does not have the reviewed filename",
    )
    try:
        payload = path.read_bytes()
    except OSError as error:
        raise ComparisonInputError(
            f"cannot read reviewed legacy report {path.name}: {error}"
        ) from error
    identity = {
        "name": path.name,
        "sha256": _sha256_bytes(payload),
        "size_bytes": len(payload),
    }
    require(
        identity["sha256"] == REVIEWED_LEGACY_REPORT_SHA256
        and identity["size_bytes"] == REVIEWED_LEGACY_REPORT_SIZE_BYTES,
        "legacy report is not the exact reviewed report",
    )

    def reject_constant(value):
        raise ValueError(f"non-finite JSON constant {value}")

    try:
        report = json.loads(
            payload.decode("utf-8"), parse_constant=reject_constant
        )
    except (UnicodeError, ValueError, json.JSONDecodeError) as error:
        raise ComparisonInputError(
            f"cannot parse reviewed legacy report {path.name}: {error}"
        ) from error
    require(
        isinstance(report, dict)
        and report.get("schema")
        == "coupfe-cardiac-step2b-publisher-comparison-v1",
        "reviewed legacy report has the wrong report schema",
    )
    curves = report.get("comparison_curves")
    require(
        isinstance(curves, dict)
        and set(curves)
        == {"meaning", "points", "sample_count", "team_count", "time_s"}
        and curves["sample_count"] == 101
        and curves["team_count"] == 10,
        "reviewed legacy report has the wrong comparison_curves schema",
    )
    times = np.asarray(curves["time_s"], dtype=float)
    require(
        times.shape == (101,)
        and np.all(np.isfinite(times))
        and np.allclose(
            times, np.arange(101, dtype=float) * 0.01, rtol=0.0, atol=5.0e-13
        ),
        "reviewed legacy report has the wrong comparison time grid",
    )
    points = curves["points"]
    require(
        isinstance(points, dict) and set(points) == {"p0", "p1"},
        "reviewed legacy report must contain exactly p0 and p1 curves",
    )
    expected_curve_fields = {
        "coupfe_m",
        "publisher_max_m",
        "publisher_mean_m",
        "publisher_min_m",
        "publisher_simula_m",
    }
    for label in ("p0", "p1"):
        require(
            isinstance(points[label], dict)
            and set(points[label]) == {"ux", "uy", "uz"},
            f"reviewed legacy report has the wrong {label} component schema",
        )
        for component in ("ux", "uy", "uz"):
            record = points[label][component]
            require(
                isinstance(record, dict) and set(record) == expected_curve_fields,
                "reviewed legacy report has the wrong "
                f"{label} {component} curve schema",
            )
            for field in expected_curve_fields:
                values = np.asarray(record[field], dtype=float)
                require(
                    values.shape == (101,) and np.all(np.isfinite(values)),
                    "reviewed legacy report has invalid "
                    f"{label} {component} {field} values",
                )
    return {"report": report, "identity": identity}


def _load_teams(
    teams_dir: Path, times: np.ndarray, hash_manifest: Path
) -> dict[str, np.ndarray]:
    reference = load_publisher_reference(teams_dir, hash_manifest)
    publisher_teams = np.asarray(reference["teams_m"], dtype=float)
    publisher_team_times = np.asarray(reference["team_times_s"], dtype=float)
    diagnostic_order = sorted(
        range(len(PUBLISHER_SELECTION)),
        key=lambda index: PUBLISHER_SELECTION[index][2],
    )
    teams = {}
    for point_index, label in enumerate(("p0", "p1")):
        team_curves = []
        for team_index in diagnostic_order:
            # Preserve this diagnostic renderer's original per-team
            # interpolation and filename-sorted reduction order. The release
            # comparator separately follows the publisher's common-grid and
            # publisher-selection order.
            team_times = publisher_team_times[team_index]
            team_curves.append(
                np.column_stack(
                    [
                        np.interp(
                            times,
                            team_times,
                            publisher_teams[team_index, point_index, component],
                        )
                        for component in range(3)
                    ]
                )
            )
        teams[label] = np.stack(team_curves)
    return teams


def compute_metrics(
    corrected_npz: Path,
    legacy_report: Path,
    teams_dir: Path,
    hash_manifest: Path = DEFAULT_HASH_MANIFEST,
) -> dict:
    corrected = load_reviewed_corrected_run(corrected_npz)
    times = corrected["times"]
    ours = corrected["ours"]
    legacy_input = load_reviewed_legacy_report(legacy_report)
    report = legacy_input["report"]
    curves = report["comparison_curves"]
    legacy_times = np.asarray(curves["time_s"], dtype=float)
    legacy = {
        lab: np.column_stack(
            [np.interp(times, legacy_times, np.asarray(curves["points"][lab][c]["coupfe_m"])) for c in ("ux", "uy", "uz")]
        )
        for lab in ("p0", "p1")
    }
    teams = _load_teams(teams_dir, times, hash_manifest)

    metrics: dict = {"versus_ten_team_mean": {}, "envelope_inside_fraction": {}, "plateau": {}}
    for lab in ("p0", "p1"):
        band = teams[lab] * 1.0e3
        mean = band.mean(axis=0)
        lo, hi = band.min(axis=0), band.max(axis=0)
        own = ours[lab] * 1.0e3
        rel = float(np.sqrt(np.sum((own - mean) ** 2) / np.sum(mean**2)))
        inside = (own >= lo) & (own <= hi)
        metrics["versus_ten_team_mean"][lab] = {"relative_l2": rel}
        metrics["envelope_inside_fraction"][lab] = {
            c: float(inside[:, i].mean()) for i, c in enumerate(("x", "y", "z"))
        }
        lo_t, hi_t = PLATEAU
        mask = (times >= lo_t) & (times <= hi_t)
        metrics["plateau"][lab] = {
            "window_s": [lo_t, hi_t],
            "corrected_z_range_mm": [
                float(own[mask, 2].min()),
                float(own[mask, 2].max()),
            ],
            "team_z_band_mm": [
                float(lo[mask, 2].min()),
                float(hi[mask, 2].max()),
            ],
            "legacy_z_range_mm": [
                float((legacy[lab][mask, 2] * 1.0e3).min()),
                float((legacy[lab][mask, 2] * 1.0e3).max()),
            ],
            "legacy_sign_matches_teams": bool(
                np.sign(mean[mask, 2]).mean()
                == np.sign((legacy[lab][mask, 2] * 1.0e3)).mean()
            ),
            "corrected_sign_matches_teams": bool(
                np.sign(mean[mask, 2]).mean() == np.sign(own[mask, 2]).mean()
            ),
        }
    metrics["inputs"] = {
        "corrected_npz": {
            "name": corrected["identity"]["name"],
            "sha256": corrected["identity"]["sha256"],
        },
        "legacy_report": {
            "name": legacy_input["identity"]["name"],
            "sha256": legacy_input["identity"]["sha256"],
        },
        "ten_team_dataset_doi": "10.5281/zenodo.14260459",
    }
    metrics["publisher_material_attribution"] = PUBLISHER_MATERIAL_ATTRIBUTION
    metrics["schema"] = SCHEMA
    metrics["claim_boundary"] = (
        "Comparison of one retained corrected-setup Step 2 Case B diagnostic "
        "with the pre-correction development run and the official ten-team "
        "envelope. Its compact record is provenance-incomplete and it is not "
        "a current-release result, reproduction, validation, or pass claim."
    )
    return metrics


def _accessible_svg(fig) -> bytes:
    import io
    import re

    buffer = io.StringIO()
    fig.savefig(buffer, format="svg")
    text = buffer.getvalue()
    text = re.sub(r"^\s*<\?xml[^?]*\?>\s*", "", text, count=1)
    text = re.sub(r"<!DOCTYPE[^>]*(\[[^]]*\])?>", "", text, count=1, flags=re.DOTALL)
    text = re.sub(r"<metadata>.*?</metadata>\s*", "", text, flags=re.DOTALL)
    title = (
        f'<title id="figure-title">{FIGURE_TITLE}</title>\n'
        f'<desc id="figure-description">{FIGURE_DESCRIPTION}</desc>'
    )
    match = re.search(r"<svg[^>]*>", text)
    if match is None:
        raise RuntimeError("matplotlib SVG root not found")
    opening = match.group(0)
    if "role=" not in opening:
        opening = opening[:-1] + ' role="img" aria-labelledby="figure-title figure-description">'
    text = text[: match.start()] + opening + "\n" + title + "\n" + text[match.end():]
    return text.encode("utf-8")


def render_figure(
    metrics: dict,
    corrected_npz: Path,
    legacy_report: Path,
    teams_dir: Path,
    hash_manifest: Path = DEFAULT_HASH_MANIFEST,
) -> bytes:
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams["svg.fonttype"] = "none"
    matplotlib.rcParams["svg.hashsalt"] = (
        "coupfe-cardiac-step2b-corrected-diagnostic-v1"
    )
    from matplotlib import pyplot as plt

    corrected = load_reviewed_corrected_run(corrected_npz)
    times = corrected["times"]
    ours = {
        "p0": corrected["ours"]["p0"] * 1.0e3,
        "p1": corrected["ours"]["p1"] * 1.0e3,
    }
    teams = _load_teams(teams_dir, times, hash_manifest)

    fig, axes = plt.subplots(2, 3, figsize=(15, 7.5), sharex=True)
    for row, lab in enumerate(("p0", "p1")):
        band = teams[lab] * 1.0e3
        lo, hi = band.min(axis=0), band.max(axis=0)
        for col, comp in enumerate(("x", "y", "z")):
            axis = axes[row, col]
            axis.fill_between(
                times, lo[:, col], hi[:, col], color="0.75", alpha=0.5,
                label="ten-team envelope" if (row, col) == (0, 0) else None,
            )
            axis.plot(
                times, ours[lab][:, col], color="tab:red", lw=1.4, ls="--",
                label="retained corrected-setup diagnostic"
                if (row, col) == (0, 0)
                else None,
            )
            axis.set_ylabel(f"{lab}  u_{comp} (mm)")
            axis.grid(alpha=0.3)
            axis.set_xlim(0.0, 1.0)
            if row == 0:
                axis.set_title(f"u_{comp}")
            if row == 1:
                axis.set_xlabel("t (s)")
    axes[0, 0].legend(loc="upper right", fontsize=9)
    app_revision = str(corrected["app_revision"])[:7]
    core_revision = str(corrected["core_revision"])[:7]
    fig.suptitle(FIGURE_TITLE)
    fig.text(
        0.01,
        0.005,
        "  |  ".join(VISIBLE_LABELS) + f"  |  app {app_revision}  |  Core {core_revision}",
        fontsize=7,
    )
    fig.tight_layout(rect=(0.0, 0.02, 1.0, 1.0))
    return _accessible_svg(fig)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corrected-npz", type=Path, required=True)
    parser.add_argument("--legacy-report", type=Path, required=True)
    parser.add_argument("--teams-dir", type=Path, required=True)
    parser.add_argument(
        "--hash-manifest",
        type=Path,
        default=DEFAULT_HASH_MANIFEST,
        help="reviewed exact ten-file publisher hash manifest",
    )
    parser.add_argument("--out-metrics", type=Path, required=True)
    parser.add_argument("--out-figure", type=Path, required=True)
    args = parser.parse_args(argv)

    if args.out_metrics.name != REPORT_NAME:
        parser.error(f"--out-metrics basename must be {REPORT_NAME}")

    metrics = compute_metrics(
        args.corrected_npz,
        args.legacy_report,
        args.teams_dir,
        args.hash_manifest,
    )
    metrics_bytes = (
        json.dumps(metrics, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    args.out_metrics.write_bytes(metrics_bytes)
    figure = render_figure(
        metrics,
        args.corrected_npz,
        args.legacy_report,
        args.teams_dir,
        args.hash_manifest,
    )
    args.out_figure.write_bytes(figure)
    print(
        json.dumps(
            {
                "metrics": {"path": str(args.out_metrics), "sha256": _sha256_bytes(metrics_bytes)},
                "figure": {"path": str(args.out_figure), "sha256": _sha256_bytes(figure), "size_bytes": len(figure)},
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
