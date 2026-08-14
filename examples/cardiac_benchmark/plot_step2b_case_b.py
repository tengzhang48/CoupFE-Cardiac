"""Render the reviewed full-cycle Step-2 Case-B publisher comparison.

The checked-in JSON report contains the derived CoupFE, ten-team mean/range,
and named Simula curves on the publisher's 101-point grid.  The raw CoupFE NPZ
and publisher pickle files remain external.  By default this script writes
``docs/figures/step2_case_b_comparison.svg``; a ``.png`` output may be selected
for local visual QA.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Dict, List, Mapping, Optional, Sequence
from xml.sax.saxutils import escape


REPORT_SCHEMA = "coupfe-cardiac-step2b-publisher-comparison-v1"
REPORT_NAME = "step2_case_b_std_kappa_2x20x17_dt0p001.report.json"
REPORT_SHA256 = "098e316daaea369a2a595cf43829d28597e53d2ff5a38cf32388e01c8dfa74aa"
RESULT_NAME = "mpi4_ga_full.npz"
RESULT_SHA256 = "23312a5e0147544eb9a4e6de004a166ada2722b70d3d39742f93aacd8a0fa0e6"
RESULT_SIZE_BYTES = 3_835_766
APP_REVISION = "e9b7d9084b24f7098170a221061eb159d0b090c1"
CORE_REVISION = "e2f42ed5772850a0a23a2ce434f430c287eae5c8"
RUNTIME_SOURCE_SHA256 = (
    "6b96395761dd3203f0e9ffab90a77d6389dca13cdad43490a1deac95073480f1"
)
REFERENCE_DOI = "10.5281/zenodo.14260459"
REFERENCE_LICENSE = "CC-BY-4.0"
POINTS = ("p0", "p1")
COMPONENTS = ("ux", "uy", "uz")
SERIES = (
    "coupfe_m",
    "publisher_mean_m",
    "publisher_min_m",
    "publisher_max_m",
    "publisher_simula_m",
)
CURVE_POINTS = 101
CANONICAL_PYTHON_VERSION = (3, 12, 3)
CANONICAL_MATPLOTLIB_VERSION = "3.11.1"
CANONICAL_FREETYPE_VERSION = "2.14.3"
CANONICAL_FONT_SHA256 = (
    "3fdf69cabf06049ea70a00b5919340e2ce1e6d02b0cc3c4b44fb6801bd1e0d22"
)


def _reject_nonstandard_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant {value!r}")


def _mapping(value: Any, description: str) -> Mapping[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{description} must be an object")
    return value


def _finite(value: Any, description: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{description} must be numeric")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{description} must be finite")
    return number


def _curve(value: Any, description: str) -> List[float]:
    if not isinstance(value, list) or len(value) != CURVE_POINTS:
        raise ValueError(f"{description} must contain {CURVE_POINTS} samples")
    return [
        _finite(sample, f"{description}[{index}]")
        for index, sample in enumerate(value)
    ]


def _validate_report(path: Path) -> Dict[str, Any]:
    payload = path.read_bytes()
    digest = hashlib.sha256(payload).hexdigest()
    if path.name != REPORT_NAME or digest != REPORT_SHA256:
        raise ValueError(
            f"report identity mismatch: {path.name} SHA-256 {digest}; "
            f"expected {REPORT_NAME} SHA-256 {REPORT_SHA256}"
        )
    report = json.loads(
        payload.decode("utf-8"), parse_constant=_reject_nonstandard_constant
    )
    root = _mapping(report, str(path))
    if root.get("schema") != REPORT_SCHEMA:
        raise ValueError(f"{path} does not use schema {REPORT_SCHEMA!r}")

    scope = _mapping(root.get("scope"), "scope")
    expected_scope = {
        "benchmark": 1,
        "step": 2,
        "case": "B",
        "configuration_id": (
            "benchmark-1-step-2-case-B-active-stress-plus-pressure"
        ),
        "loads": "active stress plus ventricular pressure",
        "comparison_grid": "101 samples, 0.00--1.00 s at 0.01 s",
        "displacement_unit": "m",
    }
    if scope != expected_scope:
        raise ValueError("report scope is not the reviewed Step-2 Case-B scope")

    inputs = _mapping(root.get("inputs"), "inputs")
    result = _mapping(inputs.get("coupfe_run"), "inputs.coupfe_run")
    expected_result = {
        "filename": RESULT_NAME,
        "sha256": RESULT_SHA256,
        "size_bytes": RESULT_SIZE_BYTES,
    }
    if result != expected_result:
        raise ValueError("report does not identify the reviewed full-cycle NPZ")
    contract = _mapping(
        inputs.get("coupfe_run_contract"), "inputs.coupfe_run_contract"
    )
    expected_contract = {
        "app_revision": APP_REVISION,
        "app_tree_state": "dirty",
        "core_revision": CORE_REVISION,
        "formulation": "hex8_standard_pointwise_kappa",
        "mpi_ranks": 4,
        "profile": "paper-source-matched-full-cycle",
        "runtime_source_sha256": RUNTIME_SOURCE_SHA256,
    }
    if contract != expected_contract:
        raise ValueError("report run contract is not the reviewed full-cycle run")

    plumbing = _mapping(root.get("plumbing"), "plumbing")
    reproduction = _mapping(root.get("reproduction"), "reproduction")
    if plumbing.get("status") != "passed":
        raise ValueError("report plumbing status is not passed")
    if reproduction.get("status") != "quantified-no-paper-acceptance-threshold":
        raise ValueError("report has an unexpected reproduction status")

    errors = _mapping(reproduction.get("trajectory_errors"), "trajectory_errors")
    full = _mapping(errors.get("full_history"), "full_history")
    aggregate = _mapping(full.get("all_components"), "full_history.all_components")
    points_error = _mapping(full.get("points"), "full_history.points")
    paper_red = _mapping(
        reproduction.get("paper_relative_discrepancy"),
        "paper_relative_discrepancy",
    )
    paper_red_points = _mapping(
        paper_red.get("points"), "paper_relative_discrepancy.points"
    )
    metrics = {
        "aggregate_relative_l2": _finite(
            aggregate.get("relative_l2"), "aggregate relative L2"
        ),
        "aggregate_rmse_mm": _finite(aggregate.get("rmse_mm"), "aggregate RMSE"),
        "aggregate_max_mm": 1000.0
        * _finite(aggregate.get("max_abs_error_m"), "aggregate max error"),
    }
    for point in POINTS:
        point_error = _mapping(points_error.get(point), f"full_history.points.{point}")
        vector = _mapping(point_error.get("vector"), f"{point}.vector")
        metrics[f"{point}_relative_l2"] = _finite(
            vector.get("relative_l2_vector"), f"{point} relative L2"
        )
        metrics[f"{point}_paper_red"] = _finite(
            paper_red_points.get(point), f"{point} paper RED"
        )

    curves = _mapping(root.get("comparison_curves"), "comparison_curves")
    if curves.get("sample_count") != CURVE_POINTS or curves.get("team_count") != 10:
        raise ValueError("comparison curve counts are not 101 samples and ten teams")
    times = _curve(curves.get("time_s"), "comparison_curves.time_s")
    for index, value in enumerate(times):
        expected = index / 100.0
        if not math.isclose(value, expected, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(f"time sample {index} is {value}; expected {expected}")

    points = _mapping(curves.get("points"), "comparison_curves.points")
    if set(points) != set(POINTS):
        raise ValueError("comparison curves must contain exactly p0 and p1")
    parsed = {}
    for point in POINTS:
        components = _mapping(points[point], f"comparison_curves.points.{point}")
        if set(components) != set(COMPONENTS):
            raise ValueError(f"{point} must contain ux, uy, and uz")
        parsed[point] = {}
        for component in COMPONENTS:
            record = _mapping(components[component], f"{point}.{component}")
            if set(record) != set(SERIES):
                raise ValueError(f"{point}.{component} has unexpected curve series")
            parsed_record = {
                name: _curve(record[name], f"{point}.{component}.{name}")
                for name in SERIES
            }
            for index in range(CURVE_POINTS):
                low = parsed_record["publisher_min_m"][index]
                mean = parsed_record["publisher_mean_m"][index]
                high = parsed_record["publisher_max_m"][index]
                simula = parsed_record["publisher_simula_m"][index]
                if not low <= mean <= high or not low <= simula <= high:
                    raise ValueError(
                        f"{point}.{component} publisher envelope fails at {index}"
                    )
            if max(abs(v) for v in parsed_record["coupfe_m"]) >= 0.20:
                raise ValueError(f"{point}.{component} CoupFE curve is implausible")
            parsed[point][component] = parsed_record
    p1_z = parsed["p1"]["uz"]
    plateau_indices = range(32, 49)
    if not all(p1_z["coupfe_m"][index] < 0.0 for index in plateau_indices):
        raise ValueError("p1 uz does not retain the reviewed negative plateau")
    if not all(
        p1_z["publisher_min_m"][index] > 0.0 for index in plateau_indices
    ):
        raise ValueError("the official p1 uz envelope is not positive on the plateau")
    return {"times": times, "points": parsed, "metrics": metrics}


def _require_canonical_renderer() -> None:
    try:
        import matplotlib as mpl
        from matplotlib import font_manager, ft2font
    except ImportError as error:
        raise RuntimeError("Matplotlib is required") from error
    font_path = Path(font_manager.findfont("DejaVu Sans", fallback_to_default=False))
    actual = {
        "Python": tuple(sys.version_info[:3]),
        "Matplotlib": mpl.__version__,
        "FreeType": ft2font.__freetype_version__,
        "DejaVu Sans SHA-256": hashlib.sha256(font_path.read_bytes()).hexdigest(),
    }
    expected = {
        "Python": CANONICAL_PYTHON_VERSION,
        "Matplotlib": CANONICAL_MATPLOTLIB_VERSION,
        "FreeType": CANONICAL_FREETYPE_VERSION,
        "DejaVu Sans SHA-256": CANONICAL_FONT_SHA256,
    }
    mismatches = [
        f"{name}: found {actual[name]!r}, expected {expected[name]!r}"
        for name in expected
        if actual[name] != expected[name]
    ]
    if mismatches:
        raise RuntimeError("canonical renderer mismatch; " + "; ".join(mismatches))


def _accessible_svg(svg: str, title: str, description: str) -> str:
    svg = re.sub(r"<!DOCTYPE svg.*?>\s*", "", svg, count=1, flags=re.DOTALL)
    svg = re.sub(r"\s*<metadata>.*?</metadata>\s*", "\n", svg, flags=re.DOTALL)
    svg = re.sub(r"\s*<title>.*?</title>\s*", "\n", svg, flags=re.DOTALL)
    marker = "<svg "
    index = svg.find(marker)
    if index < 0:
        raise RuntimeError("Matplotlib output does not contain an SVG root")
    end = svg.find(">", index)
    if end < 0:
        raise RuntimeError("Matplotlib output has an unterminated SVG root")
    root = svg[index:end]
    root += ' role="img" aria-labelledby="figure-title figure-description"'
    children = (
        '\n <title id="figure-title">'
        + escape(title)
        + '</title>\n <desc id="figure-description">'
        + escape(description)
        + "</desc>"
    )
    result = svg[:index] + root + ">" + children + svg[end + 1 :]
    forbidden = ("<!DOCTYPE", "<metadata", "rdf:", "dc:", "rdf:resource")
    if any(token in result for token in forbidden):
        raise RuntimeError("SVG sanitization left external or RDF metadata")
    return "\n".join(line.rstrip() for line in result.splitlines()) + "\n"


def _render(report: Mapping[str, Any], output: Path) -> None:
    try:
        import matplotlib as mpl

        mpl.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError("Matplotlib is required") from error

    metrics = report["metrics"]
    title = "Benchmark 1, Step 2 Case B — displacement comparison"
    description = (
        "Six line charts compare the full-cycle CoupFE-Cardiac result with the "
        "official ten-team range, all-team mean, and named Simula curve from zero "
        "to one second. Rows show p0 and p1; columns show x, y, and z displacement "
        "in millimetres. Full-history relative L2 error against the all-team mean "
        f"is {100.0 * metrics['aggregate_relative_l2']:.1f}% overall, "
        f"{100.0 * metrics['p0_relative_l2']:.1f}% at p0, and "
        f"{100.0 * metrics['p1_relative_l2']:.1f}% at p1; benchmark-paper Eq. 21 "
        f"RED is {100.0 * metrics['p0_paper_red']:.1f}% at p0 and "
        f"{100.0 * metrics['p1_paper_red']:.1f}% at p1. The p1 z plateau has "
        "the opposite sign from all ten published curves. Source: "
        f"{REPORT_NAME}, application {APP_REVISION[:7]} with an exact dirty-tree "
        f"runtime manifest, Core {CORE_REVISION[:7]}, benchmark DOI {REFERENCE_DOI}. "
        "The chart is not a validation or pass claim."
    )
    style = {
        "font.family": "DejaVu Sans",
        "font.size": 9.4,
        "axes.edgecolor": "#aeb6bf",
        "axes.labelcolor": "#29313a",
        "axes.linewidth": 0.8,
        "axes.titlecolor": "#29313a",
        "axes.titlesize": 10.5,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "grid.color": "#dfe3e6",
        "grid.linewidth": 0.65,
        "grid.alpha": 0.75,
        "svg.fonttype": "none",
        "svg.hashsalt": "coupfe-cardiac-step2b-comparison-v1",
    }
    colors = {
        "coupfe": "#2367A3",
        "mean": "#343B43",
        "simula": "#B47818",
        "envelope": "#DDE2E7",
        "warning": "#8C5A0A",
    }
    with mpl.rc_context(style):
        figure, axes = plt.subplots(
            2, 3, figsize=(12.6, 8.0), sharex=True, squeeze=False
        )
        figure.subplots_adjust(
            left=0.082,
            right=0.985,
            bottom=0.155,
            top=0.735,
            wspace=0.23,
            hspace=0.29,
        )
        handles = None
        for row, point in enumerate(POINTS):
            for column, component in enumerate(COMPONENTS):
                axis = axes[row][column]
                curves = report["points"][point][component]
                to_mm = lambda values: [1000.0 * value for value in values]
                envelope = axis.fill_between(
                    report["times"],
                    to_mm(curves["publisher_min_m"]),
                    to_mm(curves["publisher_max_m"]),
                    color=colors["envelope"],
                    alpha=0.82,
                    linewidth=0.0,
                    label="Official 10-team range",
                    zorder=1,
                )
                mean = axis.plot(
                    report["times"],
                    to_mm(curves["publisher_mean_m"]),
                    color=colors["mean"],
                    linewidth=1.55,
                    linestyle=(0, (5, 3)),
                    label="Official all-team mean",
                    zorder=3,
                )[0]
                simula = axis.plot(
                    report["times"],
                    to_mm(curves["publisher_simula_m"]),
                    color=colors["simula"],
                    linewidth=1.45,
                    linestyle=(0, (5, 2, 1.2, 2)),
                    label="Official Simula",
                    zorder=4,
                )[0]
                coupfe = axis.plot(
                    report["times"],
                    to_mm(curves["coupfe_m"]),
                    color=colors["coupfe"],
                    linewidth=2.15,
                    linestyle="-",
                    label="CoupFE-Cardiac",
                    zorder=5,
                )[0]
                handles = (coupfe, mean, simula, envelope)
                axis.axhline(0.0, color="#8f969d", linewidth=0.7, zorder=2)
                axis.grid(axis="y")
                axis.set_xlim(0.0, 1.0)
                axis.margins(y=0.08)
                axis.tick_params(colors="#4c5661", labelsize=8.7)
                axis.spines["top"].set_visible(False)
                axis.spines["right"].set_visible(False)
                axis.text(
                    0.025,
                    0.94,
                    point,
                    transform=axis.transAxes,
                    ha="left",
                    va="top",
                    fontsize=9.0,
                    fontweight="bold",
                    color="#29313a",
                    bbox={
                        "boxstyle": "square,pad=0.15",
                        "facecolor": "white",
                        "edgecolor": "none",
                        "alpha": 0.82,
                    },
                )
                if row == 0:
                    axis.set_title(f"{component[1]} displacement", pad=7)
                if row == 1:
                    axis.set_xlabel("time (s)", labelpad=4)
                if point == "p1" and component == "uz":
                    axis.text(
                        0.97,
                        0.10,
                        "0.32–0.48 s: CoupFE z < 0\nall 10 official curves z > 0",
                        transform=axis.transAxes,
                        ha="right",
                        va="bottom",
                        fontsize=8.2,
                        color=colors["warning"],
                        bbox={
                            "boxstyle": "round,pad=0.28",
                            "facecolor": "#FFF7E6",
                            "edgecolor": "#D6A756",
                            "linewidth": 0.8,
                        },
                        zorder=8,
                    )

        figure.text(
            0.5,
            0.965,
            title,
            ha="center",
            va="top",
            fontsize=16,
            fontweight="bold",
            color="#202830",
        )
        figure.text(
            0.5,
            0.925,
            (
                "Full history vs 10-team mean: "
                f"global relative L2={100.0 * metrics['aggregate_relative_l2']:.1f}% · "
                f"RMSE={metrics['aggregate_rmse_mm']:.2f} mm · "
                f"max component error={metrics['aggregate_max_mm']:.2f} mm"
            ),
            ha="center",
            va="top",
            fontsize=9.6,
            color="#3E4852",
        )
        figure.text(
            0.5,
            0.893,
            (
                f"p0/p1 vector relative L2={100.0 * metrics['p0_relative_l2']:.1f}%/"
                f"{100.0 * metrics['p1_relative_l2']:.1f}% · paper Eq. 21 RED="
                f"{100.0 * metrics['p0_paper_red']:.1f}%/"
                f"{100.0 * metrics['p1_paper_red']:.1f}%"
            ),
            ha="center",
            va="top",
            fontsize=9.3,
            color="#3E4852",
        )
        figure.text(
            0.5,
            0.861,
            (
                "Closed Hex8 2×20×17 · pointwise std-κ · 4-rank generalized-α · "
                "dt=0.001 s · active stress + ventricular pressure"
            ),
            ha="center",
            va="top",
            fontsize=9.1,
            color="#59636e",
        )
        if handles is None:
            raise RuntimeError("no comparison curves were rendered")
        figure.legend(
            handles=handles,
            labels=(
                "CoupFE-Cardiac",
                "Official all-team mean",
                "Official Simula",
                "Official 10-team range",
            ),
            loc="upper center",
            bbox_to_anchor=(0.5, 0.825),
            ncol=4,
            frameon=False,
            handlelength=3.2,
            columnspacing=1.7,
            fontsize=9.2,
        )
        figure.text(
            0.018,
            0.46,
            "displacement (mm)",
            ha="center",
            va="center",
            rotation=90,
            fontsize=10,
            color="#29313a",
        )
        figure.text(
            0.5,
            0.052,
            (
                f"Report: {REPORT_NAME} | app {APP_REVISION[:7]} dirty tree, "
                f"runtime manifest {RUNTIME_SOURCE_SHA256[:12]} | Core "
                f"{CORE_REVISION[:7]} | benchmark DOI {REFERENCE_DOI} "
                f"({REFERENCE_LICENSE})"
            ),
            ha="center",
            va="bottom",
            fontsize=8.0,
            color="#59636e",
        )
        figure.text(
            0.5,
            0.023,
            (
                "Global relative L2 and paper RED use different aggregation. No "
                "acceptance threshold is defined; the p1-z branch-sign discrepancy "
                "remains. This is not a validation or convergence claim."
            ),
            ha="center",
            va="bottom",
            fontsize=8.1,
            color="#59636e",
        )

        output.parent.mkdir(parents=True, exist_ok=True)
        if output.suffix.lower() == ".svg":
            stream = io.StringIO()
            figure.savefig(
                stream,
                format="svg",
                facecolor="white",
                metadata={"Date": None},
            )
            svg = _accessible_svg(stream.getvalue(), title, description)
            output.write_text(svg, encoding="utf-8", newline="\n")
        elif output.suffix.lower() == ".png":
            figure.savefig(output, format="png", dpi=180, facecolor="white")
        else:
            raise ValueError("output must end in .svg or .png")
        plt.close(figure)


def _arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    repository = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report",
        type=Path,
        default=repository / "examples" / "cardiac_benchmark" / "results" / REPORT_NAME,
        help="reviewed full-cycle Step-2 Case-B report",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=repository / "docs" / "figures" / "step2_case_b_comparison.svg",
        help="destination .svg or local-QA .png",
    )
    parser.add_argument(
        "--canonical",
        action="store_true",
        help="require the exact renderer stack used for the checked-in SVG",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = _arguments(argv)
    if arguments.canonical:
        _require_canonical_renderer()
    report = _validate_report(arguments.report)
    _render(report, arguments.output)
    print(arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
