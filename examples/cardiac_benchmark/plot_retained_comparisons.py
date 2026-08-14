"""Render deterministic comparisons from the selected Case A/B reports.

The reports already contain the CoupFE-Cardiac histories and the attributed
benchmark all-team means on the benchmark's 101-point canonical time grid.
Consequently, this renderer does not read the separately distributed pickle
archive.  It deliberately accepts only the two configurations named below so
that a different run cannot silently replace the evidence shown in the public
figures. Case A uses the compact retained Step 0A schema for the selected
closed-mesh run. Until the new closed Step 0B run is integrated, Case B uses
the archived generic report from the historical non-benchmark open-tip mesh.

Install the optional plotting dependency with ``pip install -e '.[reference]'``
and run this file from any directory.  By default it regenerates:

* ``docs/figures/case_a_comparison.svg``
* ``docs/figures/case_b_comparison.svg``

Use ``--canonical`` to require the renderer stack used for the checked-in,
byte-reproducible SVGs. Other supported Matplotlib environments may render the
same report data but can place text or generate internal IDs differently.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys
from typing import Any, Dict, List, Mapping, Optional, Tuple
from xml.sax.saxutils import escape


GENERIC_REPORT_SCHEMA = "coupfe-cardiac-reference-comparison-v2"
STEP0A_REPORT_SCHEMA = "coupfe-cardiac-step0a-retained-comparison-v1"
REFERENCE_DOI = "10.5281/zenodo.14260459"
REFERENCE_LICENSE = "CC-BY-4.0"
CURVE_POINTS = 101
COMPONENTS = ("x", "y", "z")
POINTS = ("p0", "p1")
CANONICAL_PYTHON_VERSION = (3, 10, 8)
CANONICAL_MATPLOTLIB_VERSION = "3.10.9"
CANONICAL_FREETYPE_VERSION = "2.6.1"
CANONICAL_FONT_SHA256 = (
    "3fdf69cabf06049ea70a00b5919340e2ce1e6d02b0cc3c4b44fb6801bd1e0d22"
)


@dataclass(frozen=True)
class FigureSpec:
    case: str
    reference_case: str
    benchmark_case_label: str
    report_schema: str
    compact_report: bool
    identity_status: str
    explicit_archive_identity_fields: bool
    report_name: str
    output_name: str
    formulation: str
    formulation_label: str
    mesh_axis_names: Tuple[str, str, str]
    mesh_axes: Tuple[int, int, int]
    mesh_topology: Optional[str]
    elements: int
    nodes: int
    degrees_of_freedom: int
    dt_s: float
    time_integration_label: str
    identity_label: str
    app_revision: str
    core_revision: str
    result_name: str
    result_sha256: str
    report_sha256: str


FIGURE_SPECS = (
    FigureSpec(
        case="A",
        reference_case="step_0A",
        benchmark_case_label="Step 0A",
        report_schema=STEP0A_REPORT_SCHEMA,
        compact_report=True,
        identity_status="legacy-inferred",
        explicit_archive_identity_fields=False,
        report_name="case_a_local_pressure_4x36x32_dt0p001.report.json",
        output_name="case_a_comparison.svg",
        formulation="hex8_local_pressure_p0_condensed_logj",
        formulation_label="Hex8 Q1/P0 condensed local pressure",
        mesh_axis_names=("n_t", "n_core", "n_radial"),
        mesh_axes=(4, 36, 32),
        mesh_topology="closed_multiblock_disk",
        elements=23616,
        nodes=29885,
        degrees_of_freedom=89655,
        dt_s=0.001,
        time_integration_label="consistent-mass generalized-alpha",
        identity_label="Step 0A legacy-inferred",
        app_revision="016a4f9eec6f2a4c74d10c734ddff3e24cf343de",
        core_revision="e2f42ed5772850a0a23a2ce434f430c287eae5c8",
        result_name="caseA_ga_local_pressure_rank8_t100.npz",
        result_sha256=(
            "ba9b31ec533398be1f39fc9a898e72f77d9587c90f9b7d9e00ce91e4d2ae6a6c"
        ),
        report_sha256=(
            "bbd26f3b30819ff2b67ffb48c9ad52cc9825c7fa0486e3984673c9e349bf82b1"
        ),
    ),
    FigureSpec(
        case="B",
        reference_case="step_0B",
        benchmark_case_label="Step 0B",
        report_schema=GENERIC_REPORT_SCHEMA,
        compact_report=False,
        identity_status="",
        explicit_archive_identity_fields=False,
        report_name="case_b_local_pressure_2x36x48_dt0p002.report.json",
        output_name="case_b_comparison.svg",
        formulation="hex8_local_pressure_p0_condensed_logj",
        formulation_label="Hex8 Q1/P0 condensed local pressure",
        mesh_axis_names=("n_t", "n_mu", "n_theta"),
        mesh_axes=(2, 36, 48),
        mesh_topology=None,
        elements=3456,
        nodes=5328,
        degrees_of_freedom=15984,
        dt_s=0.002,
        time_integration_label="",
        identity_label="",
        app_revision="e07993bcf1166bd20eb87370c0b458552753e7ee",
        core_revision="454f73ce2de284262b214a2b37bd676c6aca3c0a",
        result_name="case_b_local_pressure_2x36x48_dt0p002.npz",
        result_sha256=(
            "6bcb5e0b98a044335d452be021cb765223a7ed3c21e4a2fc89496ba8d66e911a"
        ),
        report_sha256=(
            "d409ecaac0c0abf418fce3ab2f0549979d38e02b9d73381d56982f7fc4e3bf14"
        ),
    ),
)


def _reject_nonstandard_constant(value: str) -> None:
    raise ValueError(f"non-finite JSON constant {value!r}")


def _mapping(value: Any, description: str) -> Mapping[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{description} must be an object")
    return value


def _finite_number(value: Any, description: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{description} must be numeric")
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"{description} must be finite")
    return number


def _integer(value: Any, description: str) -> int:
    number = _finite_number(value, description)
    if not number.is_integer():
        raise ValueError(f"{description} must be an integer")
    return int(number)


def _curve(
    value: Any,
    description: str,
) -> List[Tuple[float, float, float]]:
    if not isinstance(value, list) or len(value) != CURVE_POINTS:
        raise ValueError(
            f"{description} must contain exactly {CURVE_POINTS} samples"
        )
    result = []
    for index, row in enumerate(value):
        if not isinstance(row, list) or len(row) != len(COMPONENTS):
            raise ValueError(f"{description}[{index}] must contain x, y, and z")
        result.append(
            tuple(
                _finite_number(component, f"{description}[{index}][{axis}]")
                for axis, component in zip(COMPONENTS, row)
            )
        )
    return result


def _validate_report(
    path: Path,
    spec: FigureSpec,
) -> Dict[str, Any]:
    try:
        report_bytes = path.read_bytes()
        report_digest = hashlib.sha256(report_bytes).hexdigest()
        if report_digest != spec.report_sha256:
            raise ValueError(
                f"report SHA-256 is {report_digest}; expected {spec.report_sha256}"
            )
        report = json.loads(
            report_bytes.decode("utf-8"),
            parse_constant=_reject_nonstandard_constant,
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(
            f"cannot read a strict JSON report from {path}: {error}"
        ) from error

    root = _mapping(report, str(path))
    if root.get("schema") != spec.report_schema:
        raise ValueError(
            f"{path} does not use report schema {spec.report_schema!r}"
        )

    if spec.compact_report:
        benchmark_identity = _mapping(
            root.get("benchmark_identity"), f"{path}: benchmark_identity"
        )
        exact_identity = {
            "benchmark": "Benchmark 1",
            "case": spec.reference_case,
            "status": spec.identity_status,
            "explicit_archive_identity_fields": (
                spec.explicit_archive_identity_fields
            ),
        }
        for key, expected in exact_identity.items():
            if benchmark_identity.get(key) != expected:
                raise ValueError(
                    f"{path}: benchmark_identity.{key} is "
                    f"{benchmark_identity.get(key)!r}; expected {expected!r}"
                )

    result = _mapping(root.get("result"), f"{path}: result")
    reference = _mapping(root.get("reference"), f"{path}: reference")
    comparison = _mapping(root.get("comparison"), f"{path}: comparison")
    configuration = _mapping(
        result.get("configuration"), f"{path}: result.configuration"
    )
    mesh = _mapping(configuration.get("mesh"), f"{path}: configuration.mesh")

    expected_source_identity = {
        "app": {
            "revision": spec.app_revision,
            "source_kind": "git-checkout",
            "tree_state": "clean",
        },
        "core": {
            "revision": spec.core_revision,
            "source_kind": "git-checkout",
            "source_url": "https://github.com/tengzhang48/CoupFE.git",
            "tree_state": "clean",
        },
    }
    source_identity = _mapping(
        result.get("source_identity"), f"{path}: result.source_identity"
    )
    if source_identity != expected_source_identity:
        raise ValueError(
            f"{path}: source_identity does not match the selected retained run"
        )

    exact_values = {
        "result.case": (result.get("case"), spec.case),
        "result.reference_case": (
            result.get("reference_case"),
            spec.reference_case,
        ),
        "reference.case": (reference.get("case"), spec.reference_case),
        "reference.doi": (reference.get("doi"), REFERENCE_DOI),
        "reference.license": (reference.get("license"), REFERENCE_LICENSE),
        "configuration.formulation": (
            configuration.get("formulation"),
            spec.formulation,
        ),
        "result.filename": (result.get("filename"), spec.result_name),
        "result.sha256": (result.get("sha256"), spec.result_sha256),
    }
    for description, (actual, expected) in exact_values.items():
        if actual != expected:
            raise ValueError(
                f"{path}: {description} is {actual!r}; expected {expected!r}"
            )

    actual_dt = _finite_number(configuration.get("dt_s"), f"{path}: dt_s")
    if actual_dt != spec.dt_s:
        raise ValueError(f"{path}: dt_s is {actual_dt}; expected {spec.dt_s}")
    actual_t_end = _finite_number(
        configuration.get("t_end_s"), f"{path}: t_end_s"
    )
    if actual_t_end != 1.0:
        raise ValueError(f"{path}: t_end_s is {actual_t_end}; expected 1.0")

    mesh_expectations = {
        **dict(zip(spec.mesh_axis_names, spec.mesh_axes)),
        "elements": spec.elements,
        "nodes": spec.nodes,
        "degrees_of_freedom": spec.degrees_of_freedom,
    }
    for key, expected in mesh_expectations.items():
        actual = _integer(mesh.get(key), f"{path}: mesh.{key}")
        if actual != expected:
            raise ValueError(
                f"{path}: mesh.{key} is {actual}; expected {expected}"
            )
    if spec.mesh_topology is not None and mesh.get("topology") != spec.mesh_topology:
        raise ValueError(
            f"{path}: mesh.topology is {mesh.get('topology')!r}; "
            f"expected {spec.mesh_topology!r}"
        )

    if spec.compact_report:
        expected_method = {
            "integrator": "generalized-alpha",
            "mass_representation": "consistent_q1_hex8",
            "generalized_alpha": {
                "alpha_f": 0.4,
                "alpha_m": 0.2,
                "beta": 0.36,
                "gamma": 0.7,
                "stage_contract": "simula-source-matched-v1",
            },
        }
        altered = [
            key
            for key, expected in expected_method.items()
            if configuration.get(key) != expected
        ]
        if altered:
            raise ValueError(
                f"{path}: selected Step 0A time/mass method differs: {altered}"
            )

    times_raw = reference.get("canonical_grid_s")
    if not isinstance(times_raw, list) or len(times_raw) != CURVE_POINTS:
        raise ValueError(
            f"{path}: canonical_grid_s must contain exactly {CURVE_POINTS} samples"
        )
    times = [
        _finite_number(value, f"{path}: canonical_grid_s[{index}]")
        for index, value in enumerate(times_raw)
    ]
    for index, value in enumerate(times):
        expected = index / (CURVE_POINTS - 1)
        if not math.isclose(value, expected, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(
                f"{path}: canonical_grid_s[{index}] is {value}; expected {expected}"
            )

    ours = _mapping(
        comparison.get("ours_on_canonical_grid_m"),
        f"{path}: comparison.ours_on_canonical_grid_m",
    )
    means = _mapping(
        reference.get("mean_curves_m"), f"{path}: reference.mean_curves_m"
    )
    red = _mapping(comparison.get("red"), f"{path}: comparison.red")
    parsed_ours = {}
    parsed_means = {}
    parsed_red = {}
    for point in POINTS:
        parsed_ours[point] = _curve(
            ours.get(point), f"{path}: CoupFE curve {point}"
        )
        parsed_means[point] = _curve(
            means.get(point), f"{path}: benchmark mean curve {point}"
        )
        point_red = _mapping(red.get(point), f"{path}: RED {point}")
        red_value = _finite_number(
            point_red.get("ours"), f"{path}: RED {point}.ours"
        )
        if red_value < 0.0:
            raise ValueError(f"{path}: RED {point}.ours must be nonnegative")
        parsed_red[point] = red_value

    return {
        "times": times,
        "ours": parsed_ours,
        "means": parsed_means,
        "red": parsed_red,
    }


def _require_canonical_renderer() -> None:
    try:
        import matplotlib as mpl
        from matplotlib import font_manager, ft2font
    except ImportError as error:
        raise RuntimeError(
            "Matplotlib is required; install the `reference` optional dependency"
        ) from error

    font_path = Path(font_manager.findfont("DejaVu Sans", fallback_to_default=False))
    font_sha256 = hashlib.sha256(font_path.read_bytes()).hexdigest()
    actual = {
        "Python": tuple(sys.version_info[:3]),
        "Matplotlib": mpl.__version__,
        "FreeType": ft2font.__freetype_version__,
        "DejaVu Sans SHA-256": font_sha256,
    }
    expected = {
        "Python": CANONICAL_PYTHON_VERSION,
        "Matplotlib": CANONICAL_MATPLOTLIB_VERSION,
        "FreeType": CANONICAL_FREETYPE_VERSION,
        "DejaVu Sans SHA-256": CANONICAL_FONT_SHA256,
    }
    mismatches = [
        f"{key}: found {actual[key]!r}, expected {expected[key]!r}"
        for key in expected
        if actual[key] != expected[key]
    ]
    if mismatches:
        raise RuntimeError(
            "canonical SVG renderer mismatch; " + "; ".join(mismatches)
        )


def _accessible_svg(svg: str, title: str, description: str) -> str:
    # Matplotlib's SVG 1.1 prologue points at an external DTD and its default
    # RDF block carries resource-valued metadata.  This standalone release SVG
    # needs neither: retain one native SVG title/description pair instead.
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
        "\n <title id=\"figure-title\">"
        + escape(title)
        + "</title>\n <desc id=\"figure-description\">"
        + escape(description)
        + "</desc>"
    )
    result = svg[:index] + root + ">" + children + svg[end + 1 :]
    forbidden = ("<!DOCTYPE", "<metadata", "rdf:", "dc:", "rdf:resource")
    if any(token in result for token in forbidden):
        raise RuntimeError("SVG sanitization left external or RDF metadata")
    if result.count('<title id="figure-title">') != 1:
        raise RuntimeError("SVG must contain exactly one accessible title")
    if result.count('<desc id="figure-description">') != 1:
        raise RuntimeError("SVG must contain exactly one accessible description")
    # Matplotlib emits spaces at the ends of many path-data lines.  Normalize
    # them here so the generated release artifacts pass the repository's
    # whitespace gate and remain byte-reproducible across clean checkouts.
    return "\n".join(line.rstrip() for line in result.splitlines()) + "\n"


def _render(
    report: Dict[str, Any],
    spec: FigureSpec,
    output_path: Path,
) -> None:
    try:
        import matplotlib as mpl
        mpl.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise RuntimeError(
            "Matplotlib is required; install the `reference` optional dependency"
        ) from error

    title = f"Benchmark 1, Case {spec.case} — displacement comparison"
    description = (
        f"Six line charts compare CoupFE-Cardiac with the benchmark all-team "
        f"mean from zero to one second. Rows show p0 and p1; columns show x, y, "
        f"and z displacement in millimetres. Reported relative discrepancies "
        f"are {report['red']['p0']:.3f} at p0 and {report['red']['p1']:.3f} at p1. "
        f"Source: {spec.report_name}, application {spec.app_revision[:7]}, "
        f"Core {spec.core_revision[:7]}, benchmark DOI {REFERENCE_DOI}. "
        + (
            f"The {spec.benchmark_case_label} identity is explicitly labelled "
            f"{spec.identity_status}. "
            if spec.compact_report
            else ""
        )
        + "The chart is not a validation or pass claim."
    )

    style = {
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
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
        "svg.hashsalt": "coupfe-cardiac-retained-comparison-v1",
    }
    ours_color = "#2367A3"
    reference_color = "#3F454B"

    with mpl.rc_context(style):
        figure, axes = plt.subplots(
            2,
            3,
            figsize=(12.4, 7.6),
            sharex=True,
            squeeze=False,
        )
        figure.subplots_adjust(
            left=0.085,
            right=0.985,
            bottom=0.16,
            top=0.79,
            wspace=0.22,
            hspace=0.28,
        )

        reference_handle = None
        ours_handle = None
        for row, point in enumerate(POINTS):
            for column, component in enumerate(COMPONENTS):
                axis = axes[row][column]
                ours_mm = [sample[column] * 1000.0 for sample in report["ours"][point]]
                reference_mm = [
                    sample[column] * 1000.0 for sample in report["means"][point]
                ]
                reference_handle = axis.plot(
                    report["times"],
                    reference_mm,
                    color=reference_color,
                    linewidth=1.65,
                    linestyle=(0, (5, 3)),
                    label="Benchmark all-team mean",
                    zorder=2,
                )[0]
                ours_handle = axis.plot(
                    report["times"],
                    ours_mm,
                    color=ours_color,
                    linewidth=2.05,
                    linestyle="-",
                    label="CoupFE-Cardiac",
                    zorder=3,
                )[0]
                axis.axhline(0.0, color="#8f969d", linewidth=0.7, zorder=1)
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
                    fontsize=9.5,
                    fontweight="bold",
                    color="#29313a",
                )
                if row == 0:
                    axis.set_title(f"{component} displacement", pad=7)
                if row == 1:
                    axis.set_xlabel("time (s)", labelpad=4)

        figure.text(
            0.5,
            0.955,
            title,
            ha="center",
            va="top",
            fontsize=16,
            fontweight="bold",
            color="#202830",
        )
        n_t, n_mu, n_theta = spec.mesh_axes
        if spec.mesh_topology == "closed_multiblock_disk":
            configuration = (
                f"Configuration: closed t/core/radial {n_t}×{n_mu}×{n_theta} "
                f"Hex8 Q1/P0 local pressure; {spec.elements:,} elements; "
                f"{spec.nodes:,} nodes; {spec.degrees_of_freedom:,} DOF"
            )
        else:
            configuration = (
                f"Configuration: {spec.formulation_label}; "
                f"mesh {n_t}×{n_mu}×{n_theta}; "
                f"{spec.elements:,} elements; {spec.nodes:,} nodes; "
                f"{spec.degrees_of_freedom:,} DOF; dt={spec.dt_s:.3f} s"
            )
        figure.text(
            0.5,
            0.912,
            configuration,
            ha="center",
            va="top",
            fontsize=9.2,
            color="#4c5661",
        )
        discrepancy = (
            "Relative discrepancy (benchmark paper Eq. 21): "
            f"p0={report['red']['p0']:.3f}; p1={report['red']['p1']:.3f}"
        )
        if spec.time_integration_label:
            discrepancy = (
                f"{spec.time_integration_label}; dt={spec.dt_s:.3f} s | "
                + discrepancy
            )
        figure.text(
            0.5,
            0.878,
            discrepancy,
            ha="center",
            va="top",
            fontsize=9.2,
            color="#4c5661",
        )
        if ours_handle is None or reference_handle is None:
            raise RuntimeError("no comparison lines were rendered")
        figure.legend(
            handles=(ours_handle, reference_handle),
            labels=("CoupFE-Cardiac", "Benchmark all-team mean"),
            loc="upper center",
            bbox_to_anchor=(0.5, 0.852),
            ncol=2,
            frameon=False,
            handlelength=3.4,
            columnspacing=2.5,
            fontsize=9.5,
        )
        figure.text(
            0.018,
            0.475,
            "displacement (mm)",
            ha="center",
            va="center",
            rotation=90,
            fontsize=10,
            color="#29313a",
        )
        provenance = (
            f"Report: {spec.report_name} | app {spec.app_revision[:7]} | "
            f"Core {spec.core_revision[:7]} | benchmark DOI {REFERENCE_DOI} "
            f"({REFERENCE_LICENSE})"
        )
        if spec.identity_label:
            provenance += f" | {spec.identity_label}"
        figure.text(
            0.5,
            0.053,
            provenance,
            ha="center",
            va="bottom",
            fontsize=8.2,
            color="#59636e",
        )
        figure.text(
            0.5,
            0.025,
            (
                "The curves show component-level agreement and discrepancies; "
                "they are not a validation or pass claim."
            ),
            ha="center",
            va="bottom",
            fontsize=8.2,
            color="#59636e",
        )

        stream = io.StringIO()
        figure.savefig(
            stream,
            format="svg",
            facecolor="white",
            metadata={"Date": None},
        )
        plt.close(figure)

    svg = _accessible_svg(stream.getvalue(), title, description)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(svg, encoding="utf-8", newline="\n")


def _arguments() -> argparse.Namespace:
    repository = Path(__file__).resolve().parents[2]
    results = repository / "examples" / "cardiac_benchmark" / "results"
    archive = results / "archive" / "truncated_polar"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case-a-report",
        type=Path,
        default=results / FIGURE_SPECS[0].report_name,
        help="selected closed-mesh Case A comparison report",
    )
    parser.add_argument(
        "--case-b-report",
        type=Path,
        default=archive / "case_b" / FIGURE_SPECS[1].report_name,
        help="archived truncated-polar fine Q1/P0 Case B comparison report",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repository / "docs" / "figures",
        help="directory for both deterministic SVG files",
    )
    parser.add_argument(
        "--canonical",
        action="store_true",
        help="require the exact renderer stack used for checked-in SVG bytes",
    )
    return parser.parse_args()


def main() -> None:
    arguments = _arguments()
    if arguments.canonical:
        _require_canonical_renderer()
    report_paths = (arguments.case_a_report, arguments.case_b_report)
    for report_path, spec in zip(report_paths, FIGURE_SPECS):
        report = _validate_report(report_path, spec)
        output_path = arguments.output_dir / spec.output_name
        _render(report, spec, output_path)
        print(output_path)


if __name__ == "__main__":
    main()
