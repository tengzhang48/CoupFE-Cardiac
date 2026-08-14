from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

import geometry
import material
import run as cardiac_run
from coupfe.assembly.assemble import assemble_residual, assemble_tangent
from local_pressure import LocalPressureHex8Operator, PAPER_MEAN_DILATATION_LAW


ROOT = Path(__file__).resolve().parents[1]
CURRENT_CORE_REVISION = "e2f42ed5772850a0a23a2ce434f430c287eae5c8"
HISTORICAL_CASE_A_REPORT = (
    ROOT
    / "examples"
    / "cardiac_benchmark"
    / "results"
    / "archive"
    / "truncated_polar"
    / "case_a"
    / "case_a_fbar_1x2x4_dt0p002.report.json"
)
CORRECTED_CASE_A_REPORT = (
    ROOT
    / "examples"
    / "cardiac_benchmark"
    / "results"
    / "archive"
    / "truncated_polar"
    / "case_a"
    / "case_a_fbar_1x2x4_dt0p002_corrected_switch.report.json"
)
HISTORICAL_CASE_A_REPORT_SHA256 = (
    "4278cfc2f282fc7a15a4b67d58e0b65ab425ab81e746067dbb786c7ea080f7c1"
)
CORRECTED_MATERIAL_MODEL_ID = (
    "holzapfel-ogden-smooth-switch-complete-energy-derivative-v1"
)

# These are numerical reproducibility tolerances retained from the preceding
# compiled serial gate.  They are not benchmark-validation acceptance limits.
LOAD_RELATIVE_TOLERANCE = 2.0e-4
LOAD_ABSOLUTE_TOLERANCE_PA = 1.0e-1
DISPLACEMENT_RELATIVE_TOLERANCE = 2.0e-3
DISPLACEMENT_ABSOLUTE_TOLERANCE_M = 2.0e-5


@pytest.mark.slow
def test_serial_driver_reproduces_corrected_case_a_fbar_result(tmp_path):
    if shutil.which("gfortran") is None:
        pytest.fail("slow gate requested but gfortran is unavailable")

    # The predecessor is preserved byte-for-byte, but it is not a numerical
    # oracle for the corrected constitutive law introduced by app 6839c13.
    # Pinning its identity here prevents a future cleanup from silently
    # relabeling or overwriting that historical result.
    assert hashlib.sha256(HISTORICAL_CASE_A_REPORT.read_bytes()).hexdigest() == (
        HISTORICAL_CASE_A_REPORT_SHA256
    )
    historical = json.loads(
        HISTORICAL_CASE_A_REPORT.read_text(encoding="utf-8")
    )["result"]
    assert historical["source_identity"]["app"]["revision"] == (
        "62ad760d2a1731bb9668897863ac026d3768194e"
    )
    assert historical["source_identity"]["core"]["revision"] == (
        "454f73ce2de284262b214a2b37bd676c6aca3c0a"
    )

    report = json.loads(CORRECTED_CASE_A_REPORT.read_text(encoding="utf-8"))
    assert report["schema"] == "coupfe-cardiac-reference-comparison-v2"
    evidence = report["result"]
    assert evidence["result_schema"] == "coupfe-cardiac-result-v1"
    assert evidence["source_identity"]["app"] == {
        "revision": "6839c13b5bc80ec06c897684c51f503e80bd4b19",
        "source_kind": "git-checkout",
        "tree_state": "clean",
    }
    expected_core = evidence["source_identity"]["core"]
    assert expected_core == {
        "revision": CURRENT_CORE_REVISION,
        "source_kind": "git-checkout",
        "source_url": "https://github.com/tengzhang48/CoupFE.git",
        "tree_state": "clean",
    }

    configuration = evidence["configuration"]
    mesh = configuration["mesh"]
    formulations = {
        "hex8_fbar": "fbar",
        "hex8_local_pressure_p0_condensed_logj": "local-pressure",
    }
    cli_formulation = formulations[configuration["formulation"]]

    output = tmp_path / evidence["filename"]
    build_dir = tmp_path / "build"
    command = [
        sys.executable,
        str(ROOT / "examples" / "cardiac_benchmark" / "run.py"),
        "--case",
        evidence["case"],
        "--mesh-topology",
        "polar-ring",
        "--integrator",
        configuration["integrator"],
        "--nonlinear-solver",
        configuration["nonlinear_solver"],
        "--formulation",
        cli_formulation,
        "--mass",
        "lumped",
        "--nt",
        str(mesh["n_t"]),
        "--nmu",
        str(mesh["n_mu"]),
        "--ntheta",
        str(mesh["n_theta"]),
        "--apex-offset",
        str(configuration["apex_offset_rad"]),
        "--perturb",
        str(configuration["model_parameters"]["mesh_perturbation_std_m"]),
        "--dt",
        str(configuration["dt_s"]),
        "--tend",
        str(configuration["t_end_s"]),
        "--build-dir",
        str(build_dir),
        "--out",
        str(output),
    ]
    if not configuration["flip_helix"]:
        command.append("--raw-helix")
    subprocess.run(command, cwd=tmp_path, check=True, timeout=900)

    with np.load(output) as result:
        expected_histories = evidence["retained_histories"]
        expected_times = np.asarray(expected_histories["times_s"], dtype=float)
        expected_tension = np.asarray(
            expected_histories["active_tension_pa"], dtype=float
        )
        expected_pressure = np.asarray(
            expected_histories["cavity_pressure_pa"], dtype=float
        )
        expected_u0 = np.asarray(expected_histories["u0_m"], dtype=float)
        expected_u1 = np.asarray(expected_histories["u1_m"], dtype=float)
        expected_steps = len(expected_times) - 1

        assert bool(result["converged"])
        assert int(result["completed_steps"]) == expected_steps
        assert int(result["expected_steps"]) == expected_steps
        assert len(result["times"]) == expected_steps + 1
        assert np.all(np.isfinite(result["u0"]))
        assert np.all(np.isfinite(result["u1"]))
        assert str(result["result_schema"]) == evidence["result_schema"]
        assert str(result["formulation"]) == configuration["formulation"]
        assert str(result["material_model_id"]) == CORRECTED_MATERIAL_MODEL_ID
        assert configuration["model_parameters"]["material_model_id"] == (
            CORRECTED_MATERIAL_MODEL_ID
        )
        assert str(result["fiber_sampling"]) == configuration["fiber_sampling"]
        assert str(result["mass_representation"]) == "lumped_row_sum"
        assert str(result["tbar_definition"]) == "analytic_parametric"
        assert str(result["tbar_source_filename"]) == ""
        assert str(result["tbar_source_sha256"]) == ""
        assert str(result["tbar_metadata_filename"]) == ""
        assert str(result["tbar_metadata_sha256"]) == ""
        assert str(result["tbar_metadata_schema"]) == ""
        assert str(result["point_sampling"]) == configuration["point_sampling"]
        assert str(result["viscous_rate"]) == configuration["viscous_rate"]
        assert str(result["app_revision"])
        assert str(result["app_tree_state"]) in {"clean", "dirty", "unknown"}
        # The result records the actually resolved Core revision; the
        # release environment pin belongs to check_runtime_core.py, not to
        # this physics reproduction, so the expected value comes from the
        # driver's own metadata resolution rather than a hardcoded SHA.
        expected_core_revision = cardiac_run._runtime_metadata()["core_revision"]
        assert str(result["core_revision"]) == expected_core_revision
        core_source_kind = str(result["core_source_kind"])
        if core_source_kind == "git-checkout":
            assert str(result["core_tree_state"]) == "clean"
        else:
            assert core_source_kind == "pep610-vcs"
            assert str(result["core_tree_state"]) == "installed"
        assert str(result["core_source_url"]) == expected_core["source_url"]
        assert str(result["driver"]) == "examples/cardiac_benchmark/run.py"
        assert "command_argv" not in result

        assert str(result["case"]) == evidence["case"]
        assert str(result["integrator"]) == configuration["integrator"]
        assert str(result["nonlinear_solver"]) == configuration["nonlinear_solver"]
        assert str(result["element_evaluation_mode"]) == "joint"
        assert bool(result["compiled_material_residual_only_available"])
        assert int(result["n_t"]) == mesh["n_t"]
        assert int(result["n_mu"]) == mesh["n_mu"]
        assert int(result["n_theta"]) == mesh["n_theta"]
        assert result["nodes"].shape == (mesh["nodes"], 3)
        assert result["elems"].shape == (mesh["elements"], 8)
        assert result["U_peak"].size == mesh["degrees_of_freedom"]
        assert float(result["apex_offset"]) == pytest.approx(
            configuration["apex_offset_rad"]
        )
        assert float(result["dt"]) == pytest.approx(configuration["dt_s"])
        assert float(result["t_end"]) == pytest.approx(configuration["t_end_s"])
        assert bool(result["flip_helix"]) is configuration["flip_helix"]

        model = configuration["model_parameters"]
        assert float(result["density"]) == pytest.approx(model["density_kg_m3"])
        assert str(result["material_kernel_formulation"]) == model[
            "material_kernel_formulation"
        ]
        assert float(result["material_kappa_pa"]) == pytest.approx(
            model["material_kappa_pa"]
        )
        assert float(result["local_pressure_bulk_modulus_pa"]) == pytest.approx(
            model["local_pressure_bulk_modulus_pa"]
        )
        assert float(result["perturb"]) == pytest.approx(
            model["mesh_perturbation_std_m"]
        )

        solver_configuration = json.loads(str(result["solver_configuration_json"]))
        historical_solver_configuration = evidence["solver_configuration"]
        assert {
            key: solver_configuration[key]
            for key in historical_solver_configuration
        } == historical_solver_configuration
        assert set(solver_configuration) == set(historical_solver_configuration) | {
            "element_evaluation_mode",
            "compiled_material_residual_only_available",
            "epicardial_normal_mode",
        }
        assert solver_configuration["element_evaluation_mode"] == "joint"
        assert solver_configuration["compiled_material_residual_only_available"] is True
        diagnostics = json.loads(str(result["nonlinear_step_diagnostics_json"]))
        expected_diagnostics = evidence["nonlinear_step_diagnostics"]
        assert len(diagnostics) == expected_steps
        np.testing.assert_allclose(
            [record["time"] for record in diagnostics],
            expected_times[1:],
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            [record["dt"] for record in diagnostics],
            configuration["dt_s"],
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_array_equal(
            [record["nonlinear_iterations"] for record in diagnostics],
            [record["nonlinear_iterations"] for record in expected_diagnostics],
        )

        np.testing.assert_allclose(
            result["times"], expected_times, rtol=0.0, atol=1.0e-15
        )
        np.testing.assert_allclose(
            result["tau"],
            expected_tension,
            rtol=LOAD_RELATIVE_TOLERANCE,
            atol=LOAD_ABSOLUTE_TOLERANCE_PA,
        )
        np.testing.assert_allclose(
            result["pres"],
            expected_pressure,
            rtol=LOAD_RELATIVE_TOLERANCE,
            atol=LOAD_ABSOLUTE_TOLERANCE_PA,
        )
        for field, expected in (("u0", expected_u0), ("u1", expected_u1)):
            np.testing.assert_allclose(
                result[field],
                expected,
                rtol=DISPLACEMENT_RELATIVE_TOLERANCE,
                atol=DISPLACEMENT_ABSOLUTE_TOLERANCE_M,
            )

        for point_name in ("p0", "p1"):
            retained_sample = configuration["sampling_points"][point_name]
            element = int(result[f"{point_name}_sampling_element"])
            natural = np.asarray(result[f"{point_name}_sampling_natural"])
            weights = np.asarray(result[f"{point_name}_sampling_weights"])
            reconstruction_error = float(
                result[f"{point_name}_sampling_reconstruction_error_m"]
            )
            assert element == retained_sample["element"]
            assert 0 <= element < mesh["elements"]
            assert natural.shape == (3,)
            assert weights.shape == (8,)
            assert np.all(natural >= -1.0 - 1.0e-9)
            assert np.all(natural <= 1.0 + 1.0e-9)
            np.testing.assert_allclose(
                natural,
                retained_sample["natural_coordinates"],
                rtol=0.0,
                atol=2.0e-11,
            )
            np.testing.assert_allclose(
                weights,
                retained_sample["weights"],
                rtol=0.0,
                atol=2.0e-14,
            )
            assert reconstruction_error == pytest.approx(
                retained_sample["reconstruction_error_m"], abs=2.0e-14
            )
            element_nodes = np.asarray(result["elems"][element], dtype=int)
            reconstructed_point = weights @ result["nodes"][element_nodes]
            np.testing.assert_allclose(
                reconstructed_point,
                result[point_name],
                rtol=0.0,
                atol=2.0e-11,
            )
            observed_reconstruction_error = np.linalg.norm(
                reconstructed_point - result[point_name]
            )
            assert observed_reconstruction_error == pytest.approx(
                reconstruction_error,
                abs=2.0e-14,
            )

        retained_peak = evidence["peak"]
        assert retained_peak["available"] is True
        peak_step = int(result["n_peak"])
        assert peak_step == retained_peak["index"]
        assert float(result["times"][peak_step]) == pytest.approx(
            retained_peak["time_s"]
        )
        assert float(result["tau"][peak_step]) == pytest.approx(
            retained_peak["active_tension_pa"],
            rel=LOAD_RELATIVE_TOLERANCE,
            abs=LOAD_ABSOLUTE_TOLERANCE_PA,
        )
        assert float(result["pres"][peak_step]) == pytest.approx(
            retained_peak["cavity_pressure_pa"],
            rel=LOAD_RELATIVE_TOLERANCE,
            abs=LOAD_ABSOLUTE_TOLERANCE_PA,
        )
        np.testing.assert_allclose(
            result["u0"][peak_step],
            retained_peak["u0_m"],
            rtol=DISPLACEMENT_RELATIVE_TOLERANCE,
            atol=DISPLACEMENT_ABSOLUTE_TOLERANCE_M,
        )
        np.testing.assert_allclose(
            result["u1"][peak_step],
            retained_peak["u1_m"],
            rtol=DISPLACEMENT_RELATIVE_TOLERANCE,
            atol=DISPLACEMENT_ABSOLUTE_TOLERANCE_M,
        )
        peak_field = np.asarray(result["U_peak"]).reshape(-1, 3)
        assert np.all(np.isfinite(peak_field))
        assert np.max(np.linalg.norm(peak_field, axis=1)) > 0.0
        for point_name in ("p0", "p1"):
            element = int(result[f"{point_name}_sampling_element"])
            weights = np.asarray(result[f"{point_name}_sampling_weights"])
            element_nodes = np.asarray(result["elems"][element], dtype=int)
            displacement_history = result[f"u{point_name[1:]}"]
            np.testing.assert_allclose(
                weights @ peak_field[element_nodes],
                displacement_history[peak_step],
                rtol=DISPLACEMENT_RELATIVE_TOLERANCE,
                atol=DISPLACEMENT_ABSOLUTE_TOLERANCE_M,
            )

        det_f = np.asarray(result["det_f_gauss_peak"], dtype=float)
        det_f_summary = evidence["det_f_gauss_peak_summary"]
        assert det_f_summary["available"] is True
        assert det_f.shape == tuple(det_f_summary["shape"])
        assert det_f.size == det_f_summary["count"]
        assert np.all(np.isfinite(det_f))
        assert np.min(det_f) > 0.0
        for observed, key in (
            (np.min(det_f), "minimum"),
            (np.mean(det_f), "mean"),
            (np.max(det_f), "maximum"),
        ):
            assert observed == pytest.approx(
                det_f_summary[key],
                rel=DISPLACEMENT_RELATIVE_TOLERANCE,
                abs=DISPLACEMENT_ABSOLUTE_TOLERANCE_M,
            )
        assert evidence["element_pressure_peak_pa_summary"] == {"available": False}
        assert np.asarray(result["element_pressure_peak_pa"]).size == 0


@pytest.mark.slow
def test_local_pressure_uses_standard_zero_kappa_kernel_and_consistent_tangent(
    tmp_path,
):
    if shutil.which("gfortran") is None:
        pytest.fail("slow gate requested but gfortran is unavailable")

    mesh = geometry.build_mesh(
        n_t=1, n_mu=2, n_theta=4, flip_helix=True, apex_offset=0.2
    )
    cardiac_run._KERNEL.clear()
    group, element, _ta_index, local_pressure = cardiac_run.build_group(
        mesh,
        0.01,
        build_dir=tmp_path / "local-pressure-build",
        formulation="local-pressure",
    )
    property_names = list(material.CardiacHex8()._mat.props.keys())
    kappa_index = property_names.index("kappa")
    assert element.props[kappa_index] == 0.0
    assert local_pressure is not None
    assert local_pressure.bulk_modulus == pytest.approx(material.HO_PARAMS["kappa"])

    ndof = 3 * mesh.n_node
    rng = np.random.default_rng(20260801)
    direction = rng.normal(size=ndof)
    direction /= np.linalg.norm(direction)
    displacement = rng.normal(scale=2.0e-5, size=ndof)
    step = 1.0e-7
    paper_pressure = LocalPressureHex8Operator(
        mesh.nodes,
        mesh.elems,
        ndof,
        bulk_modulus=material.HO_PARAMS["kappa"],
        pressure_law=PAPER_MEAN_DILATATION_LAW,
    )
    for pressure in (local_pressure, paper_pressure):
        operators = [group, pressure]
        tangent = assemble_tangent(
            operators, displacement, None, 0.01, 0.01, ndof
        )

        def residual(values):
            assembled, _ = assemble_residual(
                operators, values, None, 0.01, 0.01, ndof
            )
            return np.asarray(assembled)

        finite_difference = (
            residual(displacement + step * direction)
            - residual(displacement - step * direction)
        ) / (2.0 * step)
        tangent_action = tangent @ direction
        relative_error = np.linalg.norm(
            tangent_action - finite_difference
        ) / max(np.linalg.norm(finite_difference), 1.0)
        assert relative_error < 2.0e-7

    fbar_group, fbar_element, _ta_index, fbar_pressure = cardiac_run.build_group(
        mesh,
        0.01,
        build_dir=tmp_path / "fbar-build",
        formulation="fbar",
    )
    assert fbar_group is not group
    assert fbar_element.props[kappa_index] == pytest.approx(material.HO_PARAMS["kappa"])
    assert fbar_pressure is None
