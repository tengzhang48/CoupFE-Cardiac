from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from activation import p_of_t
from benchmark_parameters import benchmark_configuration, benchmark_metadata
import post
import publish_step0_comparison as publisher


def _identity(path):
    payload = path.read_bytes()
    return hashlib.sha256(payload).hexdigest(), len(payload)


def _write_json(path, value):
    path.write_text(
        json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def _benchmark_record():
    selected = benchmark_configuration(0, "B")
    metadata = benchmark_metadata(
        selected,
        material_parameters=selected.material_parameters,
        activation_parameters=selected.activation_parameters,
        pressure_parameters=selected.pressure_parameters,
    )
    return {
        "benchmark": 1,
        "step": 0,
        "case": "B",
        "configuration_id": selected.identity,
        "identity_scope": metadata["benchmark_identity_scope"],
        "load_contract": selected.load_contract,
        "peak_load_definition": metadata["benchmark_peak_load_definition"],
        "active_stress_enabled": False,
        "pressure_enabled": True,
        "material_parameters": dict(selected.material_parameters),
        "activation_parameters": dict(selected.activation_parameters),
        "pressure_parameters": dict(selected.pressure_parameters),
        "runtime_source_manifest": json.loads(
            metadata["benchmark_runtime_source_manifest_json"]
        ),
        "runtime_source_sha256": metadata["benchmark_runtime_source_sha256"],
    }


def _reference_and_comparison():
    grid = np.linspace(0.0, 1.0, 101)
    mean_p0 = np.column_stack((0.01 + 0.002 * grid, grid * 0.001, -grid * 0.003))
    mean_p1 = np.column_stack((0.008 + 0.001 * grid, -grid * 0.002, grid * 0.004))
    ours_p0 = mean_p0 * 1.05
    ours_p1 = mean_p1 * 0.96
    case = "step_0B"
    names = [
        "monoventricular_nonblinded_{}_group_{}.pickle".format(case, suffix)
        for suffix in post.REFERENCE_MANIFEST_SUFFIXES
    ]
    files = []
    for name, team in zip(names, post.REFERENCE_MANIFEST_SUFFIXES):
        files.append(
            {
                "filename": name,
                "sha256": hashlib.sha256(name.encode("utf-8")).hexdigest(),
                "size_bytes": 100 + len(files),
                "source_sample_count": 101,
                "source_time_start": 0.0,
                "source_time_end": 1.0,
                "team": team,
            }
        )
    reference = {
        "canonical_grid_s": grid.tolist(),
        "case": case,
        "doi": post.REFERENCE_DOI,
        "endpoint_offset_tolerance_s": post._ENDPOINT_OFFSET_TOLERANCE,
        "license": post.REFERENCE_LICENSE,
        "mean_curves_m": {"p0": mean_p0.tolist(), "p1": mean_p1.tolist()},
        "published_archive_identity": post.REFERENCE_ARCHIVE_IDENTITY,
        "selection": {
            "excluded_aliases": [],
            "policy": post.REFERENCE_SELECTION_POLICY,
            "selected_count": 10,
            "selected_files": names,
            "upstream_figures_py_identity": post.REFERENCE_FIGURES_PY_IDENTITY,
            "upstream_manifest_variable": post.REFERENCE_MANIFEST_VARIABLES[case],
        },
        "team_files": files,
    }
    teams = {team: 0.1 for team in post.REFERENCE_MANIFEST_SUFFIXES}
    comparison = {
        "metric": "relative discrepancy (benchmark paper Eq. 21)",
        "ours_on_canonical_grid_m": {
            "p0": ours_p0.tolist(),
            "p1": ours_p1.tolist(),
        },
        "red": {
            "p0": {"ours": post.red(ours_p0, mean_p0), "teams": teams},
            "p1": {"ours": post.red(ours_p1, mean_p1), "teams": teams},
        },
    }
    return reference, comparison


def _generic_report(normalized_log):
    times = np.linspace(0.0, 1.0, 1001)
    load_times = times.copy()
    load_times[1:] -= publisher.GENERALIZED_ALPHA["alpha_f"] * publisher.DT_S
    pressure = p_of_t(load_times, t_span=(0.0, 1.0))
    u0 = np.column_stack((0.01 * times, 0.001 * times, -0.002 * times))
    u1 = np.column_stack((0.008 * times, -0.001 * times, 0.003 * times))
    peak_index = int(np.argmax(pressure))
    mesh = {
        "core_half_width": 0.36,
        "degrees_of_freedom": 81,
        "elements": 8,
        "n_core": 2,
        "n_mu": 0,
        "n_radial": 2,
        "n_side": 0,
        "n_t": 1,
        "n_theta": 0,
        "nodes": 27,
        "topology": "closed_multiblock_disk",
    }
    sampling = {
        point: {
            "element": index,
            "natural_coordinates": [0.0, 0.0, 0.0],
            "reconstruction_error_m": 0.0,
            "weights": [0.125] * 8,
        }
        for index, point in enumerate(("p0", "p1"))
    }
    configuration = {
        "apex_offset_rad": 0.0,
        "benchmark": _benchmark_record(),
        "dt_s": 0.001,
        "fiber_sampling": "gp_direct_rule",
        "fiber_sampling_option": "gp-direct",
        "flip_helix": True,
        "formulation": publisher.FORMULATION,
        "generalized_alpha": {
            **publisher.GENERALIZED_ALPHA,
            "load_evaluation_times_s": load_times.tolist(),
        },
        "integrator": "generalized-alpha",
        "isotropic": False,
        "load_horizon_origin": "recorded",
        "load_horizon_s": 1.0,
        "mass_representation": "consistent_q1_hex8",
        "material_eta_pa_s": 100.0,
        "mesh": mesh,
        "method_metadata_origin": "recorded",
        "model_parameters": publisher.MODEL_PARAMETERS,
        "nonlinear_solver": "petsc-snes-mpi",
        "parameter_variant": "benchmark_eta",
        "point_sampling": "hex8_reference_isoparametric",
        "sampling_points": sampling,
        "t_end_s": 1.0,
        "tbar": {
            "definition": "laplace_presolved",
            "metadata_filename": "tbar.meta.json",
            "metadata_schema": post.TBAR_METADATA_SCHEMA,
            "metadata_sha256": "a" * 64,
            "source_filename": "tbar.npy",
            "source_sha256": "b" * 64,
        },
        "viscous_rate": "velocity_consistent_green_lagrange_at_alpha_f_stage",
        "viscous_term_active": True,
    }
    diagnostics = [
        {
            "time": float(times[index]),
            "dt": 0.001,
            "ranks": 8,
            "snes_converged_reason": 2,
            "final_residual_norm": 1.0e-12,
            "residual_acceptance_threshold": 1.0e-10,
            "function_domain_rejections": 0,
            "ksp_residual_histories": [[1.0, 1.0e-12]],
        }
        for index in range(1, len(times))
    ]
    audit = {
        "geometry": {
            "mesh_topology": "closed_multiblock_disk",
            "require_closed": True,
            "nodes": 27,
            "elements": 8,
            "intentional_traction_free_tip_faces": 0,
            "unexpected_unclassified_exterior_faces": 0,
            "passed": True,
            "failures": [],
        },
        "pressure": {
            "pressure_surface_policy": "closed_analytic_base_projection",
            "passed": True,
            "failures": [],
        },
        "robin": {"passed": True, "failures": []},
    }
    reference, comparison = _reference_and_comparison()
    return {
        "bounded_claim": "Synthetic generic report for publisher tests.",
        "comparison": comparison,
        "reference": reference,
        "result": {
            "case": "B",
            "configuration": configuration,
            "det_f_gauss_peak_summary": {
                "available": True,
                "shape": [8, 8],
                "count": 64,
                "minimum": 0.8,
                "mean": 1.0,
                "maximum": 1.2,
            },
            "element_pressure_peak_pa_summary": {
                "available": True,
                "shape": [8],
                "count": 8,
                "minimum": -1000.0,
                "mean": 0.0,
                "maximum": 1000.0,
            },
            "filename": "synthetic-step0b.npz",
            "mpi_metadata": {
                "world_size": 8,
                "contract": publisher.MPI_CONTRACT["B"],
                "implementation": publisher.MPI_IMPLEMENTATION,
                "linear_solver_profile": publisher.LINEAR_SOLVER_PROFILE,
            },
            "nonlinear_step_diagnostics": diagnostics,
            "normalized_run_log": normalized_log,
            "peak": {
                "available": True,
                "index": peak_index,
                "time_s": float(times[peak_index]),
                "active_tension_pa": 0.0,
                "cavity_pressure_pa": float(pressure[peak_index]),
                "u0_m": u0[peak_index].tolist(),
                "u1_m": u1[peak_index].tolist(),
            },
            "peak_circumferential_ring_rotation": {
                "available": False,
                "reason": "closed mesh",
            },
            "pre_solve_audit": audit,
            "reference_case": "step_0B",
            "result_schema": post.RESULT_SCHEMA,
            "retained_histories": {
                "active_tension_pa": np.zeros_like(times).tolist(),
                "cavity_pressure_pa": pressure.tolist(),
                "times_s": times.tolist(),
                "u0_m": u0.tolist(),
                "u1_m": u1.tolist(),
            },
            "runtime_versions": {"python_version": "3.10.8"},
            "sha256": "c" * 64,
            "size_bytes": 123456,
            "solver_configuration": {
                "ranks": 8,
                "linear_solver_profile": publisher.LINEAR_SOLVER_PROFILE,
                "preconditioner": publisher.PRECONDITIONER,
                "time_integrator": "generalized-alpha",
                "mass_representation": "consistent_q1_hex8",
                "local_pressure_law": "log",
            },
            "source_identity": {
                "app": {
                    "revision": "1" * 40,
                    "source_kind": "git-checkout",
                    "tree_state": "clean",
                },
                "core": {
                    "revision": "2" * 40,
                    "source_kind": "git-checkout",
                    "source_url": post.PUBLIC_CORE_URL,
                    "tree_state": "clean",
                },
            },
        },
        "schema": publisher.GENERIC_REPORT_SCHEMA,
    }


def _evidence(tmp_path, mutate_report=None, mutate_manifest=None):
    stdout = tmp_path / "synthetic-step0b.stdout.txt"
    stdout.write_text(
        "  step 1000 t=1.000s load=0.0Pa\nelapsed 12.5s\n"
        "saved -> /machine/path/synthetic-step0b.npz\n",
        encoding="utf-8",
    )
    log_sha, log_size = _identity(stdout)
    normalized_log = {
        "filename": stdout.name,
        "normalization": "UTF-8; CRLF/CR converted to LF; final LF added when nonempty",
        "sha256": log_sha,
        "size_bytes": log_size,
    }
    report = _generic_report(normalized_log)
    if mutate_report is not None:
        mutate_report(report)
    generic = tmp_path / "synthetic-step0b.report.json"
    _write_json(generic, report)
    report_sha, report_size = _identity(generic)
    mesh = report["result"]["configuration"]["mesh"]
    manifest = {
        "schema": publisher.CAMPAIGN_MANIFEST_SCHEMA,
        "application": {"revision": "1" * 40, "tree_state": "clean"},
        "core": {"revision": "2" * 40, "tree_state": "clean"},
        "attempts": {
            "full": {
                "completed_steps": 1000,
                "expected_steps": 1000,
                "elapsed_s": 12.5,
                "log": stdout.name,
                "log_sha256": log_sha,
                "log_size_bytes": log_size,
                "mesh": {
                    key: mesh[key]
                    for key in (
                        "n_t",
                        "n_core",
                        "n_radial",
                        "elements",
                        "nodes",
                        "degrees_of_freedom",
                    )
                },
                "mpi_ranks": 8,
                "output": "synthetic-step0b.npz",
                "output_sha256": report["result"]["sha256"],
                "output_size_bytes": report["result"]["size_bytes"],
                "report": generic.name,
                "report_sha256": report_sha,
                "report_size_bytes": report_size,
                "status": "completed_and_validated",
                "t_end_s": 1.0,
            }
        },
        "configuration": {
            "benchmark": 1,
            "case": "B",
            "step": 0,
            "dt_s": 0.001,
            "integrator": "generalized-alpha",
            "mass": "consistent",
            "formulation": "local-pressure",
            "mesh_topology": "closed-multiblock",
        },
        "runtime": {
            "selected_mpi_ranks": 8,
            "thread_caps": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
        },
    }
    if mutate_manifest is not None:
        mutate_manifest(manifest)
    manifest_path = tmp_path / "manifest.json"
    _write_json(manifest_path, manifest)
    return generic, manifest_path, stdout


def test_compact_step0_publisher_retains_evidence_and_strips_bulk(tmp_path):
    generic, manifest, stdout = _evidence(tmp_path)
    report = publisher.build_compact_report(generic, manifest, stdout)
    payload = publisher.encode_compact_report(report)

    assert set(report) == {
        "benchmark_identity",
        "bounded_claim",
        "comparison",
        "reference",
        "result",
        "schema",
    }
    assert report["schema"] == publisher.COMPACT_REPORT_SCHEMA
    assert report["benchmark_identity"]["status"] == "recorded"
    assert report["benchmark_identity"]["case"] == "step_0B"
    assert report["result"]["configuration"]["mesh"]["topology"] == (
        "closed_multiblock_disk"
    )
    assert report["result"]["configuration"]["mpi"]["ranks"] == 8
    assert len(payload) < 100_000
    text = payload.decode("utf-8")
    for omitted in (
        "nonlinear_step_diagnostics",
        "retained_histories",
        "load_evaluation_times_s",
        "ksp_residual_histories",
    ):
        assert omitted not in text


def test_compact_step0_publisher_is_byte_deterministic(tmp_path):
    generic, manifest, stdout = _evidence(tmp_path)
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    publisher.publish_compact_report(generic, manifest, stdout, first)
    publisher.publish_compact_report(generic, manifest, stdout, second)
    assert first.read_bytes() == second.read_bytes()


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda report: report["result"]["configuration"]["benchmark"].__setitem__(
                "configuration_id", "wrong"
            ),
            "recorded Step 0 identity",
        ),
        (
            lambda report: report["result"]["configuration"]["mesh"].__setitem__(
                "topology", "polar_ring"
            ),
            "closed multiblock mesh",
        ),
        (
            lambda report: report["result"]["nonlinear_step_diagnostics"][20].__setitem__(
                "snes_converged_reason", -3
            ),
            "accepted-step diagnostic 21",
        ),
        (
            lambda report: report["reference"]["canonical_grid_s"].__setitem__(5, 0.2),
            "canonical grid differs",
        ),
        (
            lambda report: report["comparison"]["red"]["p0"].__setitem__(
                "ours", 99.0
            ),
            "p0 RED is not reproducible",
        ),
    ],
)
def test_compact_step0_publisher_rejects_broken_report_controls(
    tmp_path, mutation, message
):
    generic, manifest, stdout = _evidence(tmp_path, mutate_report=mutation)
    with pytest.raises(publisher.PublicationError, match=message):
        publisher.build_compact_report(generic, manifest, stdout)


def test_compact_step0_publisher_rejects_broken_manifest_identity(tmp_path):
    generic, manifest, stdout = _evidence(
        tmp_path,
        mutate_manifest=lambda value: value["attempts"]["full"].__setitem__(
            "output_sha256", "0" * 64
        ),
    )
    with pytest.raises(publisher.PublicationError, match="manifest attempt identity"):
        publisher.build_compact_report(generic, manifest, stdout)


def test_compact_step0_publisher_rejects_changed_raw_stdout(tmp_path):
    generic, manifest, stdout = _evidence(tmp_path)
    stdout.write_text(stdout.read_text(encoding="utf-8") + "changed\n", encoding="utf-8")
    with pytest.raises(publisher.PublicationError, match="manifest attempt identity"):
        publisher.build_compact_report(generic, manifest, stdout)
