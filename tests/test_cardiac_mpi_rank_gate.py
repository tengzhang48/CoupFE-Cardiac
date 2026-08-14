from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

import compare_mpi_rank_gate as gate


def _balanced_ranges(size, count):
    base, remainder = divmod(count, size)
    ranges = []
    start = 0
    for rank in range(size):
        stop = start + base + (1 if rank < remainder else 0)
        ranges.append((start, stop))
        start = stop
    return np.asarray(ranges, dtype=np.int64)


def _diagnostics(times, ranks=None):
    records = []
    for index in range(1, len(times)):
        record = {
            "time": float(times[index]),
            "dt": 0.001,
            "initial_residual_norm": 1.0e-3,
            "final_residual_norm": 1.0e-12,
            "residual_acceptance_threshold": 1.0e-10,
            "snes_converged_reason": 2,
            "ksp_converged_reason": 4,
            "nonlinear_iterations": 2,
            "linear_iterations": 2,
            "function_domain_rejections": 0,
        }
        if ranks is not None:
            record["ranks"] = ranks
        records.append(record)
    return records


def _audit(nodes, elements):
    return {
        "geometry": {
            "schema": "coupfe-cardiac-pre-solve-geometry-v1",
            "passed": True,
            "require_closed": True,
            "mesh_topology": "closed_multiblock_disk",
            "nodes": nodes,
            "elements": elements,
            "unclassified_exterior_faces": 0,
            "multiply_labeled_faces": 0,
            "nonmanifold_faces": 0,
            "nonpositive_extended_jacobians": 0,
        },
        "pressure": {
            "schema": "coupfe-cardiac-pre-solve-pressure-v1",
            "passed": True,
        },
        "robin": {
            "schema": "coupfe-cardiac-pre-solve-robin-v1",
            "passed": True,
        },
    }


def _serial_solver_configuration():
    return {
        "name": "petsc-snes",
        "snes_type": "newtonls",
        "line_search_type": "bt",
        "ksp_type": "preonly",
        "pc_type": "lu",
        "function_domain_rejection_api": "nonfinite residual for PETSc BT",
        "rtol": 1.0e-9,
        "atol": 1.0e-10,
        "stol": 1.0e-12,
        "max_it": 60,
        "dirichlet_support": "none",
        "element_evaluation_mode": "joint",
        "compiled_material_residual_only_available": True,
    }


def _mpi_solver_configuration(ranks, ranges, nnz):
    return {
        "name": "petsc-snes-mpi",
        "implementation": gate.MPI_IMPLEMENTATION,
        "communicator": "PETSc.COMM_WORLD",
        "ranks": ranks,
        "rank_local_element_assembly": True,
        "global_mesh_replicated": True,
        "collective_invalid_trial_policy": "all-owned-residual-entries-inf-for-bt",
        "commit_policy": "independent-global-residual-check-before-commit",
        "function_domain_rejection_api": "nonfinite residual for PETSc BT",
        "snes_type": "newtonls",
        "line_search_type": "bt",
        "ksp_type": "preonly",
        "pc_type": "lu",
        "rtol": 1.0e-9,
        "atol": 1.0e-10,
        "stol": 1.0e-12,
        "max_it": 60,
        "dirichlet_support": "none",
        "factor_solver_type": "superlu_dist",
        "configured_factor_solver_type": "superlu_dist",
        "mass_representation": "consistent_q1_hex8",
        "mass_partition": "owned-row-csr-all-touching-elements",
        "mass_owned_row_range": ranges[0].astype(int).tolist(),
        "mass_local_nnz": int(nnz[0]),
        "line_search_configuration_api": "SNES.getLineSearch",
        "petsc4py_version": "3.18.4",
        "petsc_version": "3.18.4",
        "element_evaluation_mode": "joint",
        "compiled_material_residual_only_available": True,
    }


def _payload(role):
    n_t = gate.EXPECTED_N_T
    n_core = gate.EXPECTED_N_CORE
    n_radial = gate.EXPECTED_N_RADIAL
    n_element = n_t * (n_core * n_core + 4 * n_core * n_radial)
    n_node = (n_t + 1) * ((n_core + 1) ** 2 + 4 * n_core * n_radial)
    ndof = 3 * n_node
    times = gate.EXPECTED_TIMES.copy()
    node_index = np.arange(n_node, dtype=float)
    nodes = np.column_stack(
        (node_index * 1.0e-5, node_index * 2.0e-5, node_index * -1.0e-5)
    )
    elems = (np.arange(n_element * 8, dtype=np.int64) % n_node).reshape(n_element, 8)
    u0 = np.column_stack((times * 0.01, times * -0.02, times * 0.03))
    u1 = np.column_stack((times * -0.015, times * 0.005, times * -0.01))
    peak = np.linspace(-0.02, 0.03, ndof)
    det_f = np.full((n_element, 8), 1.01)
    rank_delta = {"serial": 0.0, "mpi1": 2.0e-14, "mpi2": -3.0e-14, "mpi4": 4.0e-14}[role]
    u0 = u0.copy()
    u0[10, 0] += rank_delta
    u1 = u1.copy()
    u1[20, 1] -= rank_delta
    peak = peak.copy()
    peak[3] += rank_delta
    det_f = det_f.copy()
    det_f[2, 4] -= rank_delta

    payload = {
        "completed_steps": 320,
        "expected_steps": 320,
        "converged": True,
        "result_schema": gate.RESULT_SCHEMA,
        "driver": "examples/cardiac_benchmark/run.py",
        "app_revision": "a" * 40,
        "app_tree_state": "clean",
        "app_source_kind": "git-checkout",
        "core_revision": "b" * 40,
        "core_tree_state": "clean",
        "core_source_kind": "git-checkout",
        "core_source_url": gate.PUBLIC_CORE_URL,
        "coupfe_version": "0.0.1",
        "numpy_version": "1.26.4",
        "scipy_version": "1.15.2",
        "python_version": "3.10.8",
        "times": times,
        "tau": np.zeros_like(times),
        "pres": 16000.0 * np.sin(np.pi * times) ** 2,
        "u0": u0,
        "u1": u1,
        "p0": np.array([0.025, 0.030, 0.0]),
        "p1": np.array([0.000, 0.030, 0.0]),
        "U_peak": peak,
        "n_peak": 170,
        "nodes": nodes,
        "elems": elems,
        "fiber": np.tile(np.array([[1.0, 0.0, 0.0]]), (n_element, 1)),
        "facets_endo": np.arange(16 * 4, dtype=np.int64).reshape(16, 4) % n_node,
        "case": "B",
        "dt": 0.001,
        "integrator": "be",
        "formulation": "hex8_standard_pointwise_kappa",
        "density": 1000.0,
        "material_kernel_formulation": "standard",
        "material_model_id": gate.MATERIAL_MODEL_ID,
        "material_kappa_pa": 1.0e6,
        "local_pressure_bulk_modulus_pa": 0.0,
        "isotropic": False,
        "mass_representation": "consistent_q1_hex8",
        "det_f_gauss_peak": det_f,
        "element_pressure_peak_pa": np.empty(0),
        "element_evaluation_mode": "joint",
        "compiled_material_residual_only_available": True,
        "nonlinear_solver": "petsc-snes",
        "pre_solve_audit_json": json.dumps(_audit(n_node, n_element), sort_keys=True),
        "fiber_sampling": "gp_direct_rule",
        "fiber_sampling_option": "gp-direct",
        "tbar_definition": "laplace_presolved",
        "tbar_source_filename": "tbar_closed.npy",
        "tbar_source_sha256": "c" * 64,
        "tbar_metadata_filename": "tbar_closed.meta.json",
        "tbar_metadata_sha256": "d" * 64,
        "tbar_metadata_schema": gate.TBAR_SCHEMA,
        "point_sampling": "hex8_reference_isoparametric",
        "p0_sampling_element": 20,
        "p0_sampling_natural": np.array([0.1, -0.2, 0.3]),
        "p0_sampling_weights": np.arange(8, dtype=float) / 28.0,
        "p0_sampling_reconstruction_error_m": 1.0e-16,
        "p1_sampling_element": 21,
        "p1_sampling_natural": np.array([-0.1, 0.2, -0.3]),
        "p1_sampling_weights": np.arange(7, -1, -1, dtype=float) / 28.0,
        "p1_sampling_reconstruction_error_m": 2.0e-16,
        "viscous_rate": "backward_difference",
        "material_eta_pa_s": 100.0,
        "viscous_term_active": True,
        "parameter_variant": "benchmark_eta",
        "mesh_topology": "closed_multiblock_disk",
        "n_side": 0,
        "n_core": n_core,
        "n_radial": n_radial,
        "core_half_width": 0.36,
        "apex_offset": 0.0,
        "flip_helix": True,
        "perturb": 0.0,
        "t_end": 0.32,
        "load_horizon": 1.0,
        "a_top": 1.0e5,
        "b_top": 5.0e3,
        "a_epi": 1.0e8,
        "b_epi": 5.0e3,
        "n_t": n_t,
        "n_mu": 0,
        "n_theta": 0,
    }
    ranks = gate.EXPECTED_MPI_RANKS.get(role)
    if ranks is None:
        payload["solver_configuration_json"] = json.dumps(
            _serial_solver_configuration(), sort_keys=True
        )
        payload["nonlinear_step_diagnostics_json"] = json.dumps(
            _diagnostics(times), sort_keys=True
        )
    else:
        ranges = _balanced_ranges(ranks, ndof)
        counts = _balanced_ranges(ranks, n_element)
        local_counts = counts[:, 1] - counts[:, 0]
        nnz = (ranges[:, 1] - ranges[:, 0]) * 12
        touching = np.maximum(local_counts, 1)
        payload.update(
            {
                "driver": "examples/cardiac_benchmark/run_mpi.py",
                "nonlinear_solver": "petsc-snes-mpi",
                "mpi_enabled": True,
                "mpi_ranks": ranks,
                "mpi_world_size": ranks,
                "mpi_implementation": gate.MPI_IMPLEMENTATION,
                "mpi_local_element_counts": local_counts,
                "mpi_partition": "coupfe.partition_elements",
                "mpi_build_layout": "isolated-rank-directories",
                "mpi_factor_solver_type": "superlu_dist",
                "mpi_mass_partition": "owned-row-csr-all-touching-elements",
                "mpi_mass_owned_row_ranges": ranges,
                "mpi_mass_local_nnz": nnz,
                "mpi_mass_touching_element_counts": touching,
                "solver_configuration_json": json.dumps(
                    _mpi_solver_configuration(ranks, ranges, nnz), sort_keys=True
                ),
                "nonlinear_step_diagnostics_json": json.dumps(
                    _diagnostics(times, ranks=ranks), sort_keys=True
                ),
            }
        )
    return payload


def _write_archives(root: Path):
    private = root / "private-volume" / "case-b"
    private.mkdir(parents=True)
    paths = {}
    for role in gate.INPUT_ROLES:
        path = private / (role + ".npz")
        np.savez(path, **_payload(role))
        paths[role] = path
    return paths


def _rewrite(path, mutate):
    with np.load(path, allow_pickle=False) as archive:
        payload = {key: archive[key].copy() for key in archive.files}
    mutate(payload)
    np.savez(path, **payload)


def _arguments(paths, report):
    return [
        "--serial",
        str(paths["serial"]),
        "--mpi1",
        str(paths["mpi1"]),
        "--mpi2",
        str(paths["mpi2"]),
        "--mpi4",
        str(paths["mpi4"]),
        "--report",
        str(report),
    ]


def test_rank_gate_writes_path_free_atomic_pass_record(tmp_path, capsys):
    paths = _write_archives(tmp_path)
    report_path = tmp_path / "public" / "rank-gate.json"
    report = gate.main(_arguments(paths, report_path))

    assert report["status"] == "passed"
    assert report["comparison_definition"]["rtol"] == 2.0e-11
    assert report["comparison_definition"]["atol"] == 2.0e-13
    assert report["comparison_definition"]["tolerances_are_cli_configurable"] is False
    assert report["source_identity"]["app"]["revision"] == "a" * 40
    for role, ranks in gate.EXPECTED_MPI_RANKS.items():
        assert report["comparisons"][role]["passed"] is True
        assert report["runs"][role]["solver"]["mpi"]["ranks"] == ranks
    text = report_path.read_text(encoding="utf-8")
    assert str(tmp_path) not in text
    assert "private-volume" not in text
    assert json.loads(text)["schema"] == gate.REPORT_SCHEMA
    assert not list(report_path.parent.glob(".*.tmp"))
    stdout = capsys.readouterr().out
    assert str(tmp_path) not in stdout
    assert "private-volume" not in stdout
    assert "rank-gate.json" in stdout


@pytest.mark.parametrize(
    "role,mutation,match",
    [
        (
            "mpi2",
            lambda payload: payload.__setitem__(
                "U_peak", np.asarray(payload["U_peak"]) + 1.0e-8
            ),
            "fixed numerical rank-gate tolerance",
        ),
        (
            "mpi4",
            lambda payload: np.asarray(payload["pres"]).__setitem__(12, 7.0),
            "invariant field pres",
        ),
        (
            "mpi1",
            lambda payload: payload.pop("mpi_mass_touching_element_counts"),
            "missing 'mpi_mass_touching_element_counts'",
        ),
        (
            "serial",
            lambda payload: payload.__setitem__("app_tree_state", "dirty"),
            "source is not clean",
        ),
        (
            "mpi1",
            lambda payload: payload.__setitem__(
                "tbar_source_filename", "/private/field.npy"
            ),
            "not a portable basename",
        ),
    ],
)
def test_broken_controls_fail_without_report(tmp_path, role, mutation, match):
    paths = _write_archives(tmp_path)
    _rewrite(paths[role], mutation)
    report_path = tmp_path / "rank-gate.json"
    with pytest.raises(gate.RankGateInputError, match=match):
        gate.main(_arguments(paths, report_path))
    assert not report_path.exists()


def test_pickle_dependent_archive_is_rejected_without_report(tmp_path):
    paths = _write_archives(tmp_path)
    _rewrite(
        paths["mpi1"],
        lambda payload: payload.__setitem__(
            "u0", np.asarray(payload["u0"], dtype=object)
        ),
    )
    report_path = tmp_path / "rank-gate.json"
    with pytest.raises(gate.RankGateInputError, match="without pickle"):
        gate.main(_arguments(paths, report_path))
    assert not report_path.exists()


def test_explicit_mesh_mismatch_fails_closed(tmp_path):
    paths = _write_archives(tmp_path)
    _rewrite(
        paths["mpi2"],
        lambda payload: payload.__setitem__("n_core", gate.EXPECTED_N_CORE + 2),
    )
    with pytest.raises(gate.RankGateInputError, match="does not match the explicit mesh"):
        gate.main(_arguments(paths, tmp_path / "rank-gate.json"))
    assert not (tmp_path / "rank-gate.json").exists()
