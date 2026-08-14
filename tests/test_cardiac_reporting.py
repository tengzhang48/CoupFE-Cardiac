from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import pickle
from pathlib import Path
import zipfile

import numpy as np
import pytest

import diagnose
import post


CORE_REVISION = "454f73ce2de284262b214a2b37bd676c6aca3c0a"
ROOT = Path(__file__).resolve().parents[1]
RELEASE_GUARD_PATH = ROOT / ".github" / "scripts" / "check_release_artifacts.py"
_RELEASE_GUARD_SPEC = importlib.util.spec_from_file_location(
    "cardiac_release_guard", RELEASE_GUARD_PATH
)
assert _RELEASE_GUARD_SPEC is not None and _RELEASE_GUARD_SPEC.loader is not None
RELEASE_GUARD = importlib.util.module_from_spec(_RELEASE_GUARD_SPEC)
_RELEASE_GUARD_SPEC.loader.exec_module(RELEASE_GUARD)


def _release_guard_payloads(spec):
    report_name = RELEASE_GUARD._truncated_polar_archive_path(spec["report"])
    log_name = RELEASE_GUARD._truncated_polar_archive_path(spec["log"])
    return {
        report_name: (ROOT / report_name).read_bytes(),
        log_name: (ROOT / log_name).read_bytes(),
    }


def _reviewed_spec(report_basename):
    return next(
        spec
        for spec in RELEASE_GUARD.TRUNCATED_POLAR_ARCHIVE_REPORT_SPECS
        if spec["report"] == report_basename
    )


def _core_diagnostics(times):
    return [
        {
            "time": float(times[index]),
            "dt": float(times[index] - times[index - 1]),
            "nonlinear_iterations": index,
        }
        for index in range(1, len(times))
    ]


def _petsc_diagnostics(times, *, domain_rejections=False, ranks=None):
    diagnostics = [
        {
            "time": float(times[index]),
            "dt": float(times[index] - times[index - 1]),
            "initial_residual_norm": 10.0,
            "final_residual_norm": 9.876543210123456e-11,
            "residual_acceptance_threshold": 1.0e-8,
            "petsc_function_norm": 9.876543210123456e-11,
            "snes_converged_reason": 3,
            "ksp_converged_reason": 4,
            "nonlinear_iterations": 2,
            "linear_iterations": 2,
            "residual_history": [10.0, 0.25, 9.876543210123456e-11],
            "assembly_seconds": 0.12345678901234566,
            "solve_seconds": 0.23456789012345677,
        }
        for index in range(1, len(times))
    ]
    if domain_rejections:
        diagnostics[0].update(
            {
                "function_domain_rejections": 1,
                "last_function_domain_error": "invalid trial det(F)",
            }
        )
        for record in diagnostics[1:]:
            record.update(
                {
                    "function_domain_rejections": 0,
                    "last_function_domain_error": None,
                }
            )
    if ranks is not None:
        for record in diagnostics:
            record["ranks"] = ranks
    return diagnostics


def _closed_pre_solve_audit(*, n_node=8, n_elem=1, include_pressure=True):
    audit = {
        "geometry": {
            "schema": "coupfe-cardiac-pre-solve-geometry-v1",
            "mesh_topology": "closed_multiblock_disk",
            "require_closed": True,
            "nodes": n_node,
            "elements": n_elem,
            "exterior_faces": 6,
            "labeled_exterior_faces": 6,
            "unclassified_exterior_faces": 0,
            "nonexterior_labeled_faces": 0,
            "multiply_labeled_faces": 0,
            "nonmanifold_faces": 0,
            "nonpositive_gauss_jacobians": 0,
            "nonpositive_extended_jacobians": 0,
            "gauss_jacobian_min_m3": 1.0e-9,
            "extended_jacobian_min_m3": 0.9e-9,
            "extended_scaled_jacobian_min": 0.25,
            "passed": True,
            "failures": [],
        },
        "robin": {
            "schema": "coupfe-cardiac-pre-solve-robin-v1",
            "active_dofs": 12,
            "spring_symmetry_error": 1.0e-13,
            "dashpot_symmetry_error": 1.0e-18,
            "passed": True,
            "failures": [],
        },
    }
    if include_pressure:
        projected_area = 1.8e-3
        audit["pressure"] = {
            "schema": "coupfe-cardiac-pre-solve-pressure-v1",
            "unit_pressure_resultant_N": [projected_area, 0.0, 0.0],
            "expected_unit_pressure_resultant_N": [projected_area, 0.0, 0.0],
            "analytic_projected_base_area_m2": projected_area,
            "relative_magnitude_error": 0.0,
            "signed_axial_ratio": 1.0,
            "relative_signed_axial_error": 0.0,
            "relative_resultant_error": 0.0,
            "transverse_fraction": 0.0,
            "unit_pressure_moment_Nm": [0.0, 0.0, 0.0],
            "normalized_moment": 0.0,
            "passed": True,
            "failures": [],
        }
    return audit


def _mpi_linear_solver_configuration(profile):
    if profile == post.MPI_DIRECT_SUPERLU_DIST_PROFILE:
        return {
            "factor_solver_type": "superlu_dist",
            "configured_factor_solver_type": "superlu_dist",
            "linear_solver_profile": profile,
            "node_aligned_ownership": False,
            "vector_block_size": 1,
            "matrix_block_size": 1,
            "near_nullspace_kind": "none",
            "near_nullspace_mode_count": 0,
            "ksp_type": "preonly",
            "pc_type": "lu",
            "ksp_rtol": None,
            "ksp_atol": None,
            "ksp_divtol": None,
            "ksp_max_it": None,
            "gmres_restart": None,
            "ksp_norm_type": "default",
            "pc_side": "default",
            "preconditioner": "global LU with SuperLU_DIST",
        }

    configuration = {
        "factor_solver_type": "none",
        "configured_factor_solver_type": "not-applicable",
        "linear_solver_profile": profile,
        "node_aligned_ownership": False,
        "vector_block_size": 1,
        "matrix_block_size": 1,
        "near_nullspace_kind": "none",
        "near_nullspace_mode_count": 0,
        "ksp_type": "fgmres",
        "ksp_rtol": 1.0e-8,
        "ksp_atol": 1.0e-12,
        "ksp_divtol": 1.0e4,
        "ksp_max_it": 200,
        "gmres_restart": 50,
        "ksp_norm_type": "unpreconditioned",
        "pc_side": "right",
        "ksp_error_if_not_converged": False,
    }
    option_prefix = "coupfe_cardiac_deadbeef_"
    if profile in {
        post.MPI_FGMRES_GAMG_RIGID_PROFILE,
        post.MPI_FGMRES_GAMG_RIGID_REBUILD_PROFILE,
    }:
        configuration.update(
            {
                "node_aligned_ownership": True,
                "vector_block_size": 3,
                "matrix_block_size": 3,
                "near_nullspace_kind": "six-rigid-body-modes",
                "near_nullspace_mode_count": 6,
                "pc_type": "gamg",
            }
        )
        options = {
            "pc_gamg_type": "agg",
            "pc_gamg_agg_nsmooths": "1",
            "pc_gamg_threshold": "0.01",
            "pc_gamg_repartition": "false",
        }
        if profile == post.MPI_FGMRES_GAMG_RIGID_REBUILD_PROFILE:
            configuration.update(
                {
                    "ksp_max_it": 400,
                    "gmres_restart": 100,
                    "preconditioner": (
                        "PETSc GAMG aggregation with six rigid-body near-null "
                        "modes and interpolation rebuilt for changed matrices"
                    ),
                }
            )
            options["pc_gamg_reuse_interpolation"] = "false"
        else:
            configuration["preconditioner"] = (
                "PETSc GAMG aggregation with six rigid-body near-null modes"
            )
    elif profile == post.MPI_FGMRES_ASM_LU_PROFILE:
        configuration.update(
            {
                "pc_type": "asm",
                "preconditioner": (
                    "restricted additive Schwarz overlap 1 with local SuperLU"
                ),
            }
        )
        options = {
            "sub_ksp_type": "preonly",
            "sub_pc_type": "lu",
            "sub_pc_factor_mat_solver_type": "superlu",
        }
    elif profile == post.MPI_FGMRES_ASM_ILU1_PROFILE:
        configuration.update(
            {
                "pc_type": "asm",
                "preconditioner": (
                    "restricted additive Schwarz overlap 1 with local ILU(1)"
                ),
            }
        )
        options = {
            "sub_ksp_type": "preonly",
            "sub_pc_type": "ilu",
            "sub_pc_factor_levels": "1",
            "sub_pc_factor_shift_type": "nonzero",
        }
    else:  # pragma: no cover - fixture programming error
        raise AssertionError(profile)
    configuration["petsc_options"] = {
        f"{option_prefix}{name}": value for name, value in options.items()
    }
    return configuration


def _write_current_result(
    path,
    *,
    case="B",
    formulation="hex8_fbar",
    fiber_sampling="cg1_gram_schmidt",
    fiber_sampling_option=None,
    closed=False,
    mass_representation="consistent_q1_hex8",
    material_eta_pa_s=100.0,
    isotropic=False,
    tbar_definition="analytic_parametric",
    tbar_source_filename="",
    tbar_source_sha256="",
    tbar_metadata_filename="",
    tbar_metadata_sha256="",
    tbar_metadata_schema="",
    solver="core-newton",
    diagnostics=None,
    domain_rejections=False,
    mpi_ranks=1,
    linear_solver_profile=post.MPI_DIRECT_SUPERLU_DIST_PROFILE,
    history_value=0.12345678901234566,
    load_horizon=1.0,
):
    times = np.array([0.0, 0.5, 1.0])
    history = np.array(
        [
            [0.0, 0.0, 0.0],
            [history_value, -0.25, 0.5],
            [0.25, -0.5, 1.0],
        ]
    )
    if diagnostics is None:
        if solver == "petsc-snes":
            diagnostics = _petsc_diagnostics(
                times, domain_rejections=domain_rejections
            )
        elif solver == "petsc-snes-mpi":
            diagnostics = _petsc_diagnostics(
                times, domain_rejections=True, ranks=mpi_ranks
            )
        else:
            diagnostics = _core_diagnostics(times)
    if solver in {"petsc-snes", "petsc-snes-mpi"}:
        configuration = {
            "name": solver,
            "snes_type": "newtonls",
            "line_search_type": "bt",
            "ksp_type": "preonly",
            "pc_type": "lu",
            "rtol": 1.0e-9,
            "atol": 1.0e-10,
            "stol": 1.0e-12,
            "max_it": 60,
        }
        if domain_rejections or solver == "petsc-snes-mpi":
            configuration["function_domain_rejection_api"] = (
                post.PETSC_FUNCTION_DOMAIN_REJECTION_API
            )
        if solver == "petsc-snes-mpi":
            configuration.update(
                {
                    **post.MPI_COMPANION_FIXED_CONFIGURATION,
                    "name": "petsc-snes-mpi",
                    "ranks": mpi_ranks,
                    "line_search_configuration_api": "SNES.getLineSearch",
                    "element_evaluation_mode": "joint",
                    "compiled_material_residual_only_available": True,
                    "petsc4py_version": "3.22.4",
                    "petsc_version": "3.22.4",
                    **_mpi_linear_solver_configuration(
                        linear_solver_profile
                    ),
                }
            )
    else:
        configuration = {
            "name": "core-newton",
            "rtol": 1.0e-8,
            "max_it": 40,
        }
    if fiber_sampling_option is None:
        fiber_sampling_option = {
            "cg1_gram_schmidt": "cg1",
            "gp_direct_rule": "gp-direct",
        }[fiber_sampling]
    local = formulation in post.LOCAL_PRESSURE_FORMULATIONS
    standard_kappa = formulation == "hex8_standard_pointwise_kappa"
    payload = {
        "result_schema": post.RESULT_SCHEMA,
        "case": case,
        "integrator": "be",
        "formulation": formulation,
        "fiber_sampling": fiber_sampling,
        "fiber_sampling_option": fiber_sampling_option,
        "tbar_definition": tbar_definition,
        "tbar_source_filename": tbar_source_filename,
        "tbar_source_sha256": tbar_source_sha256,
        "tbar_metadata_filename": tbar_metadata_filename,
        "tbar_metadata_sha256": tbar_metadata_sha256,
        "tbar_metadata_schema": tbar_metadata_schema,
        "point_sampling": "hex8_reference_isoparametric",
        "viscous_rate": "backward_difference",
        "viscous_term_active": material_eta_pa_s > 0.0,
        "material_eta_pa_s": material_eta_pa_s,
        "parameter_variant": (
            "benchmark_eta"
            if material_eta_pa_s == 100.0
            else "eta_zero_sensitivity"
            if material_eta_pa_s == 0.0
            else "eta_sensitivity"
        ),
        "mass_representation": mass_representation,
        "isotropic": isotropic,
        "material_model_id": post.COMPLETE_SWITCH_ENERGY_MATERIAL_MODEL_ID,
        "dt": 0.5,
        "t_end": 1.0,
        "load_horizon": load_horizon,
        "apex_offset": 0.2,
        "n_t": 1,
        "n_mu": 1,
        "n_theta": 1,
        "flip_helix": True,
        "density": 1000.0,
        "a_top": 1.0e5,
        "b_top": 5.0e3,
        "a_epi": 1.0e8,
        "b_epi": 5.0e3,
        "perturb": 0.0,
        "python_version": "3.10.8",
        "numpy_version": "1.26.4",
        "scipy_version": "1.15.2",
        "coupfe_version": "0.0.1",
        "converged": True,
        "completed_steps": 2,
        "expected_steps": 2,
        "app_revision": "a" * 40,
        "app_tree_state": "clean",
        "app_source_kind": "git-checkout",
        "core_revision": CORE_REVISION,
        "core_tree_state": "clean",
        "core_source_kind": "git-checkout",
        "core_source_url": post.PUBLIC_CORE_URL,
        "times": times,
        "tau": np.zeros(3),
        "pres": np.array([0.0, 1.0, 0.0]),
        "u0": history,
        "u1": history * 0.5,
        "n_peak": 1,
        "nodes": np.array(
            [
                [-1.0, -1.0, -1.0],
                [1.0, -1.0, -1.0],
                [1.0, 1.0, -1.0],
                [-1.0, 1.0, -1.0],
                [-1.0, -1.0, 1.0],
                [1.0, -1.0, 1.0],
                [1.0, 1.0, 1.0],
                [-1.0, 1.0, 1.0],
            ]
        ),
        "elems": np.arange(8, dtype=int).reshape(1, 8),
        "U_peak": np.zeros(24),
        "nonlinear_solver": solver,
        "solver_configuration_json": json.dumps(configuration),
        "nonlinear_step_diagnostics_json": json.dumps(diagnostics),
        "det_f_gauss_peak": np.full((1, 8), 1.0123456789012346),
        "material_kernel_formulation": (
            "standard" if local or standard_kappa else "fbar_mechanics"
        ),
        "material_kappa_pa": 0.0 if local else 1.0e6,
        "local_pressure_bulk_modulus_pa": 1.0e6 if local else 0.0,
        "local_pressure_volume_law": (
            "paper-j2-of-reference-volume-weighted-geometric-mean-j-v1"
            if formulation
            == "hex8_local_pressure_p0_condensed_mean_logj_paper_j2"
            else "linear-reference-volume-mean-log-j-v1"
            if local
            else "not-applicable"
        ),
    }
    if closed:
        payload.update(
            {
                "mesh_topology": "closed_multiblock_disk",
                "n_mu": 0,
                "n_theta": 0,
                "n_side": 0,
                "n_core": 4,
                "n_radial": 1,
                "core_half_width": 0.36,
                "apex_offset": 0.0,
                "pre_solve_audit_json": json.dumps(
                    _closed_pre_solve_audit(include_pressure=case == "B"),
                    sort_keys=True,
                    separators=(",", ":"),
                ),
            }
        )
    for point in ("p0", "p1"):
        payload[f"{point}_sampling_element"] = 0
        payload[f"{point}_sampling_natural"] = np.zeros(3)
        payload[f"{point}_sampling_weights"] = np.full(8, 0.125)
        payload[f"{point}_sampling_reconstruction_error_m"] = 1.0e-16
    payload["element_pressure_peak_pa"] = (
        np.array([1234.5678901234567]) if local else np.empty(0)
    )
    if solver == "petsc-snes-mpi":
        quotient, remainder = divmod(len(payload["elems"]), mpi_ranks)
        local_counts = np.full(mpi_ranks, quotient, dtype=np.int64)
        local_counts[:remainder] += 1
        payload.update(
            {
                "driver": "examples/cardiac_benchmark/run_mpi.py",
                "mpi_enabled": True,
                "mpi_ranks": mpi_ranks,
                "mpi_world_size": mpi_ranks,
                "mpi_local_element_counts": local_counts,
                "mpi_implementation": post.MPI_COMPANION_IMPLEMENTATION,
                "mpi_partition": "coupfe.partition_elements",
                "mpi_build_layout": "isolated-rank-directories",
                "mpi_factor_solver_type": (
                    "superlu_dist"
                    if linear_solver_profile
                    == post.MPI_DIRECT_SUPERLU_DIST_PROFILE
                    else "none"
                ),
                "mpi_linear_solver_profile": linear_solver_profile,
                "element_evaluation_mode": "joint",
                "compiled_material_residual_only_available": True,
            }
        )
    np.savez(path, **payload)
    return payload


def _write_mpi_result(
    path,
    *,
    mpi_ranks=2,
    linear_solver_profile=post.MPI_DIRECT_SUPERLU_DIST_PROFILE,
):
    return _write_current_result(
        path,
        formulation="hex8_local_pressure_p0_condensed_logj",
        mass_representation="lumped_row_sum",
        solver="petsc-snes-mpi",
        mpi_ranks=mpi_ranks,
        linear_solver_profile=linear_solver_profile,
    )


def _downgrade_to_legacy_direct_profile(payload):
    payload.pop("mpi_linear_solver_profile")
    configuration = json.loads(str(payload["solver_configuration_json"]))
    for field in (
        "linear_solver_profile",
        "node_aligned_ownership",
        "vector_block_size",
        "matrix_block_size",
        "near_nullspace_kind",
        "near_nullspace_mode_count",
        "ksp_rtol",
        "ksp_atol",
        "ksp_divtol",
        "ksp_max_it",
        "gmres_restart",
        "ksp_norm_type",
        "pc_side",
        "preconditioner",
    ):
        configuration.pop(field)
    payload["solver_configuration_json"] = json.dumps(configuration)


def _write_closed_mpi_result(
    path, *, mpi_ranks=2, case="B", formulation="std-kappa"
):
    formulation_label = {
        "std-kappa": "hex8_standard_pointwise_kappa",
        "local-pressure": "hex8_local_pressure_p0_condensed_logj",
        "local-pressure-paper": (
            "hex8_local_pressure_p0_condensed_mean_logj_paper_j2"
        ),
    }[formulation]
    implementation = {
        "std-kappa": post.MPI_CLOSED_STD_KAPPA_IMPLEMENTATION,
        "local-pressure": post.MPI_CLOSED_LOCAL_PRESSURE_IMPLEMENTATION,
        "local-pressure-paper": (
            post.MPI_CLOSED_LOCAL_PRESSURE_PAPER_IMPLEMENTATION
        ),
    }[formulation]
    payload = _write_current_result(
        path,
        case=case,
        formulation=formulation_label,
        fiber_sampling="gp_direct_rule",
        closed=True,
        mass_representation="consistent_q1_hex8",
        tbar_definition="laplace_presolved",
        tbar_source_filename="closed_tbar.npy",
        tbar_source_sha256=hashlib.sha256(b"closed tbar field").hexdigest(),
        tbar_metadata_filename="closed_tbar.meta.json",
        tbar_metadata_sha256=hashlib.sha256(
            b"closed tbar metadata"
        ).hexdigest(),
        tbar_metadata_schema=post.TBAR_METADATA_SCHEMA,
        solver="petsc-snes-mpi",
        mpi_ranks=mpi_ranks,
    )
    configuration = json.loads(str(payload["solver_configuration_json"]))
    quotient, remainder = divmod(24, mpi_ranks)
    widths = np.full(mpi_ranks, quotient, dtype=np.int64)
    widths[:remainder] += 1
    stops = np.cumsum(widths)
    starts = np.concatenate((np.array([0], dtype=np.int64), stops[:-1]))
    row_ranges = np.column_stack((starts, stops))
    local_nnz = np.full(mpi_ranks, 64, dtype=np.int64)
    configuration.update(
        {
            "implementation": implementation,
            "mass_representation": "consistent_q1_hex8",
            "mass_partition": "owned-row-csr-all-touching-elements",
            "mass_owned_row_range": row_ranges[0].astype(int).tolist(),
            "mass_local_nnz": int(local_nnz[0]),
            "local_pressure_law": (
                "paper"
                if formulation == "local-pressure-paper"
                else "log"
                if formulation == "local-pressure"
                else "not-applicable"
            ),
        }
    )
    payload.update(
        {
            "solver_configuration_json": json.dumps(configuration),
            "mpi_implementation": implementation,
            "mpi_mass_partition": "owned-row-csr-all-touching-elements",
            "mpi_mass_owned_row_ranges": row_ranges,
            "mpi_mass_local_nnz": local_nnz,
            "mpi_mass_touching_element_counts": np.ones(
                mpi_ranks, dtype=np.int64
            ),
        }
    )
    np.savez(path, **payload)
    return payload


def _write_generalized_alpha_closed_mpi_result(
    path,
    *,
    mpi_ranks=2,
    formulation="local-pressure",
    case="A",
    benchmark_step=None,
):
    payload = _write_closed_mpi_result(
        path,
        mpi_ranks=mpi_ranks,
        case=case,
        formulation=formulation,
    )
    implementation = {
        "std-kappa": post.MPI_CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION,
        "local-pressure": (
            post.MPI_CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION
        ),
        "local-pressure-paper": (
            post.MPI_CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION
        ),
    }[formulation]
    configuration = json.loads(str(payload["solver_configuration_json"]))
    alpha = {
        "alpha_m": 0.2,
        "alpha_f": 0.4,
        "gamma": 0.7,
        "beta": 0.36,
        "parameter_source": "finsberg/cardiac_benchmark problem.py defaults",
        "acceleration_stage": "alpha_m*a_n + (1-alpha_m)*a_np1",
        "force_stage": "alpha_f*x_n + (1-alpha_f)*x_np1",
        "load_time": "t_np1 - alpha_f*dt",
    }
    configuration.update(
        {
            "implementation": implementation,
            "time_integrator": "generalized-alpha",
            "material_batch_time_integrator": (
                "generalized-alpha-source-matched"
            ),
            "compiled_material_dt": float(payload["dt"]),
            "acceleration_stage": "1-alpha_m",
            "force_stage": "1-alpha_f",
            "material_viscous_rate": "sym(F_stage^T*grad(v_stage))",
            "nonlinear_initial_guess": "accepted-u_n-like-simula",
            "generalized_alpha": alpha,
        }
    )
    if formulation in {"local-pressure", "local-pressure-paper"}:
        configuration["local_pressure_stage_alpha_f"] = 0.4
    times = np.asarray(payload["times"], dtype=float)
    load_times = times.copy()
    load_times[0] = 0.0
    load_times[1:] -= 0.4 * float(payload["dt"])
    source_tau = np.zeros_like(times)
    source_pressure = np.zeros_like(times)
    selected = None
    if benchmark_step is not None:
        selected = post.benchmark_configuration(benchmark_step, case)
    active_enabled = (
        case == "A" if selected is None else selected.active_stress_enabled
    )
    pressure_enabled = (
        case == "B" if selected is None else selected.pressure_enabled
    )
    if active_enabled:
        source_tau[1:] = post.tau_of_t(
            load_times[1:],
            t_span=(0.0, float(payload["load_horizon"])),
            **(
                {}
                if selected is None
                else {"p": selected.activation_parameters}
            ),
        )
    if pressure_enabled:
        source_pressure[1:] = post.p_of_t(
            load_times[1:],
            t_span=(0.0, float(payload["load_horizon"])),
            **(
                {}
                if selected is None
                else {"p": selected.pressure_parameters}
            ),
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
            "tau": source_tau,
            "pres": source_pressure,
            "solver_configuration_json": json.dumps(configuration),
            "mpi_implementation": implementation,
        }
    )
    if selected is not None:
        payload.update(
            post.benchmark_metadata(
                selected,
                material_parameters=selected.material_parameters,
                activation_parameters=selected.activation_parameters,
                pressure_parameters=selected.pressure_parameters,
            )
        )
    if formulation in {"local-pressure", "local-pressure-paper"}:
        payload["element_pressure_peak_stage"] = "alpha_f_force_stage"
    np.savez(path, **payload)
    return payload


def _replace_json_payload(payload, field, transform):
    value = json.loads(str(payload[field]))
    transform(value)
    payload[field] = json.dumps(value)


def _write_reference(path, *, start=0.0, end=1.0, malformed=False):
    if malformed:
        path.write_bytes(b"not a pickle")
        return
    times = np.linspace(start, end, 101)
    displacement = {}
    for point_factor, point in enumerate(("p0", "p1"), start=1):
        displacement[point] = {
            "ux": point_factor * times,
            "uy": -0.5 * point_factor * times,
            "uz": 0.25 * point_factor * times,
        }
    with path.open("wb") as stream:
        pickle.dump({"time": times, "displacement": displacement}, stream)


def _write_reference_manifest(
    data,
    *,
    case="step_0B",
    alias=True,
    overrides=None,
):
    overrides = overrides or {}
    paths = {}
    for suffix in post.REFERENCE_MANIFEST_SUFFIXES:
        path = data / post._reference_filename(case, suffix)
        _write_reference(path, **overrides.get(suffix, {}))
        paths[suffix] = path
    if alias:
        alias_path = data / post._reference_filename(
            case, post.REFERENCE_EXCLUDED_ALIAS_SUFFIX
        )
        alias_path.write_bytes(
            paths[post.REFERENCE_ALIAS_TARGET_SUFFIX].read_bytes()
        )
        paths[post.REFERENCE_EXCLUDED_ALIAS_SUFFIX] = alias_path
    return paths


def test_distorted_affine_hex8_centroid_gradient_is_exact():
    mapping = np.array(
        [[0.7, 0.2, -0.1], [0.1, 0.5, 0.15], [-0.05, 0.12, 0.8]]
    )
    nodes = diagnose._NAT @ mapping.T + np.array([0.3, 1.4, 1.7])
    displacement_gradient = np.array(
        [[0.08, -0.03, 0.02], [0.04, 0.05, -0.01], [-0.02, 0.06, 0.07]]
    )
    displacement = nodes @ displacement_gradient.T + np.array([0.1, -0.2, 0.3])
    np.testing.assert_allclose(
        diagnose.F_centroid(nodes, displacement),
        np.eye(3) + displacement_gradient,
        rtol=1.0e-13,
        atol=1.0e-13,
    )


@pytest.mark.parametrize(
    "formulation",
    [
        "hex8_fbar",
        "hex8_local_pressure_p0_condensed_logj",
        "hex8_local_pressure_p0_condensed_mean_logj_paper_j2",
        "hex8_standard_pointwise_kappa",
    ],
)
def test_current_formulations_and_hex8_sampling_metadata_are_accepted(
    tmp_path, formulation
):
    path = tmp_path / "current.npz"
    _write_current_result(path, formulation=formulation)
    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, path)
    assert result["formulation"] == formulation
    assert result["point_sampling"] == "hex8_reference_isoparametric"
    assert result["det_f_gauss_peak"].shape == (1, 8)
    if "local_pressure" in formulation:
        assert result["element_pressure_peak_pa"].shape == (1,)


def test_closed_multiblock_standard_gp_direct_result_is_accepted(tmp_path):
    path = tmp_path / "closed-current.npz"
    _write_current_result(
        path,
        formulation="hex8_standard_pointwise_kappa",
        fiber_sampling="gp_direct_rule",
        closed=True,
    )

    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, path)

    assert result["formulation"] == "hex8_standard_pointwise_kappa"
    assert result["fiber_sampling"] == "gp_direct_rule"
    assert result["fiber_sampling_option"] == "gp-direct"
    assert result["mass_representation"] == "consistent_q1_hex8"
    assert result["method_metadata_origin"] == "recorded"
    assert result["material_eta_pa_s"] == 100.0
    assert result["viscous_term_active"] is True
    assert result["parameter_variant"] == "benchmark_eta"
    assert result["isotropic"] is False
    assert result["model_metadata"]["material_model_id"] == (
        post.COMPLETE_SWITCH_ENERGY_MATERIAL_MODEL_ID
    )
    assert result["tbar_identity"] == {
        "definition": "analytic_parametric",
        "source_filename": "",
        "source_sha256": "",
        "metadata_filename": "",
        "metadata_sha256": "",
        "metadata_schema": "",
    }
    assert result["mesh"]["topology"] == "closed_multiblock_disk"
    assert result["mesh"]["elements"] == 1
    assert result["mesh"]["n_mu"] * result["mesh"]["n_theta"] == 0
    assert result["model_metadata"]["material_kernel_formulation"] == "standard"
    assert result["model_metadata"]["material_kappa_pa"] == 1.0e6
    assert result["model_metadata"]["local_pressure_bulk_modulus_pa"] == 0.0
    assert result["pre_solve_audit"]["geometry"][
        "unclassified_exterior_faces"
    ] == 0
    rotation = post._peak_circumferential_ring_rotation(result)
    assert rotation["available"] is False
    assert "closed multiblock" in rotation["reason"]
    zero_curve = np.zeros((len(post.CANONICAL_TIME_GRID), 3))
    report = post._build_report(
        path,
        result,
        {"synthetic": {"p0": zero_curve, "p1": zero_curve}},
        [],
        {},
        {"p0": zero_curve, "p1": zero_curve},
    )
    assert report["result"]["pre_solve_audit"]["geometry"]["passed"] is True
    configuration = report["result"]["configuration"]
    assert configuration["mass_representation"] == "consistent_q1_hex8"
    assert configuration["method_metadata_origin"] == "recorded"
    assert configuration["fiber_sampling_option"] == "gp-direct"
    assert configuration["material_eta_pa_s"] == 100.0
    assert configuration["viscous_term_active"] is True
    assert configuration["parameter_variant"] == "benchmark_eta"
    assert configuration["tbar"] == result["tbar_identity"]
    assert "retained pre-solve" in report["bounded_claim"]


def test_load_horizon_is_validated_retained_and_reported(tmp_path):
    path = tmp_path / "fixed-load-horizon.npz"
    _write_current_result(path, load_horizon=1.5)
    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, path)

    assert result["load_horizon"] == 1.5
    assert result["load_horizon_origin"] == "recorded"
    zero_curve = np.zeros((len(post.CANONICAL_TIME_GRID), 3))
    report = post._build_report(
        path,
        result,
        {"synthetic": {"p0": zero_curve, "p1": zero_curve}},
        [],
        {},
        {"p0": zero_curve, "p1": zero_curve},
    )
    configuration = report["result"]["configuration"]
    assert configuration["load_horizon_s"] == 1.5
    assert configuration["load_horizon_origin"] == "recorded"


def test_missing_load_horizon_has_explicit_legacy_backfill(tmp_path):
    path = tmp_path / "implicit-load-horizon.npz"
    payload = _write_current_result(path)
    payload.pop("load_horizon")
    np.savez(path, **payload)
    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, path)

    assert result["load_horizon"] == result["t_end"] == 1.0
    assert result["load_horizon_origin"] == "implicit_t_end_legacy"


@pytest.mark.parametrize("load_horizon", [0.5, 1.25, float("nan")])
def test_invalid_load_horizon_is_rejected(tmp_path, load_horizon):
    path = tmp_path / "invalid-load-horizon.npz"
    _write_current_result(path, load_horizon=load_horizon)
    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match="invalid load_horizon"):
            post._load_validated_result(archive, path)


def test_portable_laplace_tbar_identity_is_retained_without_a_path(tmp_path):
    path = tmp_path / "portable-laplace.npz"
    source_sha256 = hashlib.sha256(b"retained tbar field").hexdigest()
    metadata_sha256 = hashlib.sha256(b"retained tbar metadata").hexdigest()
    _write_current_result(
        path,
        tbar_definition="laplace_presolved",
        tbar_source_filename="closed_tbar.npy",
        tbar_source_sha256=source_sha256,
        tbar_metadata_filename="closed_tbar.meta.json",
        tbar_metadata_sha256=metadata_sha256,
        tbar_metadata_schema=post.TBAR_METADATA_SCHEMA,
    )

    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, path)

    assert result["tbar_identity"] == {
        "definition": "laplace_presolved",
        "source_filename": "closed_tbar.npy",
        "source_sha256": source_sha256,
        "metadata_filename": "closed_tbar.meta.json",
        "metadata_sha256": metadata_sha256,
        "metadata_schema": post.TBAR_METADATA_SCHEMA,
    }


def test_legacy_absolute_tbar_input_is_hashed_and_sanitized(tmp_path):
    tbar_path = tmp_path / "private" / "legacy_tbar.npy"
    tbar_path.parent.mkdir()
    tbar_path.write_bytes(b"legacy development field")
    result_path = tmp_path / "legacy-absolute.npz"
    payload = _write_current_result(result_path)
    payload["tbar_definition"] = f"laplace_presolved:{tbar_path.resolve()}"
    payload.pop("tbar_source_filename")
    payload.pop("tbar_source_sha256")
    payload.pop("tbar_metadata_filename")
    payload.pop("tbar_metadata_sha256")
    payload.pop("tbar_metadata_schema")
    np.savez(result_path, **payload)

    with np.load(result_path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, result_path)

    assert result["tbar_identity"] == {
        "definition": "laplace_presolved",
        "source_filename": tbar_path.name,
        "source_sha256": hashlib.sha256(tbar_path.read_bytes()).hexdigest(),
        "metadata_filename": "",
        "metadata_sha256": "",
        "metadata_schema": "",
    }
    zero_curve = np.zeros((len(post.CANONICAL_TIME_GRID), 3))
    report = post._build_report(
        result_path,
        result,
        {"synthetic": {"p0": zero_curve, "p1": zero_curve}},
        [],
        {},
        {"p0": zero_curve, "p1": zero_curve},
    )
    assert report["result"]["configuration"]["tbar"] == result["tbar_identity"]
    assert str(tmp_path) not in json.dumps(report)


def test_closed_case_a_does_not_require_a_pressure_audit(tmp_path):
    path = tmp_path / "closed-case-a.npz"
    _write_current_result(path, case="A", closed=True)

    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, path)

    assert result["case"] == "A"
    assert set(result["pre_solve_audit"]) == {"geometry", "robin"}


def test_closed_case_b_still_requires_a_pressure_audit(tmp_path):
    path = tmp_path / "closed-case-b-no-pressure.npz"
    payload = _write_current_result(path, case="B", closed=True)
    audit = json.loads(payload["pre_solve_audit_json"])
    audit.pop("pressure")
    payload["pre_solve_audit_json"] = json.dumps(audit)
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match="incomplete closed-mesh pre-solve audit"):
            post._load_validated_result(archive, path)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("mass", "unsupported mass representation"),
        ("fiber-option", "inconsistent fiber sampling option/method"),
        ("material-model", "unsupported material model identity"),
        ("eta", "invalid material_eta_pa_s"),
        ("viscous-active", "inconsistent viscous activity metadata"),
        ("variant", "inconsistent viscosity parameter variant"),
        ("isotropic", "non-boolean result field 'isotropic'"),
        ("tbar-filename", "non-portable Laplace tbar filename"),
        ("tbar-sha", "invalid Laplace tbar SHA-256"),
        ("missing-tbar-sha", "incomplete tbar provenance"),
        ("tbar-metadata-filename", "non-portable Laplace tbar metadata filename"),
        ("tbar-metadata-sha", "invalid Laplace tbar metadata SHA-256"),
        ("tbar-metadata-schema", "unsupported Laplace tbar metadata schema"),
        ("missing-tbar-metadata-sha", "incomplete tbar metadata provenance"),
    ],
)
def test_current_result_rejects_inconsistent_result_defining_metadata(
    tmp_path, mutation, message
):
    path = tmp_path / f"bad-{mutation}.npz"
    payload = _write_current_result(path)
    if mutation == "mass":
        payload["mass_representation"] = "diagonal"
    elif mutation == "fiber-option":
        payload["fiber_sampling_option"] = "gp-direct"
    elif mutation == "material-model":
        payload["material_model_id"] = "unspecified-cardiac-law"
    elif mutation == "eta":
        payload["material_eta_pa_s"] = -1.0
    elif mutation == "viscous-active":
        payload["viscous_term_active"] = False
    elif mutation == "variant":
        payload["parameter_variant"] = "eta_sensitivity"
    elif mutation == "isotropic":
        payload["isotropic"] = "false"
    elif mutation.startswith("tbar-") or mutation.startswith("missing-tbar-"):
        payload["tbar_definition"] = "laplace_presolved"
        payload["tbar_source_filename"] = "closed_tbar.npy"
        payload["tbar_source_sha256"] = "1" * 64
        payload["tbar_metadata_filename"] = "closed_tbar.meta.json"
        payload["tbar_metadata_sha256"] = "2" * 64
        payload["tbar_metadata_schema"] = post.TBAR_METADATA_SCHEMA
        if mutation == "tbar-filename":
            payload["tbar_source_filename"] = "../closed_tbar.npy"
        elif mutation == "tbar-sha":
            payload["tbar_source_sha256"] = "not-a-sha"
        elif mutation == "missing-tbar-sha":
            payload.pop("tbar_source_sha256")
        elif mutation == "tbar-metadata-filename":
            payload["tbar_metadata_filename"] = "../closed_tbar.meta.json"
        elif mutation == "tbar-metadata-sha":
            payload["tbar_metadata_sha256"] = "not-a-sha"
        elif mutation == "tbar-metadata-schema":
            payload["tbar_metadata_schema"] = "unknown"
        else:
            payload.pop("tbar_metadata_sha256")
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match=message):
            post._load_validated_result(archive, path)


def test_historical_polar_result_without_topology_field_remains_supported(tmp_path):
    path = tmp_path / "historical-polar.npz"
    _write_current_result(path)
    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, path)
    assert result["mesh"]["topology"] == "polar_ring"
    assert result["pre_solve_audit"] is None


_PRE_METHOD_METADATA_FIELDS = {
    "mass_representation",
    "fiber_sampling_option",
    "isotropic",
    "material_model_id",
    "material_eta_pa_s",
    "viscous_term_active",
    "parameter_variant",
    "tbar_definition",
    "tbar_source_filename",
    "tbar_source_sha256",
    "tbar_metadata_filename",
    "tbar_metadata_sha256",
    "tbar_metadata_schema",
}


@pytest.mark.parametrize(
    "app_revision",
    [
        "62ad760d2a1731bb9668897863ac026d3768194e",
        "6839c13b5bc80ec06c897684c51f503e80bd4b19",
    ],
)
def test_immutable_pre_metadata_result_gets_narrow_reviewed_defaults(
    tmp_path, app_revision
):
    path = tmp_path / "immutable-pre-metadata.npz"
    payload = _write_current_result(path)
    for field in _PRE_METHOD_METADATA_FIELDS:
        payload.pop(field)
    payload["app_revision"] = app_revision
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, path)

    assert result["method_metadata_origin"] == (
        "reviewed-predecessor-source-checkpoint"
    )
    assert result["mass_representation"] == "lumped_row_sum"
    assert result["fiber_sampling_option"] == "cg1"
    assert result["material_eta_pa_s"] == 100.0
    assert result["parameter_variant"] == "benchmark_eta"
    expected_model_id = (
        post.COMPLETE_SWITCH_ENERGY_MATERIAL_MODEL_ID
        if app_revision == "6839c13b5bc80ec06c897684c51f503e80bd4b19"
        else post.LEGACY_SWITCH_STRESS_MATERIAL_MODEL_ID
    )
    assert result["model_metadata"]["material_model_id"] == expected_model_id
    assert result["isotropic"] is False
    assert result["tbar_identity"]["definition"] == "analytic_parametric"


def test_unreviewed_result_cannot_omit_current_method_metadata(tmp_path):
    path = tmp_path / "unreviewed-pre-metadata.npz"
    payload = _write_current_result(path)
    for field in _PRE_METHOD_METADATA_FIELDS:
        payload.pop(field)
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match="immutable reviewed predecessor"):
            post._load_validated_result(archive, path)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("unclassified", "uncleared closed-mesh geometry audit"),
        ("pressure-failed", "failed pre-solve pressure audit"),
        ("pressure-missing-signed", "invalid pressure audit signed axial ratio"),
        (
            "pressure-tampered-signed",
            "inconsistent pre-solve pressure audit signed_axial_ratio",
        ),
        ("pressure-reversed", "out-of-tolerance pre-solve pressure audit"),
        ("robin-failed", "failed pre-solve Robin audit"),
        ("element-count", "geometry audit/mesh size disagreement"),
    ],
)
def test_closed_result_rejects_broken_pre_solve_evidence(
    tmp_path, mutation, message
):
    path = tmp_path / f"broken-{mutation}.npz"
    payload = _write_current_result(
        path,
        formulation="hex8_standard_pointwise_kappa",
        fiber_sampling="gp_direct_rule",
        closed=True,
    )
    audit = json.loads(payload["pre_solve_audit_json"])
    if mutation == "unclassified":
        audit["geometry"]["unclassified_exterior_faces"] = 1
    elif mutation == "pressure-failed":
        audit["pressure"]["passed"] = False
        audit["pressure"]["failures"] = ["pressure mismatch"]
    elif mutation == "pressure-missing-signed":
        audit["pressure"].pop("signed_axial_ratio")
    elif mutation == "pressure-tampered-signed":
        audit["pressure"]["signed_axial_ratio"] = -1.0
    elif mutation == "pressure-reversed":
        pressure = audit["pressure"]
        projected_area = pressure["analytic_projected_base_area_m2"]
        pressure.update(
            {
                "unit_pressure_resultant_N": [-projected_area, 0.0, 0.0],
                "relative_magnitude_error": 0.0,
                "signed_axial_ratio": -1.0,
                "relative_signed_axial_error": 2.0,
                "relative_resultant_error": 2.0,
                "transverse_fraction": 0.0,
            }
        )
    elif mutation == "robin-failed":
        audit["robin"]["passed"] = False
        audit["robin"]["failures"] = ["matrix mismatch"]
    elif mutation == "element-count":
        payload["elems"] = np.tile(payload["elems"], (2, 1))
    payload["pre_solve_audit_json"] = json.dumps(audit)
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match=message):
            post._load_validated_result(archive, path)


def test_closed_result_does_not_bypass_clean_source_gate(tmp_path):
    path = tmp_path / "dirty-closed.npz"
    payload = _write_current_result(
        path,
        formulation="hex8_standard_pointwise_kappa",
        fiber_sampling="gp_direct_rule",
        closed=True,
    )
    payload["app_tree_state"] = "dirty"
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match="requires a clean source tree"):
            post._load_validated_result(archive, path)


def test_malformed_matching_reference_file_fails_closed(tmp_path):
    data = tmp_path / "results_time_curves" / "data"
    data.mkdir(parents=True)
    _write_reference_manifest(
        data,
        overrides={"chimera": {"malformed": True}},
    )
    with pytest.raises(RuntimeError, match="malformed|could not read|reference"):
        post.load_reference("step_0B", tmp_path)


def test_reference_selection_is_exact_and_deduplicates_known_alias(tmp_path):
    data = tmp_path / "results_time_curves" / "data"
    data.mkdir(parents=True)
    paths = _write_reference_manifest(data)

    reference, provenance, selection = post.load_reference(
        "step_0B",
        tmp_path,
        return_provenance=True,
        return_selection=True,
    )

    assert list(reference) == list(post.REFERENCE_MANIFEST_SUFFIXES)
    assert [entry["team"] for entry in provenance] == list(
        post.REFERENCE_MANIFEST_SUFFIXES
    )
    assert selection["selected_count"] == 10
    assert selection["selected_files"] == [
        post._reference_filename("step_0B", suffix)
        for suffix in post.REFERENCE_MANIFEST_SUFFIXES
    ]
    assert selection["upstream_figures_py_identity"] == (
        post.REFERENCE_FIGURES_PY_IDENTITY
    )
    assert selection["upstream_manifest_variable"] == "TEAMS_DATASET_0B"
    alias = selection["excluded_aliases"]
    assert len(alias) == 1
    assert alias[0]["sha256"] == hashlib.sha256(
        paths[post.REFERENCE_EXCLUDED_ALIAS_SUFFIX].read_bytes()
    ).hexdigest()
    assert alias[0]["identical_to_selected_filename"].endswith(
        "_simvascular_p2.pickle"
    )
    assert "excluded" in alias[0]["reason"]


def test_reference_selection_rejects_mutated_known_alias(tmp_path):
    data = tmp_path / "results_time_curves" / "data"
    data.mkdir(parents=True)
    paths = _write_reference_manifest(data)
    paths[post.REFERENCE_EXCLUDED_ALIAS_SUFFIX].write_bytes(b"mutated alias")

    with pytest.raises(RuntimeError, match="not byte-identical"):
        post.load_reference("step_0B", tmp_path)


def test_reference_selection_rejects_selected_file_schema_change(tmp_path):
    data = tmp_path / "results_time_curves" / "data"
    data.mkdir(parents=True)
    paths = _write_reference_manifest(data, alias=False)
    with paths["cheart"].open("wb") as stream:
        pickle.dump({"time": np.linspace(0.0, 1.0, 101)}, stream)

    with pytest.raises(RuntimeError, match="missing required data"):
        post.load_reference("step_0B", tmp_path)


def test_reference_selection_rejects_missing_and_unexpected_files(tmp_path):
    data = tmp_path / "results_time_curves" / "data"
    data.mkdir(parents=True)
    paths = _write_reference_manifest(data)
    paths["comsol"].unlink()
    with pytest.raises(RuntimeError, match="missing selected files"):
        post.load_reference("step_0B", tmp_path)

    _write_reference(paths["comsol"])
    unexpected = data / "monoventricular_nonblinded_step_0B_group_extra.pickle"
    _write_reference(unexpected)
    with pytest.raises(RuntimeError, match="unexpected matching files"):
        post.load_reference("step_0B", tmp_path)


def test_reference_endpoint_offsets_are_resampled_to_canonical_grid(tmp_path):
    data = tmp_path / "results_time_curves" / "data"
    data.mkdir(parents=True)
    _write_reference_manifest(
        data,
        overrides={
            "ambit": {"start": 0.001, "end": 1.0},
            "simula": {"start": 0.0, "end": 0.999},
        },
    )
    reference, provenance = post.load_reference(
        "step_0B", tmp_path, return_provenance=True
    )
    for record in reference.values():
        np.testing.assert_array_equal(record["t"], np.linspace(0.0, 1.0, 101))
        assert record["p0"].shape == (101, 3)
        assert np.all(np.isfinite(record["p0"]))
    assert {entry["source_time_start"] for entry in provenance} == {0.0, 0.001}
    assert {entry["source_time_end"] for entry in provenance} == {0.999, 1.0}


def test_reference_resampling_rejects_missing_source_coverage():
    times = np.linspace(0.01, 1.0, 100)
    with pytest.raises(RuntimeError, match="does not cover"):
        post.resample(times, np.zeros((100, 3)), description="short team")


def test_peak_ring_rotation_uses_centered_least_squares_and_sign_convention():
    n_t, n_mu, n_theta = 2, 2, 8
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    base_angles = (10.0, 20.0, 30.0)
    apex_angles = (-35.0, -25.0, -15.0)
    nodes = []
    deformed = []
    for layer in range(n_t + 1):
        for longitudinal_ring in range(n_mu + 1):
            radius = 1.0 + 0.1 * layer + 0.2 * longitudinal_ring
            reference_yz = np.stack(
                (radius * np.cos(theta), radius * np.sin(theta)), axis=1
            )
            angle = np.deg2rad(
                base_angles[layer]
                + (apex_angles[layer] - base_angles[layer])
                * longitudinal_ring
                / n_mu
            )
            rotation = np.array(
                [[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]]
            )
            deformed_yz = 1.2 * reference_yz @ rotation.T
            reference_yz += np.array([0.4 * layer, -0.3 * longitudinal_ring])
            deformed_yz += np.array([-2.0 + layer, 3.0 - longitudinal_ring])
            x = np.full((n_theta, 1), float(longitudinal_ring))
            nodes.append(np.concatenate((x, reference_yz), axis=1))
            deformed.append(np.concatenate((x, deformed_yz), axis=1))
    nodes = np.concatenate(nodes)
    deformed = np.concatenate(deformed)
    result = {
        "case": "B",
        "mesh": {"n_t": n_t, "n_mu": n_mu, "n_theta": n_theta},
        "nodes": nodes,
        "peak_displacement": (deformed - nodes).reshape(-1),
        "times": np.array([0.0, 0.5]),
        "n_peak": 1,
    }

    profile = post._peak_circumferential_ring_rotation(result)

    assert profile["available"] is True
    assert "+y toward +z" in profile["sign_convention"]
    np.testing.assert_allclose(
        [layer["base_rotation_degrees"] for layer in profile["layers"]],
        base_angles,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        [layer["apex_rotation_degrees"] for layer in profile["layers"]],
        apex_angles,
        atol=1.0e-12,
    )
    np.testing.assert_allclose(
        [
            layer["relative_apex_minus_base_degrees"]
            for layer in profile["layers"]
        ],
        [-45.0, -45.0, -45.0],
        atol=1.0e-12,
    )


def test_peak_ring_rotation_reports_unavailable_without_complete_rings():
    result = {
        "case": "B",
        "mesh": {"n_t": 2, "n_mu": 2, "n_theta": 8},
        "nodes": np.zeros((3, 3)),
        "peak_displacement": np.zeros(9),
    }
    profile = post._peak_circumferential_ring_rotation(result)
    assert profile["available"] is False
    assert "complete ordered structured rings" in profile["reason"]


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("snes_converged_reason", -2, "SNES reason"),
        ("ksp_converged_reason", -3, "KSP reason"),
        ("final_residual_norm", 2.0e-8, "above its acceptance threshold"),
        ("initial_residual_norm", "NaN", "PETSc initial residual"),
    ],
)
def test_bad_petsc_solver_diagnostics_are_rejected(tmp_path, field, value, message):
    path = tmp_path / "bad.npz"
    times = np.array([0.0, 0.5, 1.0])
    diagnostics = _petsc_diagnostics(times)
    diagnostics[0][field] = value
    _write_current_result(
        path,
        solver="petsc-snes",
        diagnostics=diagnostics,
    )
    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match=message):
            post._load_validated_result(archive, path)


def test_historical_mpi_companion_archive_is_validated_and_retained(tmp_path):
    path = tmp_path / "mpi-valid.npz"
    payload = _write_mpi_result(path, mpi_ranks=2)
    _downgrade_to_legacy_direct_profile(payload)
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, path)

    assert result["solver_name"] == "petsc-snes-mpi"
    assert result["mpi_metadata"] == {
        "implementation": post.MPI_COMPANION_IMPLEMENTATION,
        "world_size": 2,
        "local_element_counts": [1, 0],
        "partition": "coupfe.partition_elements",
        "build_layout": "isolated-rank-directories",
        "factor_solver_type": "superlu_dist",
        "linear_solver_profile": post.MPI_DIRECT_SUPERLU_DIST_PROFILE,
        "element_evaluation_mode": "joint",
        "compiled_material_residual_only_available": True,
    }
    assert all(record["ranks"] == 2 for record in result["solver_diagnostics"])


@pytest.mark.parametrize(
    "profile",
    [
        post.MPI_FGMRES_GAMG_RIGID_PROFILE,
        post.MPI_FGMRES_GAMG_RIGID_REBUILD_PROFILE,
        post.MPI_FGMRES_ASM_LU_PROFILE,
        post.MPI_FGMRES_ASM_ILU1_PROFILE,
    ],
)
def test_reviewed_iterative_mpi_profiles_are_validated_and_retained(
    tmp_path, profile
):
    path = tmp_path / f"mpi-{profile}.npz"
    _write_mpi_result(path, mpi_ranks=2, linear_solver_profile=profile)

    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, path)

    assert result["mpi_metadata"]["linear_solver_profile"] == profile
    assert result["mpi_metadata"]["factor_solver_type"] == "none"
    assert result["solver_configuration"]["ksp_type"] == "fgmres"


def test_iterative_mpi_profile_rejects_changed_petsc_option(tmp_path):
    path = tmp_path / "mpi-iterative-option.npz"
    payload = _write_mpi_result(
        path,
        mpi_ranks=2,
        linear_solver_profile=post.MPI_FGMRES_GAMG_RIGID_REBUILD_PROFILE,
    )
    configuration = json.loads(str(payload["solver_configuration_json"]))
    option_name = next(
        name
        for name in configuration["petsc_options"]
        if name.endswith("pc_gamg_threshold")
    )
    configuration["petsc_options"][option_name] = "0.02"
    payload["solver_configuration_json"] = json.dumps(configuration)
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match="unexpected PETSc option"):
            post._load_validated_result(archive, path)


def test_iterative_mpi_profile_cannot_be_elided_as_legacy_direct(tmp_path):
    path = tmp_path / "mpi-iterative-profile-elided.npz"
    payload = _write_mpi_result(
        path,
        mpi_ranks=2,
        linear_solver_profile=post.MPI_FGMRES_ASM_LU_PROFILE,
    )
    payload.pop("mpi_linear_solver_profile")
    configuration = json.loads(str(payload["solver_configuration_json"]))
    configuration.pop("linear_solver_profile")
    payload["solver_configuration_json"] = json.dumps(configuration)
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(
            RuntimeError, match="profile configuration for 'factor_solver_type'"
        ):
            post._load_validated_result(archive, path)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("unsupported", "unsupported MPI linear-solver profile"),
        ("archive-config-disagreement", "profile archive/config disagreement"),
        ("missing-archive-profile", "profile archive/config disagreement"),
        ("missing-configuration-profile", "profile archive/config disagreement"),
        (
            "both-profiles-elided",
            "explicit-profile MPI configuration fields without duplicated",
        ),
        ("configuration", "profile configuration for 'pc_type'"),
        ("missing-null-setting", "profile configuration for 'ksp_rtol'"),
    ],
)
def test_mpi_linear_solver_profile_disagreement_is_rejected(
    tmp_path, mutation, message
):
    path = tmp_path / f"mpi-profile-{mutation}.npz"
    payload = _write_mpi_result(path, mpi_ranks=2)
    if mutation == "unsupported":
        payload["mpi_linear_solver_profile"] = "unreviewed-profile"
    elif mutation == "archive-config-disagreement":
        payload["mpi_linear_solver_profile"] = post.MPI_FGMRES_ASM_LU_PROFILE
    elif mutation == "missing-archive-profile":
        payload.pop("mpi_linear_solver_profile")
    elif mutation == "missing-configuration-profile":
        _replace_json_payload(
            payload,
            "solver_configuration_json",
            lambda configuration: configuration.pop("linear_solver_profile"),
        )
    elif mutation == "both-profiles-elided":
        payload.pop("mpi_linear_solver_profile")
        _replace_json_payload(
            payload,
            "solver_configuration_json",
            lambda configuration: configuration.pop("linear_solver_profile"),
        )
    elif mutation == "configuration":
        _replace_json_payload(
            payload,
            "solver_configuration_json",
            lambda configuration: configuration.__setitem__("pc_type", "asm"),
        )
    elif mutation == "missing-null-setting":
        _replace_json_payload(
            payload,
            "solver_configuration_json",
            lambda configuration: configuration.pop("ksp_rtol"),
        )
    else:  # pragma: no cover - parameterized helper programming error
        raise AssertionError(mutation)
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match=message):
            post._load_validated_result(archive, path)


def test_closed_case_b_mpi_archive_and_report_retain_mass_provenance(tmp_path):
    path = tmp_path / "mpi-closed-valid.npz"
    _write_closed_mpi_result(path, mpi_ranks=2)

    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, path)

    assert result["solver_name"] == "petsc-snes-mpi"
    assert result["mpi_metadata"] == {
        "implementation": post.MPI_CLOSED_STD_KAPPA_IMPLEMENTATION,
        "world_size": 2,
        "local_element_counts": [1, 0],
        "partition": "coupfe.partition_elements",
        "build_layout": "isolated-rank-directories",
        "factor_solver_type": "superlu_dist",
        "linear_solver_profile": post.MPI_DIRECT_SUPERLU_DIST_PROFILE,
        "element_evaluation_mode": "joint",
        "compiled_material_residual_only_available": True,
        "contract": "closed_case_b_std_kappa_consistent",
        "mass_partition": {
            "partition": "owned-row-csr-all-touching-elements",
            "owned_row_ranges": [[0, 12], [12, 24]],
            "local_nnz": [64, 64],
            "touching_element_counts": [1, 1],
        },
    }

    zero_curve = np.zeros((len(post.CANONICAL_TIME_GRID), 3))
    report = post._build_report(
        path,
        result,
        {"synthetic": {"p0": zero_curve, "p1": zero_curve}},
        [],
        {},
        {"p0": zero_curve, "p1": zero_curve},
    )
    assert report["result"]["mpi_metadata"] == result["mpi_metadata"]


def test_closed_case_a_mpi_archive_uses_active_step0_contract(tmp_path):
    path = tmp_path / "mpi-closed-case-a-valid.npz"
    _write_closed_mpi_result(path, mpi_ranks=2, case="A")

    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(
            archive, path, requested_case="step_0A"
        )

    assert result["case"] == "A"
    assert set(result["pre_solve_audit"]) == {"geometry", "robin"}
    assert result["mpi_metadata"]["contract"] == (
        "closed_case_a_std_kappa_consistent"
    )


def test_closed_case_a_local_pressure_archive_uses_distinct_contract(tmp_path):
    path = tmp_path / "mpi-closed-case-a-local-pressure-valid.npz"
    _write_closed_mpi_result(
        path,
        mpi_ranks=2,
        case="A",
        formulation="local-pressure",
    )

    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(
            archive, path, requested_case="step_0A"
        )

    assert result["case"] == "A"
    assert result["formulation"] == "hex8_local_pressure_p0_condensed_logj"
    assert result["mpi_metadata"]["implementation"] == (
        post.MPI_CLOSED_LOCAL_PRESSURE_IMPLEMENTATION
    )
    assert result["mpi_metadata"]["contract"] == (
        "closed_case_a_local_pressure_consistent"
    )


def test_closed_case_a_paper_local_pressure_archive_uses_distinct_contract(
    tmp_path,
):
    path = tmp_path / "mpi-closed-case-a-local-pressure-paper-valid.npz"
    _write_closed_mpi_result(
        path,
        mpi_ranks=2,
        case="A",
        formulation="local-pressure-paper",
    )

    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(
            archive, path, requested_case="step_0A"
        )

    assert result["formulation"] == (
        "hex8_local_pressure_p0_condensed_mean_logj_paper_j2"
    )
    assert result["model_metadata"]["local_pressure_volume_law"] == (
        "paper-j2-of-reference-volume-weighted-geometric-mean-j-v1"
    )
    assert result["mpi_metadata"]["contract"] == (
        "closed_case_a_local_pressure_mean_logj_paper_j2_consistent"
    )


@pytest.mark.parametrize("volume_law", [None, "linear-reference-volume-mean-log-j-v1"])
def test_paper_local_pressure_requires_exact_volume_law_metadata(
    tmp_path, volume_law
):
    path = tmp_path / "mpi-closed-case-a-local-pressure-paper-invalid.npz"
    payload = _write_closed_mpi_result(
        path, case="A", formulation="local-pressure-paper"
    )
    if volume_law is None:
        payload.pop("local_pressure_volume_law")
    else:
        payload["local_pressure_volume_law"] = volume_law
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match="volume-law metadata"):
            post._load_validated_result(archive, path)


def test_paper_local_pressure_cannot_borrow_log_implementation_identity(tmp_path):
    path = tmp_path / "mpi-closed-case-a-local-pressure-paper-wrong-id.npz"
    payload = _write_closed_mpi_result(
        path, case="A", formulation="local-pressure-paper"
    )
    payload["mpi_implementation"] = post.MPI_CLOSED_LOCAL_PRESSURE_IMPLEMENTATION
    configuration = json.loads(str(payload["solver_configuration_json"]))
    configuration["implementation"] = post.MPI_CLOSED_LOCAL_PRESSURE_IMPLEMENTATION
    payload["solver_configuration_json"] = json.dumps(configuration)
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match="discretization contract"):
            post._load_validated_result(archive, path)


@pytest.mark.parametrize(
    ("formulation", "wrong_law"),
    [("local-pressure", "paper"), ("local-pressure-paper", "log")],
)
def test_closed_local_pressure_implementation_rejects_wrong_configured_law(
    tmp_path, formulation, wrong_law
):
    path = tmp_path / f"mpi-closed-case-a-{formulation}-wrong-law.npz"
    payload = _write_closed_mpi_result(
        path, case="A", formulation=formulation
    )
    configuration = json.loads(str(payload["solver_configuration_json"]))
    configuration["local_pressure_law"] = wrong_law
    payload["solver_configuration_json"] = json.dumps(configuration)
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match="discretization contract"):
            post._load_validated_result(archive, path)


@pytest.mark.parametrize(
    ("formulation", "contract"),
    [
        ("std-kappa", "closed_case_a_std_kappa_consistent_generalized_alpha"),
        (
            "local-pressure",
            "closed_case_a_local_pressure_consistent_generalized_alpha",
        ),
        (
            "local-pressure-paper",
            "closed_case_a_local_pressure_mean_logj_paper_j2_consistent_"
            "generalized_alpha",
        ),
    ],
)
def test_closed_case_a_generalized_alpha_archive_is_fail_closed_and_reportable(
    tmp_path, formulation, contract
):
    path = tmp_path / f"mpi-closed-case-a-{formulation}-ga-valid.npz"
    _write_generalized_alpha_closed_mpi_result(
        path, mpi_ranks=2, formulation=formulation
    )

    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(
            archive, path, requested_case="step_0A"
        )

    assert result["integrator"] == "generalized-alpha"
    assert result["generalized_alpha"]["alpha_m"] == 0.2
    assert result["generalized_alpha"]["alpha_f"] == 0.4
    assert result["mpi_metadata"]["contract"] == contract


def test_closed_step0_case_b_generalized_alpha_archive_is_reportable(tmp_path):
    path = tmp_path / "mpi-closed-step0-case-b-std-kappa-ga-valid.npz"
    _write_generalized_alpha_closed_mpi_result(
        path,
        mpi_ranks=2,
        formulation="std-kappa",
        case="B",
        benchmark_step=0,
    )

    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(
            archive, path, requested_case="step_0B"
        )

    assert result["benchmark"]["step"] == 0
    assert result["benchmark"]["load_contract"] == "pressure-only"
    assert result["integrator"] == "generalized-alpha"
    assert result["mpi_metadata"]["contract"] == (
        "closed_case_b_std_kappa_consistent_generalized_alpha"
    )
    np.testing.assert_array_equal(
        result["histories"]["tau"],
        np.zeros_like(result["histories"]["tau"]),
    )
    assert np.max(result["histories"]["pres"]) > 0.0


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda payload: payload.__setitem__("benchmark_step", 2),
            "inconsistent with its step/case",
        ),
        (
            lambda payload: payload.__setitem__(
                "pres",
                post.p_of_t(np.asarray(payload["times"], dtype=float)),
            ),
            "inconsistent generalized-alpha pressure history",
        ),
    ],
)
def test_step0_case_b_generalized_alpha_identity_and_pressure_stage_fail_closed(
    tmp_path, mutation, message
):
    path = tmp_path / "mpi-closed-step0-case-b-ga-broken.npz"
    payload = _write_generalized_alpha_closed_mpi_result(
        path,
        formulation="std-kappa",
        case="B",
        benchmark_step=0,
    )
    mutation(payload)
    np.savez(path, **payload)
    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match=message):
            post._load_validated_result(
                archive, path, requested_case="step_0B"
            )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (
            lambda payload: payload.__setitem__(
                "generalized_alpha_alpha_f", 0.3
            ),
            "unexpected generalized-alpha field",
        ),
        (
            lambda payload: np.asarray(
                payload["load_evaluation_times_s"]
            ).__setitem__(1, 0.5),
            "inconsistent generalized-alpha load times",
        ),
        (
            lambda payload: payload.__setitem__(
                "element_pressure_peak_stage", "endpoint"
            ),
            "inconsistent generalized-alpha pressure stage",
        ),
        (
            lambda payload: payload.__setitem__(
                "tau", np.asarray(payload["tau"]) + np.array([0.0, 7.0, 0.0])
            ),
            "inconsistent generalized-alpha activation history",
        ),
        (
            lambda payload: payload.__setitem__(
                "material_model_id", post.LEGACY_SWITCH_STRESS_MATERIAL_MODEL_ID
            ),
            "unreviewed generalized-alpha material model",
        ),
        (
            lambda payload: payload.__setitem__("load_horizon", 2.0),
            "canonical generalized-alpha load horizon",
        ),
    ],
)
def test_generalized_alpha_archive_broken_controls_are_rejected(
    tmp_path, mutation, message
):
    path = tmp_path / "mpi-closed-case-a-ga-broken.npz"
    payload = _write_generalized_alpha_closed_mpi_result(path)
    mutation(payload)
    np.savez(path, **payload)
    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match=message):
            post._load_validated_result(
                archive, path, requested_case="step_0A"
            )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("missing-mass", "incomplete MPI mass provenance"),
        ("analytic-tbar", "closed Case B MPI companion discretization"),
        ("row-gap", "noncontiguous MPI mass row ownership"),
        ("configuration-mismatch", "mass archive/configuration disagreement"),
        ("unsupported-implementation", "unsupported MPI companion implementation"),
    ],
)
def test_closed_case_b_mpi_broken_controls_are_rejected(
    tmp_path, mutation, message
):
    path = tmp_path / f"mpi-closed-{mutation}.npz"
    payload = _write_closed_mpi_result(path, mpi_ranks=2)
    if mutation == "missing-mass":
        payload.pop("mpi_mass_local_nnz")
    elif mutation == "analytic-tbar":
        payload.update(
            {
                "tbar_definition": "analytic_parametric",
                "tbar_source_filename": "",
                "tbar_source_sha256": "",
                "tbar_metadata_filename": "",
                "tbar_metadata_sha256": "",
                "tbar_metadata_schema": "",
            }
        )
    elif mutation == "row-gap":
        payload["mpi_mass_owned_row_ranges"] = np.array(
            [[0, 11], [12, 24]], dtype=np.int64
        )
    elif mutation == "configuration-mismatch":
        _replace_json_payload(
            payload,
            "solver_configuration_json",
            lambda configuration: configuration.__setitem__(
                "mass_local_nnz", 63
            ),
        )
    elif mutation == "unsupported-implementation":
        payload["mpi_implementation"] = "unreviewed-mpi-implementation"
    else:  # pragma: no cover - parameterized helper programming error
        raise AssertionError(mutation)
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match=message):
            post._load_validated_result(archive, path)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("fixed-setting", "configuration.*'rtol'"),
        ("communicator", "configuration.*'communicator'"),
        ("archive-rank", "archive/configuration rank disagreement"),
        ("partition", "unsupported MPI partition policy"),
        ("counts", "partition/count disagreement"),
        ("factor", "factor-solver archive/config disagreement"),
        ("evaluation", "element-evaluation archive/config disagreement"),
        ("diagnostic-rank", "diagnostic/configuration rank disagreement"),
        ("missing-diagnostic-rank", "missing MPI rank metadata"),
        ("historical-semantics", "historical MPI companion discretization"),
        ("missing-archive-field", "incomplete MPI companion provenance"),
    ],
)
def test_mpi_companion_metadata_disagreement_is_rejected(
    tmp_path, mutation, message
):
    path = tmp_path / f"mpi-{mutation}.npz"
    payload = _write_mpi_result(path, mpi_ranks=2)

    if mutation == "fixed-setting":
        _replace_json_payload(
            payload,
            "solver_configuration_json",
            lambda configuration: configuration.__setitem__("rtol", 1.0e-8),
        )
    elif mutation == "communicator":
        _replace_json_payload(
            payload,
            "solver_configuration_json",
            lambda configuration: configuration.__setitem__(
                "communicator", "PETSc.COMM_SELF"
            ),
        )
    elif mutation == "archive-rank":
        payload["mpi_world_size"] = 1
    elif mutation == "partition":
        payload["mpi_partition"] = "unreviewed.partition"
    elif mutation == "counts":
        payload["mpi_local_element_counts"] = np.array([1, 1])
    elif mutation == "factor":
        payload["mpi_factor_solver_type"] = "mumps"
    elif mutation == "evaluation":
        payload["element_evaluation_mode"] = "split"
    elif mutation in {"diagnostic-rank", "missing-diagnostic-rank"}:
        diagnostics = json.loads(str(payload["nonlinear_step_diagnostics_json"]))
        if mutation == "diagnostic-rank":
            diagnostics[0]["ranks"] = 1
        else:
            diagnostics[0].pop("ranks")
        payload["nonlinear_step_diagnostics_json"] = json.dumps(diagnostics)
    elif mutation == "historical-semantics":
        payload["mass_representation"] = "consistent_q1_hex8"
    elif mutation == "missing-archive-field":
        payload.pop("mpi_build_layout")
    else:  # pragma: no cover - parameterized helper programming error
        raise AssertionError(mutation)

    np.savez(path, **payload)
    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match=message):
            post._load_validated_result(archive, path)


def test_non_mpi_solver_rejects_mpi_archive_fields(tmp_path):
    path = tmp_path / "serial-with-mpi-field.npz"
    payload = _write_current_result(path)
    payload["mpi_enabled"] = True
    np.savez(path, **payload)

    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match="MPI companion fields"):
            post._load_validated_result(archive, path)


def test_function_domain_rejection_diagnostics_are_retained(tmp_path):
    path = tmp_path / "domain-rejections.npz"
    _write_current_result(
        path,
        solver="petsc-snes",
        domain_rejections=True,
    )
    with np.load(path, allow_pickle=False) as archive:
        result = post._load_validated_result(archive, path)

    assert result["solver_configuration"]["function_domain_rejection_api"] == (
        post.PETSC_FUNCTION_DOMAIN_REJECTION_API
    )
    assert result["solver_diagnostics"][0]["function_domain_rejections"] == 1
    assert result["solver_diagnostics"][0]["last_function_domain_error"] == (
        "invalid trial det(F)"
    )
    assert result["solver_diagnostics"][1]["function_domain_rejections"] == 0
    assert result["solver_diagnostics"][1]["last_function_domain_error"] is None


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("partial", "inconsistent function-domain rejection diagnostic"),
        ("negative", "function-domain rejection count"),
        ("missing-detail", "inconsistent function-domain rejection detail"),
        ("unexpected-detail", "inconsistent function-domain rejection detail"),
    ],
)
def test_bad_function_domain_rejection_diagnostics_are_rejected(
    tmp_path, mutation, message
):
    path = tmp_path / "bad-domain-rejections.npz"
    times = np.array([0.0, 0.5, 1.0])
    diagnostics = _petsc_diagnostics(times, domain_rejections=True)
    if mutation == "partial":
        diagnostics[0].pop("last_function_domain_error")
    elif mutation == "negative":
        diagnostics[0]["function_domain_rejections"] = -1
    elif mutation == "missing-detail":
        diagnostics[0]["last_function_domain_error"] = None
    elif mutation == "unexpected-detail":
        diagnostics[1]["last_function_domain_error"] = "stale detail"
    _write_current_result(
        path,
        solver="petsc-snes",
        diagnostics=diagnostics,
        domain_rejections=True,
    )
    with np.load(path, allow_pickle=False) as archive:
        with pytest.raises(RuntimeError, match=message):
            post._load_validated_result(archive, path)


def test_json_report_is_full_precision_path_free_and_hash_identified(tmp_path):
    result_path = tmp_path / "nested" / "result.npz"
    result_path.parent.mkdir()
    precise = 0.12345678901234566
    _write_current_result(
        result_path,
        history_value=precise,
        solver="petsc-snes",
        domain_rejections=True,
    )
    data = tmp_path / "results_time_curves" / "data"
    data.mkdir(parents=True)
    reference_paths = _write_reference_manifest(data)
    reference_path = reference_paths[post.REFERENCE_MANIFEST_SUFFIXES[0]]
    run_log = tmp_path / "private" / "run.stdout"
    run_log.parent.mkdir()
    run_log.write_bytes(b"line one\r\nline two")
    json_path = tmp_path / "comparison.json"

    post.main(
        [
            str(result_path),
            "--reference-dir",
            str(tmp_path),
            "--json",
            str(json_path),
            "--run-log",
            str(run_log),
            "--supersedes-report-sha256",
            "1" * 64,
        ]
    )
    raw_report = json_path.read_text(encoding="utf-8")
    report = json.loads(raw_report)
    assert report["schema"] == post.REPORT_SCHEMA
    assert str(tmp_path) not in raw_report
    assert report["result"]["filename"] == "result.npz"
    assert report["result"]["sha256"] == hashlib.sha256(
        result_path.read_bytes()
    ).hexdigest()
    assert report["result"]["retained_histories"]["u0_m"][1][0] == precise
    assert report["result"]["solver_configuration"][
        "function_domain_rejection_api"
    ] == post.PETSC_FUNCTION_DOMAIN_REJECTION_API
    assert report["result"]["nonlinear_step_diagnostics"][0][
        "function_domain_rejections"
    ] == 1
    assert report["reference"]["canonical_grid_s"] == np.linspace(0, 1, 101).tolist()
    assert report["reference"]["published_archive_identity"] == (
        post.REFERENCE_ARCHIVE_IDENTITY
    )
    assert report["reference"]["selection"]["selected_count"] == 10
    assert len(report["reference"]["selection"]["excluded_aliases"]) == 1
    assert report["correction"] == {
        "predecessor_repository_revision": post.CORRECTION_PREDECESSOR_REVISION,
        "reason": post.CORRECTION_REASON,
        "supersedes_report_sha256": "1" * 64,
    }
    team_file = report["reference"]["team_files"][0]
    assert team_file["filename"] == reference_path.name
    assert team_file["sha256"] == hashlib.sha256(reference_path.read_bytes()).hexdigest()
    normalized = b"line one\nline two\n"
    assert report["result"]["normalized_run_log"]["sha256"] == hashlib.sha256(
        normalized
    ).hexdigest()
    assert "example-level evidence" in report["bounded_claim"]
    assert "mesh/time convergence" in report["bounded_claim"]


def test_diagnose_prefers_stored_eight_gauss_point_det_f(tmp_path):
    path = tmp_path / "spatial.npz"
    nodes = diagnose._NAT.astype(float) * 0.5 + np.array([0.0, 2.0, 2.0])
    stored = np.arange(8, dtype=float).reshape(1, 8) / 100.0 + 1.01
    np.savez(
        path,
        result_schema=post.RESULT_SCHEMA,
        converged=True,
        completed_steps=2,
        expected_steps=2,
        times=np.array([0.0, 0.5, 1.0]),
        n_peak=1,
        nodes=nodes,
        elems=np.arange(8, dtype=int).reshape(1, 8),
        fiber=np.array([[1.0, 0.0, 0.0]]),
        facets_endo=np.array([[0, 1, 2, 3]], dtype=int),
        U_peak=np.zeros(24),
        det_f_gauss_peak=stored,
        element_pressure_peak_pa=np.array([123.0]),
    )
    report = diagnose.analyze_result(path)
    assert report["det_f_sampling"] == "stored 8-Gauss-point field"
    np.testing.assert_array_equal(report["det_f"], stored)
    np.testing.assert_array_equal(report["centroid_det_f"], np.ones(1))
    np.testing.assert_array_equal(report["element_pressure_peak_pa"], [123.0])


@pytest.mark.parametrize(
    "spec",
    RELEASE_GUARD.TRUNCATED_POLAR_ARCHIVE_REPORT_SPECS,
    ids=lambda spec: spec["report"].removesuffix(".report.json"),
)
def test_archived_truncated_polar_report_passes_semantic_release_guard(spec):
    RELEASE_GUARD._validate_archived_truncated_polar_report(
        _release_guard_payloads(spec), ROOT, spec
    )


def test_reviewed_reports_preserve_historical_and_recovery_checkpoints():
    historical_reports = {
        "case_a_fbar_1x2x4_dt0p002.report.json",
        "case_b_fbar_2x24x32_dt0p002.report.json",
        "case_b_fbar_2x36x48_dt0p002.report.json",
        "case_b_local_pressure_2x12x16_dt0p002.report.json",
        "case_b_local_pressure_2x12x16_dt0p004.report.json",
        "case_b_local_pressure_2x24x32_dt0p002.report.json",
    }
    specs = RELEASE_GUARD.TRUNCATED_POLAR_ARCHIVE_REPORT_SPECS
    assert {
        spec["report"]
        for spec in specs
        if spec["app_ref"] == RELEASE_GUARD.CURRENT_RESULT_APP_REF
    } == historical_reports
    recovery_specs = [
        spec
        for spec in specs
        if spec["app_ref"] == RELEASE_GUARD.DOMAIN_RECOVERY_APP_REF
    ]
    assert recovery_specs
    assert all(
        spec.get("function_domain_diagnostics") is True
        for spec in recovery_specs
    )
    assert {spec["app_ref"] for spec in specs}.issubset(
        RELEASE_GUARD.CURRENT_RESULT_SOURCE_CHECKPOINTS
    )
    assert (
        RELEASE_GUARD.CURRENT_RESULT_SOURCE_CHECKPOINTS[
            RELEASE_GUARD.CURRENT_RESULT_APP_REF
        ]
        == RELEASE_GUARD.CURRENT_RESULT_SOURCE_SHA256
    )
    reporting_path = "examples/cardiac_benchmark/post.py"
    assert all(
        reporting_path not in hashes
        for hashes in RELEASE_GUARD.CURRENT_RESULT_SOURCE_CHECKPOINTS.values()
    )
    assert reporting_path in RELEASE_GUARD.CURRENT_REPORTING_SOURCE_SHA256


def test_configuration_equivalent_commands_pin_historical_lumped_mass():
    readme = (
        ROOT / "examples" / "cardiac_benchmark" / "results" / "README.md"
    ).read_text()
    section = readme.split("Configuration-equivalent driver commands are:", 1)[1]
    code_block = section.split("```bash", 1)[1].split("```", 1)[0]
    commands = [command for command in code_block.strip().split("\n\n") if command]

    assert len(commands) == 8
    assert all("--mass lumped" in command for command in commands)


def test_current_sources_match_reviewed_release_hashes_without_relabeling_history():
    assert {
        "examples/cardiac_benchmark/boundary_audit.py",
        "examples/cardiac_benchmark/benchmark_parameters.py",
        "examples/cardiac_benchmark/consistent_mass.py",
        "examples/cardiac_benchmark/distributed_local_pressure.py",
        "examples/cardiac_benchmark/distributed_mass.py",
        "examples/cardiac_benchmark/distributed_material.py",
        "examples/cardiac_benchmark/distributed_solver.py",
        "examples/cardiac_benchmark/generalized_alpha.py",
        "examples/cardiac_benchmark/geometry.py",
        "examples/cardiac_benchmark/run.py",
        "examples/cardiac_benchmark/run_mpi.py",
        "examples/cardiac_benchmark/structural_directions.py",
        "examples/cardiac_benchmark/tbar_laplace.py",
        "examples/cardiac_benchmark/viscous_evidence.py",
    }.issubset(RELEASE_GUARD.CURRENT_RELEASE_SOURCE_SHA256)
    payloads = {
        name: (ROOT / name).read_bytes()
        for name in RELEASE_GUARD.CURRENT_RELEASE_SOURCE_SHA256
    }
    RELEASE_GUARD._validate_current_release_source_hashes(payloads, ROOT)
    for historical_path in (
        "examples/cardiac_benchmark/geometry.py",
        "examples/cardiac_benchmark/run.py",
        "examples/cardiac_benchmark/run_mpi.py",
        "examples/cardiac_benchmark/tbar_laplace.py",
    ):
        assert (
            RELEASE_GUARD.CURRENT_RELEASE_SOURCE_SHA256[historical_path]
            != RELEASE_GUARD.STEP0B_PREFIX_DIAGNOSTIC_RUNTIME_SOURCE_MANIFEST[
                historical_path
            ]
        )
    run_path = "examples/cardiac_benchmark/run.py"
    assert (
        RELEASE_GUARD.CURRENT_RELEASE_SOURCE_SHA256[run_path]
        != RELEASE_GUARD.DOMAIN_RECOVERY_SOURCE_SHA256[run_path]
    )


def test_current_reporting_source_is_hashed_separately_from_simulation():
    reporting_paths = {
        "examples/cardiac_benchmark/post.py",
        "examples/cardiac_benchmark/compare_fenics_case_b.py",
        "examples/cardiac_benchmark/compare_mpi_rank_gate.py",
        "examples/cardiac_benchmark/compare_step2b_case_b.py",
        "examples/cardiac_benchmark/plot_step2b_case_b.py",
        "examples/cardiac_benchmark/step2b_case_b_reference_hashes.json",
        "examples/cardiac_benchmark/step2b_case_b_runtime_source_hashes.json",
    }
    assert reporting_paths.issubset(
        RELEASE_GUARD.CURRENT_REPORTING_SOURCE_SHA256
    )
    assert all(
        reporting_paths.isdisjoint(hashes)
        for hashes in RELEASE_GUARD.CURRENT_RESULT_SOURCE_CHECKPOINTS.values()
    )
    payloads = {
        name: (ROOT / name).read_bytes()
        for name in RELEASE_GUARD.CURRENT_REPORTING_SOURCE_SHA256
    }
    RELEASE_GUARD._validate_current_reporting_source_hashes(payloads, ROOT)


def test_packaged_sources_may_match_one_of_multiple_reviewed_checkpoints(
    monkeypatch,
):
    old_payload = b"old source"
    new_payload = b"new source"
    path = "examples/cardiac_benchmark/solver.py"
    checkpoints = {
        "1" * 40: {path: hashlib.sha256(old_payload).hexdigest()},
        "2" * 40: {path: hashlib.sha256(new_payload).hexdigest()},
    }
    monkeypatch.setattr(
        RELEASE_GUARD,
        "CURRENT_RESULT_SOURCE_CHECKPOINTS",
        checkpoints,
    )

    assert RELEASE_GUARD._validate_result_source_checkpoints(
        {path: new_payload}, ROOT
    ) == "2" * 40
    with pytest.raises(SystemExit, match="exactly one reviewed checkpoint"):
        RELEASE_GUARD._validate_result_source_checkpoints(
            {path: b"unreviewed source"}, ROOT
        )


def test_release_guard_accepts_reviewed_function_domain_rejection_evidence():
    spec = copy.deepcopy(
        _reviewed_spec("case_b_local_pressure_2x12x16_dt0p004.report.json")
    )
    report_name = RELEASE_GUARD._truncated_polar_archive_path(spec["report"])
    result = json.loads((ROOT / report_name).read_bytes())["result"]
    spec["function_domain_diagnostics"] = True
    result["solver_configuration"]["function_domain_rejection_api"] = (
        RELEASE_GUARD.PETSC_FUNCTION_DOMAIN_REJECTION_API
    )
    for record in result["nonlinear_step_diagnostics"]:
        record["function_domain_rejections"] = 0
        record["last_function_domain_error"] = None
    result["nonlinear_step_diagnostics"][10].update(
        {
            "function_domain_rejections": 2,
            "last_function_domain_error": "invalid trial det(F)",
        }
    )

    RELEASE_GUARD._validate_current_report_solver(
        result, spec, report_name, ROOT
    )
    result["nonlinear_step_diagnostics"][10][
        "last_function_domain_error"
    ] = None
    with pytest.raises(SystemExit, match="function-domain rejection evidence"):
        RELEASE_GUARD._validate_current_report_solver(
            result, spec, report_name, ROOT
        )


def test_fine_local_pressure_report_retains_recovered_domain_trials():
    spec = _reviewed_spec(
        "case_b_local_pressure_2x36x48_dt0p002.report.json"
    )
    report_name = RELEASE_GUARD._truncated_polar_archive_path(spec["report"])
    result = json.loads((ROOT / report_name).read_bytes())["result"]
    positive_rejections = {
        index: record["function_domain_rejections"]
        for index, record in enumerate(
            result["nonlinear_step_diagnostics"], start=1
        )
        if record["function_domain_rejections"]
    }

    assert spec["app_ref"] == RELEASE_GUARD.DOMAIN_RECOVERY_APP_REF
    assert spec["function_domain_diagnostics"] is True
    assert positive_rejections == {277: 83, 279: 85}
    assert sum(positive_rejections.values()) == 168
    rotation = result["peak_circumferential_ring_rotation"]
    assert [
        layer["relative_apex_minus_base_degrees"]
        for layer in rotation["layers"]
    ] == pytest.approx(
        [-46.410399086878016, -46.32253057989533, -48.67835724664141],
        rel=1.0e-14,
        abs=1.0e-14,
    )
    RELEASE_GUARD._validate_current_report_solver(
        result, spec, report_name, ROOT
    )


def _retained_figure_payloads():
    names = {
        RELEASE_GUARD.RETAINED_FIGURE_RENDERER,
        RELEASE_GUARD.STEP2B_FIGURE_RENDERER,
        RELEASE_GUARD.STEP2B_RAW_STDOUT,
        RELEASE_GUARD.TIP_REFINE_FIGURE_RENDERER,
        RELEASE_GUARD.STEP2B_RERUN_FIGURE_RENDERER,
        *RELEASE_GUARD.RETAINED_FIGURE_SPECS,
        *(
            spec["report"]
            for spec in RELEASE_GUARD.RETAINED_FIGURE_SPECS.values()
        ),
    }
    return {name: (ROOT / name).read_bytes() for name in names}


def _step0b_prefix_payloads():
    names = {
        RELEASE_GUARD.STEP0B_PREFIX_DIAGNOSTIC_REPORT,
        *RELEASE_GUARD.STEP0B_PREFIX_DIAGNOSTIC_RUNTIME_SOURCE_MANIFEST,
    }
    return {name: (ROOT / name).read_bytes() for name in names}


def test_retained_comparison_figures_pass_semantic_release_guard():
    RELEASE_GUARD._validate_retained_figures(
        _retained_figure_payloads(), ROOT
    )


def test_fine_case_a_compact_report_passes_semantic_release_guard():
    payload = (ROOT / RELEASE_GUARD.STEP0A_RETAINED_REPORT).read_bytes()
    assert len(payload) == RELEASE_GUARD.STEP0A_RETAINED_REPORT_SIZE_BYTES
    assert len(payload) < 100_000
    report = json.loads(payload)
    RELEASE_GUARD._validate_step0a_retained_report_semantics(
        report, RELEASE_GUARD.STEP0A_RETAINED_REPORT, ROOT
    )
    assert report["benchmark_identity"]["status"] == "legacy-inferred"
    assert report["result"]["configuration"]["mesh"]["topology"] == (
        "closed_multiblock_disk"
    )
    assert report["result"]["configuration"]["mpi"]["ranks"] == 8


def test_step0b_prefix_diagnostic_passes_semantic_release_guard():
    name = RELEASE_GUARD.STEP0B_PREFIX_DIAGNOSTIC_REPORT
    payloads = _step0b_prefix_payloads()
    payload = payloads[name]
    assert len(payload) == RELEASE_GUARD.STEP0B_PREFIX_DIAGNOSTIC_REPORT_SIZE_BYTES
    assert len(payload) < 20_000
    assert len(RELEASE_GUARD.STEP0B_PREFIX_DIAGNOSTIC_RUNTIME_SOURCE_MANIFEST) == 20
    RELEASE_GUARD._validate_step0b_prefix_diagnostic(payloads, ROOT)
    report = json.loads(payload)
    RELEASE_GUARD._validate_step0b_prefix_diagnostic_semantics(
        report, name, ROOT
    )
    assert name in RELEASE_GUARD.BENCHMARK_RESULT_SOURCES
    assert name in RELEASE_GUARD.PUBLIC_RELEASE_FILES
    assert name in RELEASE_GUARD.EXPECTED_SDIST_FILES
    assert report["decision"]["full_1s_status"] == "paused"
    assert report["completion"]["coarse_2x20x17"]["status"] == "complete"
    assert report["completion"]["wall_only_4x20x17"]["status"] == "complete"
    assert report["configuration"]["fiber_direction_reconstruction"] == (
        "toolkit-physical-coordinate-u-v-v1"
    )


def test_step0b_prefix_diagnostic_missing_file_fails_closed():
    with pytest.raises(SystemExit, match="missing the Step 0B prefix"):
        RELEASE_GUARD._validate_step0b_prefix_diagnostic({}, ROOT)


def test_step0b_prefix_diagnostic_changed_bytes_fail_closed():
    name = RELEASE_GUARD.STEP0B_PREFIX_DIAGNOSTIC_REPORT
    payload = (ROOT / name).read_bytes() + b"\n"
    with pytest.raises(SystemExit, match="differs from the reviewed compact report"):
        RELEASE_GUARD._validate_step0b_prefix_diagnostic(
            {name: payload}, ROOT
        )


def test_step0b_prefix_runtime_manifest_digest_fails_closed(monkeypatch):
    payloads = _step0b_prefix_payloads()
    manifest = dict(
        RELEASE_GUARD.STEP0B_PREFIX_DIAGNOSTIC_RUNTIME_SOURCE_MANIFEST
    )
    manifest["examples/cardiac_benchmark/run_mpi.py"] = "0" * 64
    monkeypatch.setattr(
        RELEASE_GUARD,
        "STEP0B_PREFIX_DIAGNOSTIC_RUNTIME_SOURCE_MANIFEST",
        manifest,
    )
    with pytest.raises(SystemExit, match="manifest digest differs"):
        RELEASE_GUARD._validate_step0b_prefix_diagnostic(payloads, ROOT)


def test_step0b_prefix_current_runtime_source_byte_fails_closed():
    payloads = _step0b_prefix_payloads()
    name = "examples/cardiac_benchmark/run_mpi.py"
    payloads[name] += b"\n"
    with pytest.raises(
        SystemExit, match="current Step 0B runtime paths differ"
    ):
        RELEASE_GUARD._validate_step0b_prefix_diagnostic(payloads, ROOT)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("not-object", "must contain a JSON object"),
        ("schema", "altered Step 0B schema"),
        ("benchmark-case", "altered Step 0B benchmark identity"),
        ("pressure-mode", "altered Step 0B benchmark identity"),
        ("application-revision", "altered Step 0B source identity"),
        ("core-revision", "altered Step 0B source identity"),
        ("core-tree-state", "altered Step 0B source identity"),
        ("runtime-source", "altered Step 0B source identity"),
        ("configuration", "altered Step 0B configuration"),
        ("generalized-alpha", "altered Step 0B configuration"),
        ("coarse-steps", "altered Step 0B coarse_2x20x17 completion"),
        ("coarse-ranks", "altered Step 0B coarse_2x20x17 completion"),
        ("wall-convergence", "altered Step 0B wall_only_4x20x17 completion"),
        ("wall-audit", "altered Step 0B wall_only_4x20x17 completion"),
        ("wall-domain-rejections", "altered Step 0B wall_only_4x20x17 completion"),
        ("decision", "altered Step 0B paused decision"),
    ],
)
def test_step0b_prefix_diagnostic_semantic_broken_controls_fail_closed(
    mutation, message
):
    name = RELEASE_GUARD.STEP0B_PREFIX_DIAGNOSTIC_REPORT
    report = json.loads((ROOT / name).read_bytes())
    if mutation == "not-object":
        report = []
    elif mutation == "schema":
        report["schema"] = "altered"
    elif mutation == "benchmark-case":
        report["benchmark_identity"]["case"] = "step_0A"
    elif mutation == "pressure-mode":
        report["benchmark_identity"]["pressure_enabled"] = False
    elif mutation == "application-revision":
        report["source_identity"]["application_revision"] = "0" * 40
    elif mutation == "core-revision":
        report["source_identity"]["core_revision"] = "0" * 40
    elif mutation == "core-tree-state":
        report["source_identity"]["core_tree_state"] = "dirty"
    elif mutation == "runtime-source":
        report["source_identity"]["benchmark_runtime_source_sha256"] = "0" * 64
    elif mutation == "configuration":
        report["configuration"]["t_end_s"] = 1.0
    elif mutation == "generalized-alpha":
        report["configuration"]["generalized_alpha"]["alpha_m"] = 0.1
    elif mutation == "coarse-steps":
        report["completion"]["coarse_2x20x17"]["completed_steps"] = 319
    elif mutation == "coarse-ranks":
        report["completion"]["coarse_2x20x17"]["mpi_ranks"] = 4
    elif mutation == "wall-convergence":
        report["completion"]["wall_only_4x20x17"]["converged"] = False
    elif mutation == "wall-audit":
        report["completion"]["wall_only_4x20x17"][
            "pre_solve_audits_passed"
        ] = False
    elif mutation == "wall-domain-rejections":
        report["completion"]["wall_only_4x20x17"][
            "function_domain_rejections_total"
        ] = 1
    elif mutation == "decision":
        report["decision"]["full_1s_status"] = "approved"
    else:  # pragma: no cover - test helper programming error
        raise AssertionError(mutation)

    with pytest.raises(SystemExit, match=message):
        RELEASE_GUARD._validate_step0b_prefix_diagnostic_semantics(
            report, name, ROOT
        )


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("result-hash", "altered fine Case A result identity"),
        ("source", "altered fine Case A result identity"),
        ("mesh", "altered retained configuration"),
        ("sampling", "altered retained configuration"),
        ("reference", "altered external reference identity"),
        ("curve", "altered retained p0 curves"),
        ("red", "unreproducible p0 RED"),
        ("team-red", "altered reference-team p0 RED"),
        ("benchmark-identity", "legacy-inferred Step 0A identity"),
    ],
)
def test_fine_case_a_compact_report_broken_controls_fail_closed(
    mutation, message
):
    report = json.loads(
        (ROOT / RELEASE_GUARD.STEP0A_RETAINED_REPORT).read_bytes()
    )
    if mutation == "result-hash":
        report["result"]["sha256"] = "0" * 64
    elif mutation == "source":
        report["result"]["source_identity"]["app"]["revision"] = "0" * 40
    elif mutation == "mesh":
        report["result"]["configuration"]["mesh"]["n_radial"] = 31
    elif mutation == "sampling":
        report["result"]["configuration"]["sampling_points"]["p0"][
            "element"
        ] += 1
    elif mutation == "reference":
        report["reference"]["team_files"][0]["sha256"] = "0" * 64
    elif mutation == "curve":
        report["comparison"]["ours_on_canonical_grid_m"]["p0"][48][0] += 1e-6
    elif mutation == "red":
        report["comparison"]["red"]["p0"]["ours"] += 0.01
    elif mutation == "team-red":
        report["comparison"]["red"]["p0"]["teams"]["simula"] += 0.01
    elif mutation == "benchmark-identity":
        report["benchmark_identity"]["status"] = "recorded"
    else:  # pragma: no cover - test helper programming error
        raise AssertionError(mutation)

    with pytest.raises(SystemExit, match=message):
        RELEASE_GUARD._validate_step0a_retained_report_semantics(
            report, RELEASE_GUARD.STEP0A_RETAINED_REPORT, ROOT
        )


def test_step2b_normalized_stdout_is_hash_bound():
    payloads = _retained_figure_payloads()
    payloads[RELEASE_GUARD.STEP2B_RAW_STDOUT] += b"altered\n"
    with pytest.raises(SystemExit, match="reviewed normalized stdout"):
        RELEASE_GUARD._validate_retained_figures(payloads, ROOT)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("title", "accessible SVG title"),
        ("script", "active or external SVG markup"),
        ("external-link", "external SVG link"),
        ("doctype", "active or external SVG markup"),
        ("provenance", "reviewed provenance labels"),
    ],
)
def test_retained_svg_broken_controls_fail_closed(mutation, message):
    name = "docs/figures/case_a_comparison.svg"
    spec = RELEASE_GUARD.RETAINED_FIGURE_SPECS[name]
    payload = (ROOT / name).read_bytes()
    if mutation == "title":
        payload = payload.replace(
            spec["title"].encode("utf-8"), b"altered title", 1
        )
    elif mutation == "script":
        payload = payload.replace(b"<defs>", b"<script>bad()</script><defs>", 1)
    elif mutation == "external-link":
        payload = payload.replace(
            b'xlink:href="#', b'xlink:href="https://example.invalid/', 1
        )
    elif mutation == "doctype":
        payload = payload.replace(
            b"<svg ", b'<!DOCTYPE svg SYSTEM "https://example.invalid/x">\n<svg ', 1
        )
    elif mutation == "provenance":
        source_marker = spec["visible_markers"][1].encode("utf-8")
        payload = payload.replace(source_marker, b"app altered", 1)

    with pytest.raises(SystemExit, match=message):
        RELEASE_GUARD._validate_retained_svg(payload, name, spec, ROOT)


def test_figure_packaging_inventory_is_exact_and_wheel_remains_metadata_only(
    tmp_path,
):
    figure_sources = set(RELEASE_GUARD.RETAINED_FIGURE_SPECS)
    source_only = {
        RELEASE_GUARD.RETAINED_FIGURE_RENDERER,
        RELEASE_GUARD.STEP2B_FIGURE_RENDERER,
        "docs/figures/README.md",
        *figure_sources,
    }
    assert source_only.issubset(RELEASE_GUARD.PUBLIC_RELEASE_FILES)
    assert source_only.issubset(RELEASE_GUARD.EXPECTED_SDIST_FILES)
    assert all(
        not RELEASE_GUARD._is_forbidden_path(name) for name in figure_sources
    )
    assert RELEASE_GUARD._is_forbidden_path("docs/figures/unreviewed.svg")
    assert RELEASE_GUARD._is_forbidden_path("docs/figures/unreviewed.png")

    wheel = tmp_path / "figure-leak.whl"
    with zipfile.ZipFile(wheel, mode="w") as archive:
        archive.writestr("coupfe_cardiac-0.1.0.dist-info/METADATA", "metadata")
        archive.writestr(
            "docs/figures/case_a_comparison.svg",
            (ROOT / "docs/figures/case_a_comparison.svg").read_bytes(),
        )
    with pytest.raises(SystemExit, match="metadata-only wheel policy"):
        RELEASE_GUARD._validate_wheel(
            wheel, allow_unapproved_core_ref=False
        )


def _broken_report_payload(report_basename, mutation):
    spec = copy.deepcopy(_reviewed_spec(report_basename))
    payloads = _release_guard_payloads(spec)
    report_name = RELEASE_GUARD._truncated_polar_archive_path(spec["report"])
    report = json.loads(payloads[report_name])

    if mutation == "nonfinite":
        report["result"]["det_f_gauss_peak_summary"]["mean"] = float("nan")
    elif mutation == "absolute-path":
        report["result"]["source_identity"]["core"]["source_url"] = (
            "/" + "home/private/core"
        )
    elif mutation == "nonpositive-det-f":
        report["result"]["det_f_gauss_peak_summary"]["minimum"] = 0.0
    elif mutation == "petsc-threshold":
        diagnostic = report["result"]["nonlinear_step_diagnostics"][0]
        diagnostic["residual_acceptance_threshold"] *= 2.0
    elif mutation == "red":
        report["comparison"]["red"]["p0"]["ours"] += 0.01
    elif mutation == "team-identity":
        report["reference"]["team_files"][0]["sha256"] = "0" * 64
    elif mutation == "selection":
        report["reference"]["selection"]["selected_count"] = 11
    elif mutation == "correction":
        report["correction"]["supersedes_report_sha256"] = "0" * 64
    elif mutation == "ring-rotation":
        report["result"]["peak_circumferential_ring_rotation"]["layers"][0][
            "relative_apex_minus_base_degrees"
        ] += 1.0
    else:  # pragma: no cover - test helper programming error
        raise AssertionError(mutation)

    rendered = (
        json.dumps(report, indent=2, sort_keys=True, allow_nan=True) + "\n"
    ).encode("utf-8")
    payloads[report_name] = rendered
    spec["report_sha256"] = hashlib.sha256(rendered).hexdigest()
    return payloads, spec


@pytest.mark.parametrize(
    ("report_basename", "mutation", "message"),
    [
        (
            "case_a_fbar_1x2x4_dt0p002.report.json",
            "nonfinite",
            "strict finite JSON",
        ),
        (
            "case_a_fbar_1x2x4_dt0p002.report.json",
            "absolute-path",
            "absolute filesystem path",
        ),
        (
            "case_b_local_pressure_2x12x16_dt0p004.report.json",
            "nonpositive-det-f",
            r"positive finite 8-GP det\(F\)",
        ),
        (
            "case_b_local_pressure_2x12x16_dt0p004.report.json",
            "petsc-threshold",
            "independently recomputed PETSc residual rule",
        ),
        (
            "case_b_fbar_2x24x32_dt0p002.report.json",
            "red",
            "unreproducible p0 RED",
        ),
        (
            "case_b_fbar_2x24x32_dt0p002.report.json",
            "team-identity",
            "altered external reference identity",
        ),
        (
            "case_b_fbar_2x24x32_dt0p002.report.json",
            "selection",
            "altered external reference identity",
        ),
        (
            "case_b_fbar_2x24x32_dt0p002.report.json",
            "correction",
            "altered correction lineage",
        ),
        (
            "case_b_local_pressure_2x36x48_dt0p002.report.json",
            "ring-rotation",
            "invalid ring rotation",
        ),
    ],
)
def test_archived_truncated_polar_report_release_guard_broken_controls(
    report_basename, mutation, message
):
    payloads, spec = _broken_report_payload(report_basename, mutation)
    with pytest.raises(SystemExit, match=message):
        RELEASE_GUARD._validate_archived_truncated_polar_report(
            payloads, ROOT, spec
        )
