"""Compare a completed run with the separately distributed benchmark curves.

Reference data are not vendored. Download the CC BY 4.0 dataset from
https://doi.org/10.5281/zenodo.14260459 and pass either its extracted root or
the ``results_time_curves/data`` directory with ``--reference-dir`` (or the
``CARDIAC_BENCHMARK_DATA_DIR`` environment variable).

The reference files use Python pickle. Load only the trusted Zenodo archive;
pickle is not safe for untrusted input. Every selected file is validated and
hashed. Selection follows the exact 10-file Case A/Case B manifests in the
archive's ``results_time_curves/figures.py``; missing files and unexpected
wildcard matches abort the comparison instead of quietly changing the
reference-team population. The archive's unselected base-name SimVascular
alias is accepted only when it is byte-identical to selected SimVascular P2.
"""
from __future__ import annotations

import argparse
from collections.abc import Mapping
import hashlib
import json
from numbers import Real
import os
import pickle
import re
from pathlib import Path, PureWindowsPath
import tempfile

import numpy as np

try:  # package import
    from .activation import p_of_t, tau_of_t
    from .benchmark_parameters import (
        benchmark_configuration,
        benchmark_metadata,
        validate_load_histories,
    )
except ImportError:  # direct script import
    from activation import p_of_t, tau_of_t
    from benchmark_parameters import (
        benchmark_configuration,
        benchmark_metadata,
        validate_load_histories,
    )


RESULT_SCHEMA = "coupfe-cardiac-result-v1"
REPORT_SCHEMA = "coupfe-cardiac-reference-comparison-v2"
CASE_NAMES = {"A": "step_0A", "B": "step_0B"}
BENCHMARK_ARCHIVE_FIELDS = frozenset(
    {
        "benchmark_step",
        "benchmark_configuration_id",
        "benchmark_identity_scope",
        "benchmark_load_contract",
        "benchmark_peak_load_definition",
        "benchmark_active_stress_enabled",
        "benchmark_pressure_enabled",
        "benchmark_material_parameters_json",
        "benchmark_activation_parameters_json",
        "benchmark_pressure_parameters_json",
        "benchmark_runtime_source_manifest_json",
        "benchmark_runtime_source_sha256",
    }
)
FORMULATIONS = {
    "hex8_fbar",
    "hex8_local_pressure_p0_condensed_logj",
    "hex8_local_pressure_p0_condensed_mean_logj_paper_j2",
    "hex8_standard_pointwise_kappa",
}
LOCAL_PRESSURE_FORMULATIONS = frozenset(
    {
        "hex8_local_pressure_p0_condensed_logj",
        "hex8_local_pressure_p0_condensed_mean_logj_paper_j2",
    }
)
FIBER_SAMPLING_METHODS = {
    "cg1_gram_schmidt",
    "gp_direct_rule",
}
FIBER_SAMPLING_OPTIONS = {
    "cg1": "cg1_gram_schmidt",
    "gp-direct": "gp_direct_rule",
}
MASS_REPRESENTATIONS = {
    "consistent_q1_hex8",
    "lumped_row_sum",
}
MESH_TOPOLOGIES = {
    "polar_ring",
    "closed_multiblock_disk",
}
POINT_SAMPLING_METHODS = {
    "global_delaunay_tetra",
    "hex8_reference_isoparametric",
}
CURRENT_POINT_SAMPLING = "hex8_reference_isoparametric"
PETSC_FUNCTION_DOMAIN_REJECTION_API = "nonfinite residual for PETSc BT"
TBAR_METADATA_SCHEMA = "coupfe-cardiac-laplace-tbar-v1"
COMPLETE_SWITCH_ENERGY_MATERIAL_MODEL_ID = (
    "holzapfel-ogden-smooth-switch-complete-energy-derivative-v1"
)
LEGACY_SWITCH_STRESS_MATERIAL_MODEL_ID = (
    "holzapfel-ogden-smooth-switch-stress-without-switch-derivative-v0"
)
SUPPORTED_MATERIAL_MODEL_IDS = frozenset(
    {
        COMPLETE_SWITCH_ENERGY_MATERIAL_MODEL_ID,
        LEGACY_SWITCH_STRESS_MATERIAL_MODEL_ID,
    }
)
MPI_COMPANION_IMPLEMENTATION = "cardiac-owned-distributed-q1p0-step0"
MPI_CLOSED_STD_KAPPA_IMPLEMENTATION = (
    "cardiac-owned-distributed-closed-std-kappa-step0"
)
MPI_CLOSED_LOCAL_PRESSURE_IMPLEMENTATION = (
    "cardiac-owned-distributed-closed-local-pressure-step0"
)
MPI_CLOSED_LOCAL_PRESSURE_PAPER_IMPLEMENTATION = (
    "cardiac-owned-distributed-closed-local-pressure-mean-logj-"
    "paper-j2-step0"
)
MPI_CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION = (
    "cardiac-owned-distributed-closed-std-kappa-generalized-alpha-step0"
)
MPI_CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION = (
    "cardiac-owned-distributed-closed-local-pressure-generalized-alpha-step0"
)
MPI_CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION = (
    "cardiac-owned-distributed-closed-local-pressure-mean-logj-"
    "paper-j2-generalized-alpha-step0"
)
MPI_CLOSED_GENERALIZED_ALPHA_IMPLEMENTATIONS = frozenset(
    {
        MPI_CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION,
        MPI_CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION,
        MPI_CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION,
    }
)
MPI_CLOSED_IMPLEMENTATIONS = frozenset(
    {
        MPI_CLOSED_STD_KAPPA_IMPLEMENTATION,
        MPI_CLOSED_LOCAL_PRESSURE_IMPLEMENTATION,
        MPI_CLOSED_LOCAL_PRESSURE_PAPER_IMPLEMENTATION,
        *MPI_CLOSED_GENERALIZED_ALPHA_IMPLEMENTATIONS,
    }
)
MPI_COMPANION_IMPLEMENTATIONS = frozenset(
    {MPI_COMPANION_IMPLEMENTATION, *MPI_CLOSED_IMPLEMENTATIONS}
)
MPI_COMPANION_ARCHIVE_FIELDS = frozenset(
    {
        "mpi_enabled",
        "mpi_ranks",
        "mpi_world_size",
        "mpi_local_element_counts",
        "mpi_implementation",
        "mpi_partition",
        "mpi_build_layout",
        "mpi_factor_solver_type",
    }
)
MPI_COMPANION_OPTIONAL_ARCHIVE_FIELDS = frozenset(
    {"mpi_linear_solver_profile"}
)
MPI_MASS_PARTITION_ARCHIVE_FIELDS = frozenset(
    {
        "mpi_mass_partition",
        "mpi_mass_owned_row_ranges",
        "mpi_mass_local_nnz",
        "mpi_mass_touching_element_counts",
    }
)
MPI_COMPANION_FIXED_CONFIGURATION = {
    "implementation": MPI_COMPANION_IMPLEMENTATION,
    "settings_source": "recovered 2026-06-27 Case B development adapter",
    "communicator": "PETSc.COMM_WORLD",
    "matrix_scope": "persistent one solver instance per MPI run",
    "dirichlet_support": "none",
    "rank_local_element_assembly": True,
    "global_mesh_replicated": True,
    "collective_invalid_trial_policy": "all-owned-residual-entries-inf-for-bt",
    "function_domain_rejection_api": PETSC_FUNCTION_DOMAIN_REJECTION_API,
    "commit_policy": "independent-global-residual-check-before-commit",
    "snes_type": "newtonls",
    "line_search_type": "bt",
    "rtol": 1.0e-9,
    "atol": 1.0e-10,
    "stol": 1.0e-12,
    "max_it": 60,
}
MPI_DIRECT_SUPERLU_DIST_PROFILE = "direct-superlu-dist"
MPI_FGMRES_GAMG_RIGID_PROFILE = "fgmres-gamg-rigid"
MPI_FGMRES_GAMG_RIGID_REBUILD_PROFILE = "fgmres-gamg-rigid-rebuild"
MPI_FGMRES_ASM_LU_PROFILE = "fgmres-asm-lu"
MPI_FGMRES_ASM_ILU1_PROFILE = "fgmres-asm-ilu1"
MPI_COMPANION_LEGACY_DIRECT_CONFIGURATION = {
    "factor_solver_type": "superlu_dist",
    "configured_factor_solver_type": "superlu_dist",
    "ksp_type": "preonly",
    "pc_type": "lu",
}
MPI_COMPANION_LINEAR_SOLVER_CONFIGURATIONS = {
    MPI_DIRECT_SUPERLU_DIST_PROFILE: {
        "factor_solver_type": "superlu_dist",
        "configured_factor_solver_type": "superlu_dist",
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
    },
    MPI_FGMRES_GAMG_RIGID_PROFILE: {
        "factor_solver_type": "none",
        "configured_factor_solver_type": "not-applicable",
        "node_aligned_ownership": True,
        "vector_block_size": 3,
        "matrix_block_size": 3,
        "near_nullspace_kind": "six-rigid-body-modes",
        "near_nullspace_mode_count": 6,
        "ksp_type": "fgmres",
        "pc_type": "gamg",
        "ksp_rtol": 1.0e-8,
        "ksp_atol": 1.0e-12,
        "ksp_divtol": 1.0e4,
        "ksp_max_it": 200,
        "gmres_restart": 50,
        "ksp_norm_type": "unpreconditioned",
        "pc_side": "right",
        "ksp_error_if_not_converged": False,
        "preconditioner": (
            "PETSc GAMG aggregation with six rigid-body near-null modes"
        ),
    },
    MPI_FGMRES_GAMG_RIGID_REBUILD_PROFILE: {
        "factor_solver_type": "none",
        "configured_factor_solver_type": "not-applicable",
        "node_aligned_ownership": True,
        "vector_block_size": 3,
        "matrix_block_size": 3,
        "near_nullspace_kind": "six-rigid-body-modes",
        "near_nullspace_mode_count": 6,
        "ksp_type": "fgmres",
        "pc_type": "gamg",
        "ksp_rtol": 1.0e-8,
        "ksp_atol": 1.0e-12,
        "ksp_divtol": 1.0e4,
        "ksp_max_it": 400,
        "gmres_restart": 100,
        "ksp_norm_type": "unpreconditioned",
        "pc_side": "right",
        "ksp_error_if_not_converged": False,
        "preconditioner": (
            "PETSc GAMG aggregation with six rigid-body near-null modes and "
            "interpolation rebuilt for changed matrices"
        ),
    },
    MPI_FGMRES_ASM_LU_PROFILE: {
        "factor_solver_type": "none",
        "configured_factor_solver_type": "not-applicable",
        "node_aligned_ownership": False,
        "vector_block_size": 1,
        "matrix_block_size": 1,
        "near_nullspace_kind": "none",
        "near_nullspace_mode_count": 0,
        "ksp_type": "fgmres",
        "pc_type": "asm",
        "ksp_rtol": 1.0e-8,
        "ksp_atol": 1.0e-12,
        "ksp_divtol": 1.0e4,
        "ksp_max_it": 200,
        "gmres_restart": 50,
        "ksp_norm_type": "unpreconditioned",
        "pc_side": "right",
        "ksp_error_if_not_converged": False,
        "preconditioner": (
            "restricted additive Schwarz overlap 1 with local SuperLU"
        ),
    },
    MPI_FGMRES_ASM_ILU1_PROFILE: {
        "factor_solver_type": "none",
        "configured_factor_solver_type": "not-applicable",
        "node_aligned_ownership": False,
        "vector_block_size": 1,
        "matrix_block_size": 1,
        "near_nullspace_kind": "none",
        "near_nullspace_mode_count": 0,
        "ksp_type": "fgmres",
        "pc_type": "asm",
        "ksp_rtol": 1.0e-8,
        "ksp_atol": 1.0e-12,
        "ksp_divtol": 1.0e4,
        "ksp_max_it": 200,
        "gmres_restart": 50,
        "ksp_norm_type": "unpreconditioned",
        "pc_side": "right",
        "ksp_error_if_not_converged": False,
        "preconditioner": (
            "restricted additive Schwarz overlap 1 with local ILU(1)"
        ),
    },
}
MPI_COMPANION_LINEAR_SOLVER_PETSC_OPTIONS = {
    MPI_FGMRES_GAMG_RIGID_PROFILE: {
        "pc_gamg_type": "agg",
        "pc_gamg_agg_nsmooths": "1",
        "pc_gamg_threshold": "0.01",
        "pc_gamg_repartition": "false",
    },
    MPI_FGMRES_GAMG_RIGID_REBUILD_PROFILE: {
        "pc_gamg_type": "agg",
        "pc_gamg_agg_nsmooths": "1",
        "pc_gamg_threshold": "0.01",
        "pc_gamg_repartition": "false",
        "pc_gamg_reuse_interpolation": "false",
    },
    MPI_FGMRES_ASM_LU_PROFILE: {
        "sub_ksp_type": "preonly",
        "sub_pc_type": "lu",
        "sub_pc_factor_mat_solver_type": "superlu",
    },
    MPI_FGMRES_ASM_ILU1_PROFILE: {
        "sub_ksp_type": "preonly",
        "sub_pc_type": "ilu",
        "sub_pc_factor_levels": "1",
        "sub_pc_factor_shift_type": "nonzero",
    },
}
FULL_REVISION = re.compile(r"[0-9a-fA-F]{40}")
PUBLIC_CORE_URL = "https://github.com/tengzhang48/CoupFE.git"
PUBLIC_CORE_REF = "e2f42ed5772850a0a23a2ce434f430c287eae5c8"
SUPPORTED_PUBLIC_CORE_REFS = frozenset(
    {
        # Retained public benchmark records.
        "454f73ce2de284262b214a2b37bd676c6aca3c0a",
        # Current public dependency with explicit native joint/split support.
        PUBLIC_CORE_REF,
    }
)
IMMUTABLE_PRE_METHOD_METADATA_APP_REFS = frozenset(
    {
        "62ad760d2a1731bb9668897863ac026d3768194e",
        "e07993bcf1166bd20eb87370c0b458552753e7ee",
        # Clean checkpoint that corrected the complete smooth-switch energy
        # derivative before portable mass/tbar metadata was added.  Results
        # from this revision remain distinguishable from the older 62ad760
        # constitutive law through their exact source identity.
        "6839c13b5bc80ec06c897684c51f503e80bd4b19",
    }
)
REFERENCE_DOI = "10.5281/zenodo.14260459"
REFERENCE_LICENSE = "CC-BY-4.0"
REFERENCE_ARCHIVE_IDENTITY = {
    "filename": "benchmark_article_data.zip",
    "size_bytes": 23180741494,
    "md5": "75602be4777c4ca2262c2bcfd2134b15",
    "sha256": "134951af5e38d147b0223f0a83666eb3fe1b75acb5bfa9f1b9aa30f255f8f1f5",
}
REFERENCE_FIGURES_PY_IDENTITY = {
    "filename": "results_time_curves/figures.py",
    "size_bytes": 24076,
    "sha256": "f8f519b357349341207faea4b57bfafc1a311aadccefe4968ecbdb37339c8a5b",
}
REFERENCE_SELECTION_POLICY = (
    "Select the exact 10 Case A/Case B team files explicitly loaded by upstream "
    "results_time_curves/figures.py. Reject missing selected files and unexpected "
    "matching files. Accept the unselected base-name SimVascular file only when "
    "byte-identical to selected SimVascular P2, then exclude it as a duplicate alias."
)
REFERENCE_MANIFEST_SUFFIXES = (
    "carpentry",
    "ambit",
    "4C",
    "simula",
    "chimera",
    "cheart",
    "lifex",
    "simvascular_p1p1",
    "simvascular_p2",
    "comsol",
)
REFERENCE_MANIFEST_VARIABLES = {
    "step_0A": "TEAMS_DATASETS_0A",
    "step_0B": "TEAMS_DATASET_0B",
}
REFERENCE_EXCLUDED_ALIAS_SUFFIX = "simvascular"
REFERENCE_ALIAS_TARGET_SUFFIX = "simvascular_p2"
CORRECTION_REASON = (
    "The predecessor selected 11 files with a wildcard and double-counted the "
    "byte-identical SimVascular base-name alias and SimVascular P2 curve. This "
    "successor selects the 10 files explicitly used by upstream "
    "results_time_curves/figures.py."
)
CORRECTION_PREDECESSOR_REVISION = "7f8e726cd2a79ae2ad13ebac4d9c39bca5cec8b2"
CANONICAL_TIME_GRID = np.linspace(0.0, 1.0, 101)
_ENDPOINT_OFFSET_TOLERANCE = 1.000001e-3


def _scalar(archive, key, result_path):
    if key not in archive:
        raise RuntimeError(f"{result_path} is missing required result field {key!r}")
    value = np.asarray(archive[key])
    if value.shape != ():
        raise RuntimeError(f"{result_path} has non-scalar result field {key!r}")
    return value.item()


def _boolean_scalar(archive, key, result_path):
    value = _scalar(archive, key, result_path)
    if not isinstance(value, (bool, np.bool_)):
        raise RuntimeError(f"{result_path} has non-boolean result field {key!r}")
    return bool(value)


def _text_scalar(archive, key, result_path):
    value = _scalar(archive, key, result_path)
    if not isinstance(value, (str, np.str_)):
        raise RuntimeError(f"{result_path} has non-text result field {key!r}")
    return str(value)


def _finite_array(archive, key, result_path, *, shape=None):
    if key not in archive:
        raise RuntimeError(f"{result_path} is missing required result field {key!r}")
    value = np.asarray(archive[key], dtype=float)
    if shape is not None and value.shape != shape:
        raise RuntimeError(
            f"{result_path} has invalid {key!r} shape {value.shape}; expected {shape}"
        )
    if not np.all(np.isfinite(value)):
        raise RuntimeError(f"{result_path} has non-finite values in {key!r}")
    return value


def _integer(value, description, result_path, *, minimum=0):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise RuntimeError(f"{result_path} has invalid {description}")
    numeric = float(value)
    if not np.isfinite(numeric) or not numeric.is_integer() or numeric < minimum:
        raise RuntimeError(f"{result_path} has invalid {description}")
    return int(numeric)


def _finite_number(value, description, result_path, *, minimum=None, strict=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise RuntimeError(f"{result_path} has invalid {description}")
    numeric = float(value)
    if not np.isfinite(numeric):
        raise RuntimeError(f"{result_path} has non-finite {description}")
    if minimum is not None and (
        numeric <= minimum if strict else numeric < minimum
    ):
        raise RuntimeError(f"{result_path} has invalid {description}")
    return numeric


def _finite_json_vector(value, description, result_path, *, size=3):
    if (
        not isinstance(value, list)
        or len(value) != size
        or any(
            isinstance(item, (bool, np.bool_)) or not isinstance(item, Real)
            for item in value
        )
    ):
        raise RuntimeError(f"{result_path} has invalid {description}")
    vector = np.asarray(value, dtype=float)
    if vector.shape != (size,) or not np.all(np.isfinite(vector)):
        raise RuntimeError(f"{result_path} has invalid {description}")
    return vector


def _validate_finite_json(value, description, result_path):
    """Reject JSON NaN/Infinity and unsupported nested value types."""
    if value is None or isinstance(value, (str, bool)):
        return
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        if not np.isfinite(float(value)):
            raise RuntimeError(f"{result_path} has non-finite {description}")
        return
    if isinstance(value, list):
        for item in value:
            _validate_finite_json(item, description, result_path)
        return
    if isinstance(value, dict):
        if not all(isinstance(key, str) for key in value):
            raise RuntimeError(f"{result_path} has malformed {description}")
        for item in value.values():
            _validate_finite_json(item, description, result_path)
        return
    raise RuntimeError(f"{result_path} has malformed {description}")


def _validate_no_absolute_paths(value, description="report"):
    """Ensure a durable report cannot disclose a machine-local absolute path."""
    if isinstance(value, str):
        if Path(value).is_absolute() or PureWindowsPath(value).is_absolute():
            raise RuntimeError(f"{description} contains an absolute filesystem path")
        return
    if isinstance(value, list):
        for item in value:
            _validate_no_absolute_paths(item, description)
    elif isinstance(value, dict):
        for item in value.values():
            _validate_no_absolute_paths(item, description)


def _parse_json_field(archive, key, result_path, expected_type):
    raw = _scalar(archive, key, result_path)
    if not isinstance(raw, (str, np.str_)):
        raise RuntimeError(f"{result_path} has non-text JSON field {key!r}")
    try:
        value = json.loads(str(raw))
    except (json.JSONDecodeError, TypeError) as error:
        raise RuntimeError(f"{result_path} has malformed JSON field {key!r}") from error
    if not isinstance(value, expected_type):
        raise RuntimeError(f"{result_path} has malformed JSON field {key!r}")
    _validate_finite_json(value, key, result_path)
    return value


def _validate_benchmark_configuration_metadata(archive, result_path, case):
    """Validate exact Benchmark 1 mode metadata when the archive records it."""
    present = BENCHMARK_ARCHIVE_FIELDS.intersection(archive.files)
    if not present:
        return None
    if present != BENCHMARK_ARCHIVE_FIELDS:
        raise RuntimeError(
            f"{result_path} has incomplete benchmark-configuration metadata; "
            f"missing {sorted(BENCHMARK_ARCHIVE_FIELDS - present)}"
        )
    benchmark_step = _integer(
        _scalar(archive, "benchmark_step", result_path),
        "benchmark step",
        result_path,
    )
    try:
        configuration = benchmark_configuration(benchmark_step, case)
    except ValueError as error:
        raise RuntimeError(
            f"{result_path} has an unsupported benchmark step/case identity"
        ) from error
    material = _parse_json_field(
        archive,
        "benchmark_material_parameters_json",
        result_path,
        dict,
    )
    activation = _parse_json_field(
        archive,
        "benchmark_activation_parameters_json",
        result_path,
        dict,
    )
    pressure = _parse_json_field(
        archive,
        "benchmark_pressure_parameters_json",
        result_path,
        dict,
    )
    try:
        expected = benchmark_metadata(
            configuration,
            material_parameters=material,
            activation_parameters=activation,
            pressure_parameters=pressure,
        )
    except RuntimeError as error:
        raise RuntimeError(
            f"{result_path} has benchmark parameter metadata inconsistent "
            "with its step/case identity"
        ) from error
    for field in (
        "benchmark_configuration_id",
        "benchmark_identity_scope",
        "benchmark_load_contract",
        "benchmark_peak_load_definition",
        "benchmark_runtime_source_manifest_json",
        "benchmark_runtime_source_sha256",
    ):
        if _text_scalar(archive, field, result_path) != expected[field]:
            raise RuntimeError(
                f"{result_path} has inconsistent benchmark field {field!r}"
            )
    for field in (
        "benchmark_active_stress_enabled",
        "benchmark_pressure_enabled",
    ):
        if _boolean_scalar(archive, field, result_path) is not expected[field]:
            raise RuntimeError(
                f"{result_path} has inconsistent benchmark field {field!r}"
            )
    return {
        "benchmark": 1,
        "step": benchmark_step,
        "case": case,
        "configuration_id": configuration.identity,
        "identity_scope": expected["benchmark_identity_scope"],
        "load_contract": configuration.load_contract,
        "peak_load_definition": expected["benchmark_peak_load_definition"],
        "active_stress_enabled": configuration.active_stress_enabled,
        "pressure_enabled": configuration.pressure_enabled,
        "material_parameters": material,
        "activation_parameters": activation,
        "pressure_parameters": pressure,
        "runtime_source_manifest": json.loads(
            expected["benchmark_runtime_source_manifest_json"]
        ),
        "runtime_source_sha256": expected["benchmark_runtime_source_sha256"],
        "_configuration": configuration,
    }


def _validate_solver_metadata(archive, result_path, times, *, require):
    fields = {
        "nonlinear_solver",
        "solver_configuration_json",
        "nonlinear_step_diagnostics_json",
    }
    present = {field for field in fields if field in archive}
    if not present:
        if require:
            raise RuntimeError(
                f"{result_path} is missing current nonlinear-solver provenance"
            )
        return None, None, None
    if present != fields:
        missing = sorted(fields - present)
        raise RuntimeError(
            f"{result_path} has incomplete nonlinear-solver provenance; "
            f"missing {missing}"
        )

    solver_name = str(_scalar(archive, "nonlinear_solver", result_path))
    petsc_solvers = {"petsc-snes", "petsc-snes-mpi"}
    if solver_name not in {"core-newton", *petsc_solvers}:
        raise RuntimeError(
            f"{result_path} has unsupported nonlinear solver {solver_name!r}"
        )
    configuration = _parse_json_field(
        archive, "solver_configuration_json", result_path, dict
    )
    diagnostics = _parse_json_field(
        archive, "nonlinear_step_diagnostics_json", result_path, list
    )
    if configuration.get("name") != solver_name:
        raise RuntimeError(
            f"{result_path} has solver name/configuration disagreement"
        )
    if len(diagnostics) != len(times) - 1:
        raise RuntimeError(
            f"{result_path} has {len(diagnostics)} nonlinear diagnostics for "
            f"{len(times) - 1} completed steps"
        )

    domain_rejection_api = configuration.get("function_domain_rejection_api")
    if domain_rejection_api is not None and (
        solver_name not in petsc_solvers
        or domain_rejection_api != PETSC_FUNCTION_DOMAIN_REJECTION_API
    ):
        raise RuntimeError(
            f"{result_path} has unsupported function-domain rejection metadata"
        )
    records_domain_rejections = domain_rejection_api is not None

    dt = float(times[1] - times[0])
    for index, record in enumerate(diagnostics, start=1):
        if not isinstance(record, dict):
            raise RuntimeError(
                f"{result_path} has malformed nonlinear diagnostic at step {index}"
            )
        _validate_finite_json(
            record, f"nonlinear diagnostic at step {index}", result_path
        )
        required = {"time", "dt", "nonlinear_iterations"}
        if not required.issubset(record):
            raise RuntimeError(
                f"{result_path} has incomplete nonlinear diagnostic at step {index}"
            )
        diagnostic_time = _finite_number(
            record["time"], f"diagnostic time at step {index}", result_path
        )
        diagnostic_dt = _finite_number(
            record["dt"], f"diagnostic dt at step {index}", result_path, minimum=0.0,
            strict=True,
        )
        if not np.isclose(
            diagnostic_time, times[index], rtol=1.0e-12, atol=1.0e-14
        ) or not np.isclose(
            diagnostic_dt, dt, rtol=1.0e-12, atol=1.0e-14
        ):
            raise RuntimeError(
                f"{result_path} has time-mismatched nonlinear diagnostic at step {index}"
            )
        _integer(
            record["nonlinear_iterations"],
            f"nonlinear iteration count at step {index}",
            result_path,
        )

        if solver_name not in petsc_solvers:
            continue
        petsc_fields = {
            "initial_residual_norm",
            "final_residual_norm",
            "residual_acceptance_threshold",
            "petsc_function_norm",
            "snes_converged_reason",
            "ksp_converged_reason",
            "linear_iterations",
            "residual_history",
            "assembly_seconds",
            "solve_seconds",
        }
        if not petsc_fields.issubset(record):
            raise RuntimeError(
                f"{result_path} has incomplete PETSc diagnostic at step {index}"
            )
        domain_fields = {
            "function_domain_rejections",
            "last_function_domain_error",
        }
        present_domain_fields = domain_fields.intersection(record)
        expected_domain_fields = domain_fields if records_domain_rejections else set()
        if present_domain_fields != expected_domain_fields:
            raise RuntimeError(
                f"{result_path} has inconsistent function-domain rejection "
                f"diagnostic at step {index}"
            )
        if records_domain_rejections:
            rejections = _integer(
                record["function_domain_rejections"],
                f"function-domain rejection count at step {index}",
                result_path,
            )
            last_error = record["last_function_domain_error"]
            if (
                (rejections == 0 and last_error is not None)
                or (
                    rejections > 0
                    and (not isinstance(last_error, str) or not last_error.strip())
                )
            ):
                raise RuntimeError(
                    f"{result_path} has inconsistent function-domain rejection "
                    f"detail at step {index}"
                )
        snes_reason = _integer(
            record["snes_converged_reason"],
            f"SNES reason at step {index}",
            result_path,
            minimum=1,
        )
        if snes_reason <= 0:  # explicit for a useful failure message
            raise RuntimeError(
                f"{result_path} has nonpositive SNES reason at step {index}"
            )
        ksp_reason = _integer(
            record["ksp_converged_reason"],
            f"KSP reason at step {index}",
            result_path,
        )
        linear_iterations = _integer(
            record["linear_iterations"],
            f"linear iteration count at step {index}",
            result_path,
        )
        if ksp_reason < 0 or (linear_iterations > 0 and ksp_reason == 0):
            raise RuntimeError(
                f"{result_path} has invalid KSP reason at step {index}"
            )
        initial = _finite_number(
            record["initial_residual_norm"],
            f"PETSc initial residual at step {index}", result_path, minimum=0.0,
        )
        final = _finite_number(
            record["final_residual_norm"],
            f"PETSc final residual at step {index}", result_path, minimum=0.0,
        )
        threshold = _finite_number(
            record["residual_acceptance_threshold"],
            f"PETSc residual threshold at step {index}", result_path, minimum=0.0,
            strict=True,
        )
        _finite_number(
            record["petsc_function_norm"],
            f"PETSc function norm at step {index}", result_path, minimum=0.0,
        )
        if final > threshold:
            raise RuntimeError(
                f"{result_path} reports a PETSc final residual above its "
                f"acceptance threshold at step {index}"
            )
        history = record["residual_history"]
        if not isinstance(history, list):
            raise RuntimeError(
                f"{result_path} has malformed PETSc residual history at step {index}"
            )
        for value in history:
            _finite_number(
                value, f"PETSc residual history at step {index}", result_path,
                minimum=0.0,
            )
        for key in ("assembly_seconds", "solve_seconds"):
            _finite_number(
                record[key], f"PETSc {key} at step {index}", result_path,
                minimum=0.0,
            )

    if solver_name == "petsc-snes":
        expected_strings = {
            "snes_type": "newtonls",
            "line_search_type": "bt",
            "ksp_type": "preonly",
            "pc_type": "lu",
        }
        if any(configuration.get(key) != value for key, value in expected_strings.items()):
            raise RuntimeError(f"{result_path} has unexpected PETSc solver settings")
        for key in ("rtol", "atol", "stol"):
            if key not in configuration:
                raise RuntimeError(f"{result_path} has invalid PETSc {key} setting")
            _finite_number(
                configuration[key], f"PETSc {key} setting", result_path,
                minimum=0.0, strict=True,
            )
        _integer(configuration.get("max_it"), "PETSc max_it setting", result_path, minimum=1)

    return solver_name, configuration, diagnostics


def _mpi_configuration_value_matches(actual, expected, key, result_path):
    if isinstance(expected, bool):
        return isinstance(actual, bool) and actual is expected
    if isinstance(expected, float):
        return (
            isinstance(actual, Real)
            and not isinstance(actual, bool)
            and np.isfinite(float(actual))
            and float(actual) == expected
        )
    if isinstance(expected, int):
        try:
            return (
                _integer(actual, f"MPI configuration {key}", result_path)
                == expected
            )
        except RuntimeError:
            return False
    return actual == expected


def _validate_mpi_profile_petsc_options(configuration, profile, result_path):
    if profile == MPI_DIRECT_SUPERLU_DIST_PROFILE:
        if "petsc_options" in configuration:
            raise RuntimeError(
                f"{result_path} has unexpected PETSc options for direct MPI "
                "linear-solver profile"
            )
        if "ksp_error_if_not_converged" in configuration:
            raise RuntimeError(
                f"{result_path} has unexpected iterative KSP configuration "
                "for direct MPI linear-solver profile"
            )
        return

    expected_options = MPI_COMPANION_LINEAR_SOLVER_PETSC_OPTIONS[profile]
    options = configuration.get("petsc_options")
    if not isinstance(options, dict) or len(options) != len(expected_options):
        raise RuntimeError(
            f"{result_path} has unexpected PETSc options for MPI linear-solver "
            f"profile {profile!r}"
        )

    prefixes = set()
    matched_suffixes = set()
    for full_name, value in options.items():
        if not isinstance(full_name, str) or not isinstance(value, str):
            raise RuntimeError(
                f"{result_path} has malformed PETSc options for MPI "
                f"linear-solver profile {profile!r}"
            )
        suffixes = [
            suffix for suffix in expected_options if full_name.endswith(suffix)
        ]
        if len(suffixes) != 1:
            raise RuntimeError(
                f"{result_path} has unexpected PETSc option {full_name!r} for "
                f"MPI linear-solver profile {profile!r}"
            )
        suffix = suffixes[0]
        prefix = full_name[: -len(suffix)]
        if (
            re.fullmatch(r"coupfe_cardiac_[0-9a-f]+_", prefix) is None
            or suffix in matched_suffixes
            or value != expected_options[suffix]
        ):
            raise RuntimeError(
                f"{result_path} has unexpected PETSc option {full_name!r} for "
                f"MPI linear-solver profile {profile!r}"
            )
        prefixes.add(prefix)
        matched_suffixes.add(suffix)
    if len(prefixes) != 1 or matched_suffixes != set(expected_options):
        raise RuntimeError(
            f"{result_path} has inconsistent namespaced PETSc options for MPI "
            f"linear-solver profile {profile!r}"
        )


def _validate_mpi_companion_metadata(
    archive,
    result_path,
    *,
    solver_name,
    configuration,
    diagnostics,
    case,
    n_elem,
    ndof,
    integrator,
    dt,
    formulation,
    topology,
    apex_offset,
    mass_representation,
    fiber_sampling,
    fiber_sampling_option,
    tbar_identity,
    isotropic,
    material_eta_pa_s,
    viscous_term_active,
    parameter_variant,
    benchmark_identity,
):
    """Fail closed on either declared distributed-companion contract.

    ``petsc-snes-mpi`` is not a generic label.  Its implementation identity
    must select either the historical open-polar Q1/P0 contract or the closed
    pointwise-kappa/consistent-mass contract. Duplicate archive and JSON fields
    are intentional and are cross-checked here.
    """
    present_archive_fields = (
        MPI_COMPANION_ARCHIVE_FIELDS
        | MPI_COMPANION_OPTIONAL_ARCHIVE_FIELDS
        | MPI_MASS_PARTITION_ARCHIVE_FIELDS
    ).intersection(archive.files)
    present_base_fields = MPI_COMPANION_ARCHIVE_FIELDS.intersection(
        archive.files
    )
    if solver_name != "petsc-snes-mpi":
        if present_archive_fields:
            raise RuntimeError(
                f"{result_path} has MPI companion fields with nonlinear solver "
                f"{solver_name!r}"
            )
        return None

    missing = MPI_COMPANION_ARCHIVE_FIELDS - present_base_fields
    if missing:
        raise RuntimeError(
            f"{result_path} has incomplete MPI companion provenance; "
            f"missing {sorted(missing)}"
        )

    if _text_scalar(archive, "driver", result_path) != (
        "examples/cardiac_benchmark/run_mpi.py"
    ):
        raise RuntimeError(f"{result_path} has unexpected MPI companion driver")
    if not _boolean_scalar(archive, "mpi_enabled", result_path):
        raise RuntimeError(f"{result_path} does not mark MPI execution as enabled")

    archive_ranks = _integer(
        _scalar(archive, "mpi_ranks", result_path),
        "MPI archive rank count",
        result_path,
        minimum=1,
    )
    world_size = _integer(
        _scalar(archive, "mpi_world_size", result_path),
        "MPI world size",
        result_path,
        minimum=1,
    )
    configuration_ranks = _integer(
        configuration.get("ranks"),
        "MPI configuration rank count",
        result_path,
        minimum=1,
    )
    if archive_ranks != world_size or world_size != configuration_ranks:
        raise RuntimeError(
            f"{result_path} has MPI archive/configuration rank disagreement"
        )

    counts = np.asarray(archive["mpi_local_element_counts"])
    if (
        counts.shape != (world_size,)
        or not np.issubdtype(counts.dtype, np.integer)
        or np.any(counts < 0)
    ):
        raise RuntimeError(
            f"{result_path} has invalid MPI local-element counts"
        )
    quotient, remainder = divmod(n_elem, world_size)
    expected_counts = np.full(world_size, quotient, dtype=np.int64)
    expected_counts[:remainder] += 1
    if not np.array_equal(counts.astype(np.int64, copy=False), expected_counts):
        raise RuntimeError(
            f"{result_path} has MPI partition/count disagreement"
        )

    archive_implementation = _text_scalar(
        archive, "mpi_implementation", result_path
    )
    if archive_implementation not in MPI_COMPANION_IMPLEMENTATIONS:
        raise RuntimeError(
            f"{result_path} has unsupported MPI companion implementation"
        )
    expected_configuration = dict(MPI_COMPANION_FIXED_CONFIGURATION)
    expected_configuration["implementation"] = archive_implementation
    for key, expected in expected_configuration.items():
        actual = configuration.get(key)
        if key not in configuration or not _mpi_configuration_value_matches(
            actual, expected, key, result_path
        ):
            raise RuntimeError(
                f"{result_path} has unexpected MPI companion configuration "
                f"for {key!r}"
            )

    archive_has_profile = "mpi_linear_solver_profile" in archive.files
    configuration_has_profile = "linear_solver_profile" in configuration
    if archive_has_profile != configuration_has_profile:
        raise RuntimeError(
            f"{result_path} has MPI linear-solver profile "
            "archive/config disagreement"
        )
    unexpected_legacy_fields = set()
    if archive_has_profile:
        linear_solver_profile = _text_scalar(
            archive, "mpi_linear_solver_profile", result_path
        )
        if (
            linear_solver_profile
            not in MPI_COMPANION_LINEAR_SOLVER_CONFIGURATIONS
        ):
            raise RuntimeError(
                f"{result_path} has unsupported MPI linear-solver profile "
                f"{linear_solver_profile!r}"
            )
        if configuration["linear_solver_profile"] != linear_solver_profile:
            raise RuntimeError(
                f"{result_path} has MPI linear-solver profile "
                "archive/config disagreement"
            )
        profile_configuration = MPI_COMPANION_LINEAR_SOLVER_CONFIGURATIONS[
            linear_solver_profile
        ]
    else:
        # Reviewed direct-SuperLU archives created before explicit profiles
        # remain admissible only through the complete original direct contract.
        linear_solver_profile = MPI_DIRECT_SUPERLU_DIST_PROFILE
        profile_configuration = MPI_COMPANION_LEGACY_DIRECT_CONFIGURATION
        modern_direct_fields = (
            set(
                MPI_COMPANION_LINEAR_SOLVER_CONFIGURATIONS[
                    MPI_DIRECT_SUPERLU_DIST_PROFILE
                ]
            )
            - set(MPI_COMPANION_LEGACY_DIRECT_CONFIGURATION)
        )
        unexpected_legacy_fields = modern_direct_fields.intersection(configuration)
    for key, expected in profile_configuration.items():
        actual = configuration.get(key)
        if key not in configuration or not _mpi_configuration_value_matches(
            actual, expected, key, result_path
        ):
            raise RuntimeError(
                f"{result_path} has unexpected MPI linear-solver profile "
                f"configuration for {key!r}"
            )
    if unexpected_legacy_fields:
        raise RuntimeError(
            f"{result_path} has explicit-profile MPI configuration fields "
            "without duplicated linear-solver profile provenance: "
            f"{sorted(unexpected_legacy_fields)}"
        )
    _validate_mpi_profile_petsc_options(
        configuration, linear_solver_profile, result_path
    )

    line_search_api = configuration.get("line_search_configuration_api")
    if line_search_api not in {"SNES.getLineSearch", "namespaced PETSc option"}:
        raise RuntimeError(
            f"{result_path} has unsupported MPI line-search configuration API"
        )
    for field in ("petsc4py_version", "petsc_version"):
        value = configuration.get(field)
        if not isinstance(value, str) or not value.strip():
            raise RuntimeError(
                f"{result_path} has unavailable MPI runtime field {field!r}"
            )

    if archive_implementation != configuration["implementation"]:
        raise RuntimeError(
            f"{result_path} has MPI implementation archive/config disagreement"
        )
    if _text_scalar(archive, "mpi_partition", result_path) != (
        "coupfe.partition_elements"
    ):
        raise RuntimeError(f"{result_path} has unsupported MPI partition policy")
    if _text_scalar(archive, "mpi_build_layout", result_path) != (
        "isolated-rank-directories"
    ):
        raise RuntimeError(f"{result_path} has unsupported MPI build layout")
    archive_factor = _text_scalar(
        archive, "mpi_factor_solver_type", result_path
    )
    if archive_factor != configuration["factor_solver_type"]:
        raise RuntimeError(
            f"{result_path} has MPI factor-solver archive/config disagreement"
        )

    mass_partition_metadata = None
    present_mass_fields = MPI_MASS_PARTITION_ARCHIVE_FIELDS.intersection(
        archive.files
    )
    if (
        archive_implementation in MPI_CLOSED_IMPLEMENTATIONS
        or present_mass_fields
    ):
        missing_mass = MPI_MASS_PARTITION_ARCHIVE_FIELDS - present_mass_fields
        if missing_mass:
            raise RuntimeError(
                f"{result_path} has incomplete MPI mass provenance; "
                f"missing {sorted(missing_mass)}"
            )
        mass_partition = _text_scalar(
            archive, "mpi_mass_partition", result_path
        )
        row_ranges = np.asarray(archive["mpi_mass_owned_row_ranges"])
        local_nnz = np.asarray(archive["mpi_mass_local_nnz"])
        touching = np.asarray(archive["mpi_mass_touching_element_counts"])
        if (
            row_ranges.shape != (world_size, 2)
            or local_nnz.shape != (world_size,)
            or touching.shape != (world_size,)
            or not np.issubdtype(row_ranges.dtype, np.integer)
            or not np.issubdtype(local_nnz.dtype, np.integer)
            or not np.issubdtype(touching.dtype, np.integer)
            or np.any(local_nnz < 0)
            or np.any(touching < 0)
        ):
            raise RuntimeError(f"{result_path} has invalid MPI mass arrays")
        row_ranges = row_ranges.astype(np.int64, copy=False)
        expected_start = 0
        for start, stop in row_ranges:
            if int(start) != expected_start or int(stop) < int(start):
                raise RuntimeError(
                    f"{result_path} has noncontiguous MPI mass row ownership"
                )
            expected_start = int(stop)
        if expected_start != ndof:
            raise RuntimeError(
                f"{result_path} MPI mass rows do not cover the global system"
            )

        expected_mass_partition = (
            "owned-row-csr-all-touching-elements"
            if archive_implementation in MPI_CLOSED_IMPLEMENTATIONS
            else "replicated-diagonal-owned-entry-insertion"
        )
        if mass_partition != expected_mass_partition:
            raise RuntimeError(
                f"{result_path} has unexpected MPI mass partition policy"
            )
        configured_range = configuration.get("mass_owned_row_range")
        if (
            configuration.get("mass_representation") != mass_representation
            or configuration.get("mass_partition") != mass_partition
            or not isinstance(configured_range, list)
            or configured_range != row_ranges[0].astype(int).tolist()
            or _integer(
                configuration.get("mass_local_nnz"),
                "MPI configured local mass nnz",
                result_path,
                minimum=0,
            )
            != int(local_nnz[0])
        ):
            raise RuntimeError(
                f"{result_path} has MPI mass archive/configuration disagreement"
            )
        if (
            archive_implementation in MPI_CLOSED_IMPLEMENTATIONS
            and (np.any(local_nnz <= 0) or np.any(touching <= 0))
        ):
            raise RuntimeError(
                f"{result_path} has empty closed consistent-mass rank data"
            )
        if (
            archive_implementation == MPI_COMPANION_IMPLEMENTATION
            and np.any(touching != 0)
        ):
            raise RuntimeError(
                f"{result_path} labels touching elements for diagonal mass"
            )
        mass_partition_metadata = {
            "partition": mass_partition,
            "owned_row_ranges": row_ranges.astype(int).tolist(),
            "local_nnz": local_nnz.astype(int).tolist(),
            "touching_element_counts": touching.astype(int).tolist(),
        }

    evaluation_mode = _text_scalar(
        archive, "element_evaluation_mode", result_path
    )
    if evaluation_mode not in {"joint", "split"} or (
        configuration.get("element_evaluation_mode") != evaluation_mode
    ):
        raise RuntimeError(
            f"{result_path} has element-evaluation archive/config disagreement"
        )
    residual_only_available = _boolean_scalar(
        archive, "compiled_material_residual_only_available", result_path
    )
    configured_residual_only = configuration.get(
        "compiled_material_residual_only_available"
    )
    if (
        not isinstance(configured_residual_only, bool)
        or configured_residual_only != residual_only_available
    ):
        raise RuntimeError(
            f"{result_path} has compiled-material archive/config disagreement"
        )

    for index, record in enumerate(diagnostics, start=1):
        if "ranks" not in record:
            raise RuntimeError(
                f"{result_path} is missing MPI rank metadata at step {index}"
            )
        diagnostic_ranks = _integer(
            record["ranks"],
            f"MPI diagnostic rank count at step {index}",
            result_path,
            minimum=1,
        )
        if diagnostic_ranks != world_size:
            raise RuntimeError(
                f"{result_path} has MPI diagnostic/configuration rank "
                f"disagreement at step {index}"
            )

    expected_configuration_pressure_law = {
        MPI_COMPANION_IMPLEMENTATION: "log",
        MPI_CLOSED_STD_KAPPA_IMPLEMENTATION: "not-applicable",
        MPI_CLOSED_LOCAL_PRESSURE_IMPLEMENTATION: "log",
        MPI_CLOSED_LOCAL_PRESSURE_PAPER_IMPLEMENTATION: "paper",
        MPI_CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION: "not-applicable",
        MPI_CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION: "log",
        MPI_CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION: "paper",
    }[archive_implementation]
    configured_pressure_law = configuration.get("local_pressure_law")
    pressure_law_semantics = (
        configured_pressure_law == expected_configuration_pressure_law
        or (
            configured_pressure_law is None
            and expected_configuration_pressure_law != "paper"
        )
    )
    historical_semantics = (
        archive_implementation == MPI_COMPANION_IMPLEMENTATION
        and integrator == "be"
        and formulation == "hex8_local_pressure_p0_condensed_logj"
        and topology == "polar_ring"
        and apex_offset > 0.0
        and mass_representation == "lumped_row_sum"
        and fiber_sampling == "cg1_gram_schmidt"
        and fiber_sampling_option == "cg1"
        and tbar_identity
        == {
            "definition": "analytic_parametric",
            "source_filename": "",
            "source_sha256": "",
            "metadata_filename": "",
            "metadata_sha256": "",
            "metadata_schema": "",
        }
        and isotropic is False
        and material_eta_pa_s == 100.0
        and viscous_term_active is True
        and parameter_variant == "benchmark_eta"
        and pressure_law_semantics
    )
    closed_formulation = {
        MPI_CLOSED_STD_KAPPA_IMPLEMENTATION: "hex8_standard_pointwise_kappa",
        MPI_CLOSED_LOCAL_PRESSURE_IMPLEMENTATION: (
            "hex8_local_pressure_p0_condensed_logj"
        ),
        MPI_CLOSED_LOCAL_PRESSURE_PAPER_IMPLEMENTATION: (
            "hex8_local_pressure_p0_condensed_mean_logj_paper_j2"
        ),
        MPI_CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION: (
            "hex8_standard_pointwise_kappa"
        ),
        MPI_CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION: (
            "hex8_local_pressure_p0_condensed_logj"
        ),
        MPI_CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION: (
            "hex8_local_pressure_p0_condensed_mean_logj_paper_j2"
        ),
    }.get(archive_implementation)
    generalized_alpha_semantics = (
        archive_implementation in MPI_CLOSED_GENERALIZED_ALPHA_IMPLEMENTATIONS
    )
    expected_integrator = (
        "generalized-alpha" if generalized_alpha_semantics else "be"
    )
    reviewed_generalized_alpha_case = (
        (
            case == "A"
            and (
                benchmark_identity is None
                or (
                    benchmark_identity["step"] == 0
                    and benchmark_identity["case"] == "A"
                )
            )
        )
        or (
            case == "B"
            and benchmark_identity is not None
            and benchmark_identity["step"] in {0, 2}
            and benchmark_identity["case"] == "B"
        )
    )
    closed_semantics = (
        archive_implementation in MPI_CLOSED_IMPLEMENTATIONS
        and (
            reviewed_generalized_alpha_case
            if generalized_alpha_semantics
            else case in {"A", "B"}
        )
        and integrator == expected_integrator
        and formulation == closed_formulation
        and topology == "closed_multiblock_disk"
        and apex_offset == 0.0
        and mass_representation == "consistent_q1_hex8"
        and fiber_sampling == "gp_direct_rule"
        and fiber_sampling_option == "gp-direct"
        and isinstance(tbar_identity, dict)
        and tbar_identity.get("definition") == "laplace_presolved"
        and isotropic is False
        and material_eta_pa_s == 100.0
        and viscous_term_active is True
        and parameter_variant == "benchmark_eta"
        and pressure_law_semantics
    )
    if generalized_alpha_semantics:
        alpha_metadata = configuration.get("generalized_alpha")
        expected_alpha = {
            "alpha_m": 0.2,
            "alpha_f": 0.4,
            "gamma": 0.7,
            "beta": 0.36,
            "parameter_source": (
                "finsberg/cardiac_benchmark problem.py defaults"
            ),
            "acceleration_stage": (
                "alpha_m*a_n + (1-alpha_m)*a_np1"
            ),
            "force_stage": "alpha_f*x_n + (1-alpha_f)*x_np1",
            "load_time": "t_np1 - alpha_f*dt",
        }
        generalized_alpha_semantics = (
            isinstance(alpha_metadata, dict)
            and alpha_metadata == expected_alpha
            and configuration.get("time_integrator") == "generalized-alpha"
            and configuration.get("material_batch_time_integrator")
            == "generalized-alpha-source-matched"
            and configuration.get("compiled_material_dt") == dt
            and configuration.get("acceleration_stage") == "1-alpha_m"
            and configuration.get("force_stage") == "1-alpha_f"
            and configuration.get("material_viscous_rate")
            == "sym(F_stage^T*grad(v_stage))"
            and configuration.get("nonlinear_initial_guess")
            == "accepted-u_n-like-simula"
            and (
                archive_implementation
                not in {
                    MPI_CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION,
                    MPI_CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION,
                }
                or configuration.get("local_pressure_stage_alpha_f") == 0.4
            )
        )
        closed_semantics = closed_semantics and generalized_alpha_semantics
    if (
        archive_implementation == MPI_COMPANION_IMPLEMENTATION
        and not historical_semantics
    ):
        raise RuntimeError(
            f"{result_path} does not match the historical MPI companion "
            "discretization contract"
        )
    if (
        archive_implementation in MPI_CLOSED_IMPLEMENTATIONS
        and not closed_semantics
    ):
        raise RuntimeError(
            f"{result_path} does not match the closed Case {case} MPI companion "
            "discretization contract"
        )

    metadata = {
        "implementation": archive_implementation,
        "world_size": world_size,
        "local_element_counts": counts.astype(int).tolist(),
        "partition": "coupfe.partition_elements",
        "build_layout": "isolated-rank-directories",
        "factor_solver_type": archive_factor,
        "linear_solver_profile": linear_solver_profile,
        "element_evaluation_mode": evaluation_mode,
        "compiled_material_residual_only_available": residual_only_available,
    }
    if closed_semantics:
        formulation_contract = (
            "std_kappa"
            if archive_implementation
            in {
                MPI_CLOSED_STD_KAPPA_IMPLEMENTATION,
                MPI_CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION,
            }
            else (
                "local_pressure_mean_logj_paper_j2"
                if archive_implementation
                in {
                    MPI_CLOSED_LOCAL_PRESSURE_PAPER_IMPLEMENTATION,
                    MPI_CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION,
                }
                else "local_pressure"
            )
        )
        time_contract = (
            "generalized_alpha"
            if archive_implementation
            in MPI_CLOSED_GENERALIZED_ALPHA_IMPLEMENTATIONS
            else "be"
        )
        metadata.update(
            {
                "contract": (
                    f"closed_case_{case.lower()}_{formulation_contract}_"
                    + (
                        "consistent_generalized_alpha"
                        if time_contract == "generalized_alpha"
                        else "consistent"
                    )
                ),
                "mass_partition": mass_partition_metadata,
            }
        )
    return metadata


def _validate_sampling_metadata(archive, result_path, point_sampling, n_elem):
    if point_sampling != CURRENT_POINT_SAMPLING:
        return None
    metadata = {}
    for point in ("p0", "p1"):
        element = _integer(
            _scalar(archive, f"{point}_sampling_element", result_path),
            f"{point} sampling element",
            result_path,
        )
        if element >= n_elem:
            raise RuntimeError(f"{result_path} has out-of-range {point} sampling element")
        natural = _finite_array(
            archive, f"{point}_sampling_natural", result_path, shape=(3,)
        )
        weights = _finite_array(
            archive, f"{point}_sampling_weights", result_path, shape=(8,)
        )
        error = float(
            _scalar(archive, f"{point}_sampling_reconstruction_error_m", result_path)
        )
        if (
            np.any(np.abs(natural) > 1.0 + 1.0e-8)
            or np.any(weights < -1.0e-10)
            or not np.isclose(weights.sum(), 1.0, rtol=0.0, atol=1.0e-10)
            or not np.isfinite(error)
            or error < 0.0
        ):
            raise RuntimeError(f"{result_path} has invalid {point} Hex8 sampling metadata")
        metadata[point] = {
            "element": element,
            "natural_coordinates": natural.tolist(),
            "weights": weights.tolist(),
            "reconstruction_error_m": error,
        }
    return metadata


def _validate_passed_audit_record(record, name, schema, result_path):
    if not isinstance(record, dict):
        raise RuntimeError(
            f"{result_path} has malformed pre-solve {name} audit"
        )
    if record.get("schema") != schema:
        raise RuntimeError(
            f"{result_path} has unsupported pre-solve {name} audit schema"
        )
    if record.get("passed") is not True or record.get("failures") != []:
        raise RuntimeError(
            f"{result_path} records a failed pre-solve {name} audit"
        )


def _validate_closed_pre_solve_audit(
    archive, result_path, *, case, topology, n_node, n_elem
):
    """Validate the retained setup gate for a closed benchmark mesh."""
    audit = _parse_json_field(
        archive, "pre_solve_audit_json", result_path, dict
    )
    required = {"geometry", "robin"}
    if case == "B":
        required.add("pressure")
    if not required.issubset(audit):
        raise RuntimeError(
            f"{result_path} has incomplete closed-mesh pre-solve audit"
        )

    geometry = audit["geometry"]
    _validate_passed_audit_record(
        geometry,
        "geometry",
        "coupfe-cardiac-pre-solve-geometry-v1",
        result_path,
    )
    if (
        geometry.get("mesh_topology") != topology
        or geometry.get("require_closed") is not True
    ):
        raise RuntimeError(
            f"{result_path} has inconsistent closed-mesh geometry audit"
        )
    if _integer(
        geometry.get("nodes"), "audited node count", result_path
    ) != n_node or _integer(
        geometry.get("elements"), "audited element count", result_path
    ) != n_elem:
        raise RuntimeError(
            f"{result_path} has geometry audit/mesh size disagreement"
        )
    zero_geometry_counts = (
        "unclassified_exterior_faces",
        "nonexterior_labeled_faces",
        "multiply_labeled_faces",
        "nonmanifold_faces",
        "nonpositive_gauss_jacobians",
        "nonpositive_extended_jacobians",
    )
    if any(
        _integer(geometry.get(field), f"geometry audit {field}", result_path)
        != 0
        for field in zero_geometry_counts
    ):
        raise RuntimeError(
            f"{result_path} has an uncleared closed-mesh geometry audit"
        )
    exterior = _integer(
        geometry.get("exterior_faces"),
        "geometry audit exterior face count",
        result_path,
        minimum=1,
    )
    labeled = _integer(
        geometry.get("labeled_exterior_faces"),
        "geometry audit labeled exterior face count",
        result_path,
        minimum=1,
    )
    if exterior != labeled:
        raise RuntimeError(
            f"{result_path} has unclassified closed-mesh exterior faces"
        )
    for field in (
        "gauss_jacobian_min_m3",
        "extended_jacobian_min_m3",
        "extended_scaled_jacobian_min",
    ):
        _finite_number(
            geometry.get(field),
            f"geometry audit {field}",
            result_path,
            minimum=0.0,
            strict=True,
        )

    if "pressure" in audit:
        pressure = audit["pressure"]
        _validate_passed_audit_record(
            pressure,
            "pressure",
            "coupfe-cardiac-pre-solve-pressure-v1",
            result_path,
        )
        projected_area = _finite_number(
            pressure.get("analytic_projected_base_area_m2"),
            "pressure audit analytic projected base area",
            result_path,
            minimum=0.0,
            strict=True,
        )
        resultant = _finite_json_vector(
            pressure.get("unit_pressure_resultant_N"),
            "pressure audit unit-pressure resultant",
            result_path,
        )
        expected_resultant = _finite_json_vector(
            pressure.get("expected_unit_pressure_resultant_N"),
            "pressure audit expected unit-pressure resultant",
            result_path,
        )
        moment = _finite_json_vector(
            pressure.get("unit_pressure_moment_Nm"),
            "pressure audit unit-pressure moment",
            result_path,
        )
        analytic_resultant = np.asarray([projected_area, 0.0, 0.0])
        if not np.allclose(
            expected_resultant,
            analytic_resultant,
            rtol=1.0e-12,
            atol=1.0e-15,
        ):
            raise RuntimeError(
                f"{result_path} has inconsistent pre-solve pressure audit "
                "expected resultant"
            )

        computed_metrics = {
            "relative_magnitude_error": abs(
                np.linalg.norm(resultant) / projected_area - 1.0
            ),
            "signed_axial_ratio": resultant[0] / projected_area,
            "relative_signed_axial_error": abs(
                resultant[0] / projected_area - 1.0
            ),
            "relative_resultant_error": (
                np.linalg.norm(resultant - analytic_resultant) / projected_area
            ),
            "transverse_fraction": (
                np.linalg.norm(resultant[1:]) / projected_area
            ),
            "normalized_moment": (
                np.linalg.norm(moment) / (projected_area * 0.1)
            ),
        }
        reported_metrics = {}
        for field in computed_metrics:
            reported_metrics[field] = _finite_number(
                pressure.get(field),
                f"pressure audit {field.replace('_', ' ')}",
                result_path,
                minimum=None if field == "signed_axial_ratio" else 0.0,
            )
            if not np.isclose(
                reported_metrics[field],
                computed_metrics[field],
                rtol=1.0e-12,
                atol=1.0e-15,
            ):
                raise RuntimeError(
                    f"{result_path} has inconsistent pre-solve pressure audit "
                    f"{field}"
                )

        pressure_limits = {
            "relative_magnitude_error": 5.0e-3,
            "relative_signed_axial_error": 5.0e-3,
            "relative_resultant_error": 5.0e-3,
            "transverse_fraction": 5.0e-10,
            "normalized_moment": 5.0e-10,
        }
        for field, limit in pressure_limits.items():
            if computed_metrics[field] > limit:
                raise RuntimeError(
                    f"{result_path} has an out-of-tolerance pre-solve pressure audit"
                )

    robin = audit["robin"]
    _validate_passed_audit_record(
        robin,
        "Robin",
        "coupfe-cardiac-pre-solve-robin-v1",
        result_path,
    )
    robin_limits = {
        "spring_symmetry_error": 1.0e-6,
        "dashpot_symmetry_error": 1.0e-9,
    }
    for field, limit in robin_limits.items():
        value = _finite_number(
            robin.get(field),
            f"Robin audit {field}",
            result_path,
            minimum=0.0,
        )
        if value > limit:
            raise RuntimeError(
                f"{result_path} has an out-of-tolerance pre-solve Robin audit"
            )
    _integer(
        robin.get("active_dofs"),
        "Robin audit active DOF count",
        result_path,
        minimum=1,
    )
    return audit


def _portable_tbar_identity(archive, result_path):
    """Validate tbar provenance and return a path-free identity record."""
    definition = _text_scalar(archive, "tbar_definition", result_path)
    source_fields = {"tbar_source_filename", "tbar_source_sha256"}
    present_source_fields = source_fields.intersection(archive.files)
    metadata_fields = {
        "tbar_metadata_filename",
        "tbar_metadata_sha256",
        "tbar_metadata_schema",
    }
    present_metadata_fields = metadata_fields.intersection(archive.files)

    # Development runs before portable provenance was added retained the input
    # path in tbar_definition. Accept such a file only while the referenced
    # field is available to hash, and never copy that path into a report.
    legacy_prefix = "laplace_presolved:"
    if definition.startswith(legacy_prefix):
        if present_source_fields or present_metadata_fields:
            raise RuntimeError(
                f"{result_path} mixes legacy and portable tbar provenance"
            )
        source_text = definition[len(legacy_prefix):]
        source_path = Path(source_text).expanduser()
        if not source_path.is_absolute():
            source_path = Path(result_path).parent / source_path
        if not source_path.is_file():
            raise RuntimeError(
                f"{result_path} cannot resolve its legacy Laplace tbar source"
            )
        source_sha256, _ = _sha256_file(source_path)
        return {
            "definition": "laplace_presolved",
            "source_filename": source_path.name,
            "source_sha256": source_sha256,
            "metadata_filename": "",
            "metadata_sha256": "",
            "metadata_schema": "",
        }

    if present_source_fields and present_source_fields != source_fields:
        raise RuntimeError(f"{result_path} has incomplete tbar provenance")
    if present_source_fields != source_fields:
        raise RuntimeError(f"{result_path} is missing portable tbar provenance")
    if present_metadata_fields and present_metadata_fields != metadata_fields:
        raise RuntimeError(
            f"{result_path} has incomplete tbar metadata provenance"
        )
    if not present_metadata_fields:
        raise RuntimeError(
            f"{result_path} is missing portable tbar metadata provenance"
        )

    source_filename = _text_scalar(
        archive, "tbar_source_filename", result_path
    )
    source_sha256 = _text_scalar(archive, "tbar_source_sha256", result_path)
    metadata_filename = _text_scalar(
        archive, "tbar_metadata_filename", result_path
    )
    metadata_sha256 = _text_scalar(
        archive, "tbar_metadata_sha256", result_path
    )
    metadata_schema = _text_scalar(
        archive, "tbar_metadata_schema", result_path
    )
    if definition == "analytic_parametric":
        if (
            source_filename
            or source_sha256
            or metadata_filename
            or metadata_sha256
            or metadata_schema
        ):
            raise RuntimeError(
                f"{result_path} has source fields for analytic tbar"
            )
    elif definition == "laplace_presolved":
        if (
            not source_filename
            or source_filename in {".", ".."}
            or Path(source_filename).is_absolute()
            or PureWindowsPath(source_filename).is_absolute()
            or Path(source_filename).name != source_filename
            or PureWindowsPath(source_filename).name != source_filename
        ):
            raise RuntimeError(
                f"{result_path} has a non-portable Laplace tbar filename"
            )
        if re.fullmatch(r"[0-9a-f]{64}", source_sha256) is None:
            raise RuntimeError(
                f"{result_path} has an invalid Laplace tbar SHA-256"
            )
        expected_metadata_filename = Path(source_filename).with_suffix(
            ".meta.json"
        ).name
        if (
            not metadata_filename
            or metadata_filename != expected_metadata_filename
            or Path(metadata_filename).is_absolute()
            or PureWindowsPath(metadata_filename).is_absolute()
            or Path(metadata_filename).name != metadata_filename
            or PureWindowsPath(metadata_filename).name != metadata_filename
        ):
            raise RuntimeError(
                f"{result_path} has a non-portable Laplace tbar metadata filename"
            )
        if re.fullmatch(r"[0-9a-f]{64}", metadata_sha256) is None:
            raise RuntimeError(
                f"{result_path} has an invalid Laplace tbar metadata SHA-256"
            )
        if metadata_schema != TBAR_METADATA_SCHEMA:
            raise RuntimeError(
                f"{result_path} has an unsupported Laplace tbar metadata schema"
            )
    else:
        raise RuntimeError(
            f"{result_path} has unsupported tbar definition {definition!r}"
        )
    return {
        "definition": definition,
        "source_filename": source_filename,
        "source_sha256": source_sha256,
        "metadata_filename": metadata_filename,
        "metadata_sha256": metadata_sha256,
        "metadata_schema": metadata_schema,
    }


def _current_method_metadata(
    archive, result_path, *, app_revision, fiber_sampling, topology
):
    """Validate current fields or narrowly backfill immutable old records."""
    core_fields = {
        "mass_representation",
        "fiber_sampling_option",
        "isotropic",
        "material_model_id",
        "material_eta_pa_s",
        "viscous_term_active",
        "parameter_variant",
    }
    tbar_fields = {
        "tbar_definition",
        "tbar_source_filename",
        "tbar_source_sha256",
        "tbar_metadata_filename",
        "tbar_metadata_sha256",
        "tbar_metadata_schema",
    }
    fields = core_fields | tbar_fields
    present = fields.intersection(archive.files)
    if not present:
        if (
            app_revision not in IMMUTABLE_PRE_METHOD_METADATA_APP_REFS
            or topology != "polar_ring"
            or fiber_sampling != "cg1_gram_schmidt"
        ):
            raise RuntimeError(
                f"{result_path} is missing current method metadata and does "
                "not match an immutable reviewed predecessor"
            )
        return {
            "mass_representation": "lumped_row_sum",
            "fiber_sampling_option": "cg1",
            "isotropic": False,
            "material_model_id": (
                COMPLETE_SWITCH_ENERGY_MATERIAL_MODEL_ID
                if app_revision
                == "6839c13b5bc80ec06c897684c51f503e80bd4b19"
                else LEGACY_SWITCH_STRESS_MATERIAL_MODEL_ID
            ),
            "material_eta_pa_s": 100.0,
            "viscous_term_active": True,
            "parameter_variant": "benchmark_eta",
            "tbar_identity": {
                "definition": "analytic_parametric",
                "source_filename": "",
                "source_sha256": "",
                "metadata_filename": "",
                "metadata_sha256": "",
                "metadata_schema": "",
            },
            "metadata_origin": (
                "reviewed-predecessor-source-checkpoint"
            ),
        }
    # A short-lived development format stored an absolute Laplace-field path
    # in tbar_definition before portable source filename/hash fields existed.
    # _portable_tbar_identity validates and sanitizes that format while the
    # referenced file is still locally available.  All other current records
    # require the six core fields and a tbar definition; the helper below gives
    # the more useful provenance-specific error for missing source pairs.
    required = core_fields | {"tbar_definition"}
    if not required.issubset(present):
        raise RuntimeError(
            f"{result_path} has incomplete current method metadata; missing "
            f"{sorted(required - present)}"
        )

    mass_representation = _text_scalar(
        archive, "mass_representation", result_path
    )
    if mass_representation not in MASS_REPRESENTATIONS:
        raise RuntimeError(
            f"{result_path} has unsupported mass representation "
            f"{mass_representation!r}"
        )
    fiber_sampling_option = _text_scalar(
        archive, "fiber_sampling_option", result_path
    )
    if FIBER_SAMPLING_OPTIONS.get(fiber_sampling_option) != fiber_sampling:
        raise RuntimeError(
            f"{result_path} has inconsistent fiber sampling option/method"
        )
    isotropic = _boolean_scalar(archive, "isotropic", result_path)
    material_model_id = _text_scalar(
        archive, "material_model_id", result_path
    )
    if material_model_id not in SUPPORTED_MATERIAL_MODEL_IDS:
        raise RuntimeError(
            f"{result_path} has unsupported material model identity "
            f"{material_model_id!r}"
        )
    material_eta_pa_s = _finite_number(
        _scalar(archive, "material_eta_pa_s", result_path),
        "material_eta_pa_s",
        result_path,
        minimum=0.0,
    )
    viscous_term_active = _boolean_scalar(
        archive, "viscous_term_active", result_path
    )
    if viscous_term_active != (material_eta_pa_s > 0.0):
        raise RuntimeError(
            f"{result_path} has inconsistent viscous activity metadata"
        )
    parameter_variant = _text_scalar(
        archive, "parameter_variant", result_path
    )
    expected_parameter_variant = (
        "benchmark_eta"
        if material_eta_pa_s == 100.0
        else "eta_zero_sensitivity"
        if material_eta_pa_s == 0.0
        else "eta_sensitivity"
    )
    if parameter_variant != expected_parameter_variant:
        raise RuntimeError(
            f"{result_path} has inconsistent viscosity parameter variant"
        )
    return {
        "mass_representation": mass_representation,
        "fiber_sampling_option": fiber_sampling_option,
        "isotropic": isotropic,
        "material_model_id": material_model_id,
        "material_eta_pa_s": material_eta_pa_s,
        "viscous_term_active": viscous_term_active,
        "parameter_variant": parameter_variant,
        "tbar_identity": _portable_tbar_identity(archive, result_path),
        "metadata_origin": "recorded",
    }


def _load_validated_result(archive, result_path, requested_case=None):
    """Validate an archive and copy all reporting data before it is closed."""
    schema = str(_scalar(archive, "result_schema", result_path))
    if schema != RESULT_SCHEMA:
        raise RuntimeError(
            f"{result_path} uses unsupported result schema {schema!r}; "
            f"expected {RESULT_SCHEMA!r}"
        )

    case = str(_scalar(archive, "case", result_path))
    if case not in CASE_NAMES:
        raise RuntimeError(f"{result_path} has unsupported benchmark case {case!r}")
    benchmark_identity = _validate_benchmark_configuration_metadata(
        archive, result_path, case
    )
    reference_case = (
        CASE_NAMES[case]
        if benchmark_identity is None
        else f"step_{benchmark_identity['step']}{case}"
    )
    if requested_case is not None and requested_case != reference_case:
        raise RuntimeError(
            f"{result_path} records Case {case}, which requires {reference_case}; "
            f"--case {requested_case} would compare the wrong benchmark"
        )

    if not _boolean_scalar(archive, "converged", result_path):
        raise RuntimeError(f"{result_path} is not marked as a completed solve")
    completed_steps = _integer(
        _scalar(archive, "completed_steps", result_path),
        "completed_steps",
        result_path,
    )
    expected_steps = _integer(
        _scalar(archive, "expected_steps", result_path),
        "expected_steps",
        result_path,
    )

    times = _finite_array(archive, "times", result_path)
    if times.ndim != 1 or len(times) < 2 or not np.all(np.diff(times) > 0.0):
        raise RuntimeError(f"{result_path} has an invalid time grid")
    if completed_steps != expected_steps or expected_steps != len(times) - 1:
        raise RuntimeError(
            f"{result_path} is incomplete: completed={completed_steps}, "
            f"expected={expected_steps}, time_intervals={len(times) - 1}"
        )

    histories = {}
    for field in ("u0", "u1"):
        histories[field] = _finite_array(
            archive, field, result_path, shape=(len(times), 3)
        ).copy()
    for field in ("tau", "pres"):
        if field in archive:
            histories[field] = _finite_array(
                archive, field, result_path, shape=(len(times),)
            ).copy()
    if benchmark_identity is not None:
        if not {"tau", "pres"}.issubset(histories):
            raise RuntimeError(
                f"{result_path} is missing benchmark load histories"
            )
        try:
            validate_load_histories(
                benchmark_identity["_configuration"],
                times,
                histories["tau"],
                histories["pres"],
            )
        except RuntimeError as error:
            raise RuntimeError(
                f"{result_path} has load histories inconsistent with its "
                "benchmark identity"
            ) from error

    integrator = str(_scalar(archive, "integrator", result_path))
    if integrator not in {"newmark", "be", "generalized-alpha"}:
        raise RuntimeError(f"{result_path} has unsupported integrator {integrator!r}")
    formulation = str(_scalar(archive, "formulation", result_path))
    if formulation not in FORMULATIONS:
        raise RuntimeError(f"{result_path} has unsupported formulation {formulation!r}")
    fiber_sampling = str(_scalar(archive, "fiber_sampling", result_path))
    if fiber_sampling not in FIBER_SAMPLING_METHODS:
        raise RuntimeError(
            f"{result_path} has unsupported fiber sampling {fiber_sampling!r}"
        )
    point_sampling = str(_scalar(archive, "point_sampling", result_path))
    if point_sampling not in POINT_SAMPLING_METHODS:
        raise RuntimeError(
            f"{result_path} has unsupported point sampling {point_sampling!r}"
        )
    viscous_rate = str(_scalar(archive, "viscous_rate", result_path))
    expected_viscous_rate = (
        "velocity_consistent_green_lagrange_at_alpha_f_stage"
        if integrator == "generalized-alpha"
        else "backward_difference"
    )
    if viscous_rate != expected_viscous_rate:
        raise RuntimeError(
            f"{result_path} has unsupported viscous-rate scheme {viscous_rate!r}"
        )

    dt = float(_scalar(archive, "dt", result_path))
    t_end = float(_scalar(archive, "t_end", result_path))
    for field, value in (("dt", dt), ("t_end", t_end)):
        if not np.isfinite(value) or value <= 0.0:
            raise RuntimeError(f"{result_path} has invalid positive field {field!r}")
    time_tolerance = max(1.0e-12, 1.0e-12 * abs(t_end))
    if (
        not np.isclose(times[0], 0.0, rtol=0.0, atol=time_tolerance)
        or not np.allclose(np.diff(times), dt, rtol=1.0e-12, atol=time_tolerance)
        or not np.isclose(times[-1], t_end, rtol=1.0e-12, atol=time_tolerance)
    ):
        raise RuntimeError(
            f"{result_path} time grid is inconsistent with dt={dt} and t_end={t_end}"
        )
    generalized_alpha_metadata = None
    if integrator == "generalized-alpha":
        expected_alpha = {
            "generalized_alpha_alpha_m": 0.2,
            "generalized_alpha_alpha_f": 0.4,
            "generalized_alpha_gamma": 0.7,
            "generalized_alpha_beta": 0.36,
        }
        observed_alpha = {}
        for field, expected in expected_alpha.items():
            if field not in archive:
                raise RuntimeError(
                    f"{result_path} is missing generalized-alpha field {field!r}"
                )
            value = float(_scalar(archive, field, result_path))
            if not np.isfinite(value) or value != expected:
                raise RuntimeError(
                    f"{result_path} has unexpected generalized-alpha field {field!r}"
                )
            observed_alpha[field.removeprefix("generalized_alpha_")] = value
        if _text_scalar(
            archive, "generalized_alpha_stage_contract", result_path
        ) != "simula-source-matched-v1":
            raise RuntimeError(
                f"{result_path} has unsupported generalized-alpha stage contract"
            )
        load_evaluation_times = _finite_array(
            archive,
            "load_evaluation_times_s",
            result_path,
            shape=(len(times),),
        )
        expected_load_times = times.copy()
        expected_load_times[0] = 0.0
        expected_load_times[1:] -= 0.4 * dt
        if not np.allclose(
            load_evaluation_times,
            expected_load_times,
            rtol=0.0,
            atol=max(1.0e-15, 1.0e-13 * dt),
        ):
            raise RuntimeError(
                f"{result_path} has inconsistent generalized-alpha load times"
            )
        generalized_alpha_metadata = {
            **observed_alpha,
            "stage_contract": "simula-source-matched-v1",
            "load_evaluation_times_s": load_evaluation_times.copy(),
        }
    elif "load_evaluation_times_s" in archive:
        load_evaluation_times = _finite_array(
            archive,
            "load_evaluation_times_s",
            result_path,
            shape=(len(times),),
        )
        if not np.array_equal(load_evaluation_times, times):
            raise RuntimeError(
                f"{result_path} has inconsistent endpoint load-evaluation times"
            )
    load_horizon_recorded = "load_horizon" in archive
    load_horizon = (
        float(_scalar(archive, "load_horizon", result_path))
        if load_horizon_recorded
        else t_end
    )
    if not np.isfinite(load_horizon):
        raise RuntimeError(
            f"{result_path} has an invalid load_horizon for dt={dt} and "
            f"t_end={t_end}"
        )
    load_intervals = int(round(load_horizon / dt))
    if (
        load_horizon < t_end
        or load_intervals < expected_steps
        or not np.isclose(
            load_intervals * dt,
            load_horizon,
            rtol=1.0e-12,
            atol=1.0e-15,
        )
    ):
        raise RuntimeError(
            f"{result_path} has an invalid load_horizon for dt={dt} and "
            f"t_end={t_end}"
        )
    load_horizon_origin = (
        "recorded" if load_horizon_recorded else "implicit_t_end_legacy"
    )
    if benchmark_identity is not None and not load_horizon_recorded:
        raise RuntimeError(
            f"{result_path} is missing the load horizon for its explicit "
            "benchmark identity"
        )
    if benchmark_identity is not None and integrator != "generalized-alpha":
        configuration = benchmark_identity["_configuration"]
        schedule_times = np.linspace(0.0, load_horizon, load_intervals + 1)
        expected_tau = np.zeros(len(times), dtype=float)
        expected_pressure = np.zeros(len(times), dtype=float)
        if configuration.active_stress_enabled:
            expected_tau = tau_of_t(
                schedule_times, p=configuration.activation_parameters
            )[: len(times)]
        if configuration.pressure_enabled:
            expected_pressure = p_of_t(
                schedule_times, p=configuration.pressure_parameters
            )[: len(times)]
        if not np.allclose(
            histories["tau"], expected_tau, rtol=2.0e-6, atol=1.0e-3
        ):
            raise RuntimeError(
                f"{result_path} has activation history inconsistent with its "
                "benchmark identity"
            )
        if not np.allclose(
            histories["pres"], expected_pressure, rtol=2.0e-6, atol=1.0e-3
        ):
            raise RuntimeError(
                f"{result_path} has pressure history inconsistent with its "
                "benchmark identity"
            )
    if integrator == "generalized-alpha":
        if load_horizon != 1.0:
            raise RuntimeError(
                f"{result_path} lacks the canonical generalized-alpha load horizon"
            )
        generalized_alpha_configuration = (
            benchmark_configuration(0, "A")
            if benchmark_identity is None and case == "A"
            else (
                benchmark_identity["_configuration"]
                if benchmark_identity is not None
                else None
            )
        )
        if generalized_alpha_configuration is None or (
            generalized_alpha_configuration.benchmark_step,
            generalized_alpha_configuration.case,
        ) not in {(0, "A"), (0, "B"), (2, "B")}:
            raise RuntimeError(
                f"{result_path} uses generalized-alpha outside reviewed "
                "Step 0 Cases A/B / Step 2 Case B modes"
            )
        if not {"tau", "pres"}.issubset(histories):
            raise RuntimeError(
                f"{result_path} is missing generalized-alpha load histories"
            )
        expected_tau = np.zeros(len(times), dtype=float)
        expected_pressure = np.zeros(len(times), dtype=float)
        if generalized_alpha_configuration.active_stress_enabled:
            expected_tau[1:] = tau_of_t(
                load_evaluation_times[1:],
                p=generalized_alpha_configuration.activation_parameters,
                t_span=(0.0, load_horizon),
            )
        if generalized_alpha_configuration.pressure_enabled:
            expected_pressure[1:] = p_of_t(
                load_evaluation_times[1:],
                p=generalized_alpha_configuration.pressure_parameters,
                t_span=(0.0, load_horizon),
            )
        if not np.allclose(
            histories["tau"],
            expected_tau,
            rtol=2.0e-6,
            atol=1.0e-3,
        ):
            raise RuntimeError(
                f"{result_path} has inconsistent generalized-alpha activation history"
            )
        if not np.allclose(
            histories["pres"],
            expected_pressure,
            rtol=2.0e-6,
            atol=1.0e-3,
        ):
            raise RuntimeError(
                f"{result_path} has inconsistent generalized-alpha pressure history"
            )
    apex_offset = float(_scalar(archive, "apex_offset", result_path))
    if not np.isfinite(apex_offset) or apex_offset < 0.0:
        raise RuntimeError(f"{result_path} has invalid apex_offset")
    topology_recorded = "mesh_topology" in archive
    topology = (
        str(_scalar(archive, "mesh_topology", result_path))
        if topology_recorded
        else "polar_ring"
    )
    if topology not in MESH_TOPOLOGIES:
        raise RuntimeError(
            f"{result_path} has unsupported mesh topology {topology!r}"
        )
    mesh = {
        "topology": topology,
        "n_t": _integer(
            _scalar(archive, "n_t", result_path),
            "mesh field 'n_t'",
            result_path,
            minimum=1,
        ),
    }
    if topology == "polar_ring":
        for field in ("n_mu", "n_theta"):
            mesh[field] = _integer(
                _scalar(archive, field, result_path),
                f"mesh field {field!r}",
                result_path,
                minimum=1,
            )
        auxiliary = {"n_side", "n_core", "n_radial", "core_half_width"}
        present_auxiliary = auxiliary.intersection(archive.files)
        if present_auxiliary and present_auxiliary != auxiliary:
            raise RuntimeError(
                f"{result_path} has incomplete polar mesh topology fields"
            )
        if present_auxiliary:
            for field in ("n_side", "n_core", "n_radial"):
                value = _integer(
                    _scalar(archive, field, result_path),
                    f"mesh field {field!r}",
                    result_path,
                )
                if value != 0:
                    raise RuntimeError(
                        f"{result_path} has inconsistent polar mesh field {field!r}"
                    )
                mesh[field] = value
            core_half_width = _finite_number(
                _scalar(archive, "core_half_width", result_path),
                "mesh field 'core_half_width'",
                result_path,
                minimum=0.0,
            )
            if core_half_width != 0.0:
                raise RuntimeError(
                    f"{result_path} has inconsistent polar mesh field "
                    "'core_half_width'"
                )
            mesh["core_half_width"] = core_half_width
    else:
        if not topology_recorded:
            raise RuntimeError(
                f"{result_path} cannot infer closed mesh topology"
            )
        for field in ("n_mu", "n_theta", "n_side"):
            value = _integer(
                _scalar(archive, field, result_path),
                f"mesh field {field!r}",
                result_path,
            )
            if value != 0:
                raise RuntimeError(
                    f"{result_path} has inconsistent closed mesh field {field!r}"
                )
            mesh[field] = value
        n_core = _integer(
            _scalar(archive, "n_core", result_path),
            "mesh field 'n_core'",
            result_path,
            minimum=4,
        )
        if n_core % 2:
            raise RuntimeError(f"{result_path} has invalid closed mesh n_core")
        mesh["n_core"] = n_core
        mesh["n_radial"] = _integer(
            _scalar(archive, "n_radial", result_path),
            "mesh field 'n_radial'",
            result_path,
            minimum=1,
        )
        mesh["core_half_width"] = _finite_number(
            _scalar(archive, "core_half_width", result_path),
            "mesh field 'core_half_width'",
            result_path,
            minimum=0.1,
            strict=True,
        )
        if mesh["core_half_width"] >= 1.0 / np.sqrt(2.0):
            raise RuntimeError(
                f"{result_path} has invalid closed mesh core_half_width"
            )
        if apex_offset != 0.0:
            raise RuntimeError(
                f"{result_path} has nonzero apex_offset for a closed mesh"
            )
    flip_helix = _boolean_scalar(archive, "flip_helix", result_path)

    source = {}
    for component in ("app", "core"):
        revision = str(_scalar(archive, f"{component}_revision", result_path))
        tree_state = str(_scalar(archive, f"{component}_tree_state", result_path))
        source_kind = str(_scalar(archive, f"{component}_source_kind", result_path))
        if FULL_REVISION.fullmatch(revision) is None:
            raise RuntimeError(
                f"{result_path} has untraceable {component} revision {revision!r}; "
                "a full 40-hex Git revision is required"
            )
        if (
            component == "core"
            and revision.casefold() not in SUPPORTED_PUBLIC_CORE_REFS
        ):
            raise RuntimeError(
                f"{result_path} records Core revision {revision!r}; this result "
                "schema requires one of the supported public revisions "
                f"{sorted(SUPPORTED_PUBLIC_CORE_REFS)}"
            )
        content_anchored_step2_worktree = (
            component == "app"
            and source_kind == "git-checkout"
            and tree_state == "dirty"
            and benchmark_identity is not None
            and benchmark_identity.get("runtime_source_sha256") is not None
        )
        if (
            source_kind in {"git-checkout", "asserted"}
            and tree_state != "clean"
            and not content_anchored_step2_worktree
        ):
            raise RuntimeError(
                f"{result_path} records {component}_tree_state={tree_state!r}; "
                "benchmark comparison requires a clean source tree"
            )
        if component == "app" and source_kind not in {"git-checkout", "asserted"}:
            raise RuntimeError(
                f"{result_path} has unsupported app source kind {source_kind!r}"
            )
        identity = {
            "revision": revision.lower(),
            "tree_state": tree_state,
            "source_kind": source_kind,
        }
        if component == "core":
            if source_kind == "pep610-vcs":
                source_url = str(_scalar(archive, "core_source_url", result_path))
                if tree_state != "installed" or source_url != PUBLIC_CORE_URL:
                    raise RuntimeError(
                        f"{result_path} has an unqualified PEP 610 Core identity"
                    )
            elif source_kind not in {"git-checkout", "asserted"}:
                raise RuntimeError(
                    f"{result_path} has unsupported core source kind {source_kind!r}"
                )
            source_url = str(_scalar(archive, "core_source_url", result_path))
            if source_url != PUBLIC_CORE_URL:
                raise RuntimeError(f"{result_path} has an unexpected Core source URL")
            identity["source_url"] = source_url
        source[component] = identity

    current = point_sampling == CURRENT_POINT_SAMPLING
    if current and not {"tau", "pres"}.issubset(histories):
        raise RuntimeError(f"{result_path} is missing current retained load histories")
    model_metadata = None
    runtime_versions = None
    mass_representation = None
    fiber_sampling_option = None
    tbar_identity = None
    parameter_variant = None
    material_model_id = None
    isotropic = None
    material_eta_pa_s = None
    viscous_term_active = None
    method_metadata_origin = None
    nodes = None
    elements = None
    peak_displacement = None
    if topology == "closed_multiblock_disk" and not current:
        raise RuntimeError(
            f"{result_path} closed mesh requires current Hex8 connectivity"
        )
    if current:
        nodes = _finite_array(archive, "nodes", result_path)
        if nodes.ndim != 2 or nodes.shape[1] != 3 or len(nodes) < 8:
            raise RuntimeError(f"{result_path} has invalid current node coordinates")
        if "elems" not in archive:
            raise RuntimeError(f"{result_path} is missing current Hex8 connectivity")
        elements = np.asarray(archive["elems"])
        if (
            elements.ndim != 2
            or elements.shape[1] != 8
            or len(elements) < 1
            or not np.issubdtype(elements.dtype, np.integer)
            or np.any(elements < 0)
            or np.any(elements >= len(nodes))
        ):
            raise RuntimeError(f"{result_path} has invalid current Hex8 connectivity")
        n_elem = int(elements.shape[0])
        if topology == "polar_ring":
            expected_elements = mesh["n_t"] * mesh["n_mu"] * mesh["n_theta"]
            if n_elem != expected_elements:
                raise RuntimeError(
                    f"{result_path} has polar mesh/connectivity size disagreement"
                )
        peak_displacement = _finite_array(archive, "U_peak", result_path)
        if peak_displacement.shape != (3 * len(nodes),):
            raise RuntimeError(f"{result_path} has invalid current peak displacement")
        mesh.update(
            {
                "nodes": int(len(nodes)),
                "elements": int(n_elem),
                "degrees_of_freedom": int(len(peak_displacement)),
            }
        )
        physical_fields = {
            "density_kg_m3": ("density", True),
            "base_robin_stiffness_pa_m": ("a_top", False),
            "base_robin_damping_pa_s_m": ("b_top", False),
            "epicardial_robin_stiffness_pa_m": ("a_epi", False),
            "epicardial_robin_damping_pa_s_m": ("b_epi", False),
            "mesh_perturbation_std_m": ("perturb", False),
        }
        model_metadata = {}
        for label, (field, positive) in physical_fields.items():
            model_metadata[label] = _finite_number(
                _scalar(archive, field, result_path),
                field,
                result_path,
                minimum=0.0,
                strict=positive,
            )
        method_metadata = _current_method_metadata(
            archive,
            result_path,
            app_revision=source["app"]["revision"],
            fiber_sampling=fiber_sampling,
            topology=topology,
        )
        mass_representation = method_metadata["mass_representation"]
        fiber_sampling_option = method_metadata["fiber_sampling_option"]
        isotropic = method_metadata["isotropic"]
        material_model_id = method_metadata["material_model_id"]
        if (
            integrator == "generalized-alpha"
            and material_model_id
            != COMPLETE_SWITCH_ENERGY_MATERIAL_MODEL_ID
        ):
            raise RuntimeError(
                f"{result_path} has an unreviewed generalized-alpha material model"
            )
        material_eta_pa_s = method_metadata["material_eta_pa_s"]
        viscous_term_active = method_metadata["viscous_term_active"]
        parameter_variant = method_metadata["parameter_variant"]
        tbar_identity = method_metadata["tbar_identity"]
        method_metadata_origin = method_metadata["metadata_origin"]
        runtime_versions = {}
        for field in (
            "python_version",
            "numpy_version",
            "scipy_version",
            "coupfe_version",
        ):
            value = str(_scalar(archive, field, result_path))
            if not value or value == "unknown":
                raise RuntimeError(
                    f"{result_path} has unavailable current runtime version {field!r}"
                )
            runtime_versions[field] = value
    else:
        n_elem = mesh["n_t"] * mesh["n_mu"] * mesh["n_theta"]
    sampling_metadata = _validate_sampling_metadata(
        archive, result_path, point_sampling, n_elem
    )
    pre_solve_audit = None
    if topology == "closed_multiblock_disk":
        pre_solve_audit = _validate_closed_pre_solve_audit(
            archive,
            result_path,
            case=case,
            topology=topology,
            n_node=len(nodes),
            n_elem=n_elem,
        )
    solver_name, solver_configuration, solver_diagnostics = _validate_solver_metadata(
        archive, result_path, times, require=current
    )
    if integrator == "generalized-alpha" and solver_name != "petsc-snes-mpi":
        raise RuntimeError(
            f"{result_path} uses generalized-alpha outside its reviewed MPI path"
        )
    mpi_metadata = _validate_mpi_companion_metadata(
        archive,
        result_path,
        solver_name=solver_name,
        configuration=solver_configuration,
        diagnostics=solver_diagnostics,
        case=case,
        n_elem=n_elem,
        ndof=(int(len(peak_displacement)) if peak_displacement is not None else 0),
        integrator=integrator,
        dt=dt,
        formulation=formulation,
        topology=topology,
        apex_offset=apex_offset,
        mass_representation=mass_representation,
        fiber_sampling=fiber_sampling,
        fiber_sampling_option=fiber_sampling_option,
        tbar_identity=tbar_identity,
        isotropic=isotropic,
        material_eta_pa_s=material_eta_pa_s,
        viscous_term_active=viscous_term_active,
        parameter_variant=parameter_variant,
        benchmark_identity=benchmark_identity,
    )

    det_f = None
    pressure = None
    if "det_f_gauss_peak" in archive:
        det_f = _finite_array(
            archive, "det_f_gauss_peak", result_path, shape=(n_elem, 8)
        ).copy()
        if np.any(det_f <= 0.0):
            raise RuntimeError(f"{result_path} has nonpositive det_f_gauss_peak")
    elif current:
        raise RuntimeError(f"{result_path} is missing current 8-GP det(F) evidence")

    if formulation in LOCAL_PRESSURE_FORMULATIONS:
        if point_sampling != CURRENT_POINT_SAMPLING:
            raise RuntimeError(
                f"{result_path} pairs local pressure with legacy point sampling"
            )
        pressure = _finite_array(
            archive, "element_pressure_peak_pa", result_path, shape=(n_elem,)
        ).copy()
        if (
            str(_scalar(archive, "material_kernel_formulation", result_path))
            != "standard"
        ):
            raise RuntimeError(
                f"{result_path} has inconsistent local-pressure kernel metadata"
            )
        material_kappa = float(_scalar(archive, "material_kappa_pa", result_path))
        pressure_bulk = float(
            _scalar(archive, "local_pressure_bulk_modulus_pa", result_path)
        )
        if (
            material_kappa != 0.0
            or not np.isfinite(pressure_bulk)
            or pressure_bulk <= 0.0
        ):
            raise RuntimeError(
                f"{result_path} has inconsistent local-pressure material metadata"
            )
        expected_pressure_law = {
            "hex8_local_pressure_p0_condensed_logj": (
                "linear-reference-volume-mean-log-j-v1"
            ),
            "hex8_local_pressure_p0_condensed_mean_logj_paper_j2": (
                "paper-j2-of-reference-volume-weighted-geometric-mean-j-v1"
            ),
        }[formulation]
        if "local_pressure_volume_law" in archive:
            pressure_law = _text_scalar(
                archive, "local_pressure_volume_law", result_path
            )
            if pressure_law != expected_pressure_law:
                raise RuntimeError(
                    f"{result_path} has inconsistent local-pressure volume-law metadata"
                )
        elif formulation == (
            "hex8_local_pressure_p0_condensed_mean_logj_paper_j2"
        ):
            raise RuntimeError(
                f"{result_path} is missing local-pressure volume-law metadata"
            )
        else:
            pressure_law = "linear-reference-volume-mean-log-j-v1"
        if integrator == "generalized-alpha" and _text_scalar(
            archive, "element_pressure_peak_stage", result_path
        ) != "alpha_f_force_stage":
            raise RuntimeError(
                f"{result_path} has inconsistent generalized-alpha pressure stage"
            )
    elif current:
        expected_kernel = (
            "standard"
            if formulation == "hex8_standard_pointwise_kappa"
            else "fbar_mechanics"
        )
        if str(
            _scalar(archive, "material_kernel_formulation", result_path)
        ) != expected_kernel:
            raise RuntimeError(
                f"{result_path} has inconsistent material kernel metadata"
            )
        material_kappa = float(_scalar(archive, "material_kappa_pa", result_path))
        pressure_bulk = float(
            _scalar(archive, "local_pressure_bulk_modulus_pa", result_path)
        )
        if not np.isfinite(material_kappa) or material_kappa <= 0.0 or pressure_bulk != 0.0:
            raise RuntimeError(
                f"{result_path} has inconsistent penalty material metadata"
            )

    if current:
        model_metadata.update(
            {
                "material_kernel_formulation": str(
                    _scalar(archive, "material_kernel_formulation", result_path)
                ),
                "material_model_id": material_model_id,
                "material_kappa_pa": material_kappa,
                "local_pressure_bulk_modulus_pa": pressure_bulk,
                "local_pressure_volume_law": (
                    pressure_law
                    if formulation in LOCAL_PRESSURE_FORMULATIONS
                    else "not-applicable"
                ),
            }
        )

    n_peak = None
    if "n_peak" in archive:
        n_peak = _integer(_scalar(archive, "n_peak", result_path), "n_peak", result_path)
        if n_peak >= len(times):
            raise RuntimeError(f"{result_path} has out-of-range n_peak")

    return {
        "schema": schema,
        "case": case,
        "reference_case": reference_case,
        "benchmark": (
            None
            if benchmark_identity is None
            else {
                key: value
                for key, value in benchmark_identity.items()
                if key != "_configuration"
            }
        ),
        "times": times.copy(),
        "histories": histories,
        "integrator": integrator,
        "generalized_alpha": generalized_alpha_metadata,
        "formulation": formulation,
        "fiber_sampling": fiber_sampling,
        "point_sampling": point_sampling,
        "viscous_rate": viscous_rate,
        "viscous_term_active": viscous_term_active,
        "material_eta_pa_s": material_eta_pa_s,
        "parameter_variant": parameter_variant,
        "mass_representation": mass_representation,
        "isotropic": isotropic,
        "fiber_sampling_option": fiber_sampling_option,
        "tbar_identity": tbar_identity,
        "method_metadata_origin": method_metadata_origin,
        "dt": dt,
        "t_end": t_end,
        "load_horizon": load_horizon,
        "load_horizon_origin": load_horizon_origin,
        "apex_offset": apex_offset,
        "flip_helix": flip_helix,
        "mesh": mesh,
        "model_metadata": model_metadata,
        "sampling_metadata": sampling_metadata,
        "pre_solve_audit": pre_solve_audit,
        "runtime_versions": runtime_versions,
        "source": source,
        "solver_name": solver_name,
        "solver_configuration": solver_configuration,
        "solver_diagnostics": solver_diagnostics,
        "mpi_metadata": mpi_metadata,
        "det_f_gauss_peak": det_f,
        "element_pressure_peak_pa": pressure,
        "n_peak": n_peak,
        "nodes": None if nodes is None else nodes.copy(),
        "peak_displacement": (
            None if peak_displacement is None else peak_displacement.copy()
        ),
    }


def _validate_completed_result(archive, result_path, requested_case=None):
    """Backward-compatible validation entry point used by existing callers."""
    result = _load_validated_result(archive, result_path, requested_case)
    return result["reference_case"], result["times"]


def resolve_reference_dir(value=None):
    """Resolve a direct data directory or an extracted Zenodo archive root."""
    raw = value or os.environ.get("CARDIAC_BENCHMARK_DATA_DIR")
    if not raw:
        raise FileNotFoundError(
            "reference data are external; pass --reference-dir or set "
            "CARDIAC_BENCHMARK_DATA_DIR (Zenodo DOI 10.5281/zenodo.14260459)"
        )
    root = Path(raw).expanduser().resolve()
    candidates = (
        root,
        root / "benchmark_article_data" / "results_time_curves" / "data",
        root / "results_time_curves" / "data",
        root / "extracted" / "benchmark_article_data" / "results_time_curves" / "data",
    )
    for candidate in candidates:
        if candidate.is_dir() and any(candidate.glob("*nonblinded_step_0*_group_*")):
            return candidate
    raise FileNotFoundError(
        f"no benchmark curve files found below {root}; expected "
        "benchmark_article_data/results_time_curves/data"
    )


def _validate_resampling_source(t_src, u_src, t_dst, description):
    t_src = np.asarray(t_src, dtype=float)
    u_src = np.asarray(u_src, dtype=float)
    t_dst = np.asarray(t_dst, dtype=float)
    if (
        t_src.ndim != 1
        or len(t_src) < 2
        or not np.all(np.isfinite(t_src))
        or not np.all(np.diff(t_src) > 0.0)
    ):
        raise RuntimeError(f"{description} has a non-finite or non-monotone time grid")
    if u_src.shape != (len(t_src), 3) or not np.all(np.isfinite(u_src)):
        raise RuntimeError(f"{description} has malformed displacement history")
    if (
        t_dst.ndim != 1
        or len(t_dst) < 2
        or not np.all(np.isfinite(t_dst))
        or not np.all(np.diff(t_dst) > 0.0)
    ):
        raise RuntimeError("canonical comparison grid is invalid")
    if (
        t_src[0] - t_dst[0] > _ENDPOINT_OFFSET_TOLERANCE
        or t_dst[-1] - t_src[-1] > _ENDPOINT_OFFSET_TOLERANCE
        or t_src[0] < t_dst[0] - _ENDPOINT_OFFSET_TOLERANCE
        or t_src[-1] > t_dst[-1] + _ENDPOINT_OFFSET_TOLERANCE
    ):
        raise RuntimeError(
            f"{description} does not cover the canonical [0, 1] comparison interval"
        )
    return t_src, u_src, t_dst


def resample(t_src, u_src, t_dst=CANONICAL_TIME_GRID, *, description="history"):
    """Resample onto a checked grid, allowing the dataset's 0.001 s offsets."""
    t_src, u_src, t_dst = _validate_resampling_source(
        t_src, u_src, t_dst, description
    )
    return np.stack(
        [np.interp(t_dst, t_src, u_src[:, component]) for component in range(3)],
        axis=1,
    )


def _parse_reference_record(record, filename):
    if not isinstance(record, Mapping):
        raise RuntimeError(f"reference file {filename} does not contain a mapping")
    if "time" not in record or "displacement" not in record:
        raise RuntimeError(f"reference file {filename} is missing required data")
    t = np.asarray(record["time"], dtype=float)
    displacement_root = record["displacement"]
    if not isinstance(displacement_root, Mapping):
        raise RuntimeError(f"reference file {filename} has malformed displacement data")
    parsed = {"source_t": t.copy(), "t": CANONICAL_TIME_GRID.copy()}
    for point in ("p0", "p1"):
        if point not in displacement_root:
            raise RuntimeError(f"reference file {filename} is missing {point}")
        displacement = displacement_root[point]
        if isinstance(displacement, np.ndarray) and displacement.shape == ():
            displacement = displacement.item()
        if not isinstance(displacement, Mapping):
            raise RuntimeError(f"reference file {filename} has malformed {point}")
        try:
            raw = np.stack(
                [
                    np.asarray(displacement[component], dtype=float)
                    for component in ("ux", "uy", "uz")
                ],
                axis=1,
            )
        except (KeyError, TypeError, ValueError) as error:
            raise RuntimeError(f"reference file {filename} has malformed {point}") from error
        parsed[point] = resample(
            t,
            raw,
            CANONICAL_TIME_GRID,
            description=f"reference file {filename} {point}",
        )
    return parsed


def _reference_filename(case, suffix):
    return f"monoventricular_nonblinded_{case}_group_{suffix}.pickle"


def load_reference(
    case="step_0A",
    reference_dir=None,
    *,
    return_provenance=False,
    return_selection=False,
):
    """Return the exact upstream 10-team manifest on the canonical grid."""
    if case == "step_2B":
        raise ValueError(
            "Step 2 Case B uses the blinded publisher selection; run "
            "compare_step2b_case_b.py for its hash-pinned comparison"
        )
    if case not in CASE_NAMES.values():
        raise ValueError(f"unsupported reference case {case!r}")
    data_dir = resolve_reference_dir(reference_dir)
    expected_names = [
        _reference_filename(case, suffix)
        for suffix in REFERENCE_MANIFEST_SUFFIXES
    ]
    alias_name = _reference_filename(case, REFERENCE_EXCLUDED_ALIAS_SUFFIX)
    matching = {
        path.name: path
        for path in data_dir.glob(f"*nonblinded_{case}_group_*")
    }
    missing = sorted(set(expected_names) - set(matching))
    if missing:
        raise RuntimeError(
            f"reference manifest for {case} is missing selected files: {missing}"
        )
    unexpected = sorted(set(matching) - set(expected_names) - {alias_name})
    if unexpected:
        raise RuntimeError(
            f"reference manifest for {case} has unexpected matching files: "
            f"{unexpected}"
        )

    out = {}
    provenance = []
    marker = f"{case}_group_"
    selected_payloads = {}
    for name in expected_names:
        path = matching[name]
        if not path.is_file():
            raise RuntimeError(f"selected reference entry is not a file: {path.name}")
        label = name.split(marker, 1)[1].rsplit(".", 1)[0]
        if not label or label in out:
            raise RuntimeError(f"duplicate or empty reference team label in {name}")
        raw = path.read_bytes()
        selected_payloads[name] = raw
        try:
            record = pickle.loads(raw)  # trusted Zenodo dataset only
            parsed = _parse_reference_record(record, name)
        except RuntimeError:
            raise
        except Exception as error:
            raise RuntimeError(f"could not read trusted reference file {name}") from error
        out[label] = parsed
        provenance.append(
            {
                "team": label,
                "filename": name,
                "sha256": hashlib.sha256(raw).hexdigest(),
                "size_bytes": len(raw),
                "source_sample_count": int(len(parsed["source_t"])),
                "source_time_start": float(parsed["source_t"][0]),
                "source_time_end": float(parsed["source_t"][-1]),
            }
        )

    excluded_aliases = []
    if alias_name in matching:
        alias_path = matching[alias_name]
        if not alias_path.is_file():
            raise RuntimeError(
                f"excluded reference alias is not a file: {alias_path.name}"
            )
        alias_raw = alias_path.read_bytes()
        target_name = _reference_filename(case, REFERENCE_ALIAS_TARGET_SUFFIX)
        if alias_raw != selected_payloads[target_name]:
            raise RuntimeError(
                f"known reference alias {alias_name} is not byte-identical to "
                f"selected {target_name}"
            )
        excluded_aliases.append(
            {
                "filename": alias_name,
                "sha256": hashlib.sha256(alias_raw).hexdigest(),
                "size_bytes": len(alias_raw),
                "identical_to_selected_filename": target_name,
                "reason": (
                    "Byte-identical base-name alias of selected SimVascular P2; "
                    "excluded because upstream figures.py does not select it."
                ),
            }
        )

    selection = {
        "policy": REFERENCE_SELECTION_POLICY,
        "upstream_figures_py_identity": REFERENCE_FIGURES_PY_IDENTITY,
        "upstream_manifest_variable": REFERENCE_MANIFEST_VARIABLES[case],
        "selected_count": len(expected_names),
        "selected_files": expected_names,
        "excluded_aliases": excluded_aliases,
    }
    if return_provenance and return_selection:
        return out, provenance, selection
    if return_provenance:
        return out, provenance
    if return_selection:
        return out, selection
    return out


def red(u, mean):
    """Relative discrepancy (benchmark paper Eq. 21)."""
    u = np.asarray(u, dtype=float)
    mean = np.asarray(mean, dtype=float)
    if u.shape != mean.shape or u.ndim != 2 or u.shape[1] != 3:
        raise ValueError("RED inputs must be matching (n, 3) histories")
    if not np.all(np.isfinite(u)) or not np.all(np.isfinite(mean)):
        raise ValueError("RED inputs must be finite")
    numerator = np.linalg.norm(u - mean, axis=1)
    denominator = np.linalg.norm(mean, axis=1) + 1.0e-30
    return float(np.mean(numerator / denominator))


def _sha256_file(path):
    digest = hashlib.sha256()
    size = 0
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
            size += len(chunk)
    return digest.hexdigest(), size


def _normalized_run_log(path):
    raw = Path(path).read_bytes()
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as error:
        raise RuntimeError("--run-log must be UTF-8 text") from error
    if "\x00" in text:
        raise RuntimeError("--run-log contains a NUL byte")
    normalized = text.replace("\r\n", "\n").replace("\r", "\n")
    if normalized and not normalized.endswith("\n"):
        normalized += "\n"
    payload = normalized.encode("utf-8")
    return {
        "filename": Path(path).name,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size_bytes": len(payload),
        "normalization": "UTF-8; CRLF/CR converted to LF; final LF added when nonempty",
    }


def _summary(values):
    if values is None:
        return {"available": False}
    flat = np.asarray(values, dtype=float).ravel()
    return {
        "available": True,
        "shape": list(np.shape(values)),
        "count": int(flat.size),
        "minimum": float(np.min(flat)),
        "mean": float(np.mean(flat)),
        "maximum": float(np.max(flat)),
    }


def _peak_circumferential_ring_rotation(result):
    """Measure Case B base/apex ring rotation about the global long axis."""
    if result["case"] != "B":
        return {
            "available": False,
            "reason": "Reported only for Case B results.",
        }
    nodes = result.get("nodes")
    peak_displacement = result.get("peak_displacement")
    mesh = result["mesh"]
    if mesh.get("topology", "polar_ring") != "polar_ring":
        return {
            "available": False,
            "reason": (
                "Metric is defined only for the ordered rings of the historical "
                "polar mesh, not the closed multiblock topology."
            ),
        }
    n_t = mesh["n_t"]
    n_mu = mesh["n_mu"]
    n_theta = mesh["n_theta"]
    expected_nodes = (n_t + 1) * (n_mu + 1) * n_theta
    if (
        nodes is None
        or peak_displacement is None
        or nodes.shape != (expected_nodes, 3)
        or peak_displacement.shape != (3 * expected_nodes,)
    ):
        return {
            "available": False,
            "reason": (
                "Result does not contain the complete ordered structured rings "
                "required by this metric."
            ),
        }

    deformed = nodes + peak_displacement.reshape(expected_nodes, 3)
    layers = []
    for layer in range(n_t + 1):
        angles = []
        for longitudinal_ring in (0, n_mu):
            first = ((layer * (n_mu + 1)) + longitudinal_ring) * n_theta
            indices = np.arange(first, first + n_theta)
            reference_yz = nodes[indices, 1:3]
            deformed_yz = deformed[indices, 1:3]
            reference_yz = reference_yz - reference_yz.mean(axis=0)
            deformed_yz = deformed_yz - deformed_yz.mean(axis=0)
            cross = np.sum(
                reference_yz[:, 0] * deformed_yz[:, 1]
                - reference_yz[:, 1] * deformed_yz[:, 0]
            )
            dot = np.sum(reference_yz * deformed_yz)
            if not np.isfinite(cross) or not np.isfinite(dot) or (
                cross == 0.0 and dot == 0.0
            ):
                raise RuntimeError(
                    "peak circumferential ring rotation is undefined for a "
                    "zero-radius or non-finite ring"
                )
            angles.append(float(np.degrees(np.arctan2(cross, dot))))
        relative = (angles[1] - angles[0] + 180.0) % 360.0 - 180.0
        layers.append(
            {
                "transmural_node_layer_index": layer,
                "transmural_coordinate": layer / n_t,
                "base_rotation_degrees": angles[0],
                "apex_rotation_degrees": angles[1],
                "relative_apex_minus_base_degrees": relative,
            }
        )

    return {
        "available": True,
        "method": "least-squares rotation of centered circumferential rings",
        "coordinate_system": (
            "global Cartesian; long axis +x; circumferential plane yz"
        ),
        "centering": (
            "subtract separate reference and deformed yz centroids per ring"
        ),
        "angle_definition": (
            "atan2(sum(Y_y*Yd_z - Y_z*Yd_y), "
            "sum(Y_y*Yd_y + Y_z*Yd_z))"
        ),
        "sign_convention": "positive is right-handed about global +x (+y toward +z)",
        "relative_angle_definition": (
            "apex minus base, wrapped to [-180, 180) degrees"
        ),
        "peak_time_s": float(result["times"][result["n_peak"]]),
        "base_longitudinal_ring_index": 0,
        "apex_longitudinal_ring_index": n_mu,
        "circumferential_nodes_per_ring": n_theta,
        "layers": layers,
    }


def _build_report(
    result_path,
    result,
    ref,
    reference_files,
    reference_selection,
    mine_on_grid,
    *,
    run_log=None,
    supersedes_report_sha256=None,
):
    labels = sorted(ref)
    comparisons = {}
    means = {}
    for point in ("p0", "p1"):
        stack = np.stack([ref[label][point] for label in labels], axis=0)
        mean = stack.mean(axis=0)
        means[point] = mean
        comparisons[point] = {
            "ours": red(mine_on_grid[point], mean),
            "teams": {label: red(ref[label][point], mean) for label in labels},
        }

    result_sha, result_size = _sha256_file(result_path)
    histories = {
        "times_s": result["times"].tolist(),
        "u0_m": result["histories"]["u0"].tolist(),
        "u1_m": result["histories"]["u1"].tolist(),
    }
    if "tau" in result["histories"]:
        histories["active_tension_pa"] = result["histories"]["tau"].tolist()
    if "pres" in result["histories"]:
        histories["cavity_pressure_pa"] = result["histories"]["pres"].tolist()

    generalized_alpha = result.get("generalized_alpha")
    if generalized_alpha is not None:
        generalized_alpha = dict(generalized_alpha)
        generalized_alpha["load_evaluation_times_s"] = np.asarray(
            generalized_alpha["load_evaluation_times_s"], dtype=float
        ).tolist()

    peak = {"available": result["n_peak"] is not None}
    if result["n_peak"] is not None:
        n_peak = result["n_peak"]
        peak.update(
            {
                "index": n_peak,
                "time_s": float(result["times"][n_peak]),
                "u0_m": result["histories"]["u0"][n_peak].tolist(),
                "u1_m": result["histories"]["u1"][n_peak].tolist(),
            }
        )
        for key, label in (("tau", "active_tension_pa"), ("pres", "cavity_pressure_pa")):
            if key in result["histories"]:
                peak[label] = float(result["histories"][key][n_peak])

    report = {
        "schema": REPORT_SCHEMA,
        "bounded_claim": (
            "This report compares one completed, source-identified CoupFE-Cardiac "
            "run with the separately distributed benchmark time curves on a common "
            "grid. It is example-level evidence, not clinical, device, or broad "
            "solver validation. It does not by itself establish mesh/time "
            "convergence, equivalence to the reference P2 discretization, or a "
            "unique twist direction. A closed-topology statement, when present, "
            "is bounded to the retained pre-solve geometry and boundary audits."
        ),
        "result": {
            "filename": Path(result_path).name,
            "sha256": result_sha,
            "size_bytes": result_size,
            "result_schema": result["schema"],
            "case": result["case"],
            "reference_case": result["reference_case"],
            "configuration": {
                "integrator": result["integrator"],
                "generalized_alpha": generalized_alpha,
                "formulation": result["formulation"],
                "mass_representation": result["mass_representation"],
                "method_metadata_origin": result["method_metadata_origin"],
                "isotropic": result["isotropic"],
                "fiber_sampling": result["fiber_sampling"],
                "fiber_sampling_option": result["fiber_sampling_option"],
                "tbar": result["tbar_identity"],
                "point_sampling": result["point_sampling"],
                "viscous_rate": result["viscous_rate"],
                "material_eta_pa_s": result["material_eta_pa_s"],
                "viscous_term_active": result["viscous_term_active"],
                "parameter_variant": result["parameter_variant"],
                "dt_s": result["dt"],
                "t_end_s": result["t_end"],
                "load_horizon_s": result["load_horizon"],
                "load_horizon_origin": result["load_horizon_origin"],
                "apex_offset_rad": result["apex_offset"],
                "flip_helix": result["flip_helix"],
                "mesh": result["mesh"],
                "nonlinear_solver": result["solver_name"],
                "model_parameters": result["model_metadata"],
                "sampling_points": result["sampling_metadata"],
            },
            "runtime_versions": result["runtime_versions"],
            "source_identity": result["source"],
            "solver_configuration": result["solver_configuration"],
            "nonlinear_step_diagnostics": result["solver_diagnostics"],
            "mpi_metadata": result["mpi_metadata"],
            "pre_solve_audit": result["pre_solve_audit"],
            "retained_histories": histories,
            "peak": peak,
            "peak_circumferential_ring_rotation": (
                _peak_circumferential_ring_rotation(result)
            ),
            "det_f_gauss_peak_summary": _summary(result["det_f_gauss_peak"]),
            "element_pressure_peak_pa_summary": _summary(
                result["element_pressure_peak_pa"]
            ),
        },
        "reference": {
            "doi": REFERENCE_DOI,
            "license": REFERENCE_LICENSE,
            "published_archive_identity": REFERENCE_ARCHIVE_IDENTITY,
            "case": result["reference_case"],
            "canonical_grid_s": CANONICAL_TIME_GRID.tolist(),
            "endpoint_offset_tolerance_s": _ENDPOINT_OFFSET_TOLERANCE,
            "selection": reference_selection,
            "team_files": reference_files,
            "mean_curves_m": {
                point: means[point].tolist() for point in ("p0", "p1")
            },
        },
        "comparison": {
            "metric": "relative discrepancy (benchmark paper Eq. 21)",
            "red": comparisons,
            "ours_on_canonical_grid_m": {
                point: mine_on_grid[point].tolist() for point in ("p0", "p1")
            },
        },
    }
    if result.get("benchmark") is not None:
        report["result"]["configuration"]["benchmark"] = result["benchmark"]
    if supersedes_report_sha256 is not None:
        if (
            not isinstance(supersedes_report_sha256, str)
            or re.fullmatch(r"[0-9a-f]{64}", supersedes_report_sha256) is None
        ):
            raise RuntimeError(
                "superseded report SHA-256 must be 64 lowercase hexadecimal digits"
            )
        report["correction"] = {
            "supersedes_report_sha256": supersedes_report_sha256,
            "predecessor_repository_revision": CORRECTION_PREDECESSOR_REVISION,
            "reason": CORRECTION_REASON,
        }
    if run_log is not None:
        report["result"]["normalized_run_log"] = _normalized_run_log(run_log)
    _validate_finite_json(report, "generated report", Path(result_path).name)
    _validate_no_absolute_paths(report, "generated report")
    return report


def _write_json_atomic(path, value):
    destination = Path(path).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
            stream.write(payload)
        os.replace(temporary, destination)
    except BaseException:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result", type=Path)
    parser.add_argument(
        "--case",
        choices=("step_0A", "step_0B", "step_2B"),
        help="optional cross-check; otherwise derive the reference case from the result",
    )
    parser.add_argument("--reference-dir")
    parser.add_argument("--plot", type=Path, help="optional output PNG")
    parser.add_argument("--json", type=Path, dest="json_output", help="machine-readable report")
    parser.add_argument(
        "--run-log",
        type=Path,
        help="optional UTF-8 stdout transcript to normalize, hash, and identify",
    )
    parser.add_argument(
        "--supersedes-report-sha256",
        help=(
            "optional exact SHA-256 of a predecessor report corrected by this "
            "generated report"
        ),
    )
    args = parser.parse_args(argv)

    with np.load(args.result, allow_pickle=False) as archive:
        result = _load_validated_result(archive, args.result, requested_case=args.case)

    if result["reference_case"] == "step_2B":
        parser.error(
            "Step 2 Case B uses the blinded publisher selection; run "
            "compare_step2b_case_b.py for its hash-pinned comparison"
        )

    ref, reference_files, reference_selection = load_reference(
        result["reference_case"],
        args.reference_dir,
        return_provenance=True,
        return_selection=True,
    )
    labels = sorted(ref)
    print(
        f"case {result['reference_case']}; reference teams ({len(labels)}): {labels}"
    )
    mine_on_grid = {
        point: resample(
            result["times"],
            result["histories"][point],
            CANONICAL_TIME_GRID,
            description=f"CoupFE result {point}",
        )
        for point in ("u0", "u1")
    }
    mine_on_grid = {"p0": mine_on_grid["u0"], "p1": mine_on_grid["u1"]}

    report = _build_report(
        args.result,
        result,
        ref,
        reference_files,
        reference_selection,
        mine_on_grid,
        run_log=args.run_log,
        supersedes_report_sha256=args.supersedes_report_sha256,
    )
    print("\n=== Relative Discrepancy (RED) vs all-team mean ===")
    for point in ("p0", "p1"):
        values = report["comparison"]["red"][point]
        team_red = sorted(values["teams"].values())
        print(
            f"  {point}: ours={values['ours']:.3f}; teams=[{team_red[0]:.3f}, "
            f"{team_red[-1]:.3f}], median={np.median(team_red):.3f}"
        )

    if args.json_output:
        _write_json_atomic(args.json_output, report)
        print(f"saved JSON report -> {args.json_output}")

    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
        for column, point in enumerate(("p0", "p1")):
            for label in labels:
                for component in range(3):
                    axes[component, column].plot(
                        CANONICAL_TIME_GRID,
                        ref[label][point][:, component] * 1.0e3,
                        color="0.6",
                        lw=0.8,
                        alpha=0.7,
                    )
            for component, name in enumerate(("x", "y", "z")):
                axes[component, column].plot(
                    CANONICAL_TIME_GRID,
                    mine_on_grid[point][:, component] * 1.0e3,
                    "r-",
                )
                axes[component, column].set_ylabel(f"u_{name} [mm]")
                axes[component, column].grid(True, alpha=0.3)
            axes[0, column].set_title(point)
            axes[2, column].set_xlabel("time [s]")
        fig.tight_layout()
        args.plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.plot, dpi=120)
        print(f"saved plot -> {args.plot}")

    return report


if __name__ == "__main__":
    main()
