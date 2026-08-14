#!/usr/bin/env python3
"""Compare one completed CoupFE Case B trajectory with explicit FEniCS arrays.

This utility reads only five caller-selected inputs: one CoupFE ``.npz``
archive, one FEniCS ``parameters.json``, one FEniCS time array, and the two
FEniCS landmark-displacement arrays.  It does not read pickle, HDF5, stress,
or caller-directory defaults.  Every input is identified by basename, byte
count, and SHA-256 in the report; resolved filesystem paths are never copied
to public output.

The comparison is deliberately narrow.  It checks the fixed Case B physical
setup, maps CoupFE p0/p1 histories to the retained FEniCS 0.001--0.999 s grid,
reports declared displacement metrics, and produces a six-panel overlay.
It is not a stress comparison, a convergence study, or broad validation.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PureWindowsPath
import re
import tempfile
from typing import Dict, Iterable, Mapping, Optional, Tuple

import numpy as np

try:  # package import
    from .activation import p_of_t
    from .benchmark_parameters import (
        RUNTIME_SOURCE_FILES,
        benchmark_configuration,
        benchmark_metadata,
    )
except ImportError:  # direct script import
    from activation import p_of_t
    from benchmark_parameters import (
        RUNTIME_SOURCE_FILES,
        benchmark_configuration,
        benchmark_metadata,
    )


REPORT_SCHEMA = "coupfe-cardiac-fenics-case-b-comparison-v1"
RESULT_SCHEMA = "coupfe-cardiac-result-v1"
HASH_MANIFEST_SCHEMA = "coupfe-cardiac-comparison-input-hashes-v1"
MATERIAL_MODEL_ID = "holzapfel-ogden-smooth-switch-complete-energy-derivative-v1"
PUBLIC_CORE_URL = "https://github.com/tengzhang48/CoupFE.git"

EXPECTED_DT_S = 1.0e-3
EXPECTED_STEPS = 1000
FENICS_SAMPLES = 999
FENICS_TIME = np.arange(1, 1000, dtype=float) * EXPECTED_DT_S
COUPFE_TIME = np.arange(0, 1001, dtype=float) * EXPECTED_DT_S
SNAP_START_S = 0.20
SNAP_END_S = 0.32
SNAP_ONSET_THRESHOLD_M = -5.0e-3

POINTS = ("p0", "p1")
COMPONENTS = ("x", "y", "z")
INPUT_ROLES = (
    "coupfe-run",
    "fenics-parameters",
    "fenics-times",
    "fenics-p0",
    "fenics-p1",
)
FULL_SHA256 = re.compile(r"[0-9a-fA-F]{64}")
FULL_GIT_REVISION = re.compile(r"[0-9a-fA-F]{40}")
NAIVE_ISO_TIMESTAMP = re.compile(
    r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}(?:\.\d{1,6})?"
)

MPI_IMPLEMENTATION = "cardiac-owned-distributed-closed-std-kappa-step0"
MPI_GENERALIZED_ALPHA_IMPLEMENTATION = (
    "cardiac-owned-distributed-closed-std-kappa-generalized-alpha-step0"
)
MPI_IMPLEMENTATIONS = {
    "be": MPI_IMPLEMENTATION,
    "generalized-alpha": MPI_GENERALIZED_ALPHA_IMPLEMENTATION,
}
BENCHMARK_IDENTITY_FIELDS = frozenset(
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
GENERALIZED_ALPHA_PARAMETERS = {
    "alpha_m": 0.2,
    "alpha_f": 0.4,
    "gamma": 0.7,
    "beta": 0.36,
}
GENERALIZED_ALPHA_CONFIGURATION = {
    **GENERALIZED_ALPHA_PARAMETERS,
    "parameter_source": "finsberg/cardiac_benchmark problem.py defaults",
    "acceleration_stage": "alpha_m*a_n + (1-alpha_m)*a_np1",
    "force_stage": "alpha_f*x_n + (1-alpha_f)*x_np1",
    "load_time": "t_np1 - alpha_f*dt",
}
MPI_BASE_FIELDS = frozenset(
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
MPI_MASS_FIELDS = frozenset(
    {
        "mpi_mass_partition",
        "mpi_mass_owned_row_ranges",
        "mpi_mass_local_nnz",
        "mpi_mass_touching_element_counts",
    }
)
MPI_FIXED_CONFIGURATION = {
    "implementation": MPI_IMPLEMENTATION,
    "settings_source": "recovered 2026-06-27 Case B development adapter",
    "communicator": "PETSc.COMM_WORLD",
    "matrix_scope": "persistent one solver instance per MPI run",
    "dirichlet_support": "none",
    "rank_local_element_assembly": True,
    "global_mesh_replicated": True,
    "collective_invalid_trial_policy": "all-owned-residual-entries-inf-for-bt",
    "function_domain_rejection_api": "nonfinite residual for PETSc BT",
    "commit_policy": "independent-global-residual-check-before-commit",
    "factor_solver_type": "superlu_dist",
    "configured_factor_solver_type": "superlu_dist",
    "snes_type": "newtonls",
    "line_search_type": "bt",
    "ksp_type": "preonly",
    "pc_type": "lu",
    "rtol": 1.0e-9,
    "atol": 1.0e-10,
    "stol": 1.0e-12,
    "max_it": 60,
}

EXPECTED_LANDMARKS_M = {
    "p0": np.array([0.025, 0.030, 0.0]),
    "p1": np.array([0.000, 0.030, 0.0]),
}

# Values retained in the local FEniCS parameters.json.  Keys not listed here
# (notably outdir/outpath) are neither trusted nor copied to the report.
EXPECTED_FENICS_NUMERIC = {
    ("problem_parameters", "alpha_epi"): 1.0e8,
    ("problem_parameters", "alpha_f"): 0.4,
    ("problem_parameters", "alpha_m"): 0.2,
    ("problem_parameters", "alpha_top"): 1.0e5,
    ("problem_parameters", "beta_epi"): 5.0e3,
    ("problem_parameters", "beta_top"): 5.0e3,
    ("problem_parameters", "dt"): 1.0e-3,
    ("problem_parameters", "p"): 0.0,
    ("problem_parameters", "rho"): 1000.0,
    ("fiber_parameters", "alpha_endo"): -60.0,
    ("fiber_parameters", "alpha_epi"): 60.0,
    ("material_parameters", "a"): 59.0,
    ("material_parameters", "b"): 8.023,
    ("material_parameters", "a_f"): 18472.0,
    ("material_parameters", "b_f"): 16.026,
    ("material_parameters", "a_s"): 2481.0,
    ("material_parameters", "b_s"): 11.12,
    ("material_parameters", "a_fs"): 216.0,
    ("material_parameters", "b_fs"): 11.436,
    ("material_parameters", "eta"): 100.0,
    ("material_parameters", "k"): 100.0,
    ("material_parameters", "kappa"): 1.0e6,
    ("mesh_parameters", "mesh_size_factor"): 1.0,
    ("mesh_parameters", "mu_apex_endo"): -np.pi,
    ("mesh_parameters", "mu_apex_epi"): -np.pi,
    ("mesh_parameters", "mu_base_endo"): -1.2722641256100204,
    ("mesh_parameters", "mu_base_epi"): -1.318116071652818,
    ("mesh_parameters", "psize_ref"): 0.005,
    ("mesh_parameters", "r_long_endo"): 0.09,
    ("mesh_parameters", "r_long_epi"): 0.097,
    ("mesh_parameters", "r_short_endo"): 0.025,
    ("mesh_parameters", "r_short_epi"): 0.035,
    ("activation_parameters", "a_max"): 5.0,
    ("activation_parameters", "a_min"): -30.0,
    ("activation_parameters", "gamma"): 0.005,
    ("activation_parameters", "sigma_0"): 150000.0,
    ("activation_parameters", "t_dias"): 0.484,
    ("activation_parameters", "t_sys"): 0.16,
    ("pressure_parameters", "a_max"): 5.0,
    ("pressure_parameters", "a_min"): -30.0,
    ("pressure_parameters", "alpha_mid"): 1.0,
    ("pressure_parameters", "alpha_pre"): 5.0,
    ("pressure_parameters", "gamma"): 0.005,
    ("pressure_parameters", "sigma_mid"): 16000.0,
    ("pressure_parameters", "sigma_pre"): 7000.0,
    ("pressure_parameters", "t_dias_pre"): 0.484,
    ("pressure_parameters", "t_sys_pre"): 0.17,
}


class ComparisonInputError(RuntimeError):
    """An input is unsafe, incomplete, or outside the declared comparison."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ComparisonInputError(message)


def sha256_file(path: Path) -> Tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
                size += len(chunk)
    except OSError as error:
        raise ComparisonInputError(
            "cannot read input file {0}: {1}".format(path.name, error)
        ) from error
    return digest.hexdigest(), size


def _identity(path: Path) -> Dict[str, object]:
    path = path.expanduser().resolve()
    require(path.is_file(), "input file does not exist: {0}".format(path))
    digest, size = sha256_file(path)
    filename = _portable_basename(path.name, "input filename")
    return {"filename": filename, "sha256": digest, "size_bytes": size}


def _portable_basename(value: str, description: str) -> str:
    """Require one platform-neutral basename with no directory component."""
    value = str(value)
    posix_name = Path(value).name
    windows_name = PureWindowsPath(value).name
    require(
        bool(value)
        and value == posix_name == windows_name
        and not Path(value).is_absolute()
        and not PureWindowsPath(value).is_absolute(),
        "{0} is not a portable basename".format(description),
    )
    return value


def _sanitized_basename(value: str) -> str:
    """Return only the final component of either a POSIX or Windows path."""
    value = str(value)
    windows_name = PureWindowsPath(value).name
    return Path(windows_name).name


def _finite_json(value, location: str = "JSON") -> None:
    if value is None or isinstance(value, (str, bool)):
        return
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        require(np.isfinite(float(value)), "non-finite value in {0}".format(location))
        return
    if isinstance(value, list):
        for index, child in enumerate(value):
            _finite_json(child, "{0}[{1}]".format(location, index))
        return
    if isinstance(value, dict):
        require(
            all(isinstance(key, str) for key in value),
            "non-string key in {0}".format(location),
        )
        for key, child in value.items():
            _finite_json(child, "{0}.{1}".format(location, key))
        return
    raise ComparisonInputError(
        "unsupported value {0} in {1}".format(type(value).__name__, location)
    )


def _load_json(path: Path, description: str):
    def reject_constant(value):
        raise ValueError("non-finite JSON constant {0}".format(value))

    try:
        value = json.loads(
            path.read_text(encoding="utf-8"), parse_constant=reject_constant
        )
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError) as error:
        raise ComparisonInputError(
            "cannot parse {0} {1}: {2}".format(description, path.name, error)
        ) from error
    _finite_json(value, description)
    return value


def _scalar(archive, key: str, path: Path):
    require(key in archive, "{0} is missing scalar {1!r}".format(path.name, key))
    value = np.asarray(archive[key])
    require(
        value.shape == (),
        "{0} field {1!r} is not scalar".format(path.name, key),
    )
    return value.item()


def _boolean_scalar(archive, key: str, path: Path) -> bool:
    value = _scalar(archive, key, path)
    require(
        isinstance(value, (bool, np.bool_)),
        "{0} field {1!r} is not boolean".format(path.name, key),
    )
    return bool(value)


def _integer(value, description: str, minimum: Optional[int] = None) -> int:
    require(
        not isinstance(value, (bool, np.bool_)),
        "{0} is boolean, not an integer".format(description),
    )
    try:
        numeric = float(value)
    except (TypeError, ValueError) as error:
        raise ComparisonInputError("{0} is not an integer".format(description)) from error
    require(
        np.isfinite(numeric) and numeric.is_integer(),
        "{0} is not an integer".format(description),
    )
    result = int(numeric)
    if minimum is not None:
        require(result >= minimum, "{0} is below {1}".format(description, minimum))
    return result


def _finite_array(archive, key: str, path: Path, shape=None) -> np.ndarray:
    require(key in archive, "{0} is missing array {1!r}".format(path.name, key))
    try:
        value = np.asarray(archive[key], dtype=float)
    except (TypeError, ValueError) as error:
        raise ComparisonInputError(
            "{0} field {1!r} is not numeric".format(path.name, key)
        ) from error
    if shape is not None:
        require(
            value.shape == shape,
            "{0} field {1!r} has shape {2}; expected {3}".format(
                path.name, key, value.shape, shape
            ),
        )
    require(
        np.all(np.isfinite(value)),
        "{0} field {1!r} contains non-finite values".format(path.name, key),
    )
    return value.copy()


def _exact_float(observed, expected: float, description: str) -> float:
    require(
        not isinstance(observed, (bool, np.bool_)),
        "{0} is boolean, not numeric".format(description),
    )
    try:
        value = float(observed)
    except (TypeError, ValueError) as error:
        raise ComparisonInputError("{0} is not numeric".format(description)) from error
    tolerance = 2.0e-13 * max(1.0, abs(float(expected)))
    require(
        np.isfinite(value) and abs(value - expected) <= tolerance,
        "{0}={1!r}; expected {2!r}".format(description, value, expected),
    )
    return value


def _nested(mapping: Mapping, keys: Iterable[str], description: str):
    value = mapping
    walked = []
    for key in keys:
        walked.append(key)
        require(
            isinstance(value, Mapping) and key in value,
            "{0} is missing {1}".format(description, ".".join(walked)),
        )
        value = value[key]
    return value


def _validate_expected_hashes(
    identities: Mapping[str, Mapping[str, object]],
    expected: Mapping[str, Mapping[str, object]],
) -> None:
    for role, record in expected.items():
        require(role in INPUT_ROLES, "unknown expected-hash role {0!r}".format(role))
        actual = identities[role]
        digest = str(record.get("sha256", ""))
        require(
            FULL_SHA256.fullmatch(digest) is not None,
            "expected hash for {0} is not SHA-256".format(role),
        )
        require(
            digest.lower() == actual["sha256"],
            "SHA-256 mismatch for {0}: expected {1}, observed {2}".format(
                role, digest.lower(), actual["sha256"]
            ),
        )
        if "filename" in record:
            require(
                str(record["filename"]) == actual["filename"],
                "filename mismatch for expected-hash role {0}".format(role),
            )
        if "size_bytes" in record:
            size = record["size_bytes"]
            require(
                isinstance(size, int)
                and not isinstance(size, bool)
                and size == actual["size_bytes"],
                "byte-count mismatch for expected-hash role {0}".format(role),
            )


def _hash_verification_record(
    identities: Mapping[str, Mapping[str, object]],
    expected: Mapping[str, Mapping[str, object]],
    *,
    retained: bool,
) -> Dict[str, object]:
    """Describe the expected-hash gate without exposing caller paths."""
    missing = [role for role in INPUT_ROLES if role not in expected]
    if retained:
        require(
            not missing,
            "retained mode requires expected SHA-256 for all five input roles; "
            "missing: {0}".format(", ".join(missing)),
        )
    verified = {
        role: {
            "filename": str(identities[role]["filename"]),
            "sha256": str(identities[role]["sha256"]),
            "size_bytes": int(identities[role]["size_bytes"]),
        }
        for role in INPUT_ROLES
        if role in expected
    }
    return {
        "mode": "retained" if retained else "development",
        "all_input_roles_required": retained,
        "all_input_roles_verified": not missing,
        "verified_inputs": verified,
        "unverified_roles": missing,
        "interpretation": (
            "All five caller-selected inputs passed their expected SHA-256 gates."
            if retained
            else (
                "Development comparison: expected hashes were checked only for the "
                "listed roles; this report is not a retained public result."
            )
        ),
    }


def load_expected_hashes(
    manifest_path: Optional[Path], entries: Iterable[str]
) -> Tuple[Dict[str, Dict[str, object]], Optional[Dict[str, object]]]:
    expected: Dict[str, Dict[str, object]] = {}
    manifest_identity = None
    if manifest_path is not None:
        manifest_path = manifest_path.expanduser().resolve()
        require(manifest_path.is_file(), "hash manifest does not exist")
        manifest_identity = _identity(manifest_path)
        manifest = _load_json(manifest_path, "expected-hash manifest")
        require(isinstance(manifest, dict), "expected-hash manifest is not an object")
        require(
            manifest.get("schema") == HASH_MANIFEST_SCHEMA,
            "expected-hash manifest has unsupported schema",
        )
        records = manifest.get("files")
        require(isinstance(records, dict) and records, "hash manifest has no files")
        for role, record in records.items():
            require(role in INPUT_ROLES, "unknown manifest role {0!r}".format(role))
            require(isinstance(record, dict), "manifest role {0} is not an object".format(role))
            expected[role] = dict(record)

    for entry in entries:
        role, separator, digest = entry.partition("=")
        require(separator == "=" and role and digest, "expected hash must be ROLE=SHA256")
        require(role in INPUT_ROLES, "unknown expected-hash role {0!r}".format(role))
        require(FULL_SHA256.fullmatch(digest) is not None, "expected hash is not SHA-256")
        digest = digest.lower()
        if role in expected:
            require(
                str(expected[role].get("sha256", "")).lower() == digest,
                "conflicting expected hashes for {0}".format(role),
            )
        else:
            expected[role] = {"sha256": digest}
    return expected, manifest_identity


def _validate_source_identity(archive, path: Path) -> Dict[str, Dict[str, str]]:
    result = {}
    for component in ("app", "core"):
        revision = str(_scalar(archive, component + "_revision", path)).lower()
        tree_state = str(_scalar(archive, component + "_tree_state", path))
        source_kind = str(_scalar(archive, component + "_source_kind", path))
        require(
            FULL_GIT_REVISION.fullmatch(revision) is not None,
            "{0} has untraceable {1} revision".format(path.name, component),
        )
        require(
            tree_state == "clean",
            "{0} {1} tree_state is {2!r}; public comparison requires clean source".format(
                path.name, component, tree_state
            ),
        )
        require(
            source_kind == "git-checkout",
            "{0} {1} source_kind must be 'git-checkout'".format(path.name, component),
        )
        result[component] = {
            "revision": revision,
            "tree_state": tree_state,
            "source_kind": source_kind,
        }
    source_url = str(_scalar(archive, "core_source_url", path))
    require(
        source_url == PUBLIC_CORE_URL,
        "{0} has an unexpected Core source URL".format(path.name),
    )
    result["core"]["source_url"] = source_url
    return result


def _exact_parameter_mapping(
    observed, expected, description: str
) -> Dict[str, float]:
    require(
        isinstance(observed, Mapping),
        "{0} is not an object".format(description),
    )
    require(
        set(observed) == set(expected),
        "{0} parameter names disagree with Step 0 Case B".format(description),
    )
    validated = {}
    for name, expected_value in expected.items():
        validated[name] = _exact_float(
            observed[name], expected_value, "{0}.{1}".format(description, name)
        )
    return validated


def _validate_step0b_identity(archive, path: Path, *, required: bool):
    """Validate explicit Step 0B identity when present or required by GA."""
    present = BENCHMARK_IDENTITY_FIELDS.intersection(archive.files)
    if not present:
        require(
            not required,
            "source-matched generalized-alpha requires explicit Step 0 Case B "
            "benchmark identity metadata",
        )
        return None
    missing = BENCHMARK_IDENTITY_FIELDS - present
    require(
        not missing,
        "CoupFE benchmark identity metadata is incomplete; missing {0}".format(
            sorted(missing)
        ),
    )
    selected = benchmark_configuration(0, "B")
    require(
        _integer(_scalar(archive, "benchmark_step", path), "benchmark step") == 0,
        "CoupFE archive is not Benchmark 1 Step 0 Case B",
    )
    expected_text = {
        "benchmark_configuration_id": selected.identity,
        "benchmark_identity_scope": (
            "physical-mode-defined-by-benchmark-step-and-configuration-id;"
            "mpi-implementation-is-numerical-provenance"
        ),
        "benchmark_load_contract": "pressure-only",
        "benchmark_peak_load_definition": (
            "argmax(abs(active_tension_pa+pressure_pa))"
        ),
    }
    for field, expected in expected_text.items():
        require(
            str(_scalar(archive, field, path)) == expected,
            "CoupFE archive has inconsistent Step 0 Case B field {0!r}".format(
                field
            ),
        )
    require(
        not _boolean_scalar(archive, "benchmark_active_stress_enabled", path)
        and _boolean_scalar(archive, "benchmark_pressure_enabled", path),
        "CoupFE archive has inconsistent Step 0 Case B load switches",
    )
    material = _exact_parameter_mapping(
        _load_embedded_json(
            _scalar(archive, "benchmark_material_parameters_json", path),
            "CoupFE benchmark material parameters",
        ),
        selected.material_parameters,
        "CoupFE benchmark material parameters",
    )
    activation = _exact_parameter_mapping(
        _load_embedded_json(
            _scalar(archive, "benchmark_activation_parameters_json", path),
            "CoupFE benchmark activation parameters",
        ),
        selected.activation_parameters,
        "CoupFE benchmark activation parameters",
    )
    pressure = _exact_parameter_mapping(
        _load_embedded_json(
            _scalar(archive, "benchmark_pressure_parameters_json", path),
            "CoupFE benchmark pressure parameters",
        ),
        selected.pressure_parameters,
        "CoupFE benchmark pressure parameters",
    )
    runtime_manifest = _load_embedded_json(
        _scalar(archive, "benchmark_runtime_source_manifest_json", path),
        "CoupFE benchmark runtime-source manifest",
    )
    require(
        isinstance(runtime_manifest, Mapping)
        and runtime_manifest
        and set(runtime_manifest) == set(RUNTIME_SOURCE_FILES)
        and all(
            isinstance(name, str)
            and not Path(name).is_absolute()
            and not PureWindowsPath(name).is_absolute()
            and FULL_SHA256.fullmatch(str(digest)) is not None
            for name, digest in runtime_manifest.items()
        ),
        "CoupFE benchmark runtime-source manifest is malformed",
    )
    aggregate = str(
        _scalar(archive, "benchmark_runtime_source_sha256", path)
    ).lower()
    encoded = json.dumps(
        dict(runtime_manifest),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    require(
        FULL_SHA256.fullmatch(aggregate) is not None
        and hashlib.sha256(encoded).hexdigest() == aggregate,
        "CoupFE benchmark runtime-source identity is inconsistent",
    )
    return {
        "benchmark": 1,
        "step": 0,
        "case": "B",
        "configuration_id": selected.identity,
        "load_contract": selected.load_contract,
        "active_stress_enabled": False,
        "pressure_enabled": True,
        "material_parameters": material,
        "activation_parameters": activation,
        "pressure_parameters": pressure,
        "runtime_source_sha256": aggregate,
    }


def _validate_coupfe_configuration(archive, path: Path) -> Dict[str, object]:
    integrator = str(_scalar(archive, "integrator", path))
    require(
        integrator in MPI_IMPLEMENTATIONS,
        "{0} has unsupported time integrator {1!r}".format(path.name, integrator),
    )
    generalized_alpha = integrator == "generalized-alpha"
    benchmark_identity = _validate_step0b_identity(
        archive, path, required=generalized_alpha
    )
    expected_text = {
        "case": "B",
        "formulation": "hex8_standard_pointwise_kappa",
        "material_kernel_formulation": "standard",
        "material_model_id": MATERIAL_MODEL_ID,
        "mass_representation": "consistent_q1_hex8",
        "fiber_sampling": "gp_direct_rule",
        "fiber_sampling_option": "gp-direct",
        "mesh_topology": "closed_multiblock_disk",
        "point_sampling": "hex8_reference_isoparametric",
        "viscous_rate": (
            "velocity_consistent_green_lagrange_at_alpha_f_stage"
            if generalized_alpha
            else "backward_difference"
        ),
        "parameter_variant": "benchmark_eta",
        "tbar_definition": "laplace_presolved",
    }
    for key, expected in expected_text.items():
        observed = str(_scalar(archive, key, path))
        require(
            observed == expected,
            "{0} {1}={2!r}; expected {3!r}".format(path.name, key, observed, expected),
        )

    nonlinear_solver = str(_scalar(archive, "nonlinear_solver", path))
    require(
        nonlinear_solver in {"petsc-snes", "petsc-snes-mpi"},
        "{0} has unsupported nonlinear solver {1!r}".format(
            path.name, nonlinear_solver
        ),
    )

    expected_numeric = {
        "dt": EXPECTED_DT_S,
        "t_end": 1.0,
        "load_horizon": 1.0,
        "density": 1000.0,
        "material_eta_pa_s": 100.0,
        "material_kappa_pa": 1.0e6,
        "local_pressure_bulk_modulus_pa": 0.0,
        "apex_offset": 0.0,
        "perturb": 0.0,
        "a_top": 1.0e5,
        "b_top": 5.0e3,
        "a_epi": 1.0e8,
        "b_epi": 5.0e3,
    }
    for key, expected in expected_numeric.items():
        _exact_float(_scalar(archive, key, path), expected, "CoupFE " + key)

    expected_bool = {
        "isotropic": False,
        "viscous_term_active": True,
        "flip_helix": True,
    }
    for key, expected in expected_bool.items():
        observed = _scalar(archive, key, path)
        require(
            isinstance(observed, (bool, np.bool_)) and bool(observed) is expected,
            "{0} {1} is not {2}".format(path.name, key, expected),
        )

    generalized_alpha_metadata = None
    if generalized_alpha:
        for name, expected in GENERALIZED_ALPHA_PARAMETERS.items():
            _exact_float(
                _scalar(archive, "generalized_alpha_" + name, path),
                expected,
                "CoupFE generalized-alpha " + name,
            )
        require(
            str(_scalar(archive, "generalized_alpha_stage_contract", path))
            == "simula-source-matched-v1",
            "CoupFE generalized-alpha stage contract is unsupported",
        )
        load_evaluation_times = _finite_array(
            archive,
            "load_evaluation_times_s",
            path,
            shape=(EXPECTED_STEPS + 1,),
        )
        expected_load_times = COUPFE_TIME.copy()
        expected_load_times[0] = 0.0
        expected_load_times[1:] -= 0.4 * EXPECTED_DT_S
        require(
            np.allclose(
                load_evaluation_times,
                expected_load_times,
                rtol=0.0,
                atol=2.0e-15,
            ),
            "CoupFE generalized-alpha load-evaluation times are inconsistent",
        )
        generalized_alpha_metadata = {
            **GENERALIZED_ALPHA_PARAMETERS,
            "stage_contract": "simula-source-matched-v1",
            "load_time": "t_np1 - alpha_f*dt",
        }

    tbar_fields = {
        "field_filename": str(_scalar(archive, "tbar_source_filename", path)),
        "field_sha256": str(_scalar(archive, "tbar_source_sha256", path)).lower(),
        "metadata_filename": str(_scalar(archive, "tbar_metadata_filename", path)),
        "metadata_sha256": str(_scalar(archive, "tbar_metadata_sha256", path)).lower(),
        "metadata_schema": str(_scalar(archive, "tbar_metadata_schema", path)),
    }
    for key in ("field_filename", "metadata_filename"):
        tbar_fields[key] = _portable_basename(
            tbar_fields[key], "CoupFE tbar {0}".format(key)
        )
    for key in ("field_sha256", "metadata_sha256"):
        require(FULL_SHA256.fullmatch(tbar_fields[key]) is not None, "CoupFE tbar hash is invalid")
    require(
        tbar_fields["metadata_schema"] == "coupfe-cardiac-laplace-tbar-v1",
        "CoupFE tbar metadata schema is unsupported",
    )

    return {
        "case": "B",
        "benchmark_identity": benchmark_identity,
        "integrator": integrator,
        "time_integrator": (
            "source-matched generalized-alpha"
            if generalized_alpha
            else "backward Euler"
        ),
        "generalized_alpha": generalized_alpha_metadata,
        "spatial_discretization": "Q1 Hex8",
        "dt_s": EXPECTED_DT_S,
        "load_horizon_s": 1.0,
        "density_kg_m3": 1000.0,
        "material_eta_pa_s": 100.0,
        "material_kappa_pa": 1.0e6,
        "material_model_id": MATERIAL_MODEL_ID,
        "mass_representation": expected_text["mass_representation"],
        "fiber_sampling": expected_text["fiber_sampling"],
        "mesh_topology": expected_text["mesh_topology"],
        "formulation": expected_text["formulation"],
        "nonlinear_solver": nonlinear_solver,
        "tbar": tbar_fields,
    }


def _configuration_value_matches(observed, expected) -> bool:
    if isinstance(expected, bool):
        return isinstance(observed, bool) and observed is expected
    if isinstance(expected, int):
        return (
            not isinstance(observed, bool)
            and isinstance(observed, (int, float))
            and np.isfinite(float(observed))
            and float(observed).is_integer()
            and int(observed) == expected
        )
    if isinstance(expected, float):
        return (
            not isinstance(observed, bool)
            and isinstance(observed, (int, float))
            and np.isfinite(float(observed))
            and float(observed) == expected
        )
    return observed == expected


def _validate_solver_and_mpi_metadata(
    archive,
    path: Path,
    solver_name: str,
    solver_configuration: Mapping,
    diagnostics,
    *,
    integrator: str,
    n_element: int,
    degrees_of_freedom: int,
) -> Optional[Dict[str, object]]:
    require(
        isinstance(solver_configuration, Mapping),
        "CoupFE solver configuration is not an object",
    )
    require(
        solver_configuration.get("name") == solver_name,
        "CoupFE solver archive/configuration names disagree",
    )
    present_mpi = (MPI_BASE_FIELDS | MPI_MASS_FIELDS).intersection(archive.files)
    if solver_name == "petsc-snes":
        require(
            integrator == "be",
            "source-matched generalized-alpha requires the reviewed MPI solver path",
        )
        require(
            not present_mpi,
            "serial CoupFE solver archive contains MPI-only metadata",
        )
        common = {
            "snes_type": "newtonls",
            "line_search_type": "bt",
            "ksp_type": "preonly",
            "pc_type": "lu",
            "function_domain_rejection_api": "nonfinite residual for PETSc BT",
        }
        for key, expected in common.items():
            require(
                solver_configuration.get(key) == expected,
                "serial CoupFE solver configuration has unexpected {0}".format(key),
            )
        return None

    require(solver_name == "petsc-snes-mpi", "unsupported CoupFE solver")
    missing = (MPI_BASE_FIELDS | MPI_MASS_FIELDS) - present_mpi
    require(
        not missing,
        "MPI CoupFE result has incomplete MPI provenance; missing {0}".format(
            sorted(missing)
        ),
    )
    require(
        str(_scalar(archive, "driver", path))
        == "examples/cardiac_benchmark/run_mpi.py",
        "MPI CoupFE result has unexpected driver",
    )
    require(_boolean_scalar(archive, "mpi_enabled", path), "MPI execution is not enabled")
    ranks = _integer(_scalar(archive, "mpi_ranks", path), "MPI rank count", 1)
    world_size = _integer(
        _scalar(archive, "mpi_world_size", path), "MPI world size", 1
    )
    configured_ranks = _integer(
        solver_configuration.get("ranks"), "configured MPI rank count", 1
    )
    require(
        ranks == world_size == configured_ranks,
        "MPI archive/configuration rank counts disagree",
    )

    counts = np.asarray(archive["mpi_local_element_counts"])
    require(
        counts.shape == (world_size,)
        and np.issubdtype(counts.dtype, np.integer)
        and np.all(counts >= 0),
        "MPI local-element counts are invalid",
    )
    quotient, remainder = divmod(n_element, world_size)
    expected_counts = np.full(world_size, quotient, dtype=np.int64)
    expected_counts[:remainder] += 1
    require(
        np.array_equal(counts.astype(np.int64, copy=False), expected_counts),
        "MPI partition and local-element counts disagree",
    )
    expected_implementation = MPI_IMPLEMENTATIONS[integrator]
    require(
        str(_scalar(archive, "mpi_implementation", path))
        == expected_implementation,
        "MPI implementation does not match the Case B time-integration contract",
    )
    expected_fixed_configuration = dict(MPI_FIXED_CONFIGURATION)
    expected_fixed_configuration["implementation"] = expected_implementation
    for key, expected in expected_fixed_configuration.items():
        require(
            _configuration_value_matches(solver_configuration.get(key), expected),
            "MPI solver configuration has unexpected {0}".format(key),
        )
    if integrator == "generalized-alpha":
        generalized_alpha_configuration = {
            "time_integrator": "generalized-alpha",
            "material_batch_time_integrator": (
                "generalized-alpha-source-matched"
            ),
            "compiled_material_dt": EXPECTED_DT_S,
            "acceleration_stage": "1-alpha_m",
            "force_stage": "1-alpha_f",
            "material_viscous_rate": "sym(F_stage^T*grad(v_stage))",
            "nonlinear_initial_guess": "accepted-u_n-like-simula",
            "generalized_alpha": GENERALIZED_ALPHA_CONFIGURATION,
        }
        for key, expected in generalized_alpha_configuration.items():
            require(
                _configuration_value_matches(
                    solver_configuration.get(key), expected
                ),
                "MPI generalized-alpha solver configuration has unexpected "
                "{0}".format(key),
            )
    require(
        solver_configuration.get("line_search_configuration_api")
        in {"SNES.getLineSearch", "namespaced PETSc option"},
        "MPI line-search configuration API is unsupported",
    )
    for field in ("petsc4py_version", "petsc_version"):
        value = solver_configuration.get(field)
        require(
            isinstance(value, str) and bool(value.strip()),
            "MPI runtime field {0} is unavailable".format(field),
        )

    require(
        str(_scalar(archive, "mpi_partition", path)) == "coupfe.partition_elements",
        "MPI partition policy is unsupported",
    )
    require(
        str(_scalar(archive, "mpi_build_layout", path)) == "isolated-rank-directories",
        "MPI build layout is unsupported",
    )
    factor = str(_scalar(archive, "mpi_factor_solver_type", path))
    require(
        factor == "superlu_dist"
        and factor == solver_configuration.get("factor_solver_type"),
        "MPI factor-solver metadata disagrees",
    )

    require(
        str(_scalar(archive, "mpi_mass_partition", path))
        == "owned-row-csr-all-touching-elements",
        "MPI consistent-mass partition policy is unsupported",
    )
    row_ranges = np.asarray(archive["mpi_mass_owned_row_ranges"])
    local_nnz = np.asarray(archive["mpi_mass_local_nnz"])
    touching = np.asarray(archive["mpi_mass_touching_element_counts"])
    require(
        row_ranges.shape == (world_size, 2)
        and local_nnz.shape == (world_size,)
        and touching.shape == (world_size,)
        and np.issubdtype(row_ranges.dtype, np.integer)
        and np.issubdtype(local_nnz.dtype, np.integer)
        and np.issubdtype(touching.dtype, np.integer)
        and np.all(local_nnz > 0)
        and np.all(touching > 0),
        "MPI consistent-mass ownership arrays are invalid",
    )
    row_ranges = row_ranges.astype(np.int64, copy=False)
    expected_start = 0
    for start, stop in row_ranges:
        require(
            int(start) == expected_start and int(stop) >= int(start),
            "MPI mass row ranges are not contiguous",
        )
        expected_start = int(stop)
    require(
        expected_start == degrees_of_freedom,
        "MPI mass rows do not cover the global system",
    )
    require(
        solver_configuration.get("mass_representation") == "consistent_q1_hex8"
        and solver_configuration.get("mass_partition")
        == "owned-row-csr-all-touching-elements"
        and solver_configuration.get("mass_owned_row_range")
        == row_ranges[0].astype(int).tolist()
        and _integer(
            solver_configuration.get("mass_local_nnz"),
            "configured local mass nnz",
            0,
        )
        == int(local_nnz[0]),
        "MPI mass archive/configuration metadata disagrees",
    )

    evaluation_mode = str(_scalar(archive, "element_evaluation_mode", path))
    residual_only = _boolean_scalar(
        archive, "compiled_material_residual_only_available", path
    )
    require(
        evaluation_mode in {"joint", "split"}
        and solver_configuration.get("element_evaluation_mode") == evaluation_mode,
        "MPI element-evaluation metadata disagrees",
    )
    require(
        isinstance(
            solver_configuration.get("compiled_material_residual_only_available"),
            bool,
        )
        and solver_configuration["compiled_material_residual_only_available"]
        is residual_only,
        "MPI residual-only capability metadata disagrees",
    )
    for index, record in enumerate(diagnostics, start=1):
        require(
            "ranks" in record
            and _integer(record["ranks"], "diagnostic MPI ranks", 1) == world_size,
            "MPI diagnostic rank count disagrees at step {0}".format(index),
        )

    return {
        "implementation": expected_implementation,
        "world_size": world_size,
        "local_element_counts": counts.astype(int).tolist(),
        "partition": "coupfe.partition_elements",
        "build_layout": "isolated-rank-directories",
        "factor_solver_type": factor,
        "mass_partition": {
            "policy": "owned-row-csr-all-touching-elements",
            "owned_row_ranges": row_ranges.astype(int).tolist(),
            "local_nnz": local_nnz.astype(int).tolist(),
            "touching_element_counts": touching.astype(int).tolist(),
        },
        "element_evaluation_mode": evaluation_mode,
        "compiled_material_residual_only_available": residual_only,
    }


def load_coupfe_run(path: Path, identity: Mapping[str, object]) -> Dict[str, object]:
    path = path.expanduser().resolve()
    try:
        archive_context = np.load(path, allow_pickle=False)
    except (OSError, ValueError) as error:
        raise ComparisonInputError(
            "cannot load CoupFE archive {0} without pickle: {1}".format(path.name, error)
        ) from error

    with archive_context as archive:
        require(
            str(_scalar(archive, "result_schema", path)) == RESULT_SCHEMA,
            "CoupFE archive has unsupported result schema",
        )
        converged = _scalar(archive, "converged", path)
        require(
            isinstance(converged, (bool, np.bool_)) and bool(converged),
            "CoupFE archive is not marked converged",
        )
        completed = int(_scalar(archive, "completed_steps", path))
        expected = int(_scalar(archive, "expected_steps", path))
        require(
            completed == expected == EXPECTED_STEPS,
            "CoupFE archive is incomplete: {0}/{1}; required {2}/{2}".format(
                completed, expected, EXPECTED_STEPS
            ),
        )
        times = _finite_array(archive, "times", path, shape=(EXPECTED_STEPS + 1,))
        require(
            np.allclose(times, COUPFE_TIME, rtol=0.0, atol=2.0e-13),
            "CoupFE time grid is not the required 0--1 s, 1 ms grid",
        )
        configuration = _validate_coupfe_configuration(archive, path)
        integrator = configuration["integrator"]
        histories = {
            "p0": _finite_array(archive, "u0", path, shape=(len(times), 3)),
            "p1": _finite_array(archive, "u1", path, shape=(len(times), 3)),
        }
        pressure = _finite_array(archive, "pres", path, shape=(len(times),))
        activation = _finite_array(archive, "tau", path, shape=(len(times),))
        if integrator == "generalized-alpha":
            load_times = COUPFE_TIME.copy()
            load_times[0] = 0.0
            load_times[1:] -= 0.4 * EXPECTED_DT_S
            expected_pressure = np.zeros_like(COUPFE_TIME)
            expected_pressure[1:] = p_of_t(
                load_times[1:], t_span=(0.0, 1.0)
            )
        else:
            expected_pressure = np.asarray(p_of_t(COUPFE_TIME), dtype=float)
        require(
            np.array_equal(pressure, expected_pressure),
            "CoupFE Case B pressure is not the fixed benchmark schedule at "
            "the selected load stage",
        )
        require(np.array_equal(activation, np.zeros_like(activation)), "CoupFE Case B activation is not zero")

        for point in POINTS:
            landmark = _finite_array(archive, point, path, shape=(3,))
            require(
                np.allclose(landmark, EXPECTED_LANDMARKS_M[point], rtol=0.0, atol=1.0e-15),
                "CoupFE {0} landmark coordinate differs from the benchmark".format(point),
            )
            error = _exact_float(
                _scalar(archive, point + "_sampling_reconstruction_error_m", path),
                float(_scalar(archive, point + "_sampling_reconstruction_error_m", path)),
                "CoupFE " + point + " reconstruction error",
            )
            require(error <= 1.0e-9, "CoupFE landmark reconstruction error is too large")

        source_identity = _validate_source_identity(archive, path)

        nodes = _finite_array(archive, "nodes", path)
        require(nodes.ndim == 2 and nodes.shape[1] == 3, "CoupFE nodes are malformed")
        require("elems" in archive, "CoupFE archive is missing elems")
        elements = np.asarray(archive["elems"])
        require(
            elements.ndim == 2
            and elements.shape[1] == 8
            and np.issubdtype(elements.dtype, np.integer)
            and np.all(elements >= 0)
            and np.all(elements < len(nodes)),
            "CoupFE Hex8 connectivity is malformed",
        )
        det_f = _finite_array(
            archive, "det_f_gauss_peak", path, shape=(len(elements), 8)
        )
        require(np.all(det_f > 0.0), "CoupFE peak field has nonpositive det(F)")

        pre_solve = _load_embedded_json(
            _scalar(archive, "pre_solve_audit_json", path), "CoupFE pre-solve audit"
        )
        require(isinstance(pre_solve, dict), "CoupFE pre-solve audit is not an object")
        for name in ("geometry", "pressure", "robin"):
            record = pre_solve.get(name)
            require(
                isinstance(record, dict) and record.get("passed") is True,
                "CoupFE {0} pre-solve audit did not pass".format(name),
            )
        require(
            pre_solve["geometry"].get("mesh_topology") == "closed_multiblock_disk"
            and pre_solve["geometry"].get("unclassified_exterior_faces") == 0
            and pre_solve["geometry"].get("nonpositive_extended_jacobians") == 0,
            "CoupFE closed-geometry audit metadata is inconsistent",
        )

        diagnostics = _load_embedded_json(
            _scalar(archive, "nonlinear_step_diagnostics_json", path),
            "CoupFE nonlinear diagnostics",
        )
        solver_configuration = _load_embedded_json(
            _scalar(archive, "solver_configuration_json", path),
            "CoupFE solver configuration",
        )
        require(
            isinstance(diagnostics, list) and len(diagnostics) == EXPECTED_STEPS,
            "CoupFE archive does not retain 1000 step diagnostics",
        )
        iteration_counts = []
        residual_ratios = []
        for index, record in enumerate(diagnostics, start=1):
            require(isinstance(record, dict), "malformed diagnostic at step {0}".format(index))
            _exact_float(record.get("time"), times[index], "diagnostic time")
            require(int(record.get("snes_converged_reason", 0)) > 0, "non-converged diagnostic")
            final_norm = float(record.get("final_residual_norm", np.nan))
            threshold = float(record.get("residual_acceptance_threshold", np.nan))
            require(
                np.isfinite(final_norm)
                and np.isfinite(threshold)
                and threshold > 0.0
                and final_norm <= threshold * (1.0 + 1.0e-9) + 1.0e-14,
                "unaccepted residual at diagnostic step {0}".format(index),
            )
            iteration_counts.append(int(record.get("nonlinear_iterations", -1)))
            residual_ratios.append(final_norm / threshold)
        require(min(iteration_counts) >= 0, "negative nonlinear iteration count")
        mpi_metadata = _validate_solver_and_mpi_metadata(
            archive,
            path,
            configuration["nonlinear_solver"],
            solver_configuration,
            diagnostics,
            integrator=integrator,
            n_element=len(elements),
            degrees_of_freedom=3 * len(nodes),
        )

        configuration.update(
            {
                "mesh": {
                    "nodes": int(len(nodes)),
                    "elements": int(len(elements)),
                    "degrees_of_freedom": int(3 * len(nodes)),
                    "n_t": int(_scalar(archive, "n_t", path)),
                    "n_core": int(_scalar(archive, "n_core", path)),
                    "n_radial": int(_scalar(archive, "n_radial", path)),
                    "core_half_width": float(_scalar(archive, "core_half_width", path)),
                },
                "solver_evidence": {
                    "name": configuration["nonlinear_solver"],
                    "diagnostic_steps": len(diagnostics),
                    "nonlinear_iterations_min": int(min(iteration_counts)),
                    "nonlinear_iterations_max": int(max(iteration_counts)),
                    "maximum_final_residual_fraction": float(max(residual_ratios)),
                    "peak_det_f_minimum": float(np.min(det_f)),
                    "peak_det_f_maximum": float(np.max(det_f)),
                },
                "mpi": mpi_metadata,
            }
        )

    return {
        "identity": dict(identity),
        "times": times,
        "histories": histories,
        "configuration": configuration,
        "source_identity": source_identity,
    }


def _load_embedded_json(raw, description: str):
    try:
        value = json.loads(str(raw), parse_constant=lambda item: (_ for _ in ()).throw(ValueError(item)))
    except (TypeError, ValueError, json.JSONDecodeError) as error:
        raise ComparisonInputError("malformed {0}".format(description)) from error
    _finite_json(value, description)
    return value


def _validate_fenics_parameters(parameters: Mapping) -> Dict[str, object]:
    require(isinstance(parameters, Mapping), "FEniCS parameters are not an object")
    require(parameters.get("benchmark") == 1, "FEniCS output is not benchmark 1")
    require(parameters.get("step") == 0, "FEniCS output is not step 0")
    require(str(parameters.get("case", "")).lower() == "b", "FEniCS output is not Case B")
    require(parameters.get("zero_activation") is True, "FEniCS activation is not zero")
    require(parameters.get("zero_pressure") is False, "FEniCS pressure is disabled")
    require(
        _nested(parameters, ("problem_parameters", "function_space"), "FEniCS parameters")
        == "P_2",
        "FEniCS displacement space is not P_2",
    )
    require(
        _nested(parameters, ("fiber_parameters", "function_space"), "FEniCS parameters")
        == "P_2",
        "FEniCS fiber space is not P_2",
    )
    for keys, expected in EXPECTED_FENICS_NUMERIC.items():
        value = _nested(parameters, keys, "FEniCS parameters")
        _exact_float(value, expected, "FEniCS " + ".".join(keys))

    geometry_filename = _sanitized_basename(
        str(parameters.get("geometry_path", ""))
    )
    require(bool(geometry_filename), "FEniCS geometry_path has no basename")

    producer_timestamp = str(parameters.get("timestamp", ""))
    require(
        NAIVE_ISO_TIMESTAMP.fullmatch(producer_timestamp) is not None,
        "FEniCS timestamp is not the recorded naive ISO-8601 form",
    )
    try:
        parsed_timestamp = datetime.fromisoformat(producer_timestamp)
    except ValueError as error:
        raise ComparisonInputError(
            "FEniCS timestamp is not a valid ISO-8601 value"
        ) from error
    require(
        parsed_timestamp.tzinfo is None,
        "FEniCS timestamp unexpectedly records a timezone",
    )

    # This is an allowlisted, path-sanitized record. In particular, caller
    # outdir/outpath strings are intentionally absent.
    return {
        "benchmark": 1,
        "step": 0,
        "case": "B",
        "zero_activation": True,
        "zero_pressure": False,
        "time_integrator": "generalized-alpha",
        "generalized_alpha": {"alpha_f": 0.4, "alpha_m": 0.2},
        "spatial_discretization": "P2 tetrahedra",
        "fiber_function_space": "P_2",
        "dt_s": EXPECTED_DT_S,
        "density_kg_m3": 1000.0,
        "material_eta_pa_s": 100.0,
        "material_kappa_pa": 1.0e6,
        "geometry_filename": geometry_filename,
        "producer_recorded_timestamp": producer_timestamp,
        "producer_timestamp_timezone": "not-recorded-in-source",
    }


def _load_npy(path: Path, description: str, shape) -> np.ndarray:
    try:
        value = np.load(path, allow_pickle=False)
    except (OSError, ValueError) as error:
        raise ComparisonInputError(
            "cannot load {0} {1} without pickle: {2}".format(description, path.name, error)
        ) from error
    try:
        value = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as error:
        raise ComparisonInputError("{0} is not numeric".format(description)) from error
    require(value.shape == shape, "{0} has shape {1}; expected {2}".format(description, value.shape, shape))
    require(np.all(np.isfinite(value)), "{0} contains non-finite values".format(description))
    return value.copy()


def load_fenics_reference(
    parameters_path: Path,
    times_path: Path,
    p0_path: Path,
    p1_path: Path,
    identities: Mapping[str, Mapping[str, object]],
) -> Dict[str, object]:
    parameters_path = parameters_path.expanduser().resolve()
    times_path = times_path.expanduser().resolve()
    p0_path = p0_path.expanduser().resolve()
    p1_path = p1_path.expanduser().resolve()
    parameters = _load_json(parameters_path, "FEniCS parameters")
    configuration = _validate_fenics_parameters(parameters)
    times = _load_npy(times_path, "FEniCS time array", (FENICS_SAMPLES,))
    p0 = _load_npy(p0_path, "FEniCS p0 displacement", (FENICS_SAMPLES, 3))
    p1 = _load_npy(p1_path, "FEniCS p1 displacement", (FENICS_SAMPLES, 3))
    require(
        np.allclose(times, FENICS_TIME, rtol=0.0, atol=2.0e-13),
        "FEniCS time grid is not the retained 0.001--0.999 s convention",
    )
    return {
        "identities": {
            role: dict(identities[role])
            for role in ("fenics-parameters", "fenics-times", "fenics-p0", "fenics-p1")
        },
        "times": times,
        "histories": {"p0": p0, "p1": p1},
        "configuration": configuration,
    }


def map_to_common_grid(
    source_times: np.ndarray,
    source_history: np.ndarray,
    target_times: np.ndarray,
) -> Tuple[np.ndarray, bool, str]:
    source_times = np.asarray(source_times, dtype=float)
    source_history = np.asarray(source_history, dtype=float)
    target_times = np.asarray(target_times, dtype=float)
    require(
        source_times.ndim == target_times.ndim == 1
        and len(source_times) >= 2
        and np.all(np.diff(source_times) > 0.0)
        and np.all(np.diff(target_times) > 0.0),
        "invalid comparison time grid",
    )
    require(
        source_history.shape == (len(source_times), 3)
        and np.all(np.isfinite(source_history)),
        "invalid comparison displacement history",
    )
    require(
        target_times[0] >= source_times[0] and target_times[-1] <= source_times[-1],
        "target grid requires endpoint extrapolation",
    )
    insertion = np.searchsorted(source_times, target_times)
    right = np.clip(insertion, 0, len(source_times) - 1)
    left = np.clip(insertion - 1, 0, len(source_times) - 1)
    choose_left = (
        np.abs(source_times[left] - target_times)
        <= np.abs(source_times[right] - target_times)
    )
    positions = np.where(choose_left, left, right)
    identity = bool(
        np.all(np.diff(positions) > 0)
        and np.allclose(
            source_times[positions], target_times, rtol=0.0, atol=2.0e-13
        )
    )
    if identity:
        return source_history[positions].copy(), True, "identity_index_selection"
    mapped = np.column_stack(
        [np.interp(target_times, source_times, source_history[:, component]) for component in range(3)]
    )
    return mapped, False, "linear_interpolation_without_extrapolation"


def onset_time(times: np.ndarray, history: np.ndarray) -> Optional[float]:
    """First downward p0/p1 u_z=-5 mm crossing, linearly interpolated."""
    times = np.asarray(times, dtype=float)
    history = np.asarray(history, dtype=float)
    require(history.shape == (len(times), 3), "onset history shape mismatch")
    uz = history[:, 2]
    indices = np.flatnonzero(uz <= SNAP_ONSET_THRESHOLD_M)
    if not len(indices):
        return None
    index = int(indices[0])
    if index == 0 or uz[index] == uz[index - 1]:
        return float(times[index])
    fraction = (SNAP_ONSET_THRESHOLD_M - uz[index - 1]) / (uz[index] - uz[index - 1])
    return float(times[index - 1] + fraction * (times[index] - times[index - 1]))


def _phase_masks(times: np.ndarray) -> Dict[str, np.ndarray]:
    return {
        "full_0p001_to_0p999_s": np.ones(len(times), dtype=bool),
        "before_snap_t_lt_0p20_s": times < SNAP_START_S,
        "snap_0p20_to_0p32_s_inclusive": (times >= SNAP_START_S) & (times <= SNAP_END_S),
        "after_snap_t_gt_0p32_s": times > SNAP_END_S,
    }


def vector_metrics(
    times: np.ndarray, candidate: np.ndarray, reference: np.ndarray
) -> Dict[str, object]:
    times = np.asarray(times, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    reference = np.asarray(reference, dtype=float)
    require(candidate.shape == reference.shape == (len(times), 3), "metric shape mismatch")
    error = candidate - reference
    result: Dict[str, object] = {"phases": {}}
    for phase, mask in _phase_masks(times).items():
        require(np.any(mask), "empty metric phase {0}".format(phase))
        phase_error = error[mask]
        phase_reference = reference[mask]
        vector_error = np.linalg.norm(phase_error, axis=1)
        denominator = float(np.linalg.norm(phase_reference))
        require(denominator > 0.0, "reference norm is zero in phase {0}".format(phase))
        maximum_index = int(np.argmax(vector_error))
        result["phases"][phase] = {
            "samples": int(np.count_nonzero(mask)),
            "vector_rmse_mm": float(1.0e3 * np.sqrt(np.mean(vector_error ** 2))),
            "relative_l2": float(np.linalg.norm(phase_error) / denominator),
            "component_rmse_mm": {
                component: float(1.0e3 * np.sqrt(np.mean(phase_error[:, index] ** 2)))
                for index, component in enumerate(COMPONENTS)
            },
            "maximum_vector_error_mm": float(1.0e3 * vector_error[maximum_index]),
            "maximum_vector_error_time_s": float(times[mask][maximum_index]),
        }
    candidate_onset = onset_time(times, candidate)
    reference_onset = onset_time(times, reference)
    result["snap_onset"] = {
        "definition": (
            "first downward crossing of landmark u_z=-5 mm, with linear "
            "interpolation between adjacent retained samples"
        ),
        "threshold_mm": -5.0,
        "candidate_time_s": candidate_onset,
        "reference_time_s": reference_onset,
        "candidate_minus_reference_s": (
            None
            if candidate_onset is None or reference_onset is None
            else float(candidate_onset - reference_onset)
        ),
    }
    return result


def build_report(
    coupfe: Mapping[str, object],
    fenics: Mapping[str, object],
    expected_hash_manifest: Optional[Mapping[str, object]],
    hash_verification: Optional[Mapping[str, object]] = None,
) -> Tuple[Dict[str, object], Dict[str, object]]:
    if hash_verification is None:
        hash_verification = {
            "mode": "development",
            "all_input_roles_required": False,
            "all_input_roles_verified": False,
            "verified_inputs": {},
            "unverified_roles": list(INPUT_ROLES),
            "interpretation": (
                "Development comparison: no caller-expected hashes were supplied; "
                "this report is not a retained public result."
            ),
        }
    retained = hash_verification.get("mode") == "retained"
    require(
        hash_verification.get("mode") in {"development", "retained"},
        "comparison hash-verification mode is invalid",
    )
    require(
        not retained or hash_verification.get("all_input_roles_verified") is True,
        "retained report requires all five input hashes to be verified",
    )
    coupfe_integrator = str(coupfe["configuration"]["integrator"])
    coupfe_method_label = (
        "CoupFE Q1-Hex8, source-matched generalized-alpha"
        if coupfe_integrator == "generalized-alpha"
        else "CoupFE Q1-Hex8, backward Euler"
    )
    common_times = np.asarray(fenics["times"], dtype=float)
    mapped = {}
    mappings = {}
    comparisons = {}
    for point in POINTS:
        mapped_history, identity, method = map_to_common_grid(
            np.asarray(coupfe["times"]),
            np.asarray(coupfe["histories"][point]),
            common_times,
        )
        mapped[point] = mapped_history
        mappings[point] = {"identity_sampling": identity, "method": method}
        comparisons[point] = vector_metrics(
            common_times, mapped_history, np.asarray(fenics["histories"][point])
        )

    report = {
        "schema": REPORT_SCHEMA,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "evidence_status": {
            "classification": (
                "clean-source-hash-gated-landmark-comparison"
                if retained
                else "development-landmark-comparison"
            ),
            "public_retained_candidate": retained,
            "interpretation": (
                (
                    "Both CoupFE source identities are clean and traceable, and all five "
                    "input files passed caller-expected SHA-256 gates. "
                )
                if retained
                else (
                    "Both CoupFE source identities are clean and traceable, but this "
                    "development report did not require expected hashes for all inputs. "
                )
            )
            + (
                "The FEniCS producer revision was not retained contemporaneously."
            ),
        },
        "bounded_claim": (
            "This report compares p0/p1 displacement trajectories for one completed "
            "CoupFE Case B run with four hash-identified local FEniCS output files. "
            "It does not establish stress agreement, mesh/time convergence, equivalence "
            "of the named Q1-Hex8 and P2-tetrahedron discretizations, clinical "
            "validity, or broad solver validation. The "
            "producing FEniCS source revision was not retained contemporaneously; a "
            "nearby clone can provide context but is not producer identity."
        ),
        "excluded_inputs": (
            "No pressure_model.npy, pickle, HDF5, stress field, or solver executable "
            "is read by this comparison."
        ),
        "metric_definitions": {
            "vector_rmse": "sqrt(mean_t(||u_CoupFE(t)-u_FEniCS(t)||_2^2))",
            "relative_l2": "||u_CoupFE-u_FEniCS||_F / ||u_FEniCS||_F",
            "snap_window_s": [SNAP_START_S, SNAP_END_S],
            "snap_onset": (
                "first downward landmark u_z=-5 mm crossing, linearly interpolated "
                "between adjacent samples; threshold is fixed by the method"
            ),
        },
        "inputs": {
            "coupfe_run": dict(coupfe["identity"]),
            "fenics": dict(fenics["identities"]),
            "expected_hash_manifest": (
                None if expected_hash_manifest is None else dict(expected_hash_manifest)
            ),
            "hash_verification": dict(hash_verification),
            "path_policy": (
                "Only basenames, byte counts, and SHA-256 identities are retained; "
                "resolved caller paths and parameters.json outpath/outdir values are omitted."
            ),
        },
        "configuration": {
            "coupfe": dict(coupfe["configuration"]),
            "fenics": dict(fenics["configuration"]),
            "coupfe_source_identity": dict(coupfe["source_identity"]),
            "fenics_source_identity": {
                "producing_revision": "not-retained-contemporaneously",
                "interpretation": (
                    "The parameters and output arrays have exact file hashes. A nearby "
                    "clean source clone is contextual and is not claimed as their producer."
                ),
            },
        },
        "common_grid": {
            "convention": (
                "FEniCS retained-output convention: 999 samples from 0.001 through "
                "0.999 s at 0.001 s intervals"
            ),
            "samples": int(len(common_times)),
            "start_s": float(common_times[0]),
            "end_s": float(common_times[-1]),
            "increment_s": EXPECTED_DT_S,
            "coupfe_mapping": mappings,
            "interpretation": (
                "For the required CoupFE dt=0.001 s archive, all common-grid samples "
                "coincide exactly; no interpolation is used."
            ),
        },
        "comparison": comparisons,
    }
    _finite_json(report, "comparison report")
    figure_data = {
        "times": common_times,
        "coupfe": mapped,
        "fenics": fenics["histories"],
        "comparison": comparisons,
        "title": (
            "Case B landmark displacement comparison"
            if retained
            else "DEVELOPMENT — Case B landmark displacement comparison"
        ),
        "retention_note": (
            "All five input identities passed caller-expected SHA-256 gates."
            if retained
            else "Development output—not a retained public result."
        ),
        "source_labels": {
            "coupfe": coupfe_method_label,
            "fenics": "FEniCS P2-tet, generalized-alpha",
        },
    }
    return report, figure_data


def _write_json_temporary(path: Path, report: Mapping) -> Path:
    payload = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    stream = tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        prefix="." + path.name + ".",
        suffix=".tmp",
        dir=str(path.parent),
        delete=False,
    )
    temporary = Path(stream.name)
    try:
        with stream:
            stream.write(payload)
        return temporary
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _write_figure_temporary(path: Path, data: Mapping[str, object]) -> Path:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise ComparisonInputError(
            "figure generation requires matplotlib; install .[reference]"
        ) from error

    path.parent.mkdir(parents=True, exist_ok=True)
    colors = {
        "fenics": "#173F5F",
        "coupfe": "#D97706",
        "snap": "#F2C94C",
        "grid": "#D1D5DB",
    }
    fig, axes = plt.subplots(2, 3, figsize=(13.8, 7.6), sharex=True)
    for row, point in enumerate(POINTS):
        for column, component in enumerate(COMPONENTS):
            axis = axes[row, column]
            axis.axvspan(
                SNAP_START_S,
                SNAP_END_S,
                color=colors["snap"],
                alpha=0.13,
                linewidth=0.0,
                label="Declared snap window" if (row, column) == (0, 0) else None,
            )
            axis.plot(
                data["times"],
                1.0e3 * np.asarray(data["fenics"][point])[:, column],
                color=colors["fenics"],
                linewidth=2.1,
                label=data["source_labels"]["fenics"] if (row, column) == (0, 0) else None,
            )
            axis.plot(
                data["times"],
                1.0e3 * np.asarray(data["coupfe"][point])[:, column],
                color=colors["coupfe"],
                linewidth=1.9,
                linestyle="--",
                label=data["source_labels"]["coupfe"] if (row, column) == (0, 0) else None,
            )
            axis.axhline(0.0, color="#9CA3AF", linewidth=0.65, zorder=0)
            axis.set_title("{0}: {1}-component".format(point, component), fontsize=10.5)
            axis.set_ylabel("Displacement (mm)")
            if row == 1:
                axis.set_xlabel("Time (s)")
            axis.grid(color=colors["grid"], alpha=0.55, linewidth=0.55)
            axis.tick_params(labelsize=8.5)
            if column == 2:
                full = data["comparison"][point]["phases"][
                    "full_0p001_to_0p999_s"
                ]
                axis.text(
                    0.975,
                    0.055,
                    "vector RMSE {0:.2f} mm\nrelative L2 {1:.1%}".format(
                        full["vector_rmse_mm"], full["relative_l2"]
                    ),
                    transform=axis.transAxes,
                    ha="right",
                    va="bottom",
                    fontsize=8.5,
                    color="#374151",
                    bbox={
                        "boxstyle": "square,pad=0.35",
                        "facecolor": "white",
                        "edgecolor": "#D1D5DB",
                        "alpha": 0.88,
                    },
                )

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.935),
        ncol=3,
        frameon=False,
        fontsize=9,
    )
    fig.suptitle(
        data["title"],
        fontsize=15.5,
        y=0.992,
    )
    fig.text(
        0.5,
        0.018,
        (
            "Common grid: retained FEniCS 0.001–0.999 s samples; CoupFE samples "
            "coincide at dt=1 ms. Landmark comparison only—no stress or convergence claim.\n"
            + data["retention_note"]
        ),
        ha="center",
        fontsize=8.6,
        color="#4B5563",
    )
    fig.tight_layout(rect=(0.035, 0.072, 0.995, 0.875))

    suffix = path.suffix.lower().lstrip(".") or "png"
    stream = tempfile.NamedTemporaryFile(
        prefix="." + path.name + ".",
        suffix=".tmp",
        dir=str(path.parent),
        delete=False,
    )
    temporary = Path(stream.name)
    stream.close()
    try:
        fig.savefig(
            temporary,
            format=suffix,
            dpi=190,
            bbox_inches="tight",
            metadata={
                "Title": data["title"],
                "Description": (
                    "Hash-identified FEniCS and CoupFE p0/p1 trajectories; no filesystem paths"
                ),
            },
        )
        return temporary
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    finally:
        plt.close(fig)


def write_outputs_atomic(
    report_path: Path,
    figure_path: Path,
    report: Mapping[str, object],
    figure_data: Mapping[str, object],
) -> None:
    report_path = report_path.expanduser().resolve()
    figure_path = figure_path.expanduser().resolve()
    require(report_path != figure_path, "report and figure paths must differ")
    _finite_json(report, "comparison report")
    report_temp = None
    figure_temp = None
    try:
        report_temp = _write_json_temporary(report_path, report)
        figure_temp = _write_figure_temporary(figure_path, figure_data)
        os.replace(figure_temp, figure_path)
        figure_temp = None
        os.replace(report_temp, report_path)
        report_temp = None
    finally:
        if report_temp is not None:
            report_temp.unlink(missing_ok=True)
        if figure_temp is not None:
            figure_temp.unlink(missing_ok=True)


def _input_paths(args) -> Dict[str, Path]:
    return {
        "coupfe-run": args.coupfe_run.expanduser().resolve(),
        "fenics-parameters": args.fenics_parameters.expanduser().resolve(),
        "fenics-times": args.fenics_times.expanduser().resolve(),
        "fenics-p0": args.fenics_p0.expanduser().resolve(),
        "fenics-p1": args.fenics_p1.expanduser().resolve(),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coupfe-run", type=Path, required=True)
    parser.add_argument("--fenics-parameters", type=Path, required=True)
    parser.add_argument("--fenics-times", type=Path, required=True)
    parser.add_argument("--fenics-p0", type=Path, required=True)
    parser.add_argument("--fenics-p1", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--figure", type=Path, required=True)
    parser.add_argument(
        "--expected-hashes",
        type=Path,
        help="optional JSON manifest of expected role/basename/SHA-256 identities",
    )
    parser.add_argument(
        "--expect-sha256",
        action="append",
        default=[],
        metavar="ROLE=SHA256",
        help="optional expected input hash; may be repeated",
    )
    parser.add_argument(
        "--retained",
        action="store_true",
        help=(
            "create a retained-result candidate; requires expected SHA-256 values "
            "for all five input roles"
        ),
    )
    args = parser.parse_args(argv)

    paths = _input_paths(args)
    identities = {role: _identity(path) for role, path in paths.items()}
    expected, manifest_identity = load_expected_hashes(
        args.expected_hashes, args.expect_sha256
    )
    _validate_expected_hashes(identities, expected)
    hash_verification = _hash_verification_record(
        identities, expected, retained=args.retained
    )
    for role, identity in identities.items():
        identity["expected_sha256_checked"] = role in expected

    coupfe = load_coupfe_run(paths["coupfe-run"], identities["coupfe-run"])
    fenics = load_fenics_reference(
        paths["fenics-parameters"],
        paths["fenics-times"],
        paths["fenics-p0"],
        paths["fenics-p1"],
        identities,
    )
    report, figure_data = build_report(
        coupfe, fenics, manifest_identity, hash_verification
    )
    write_outputs_atomic(args.report, args.figure, report, figure_data)

    for point in POINTS:
        full = report["comparison"][point]["phases"]["full_0p001_to_0p999_s"]
        onset = report["comparison"][point]["snap_onset"]
        print(
            "{0}: RMSE={1:.6f} mm; relative-L2={2:.6g}; "
            "onset CoupFE/FEniCS={3}/{4} s".format(
                point,
                full["vector_rmse_mm"],
                full["relative_l2"],
                onset["candidate_time_s"],
                onset["reference_time_s"],
            )
        )
    print("saved report -> {0}".format(args.report.name))
    print("saved figure -> {0}".format(args.figure.name))
    return report


if __name__ == "__main__":
    main()
