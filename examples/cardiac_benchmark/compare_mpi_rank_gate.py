#!/usr/bin/env python3
"""Create a fail-closed serial/MPI rank-equivalence record for closed Case B.

The utility compares four caller-selected completed archives: the serial
``petsc-snes`` result and the matching 1-, 2-, and 4-rank
``petsc-snes-mpi`` results.  It accepts only the fixed closed Case B method
used by the public benchmark and writes a report only when every provenance,
configuration, completion, audit, exact-invariant, and numerical-equivalence
gate passes.

Filesystem paths are never copied into the report.  Inputs are identified by
portable basename, byte count, and SHA-256 only.  NumPy archives are always
opened with ``allow_pickle=False``.
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
from typing import Dict, Mapping, Optional, Sequence, Tuple

import numpy as np


REPORT_SCHEMA = "coupfe-cardiac-case-b-mpi-rank-gate-v1"
RESULT_SCHEMA = "coupfe-cardiac-result-v1"
MATERIAL_MODEL_ID = (
    "holzapfel-ogden-smooth-switch-complete-energy-derivative-v1"
)
MPI_IMPLEMENTATION = "cardiac-owned-distributed-closed-std-kappa-step0"
PUBLIC_CORE_URL = "https://github.com/tengzhang48/CoupFE.git"
TBAR_SCHEMA = "coupfe-cardiac-laplace-tbar-v1"

EXPECTED_DT_S = 1.0e-3
EXPECTED_T_END_S = 0.32
EXPECTED_LOAD_HORIZON_S = 1.0
EXPECTED_STEPS = 320
EXPECTED_TIMES = np.linspace(0.0, EXPECTED_T_END_S, EXPECTED_STEPS + 1)
EXPECTED_N_T = 2
EXPECTED_N_CORE = 20
EXPECTED_N_RADIAL = 17
EXPECTED_CORE_HALF_WIDTH = 0.36
EXPECTED_MPI_RANKS = {"mpi1": 1, "mpi2": 2, "mpi4": 4}
INPUT_ROLES = ("serial", "mpi1", "mpi2", "mpi4")

# This is the precedent used by the historical serial/MPI gate.  These values
# are deliberately not command-line options: a failed gate cannot be made to
# pass by retuning its tolerance.
NUMERICAL_RTOL = 2.0e-11
NUMERICAL_ATOL = 2.0e-13
TOLERANT_FIELDS = ("u0", "u1", "U_peak", "det_f_gauss_peak")
EXACT_FLOAT_ARRAY_FIELDS = (
    "times",
    "tau",
    "pres",
    "nodes",
    "p0",
    "p1",
    "fiber",
    "p0_sampling_natural",
    "p0_sampling_weights",
    "p1_sampling_natural",
    "p1_sampling_weights",
)
EXACT_INTEGER_ARRAY_FIELDS = ("elems", "facets_endo")

# These common text fields describe intentionally different execution paths.
# Every other common integer, boolean, or text field is an invariant and is
# compared exactly.
RUN_DEPENDENT_COMMON_FIELDS = frozenset(
    {
        "driver",
        "nonlinear_solver",
        "solver_configuration_json",
        "nonlinear_step_diagnostics_json",
    }
)

FULL_SHA256 = re.compile(r"[0-9a-fA-F]{64}")
FULL_GIT_REVISION = re.compile(r"[0-9a-fA-F]{40}")


class RankGateInputError(RuntimeError):
    """An input is unsafe, incomplete, or outside the declared rank gate."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RankGateInputError(message)


def _portable_basename(value: str, description: str) -> str:
    value = str(value)
    require(
        bool(value)
        and value == Path(value).name == PureWindowsPath(value).name
        and not Path(value).is_absolute()
        and not PureWindowsPath(value).is_absolute(),
        "{0} is not a portable basename".format(description),
    )
    return value


def _sha256_file(path: Path) -> Tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    try:
        with path.open("rb") as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(chunk)
                size += len(chunk)
    except OSError as error:
        raise RankGateInputError(
            "cannot read input {0}: {1}".format(path.name, error)
        ) from error
    return digest.hexdigest(), size


def _identity(path: Path) -> Dict[str, object]:
    resolved = path.expanduser().resolve()
    require(resolved.is_file(), "input file does not exist: {0}".format(path.name))
    filename = _portable_basename(resolved.name, "input filename")
    digest, size = _sha256_file(resolved)
    return {"filename": filename, "sha256": digest, "size_bytes": size}


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
    raise RankGateInputError(
        "unsupported {0} value in {1}".format(type(value).__name__, location)
    )


def _embedded_json(raw, description: str):
    def reject_constant(value):
        raise ValueError("non-finite JSON constant {0}".format(value))

    try:
        value = json.loads(str(raw), parse_constant=reject_constant)
    except (TypeError, ValueError, json.JSONDecodeError) as error:
        raise RankGateInputError("malformed {0}".format(description)) from error
    _finite_json(value, description)
    return value


def _scalar(archive, key: str, filename: str):
    require(key in archive, "{0} is missing {1!r}".format(filename, key))
    try:
        value = np.asarray(archive[key])
    except (TypeError, ValueError) as error:
        raise RankGateInputError(
            "{0} field {1!r} cannot be loaded without pickle".format(filename, key)
        ) from error
    require(value.shape == (), "{0} field {1!r} is not scalar".format(filename, key))
    return value.item()


def _array(archive, key: str, filename: str, *, kind: Optional[str] = None):
    require(key in archive, "{0} is missing {1!r}".format(filename, key))
    try:
        value = np.asarray(archive[key])
    except (TypeError, ValueError) as error:
        raise RankGateInputError(
            "{0} field {1!r} cannot be loaded without pickle".format(filename, key)
        ) from error
    require(value.dtype.kind != "O", "{0} field {1!r} is an object array".format(filename, key))
    if kind == "float":
        require(
            value.dtype.kind in "fiu",
            "{0} field {1!r} is not numeric".format(filename, key),
        )
        value = value.astype(float, copy=False)
        require(
            np.all(np.isfinite(value)),
            "{0} field {1!r} contains non-finite values".format(filename, key),
        )
    elif kind == "integer":
        require(
            value.dtype.kind in "iu",
            "{0} field {1!r} is not integer".format(filename, key),
        )
    return value.copy()


def _integer(value, description: str, minimum: Optional[int] = None) -> int:
    require(not isinstance(value, (bool, np.bool_)), "{0} is boolean".format(description))
    try:
        numeric = float(value)
    except (TypeError, ValueError) as error:
        raise RankGateInputError("{0} is not an integer".format(description)) from error
    require(
        np.isfinite(numeric) and numeric.is_integer(),
        "{0} is not an integer".format(description),
    )
    result = int(numeric)
    if minimum is not None:
        require(result >= minimum, "{0} is below {1}".format(description, minimum))
    return result


def _exact_number(observed, expected: float, description: str) -> float:
    require(not isinstance(observed, (bool, np.bool_)), "{0} is boolean".format(description))
    try:
        value = float(observed)
    except (TypeError, ValueError) as error:
        raise RankGateInputError("{0} is not numeric".format(description)) from error
    require(
        np.isfinite(value) and value == float(expected),
        "{0}={1!r}; expected {2!r}".format(description, value, expected),
    )
    return value


def _source_identity(archive, filename: str) -> Dict[str, Dict[str, str]]:
    identity: Dict[str, Dict[str, str]] = {}
    for component in ("app", "core"):
        revision = str(_scalar(archive, component + "_revision", filename)).lower()
        tree_state = str(_scalar(archive, component + "_tree_state", filename))
        source_kind = str(_scalar(archive, component + "_source_kind", filename))
        require(
            FULL_GIT_REVISION.fullmatch(revision) is not None,
            "{0} has untraceable {1} revision".format(filename, component),
        )
        require(
            tree_state == "clean",
            "{0} {1} source is not clean".format(filename, component),
        )
        require(
            source_kind == "git-checkout",
            "{0} {1} source is not a public Git checkout".format(filename, component),
        )
        identity[component] = {
            "revision": revision,
            "tree_state": tree_state,
            "source_kind": source_kind,
        }
    core_url = str(_scalar(archive, "core_source_url", filename))
    require(core_url == PUBLIC_CORE_URL, "{0} has unexpected Core URL".format(filename))
    identity["core"]["source_url"] = core_url
    return identity


def _validate_fixed_configuration(archive, filename: str, mesh: Mapping[str, object]) -> Dict[str, object]:
    expected_text = {
        "result_schema": RESULT_SCHEMA,
        "case": "B",
        "integrator": "be",
        "formulation": "hex8_standard_pointwise_kappa",
        "material_kernel_formulation": "standard",
        "material_model_id": MATERIAL_MODEL_ID,
        "mass_representation": "consistent_q1_hex8",
        "fiber_sampling": "gp_direct_rule",
        "fiber_sampling_option": "gp-direct",
        "tbar_definition": "laplace_presolved",
        "tbar_metadata_schema": TBAR_SCHEMA,
        "point_sampling": "hex8_reference_isoparametric",
        "viscous_rate": "backward_difference",
        "parameter_variant": "benchmark_eta",
        "mesh_topology": "closed_multiblock_disk",
    }
    for key, expected in expected_text.items():
        observed = str(_scalar(archive, key, filename))
        require(
            observed == expected,
            "{0} {1}={2!r}; expected {3!r}".format(filename, key, observed, expected),
        )

    expected_numbers = {
        "dt": EXPECTED_DT_S,
        "t_end": EXPECTED_T_END_S,
        "load_horizon": EXPECTED_LOAD_HORIZON_S,
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
        "core_half_width": float(mesh["core_half_width"]),
    }
    for key, expected in expected_numbers.items():
        _exact_number(_scalar(archive, key, filename), expected, filename + " " + key)

    for key, expected in {
        "isotropic": False,
        "viscous_term_active": True,
        "flip_helix": True,
    }.items():
        observed = _scalar(archive, key, filename)
        require(
            isinstance(observed, (bool, np.bool_)) and bool(observed) is expected,
            "{0} {1} is not {2}".format(filename, key, expected),
        )

    for key, expected in (
        ("n_t", mesh["n_t"]),
        ("n_core", mesh["n_core"]),
        ("n_radial", mesh["n_radial"]),
        ("n_side", 0),
        ("n_mu", 0),
        ("n_theta", 0),
    ):
        require(
            _integer(_scalar(archive, key, filename), filename + " " + key) == int(expected),
            "{0} {1} does not match the explicit mesh".format(filename, key),
        )

    for field in ("tbar_source_filename", "tbar_metadata_filename"):
        _portable_basename(str(_scalar(archive, field, filename)), filename + " " + field)
    for field in ("tbar_source_sha256", "tbar_metadata_sha256"):
        digest = str(_scalar(archive, field, filename)).lower()
        require(FULL_SHA256.fullmatch(digest) is not None, "{0} has invalid {1}".format(filename, field))

    return {
        "material_model_id": MATERIAL_MODEL_ID,
        "dt_s": EXPECTED_DT_S,
        "t_end_s": EXPECTED_T_END_S,
        "load_horizon_s": EXPECTED_LOAD_HORIZON_S,
        "material_eta_pa_s": 100.0,
        "material_kappa_pa": 1.0e6,
        "mass_representation": "consistent_q1_hex8",
        "fiber_sampling": "gp_direct_rule",
        "tbar_definition": "laplace_presolved",
    }


def _validate_audit(archive, filename: str, *, expected_nodes: int, expected_elements: int):
    audit = _embedded_json(
        _scalar(archive, "pre_solve_audit_json", filename),
        filename + " pre-solve audit",
    )
    require(isinstance(audit, dict), "{0} pre-solve audit is not an object".format(filename))
    for name in ("geometry", "pressure", "robin"):
        record = audit.get(name)
        require(
            isinstance(record, dict) and record.get("passed") is True,
            "{0} {1} audit did not pass".format(filename, name),
        )
    geometry = audit["geometry"]
    require(
        geometry.get("schema") == "coupfe-cardiac-pre-solve-geometry-v1"
        and geometry.get("mesh_topology") == "closed_multiblock_disk"
        and geometry.get("require_closed") is True
        and geometry.get("nodes") == expected_nodes
        and geometry.get("elements") == expected_elements
        and geometry.get("unclassified_exterior_faces") == 0
        and geometry.get("multiply_labeled_faces") == 0
        and geometry.get("nonmanifold_faces") == 0
        and geometry.get("nonpositive_extended_jacobians") == 0,
        "{0} closed-geometry audit is incomplete or inconsistent".format(filename),
    )
    require(
        audit["pressure"].get("schema") == "coupfe-cardiac-pre-solve-pressure-v1"
        and audit["robin"].get("schema") == "coupfe-cardiac-pre-solve-robin-v1",
        "{0} boundary-audit schemas are unsupported".format(filename),
    )
    canonical = json.dumps(audit, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return audit, hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def _diagnostic_summary(archive, filename: str, times: np.ndarray, ranks: Optional[int]):
    diagnostics = _embedded_json(
        _scalar(archive, "nonlinear_step_diagnostics_json", filename),
        filename + " nonlinear diagnostics",
    )
    require(
        isinstance(diagnostics, list) and len(diagnostics) == EXPECTED_STEPS,
        "{0} does not contain all {1} step diagnostics".format(filename, EXPECTED_STEPS),
    )
    nonlinear = []
    linear = []
    ratios = []
    rejections = []
    snes_reasons: Dict[str, int] = {}
    ksp_reasons: Dict[str, int] = {}
    for index, record in enumerate(diagnostics, start=1):
        require(isinstance(record, dict), "{0} diagnostic {1} is malformed".format(filename, index))
        _finite_json(record, "{0} diagnostic {1}".format(filename, index))
        _exact_number(record.get("time"), float(times[index]), "diagnostic time")
        _exact_number(record.get("dt"), EXPECTED_DT_S, "diagnostic dt")
        snes = _integer(record.get("snes_converged_reason"), "SNES reason")
        ksp = _integer(record.get("ksp_converged_reason"), "KSP reason")
        require(snes > 0 and ksp > 0, "{0} diagnostic {1} did not converge".format(filename, index))
        final = float(record.get("final_residual_norm", np.nan))
        threshold = float(record.get("residual_acceptance_threshold", np.nan))
        require(
            np.isfinite(final)
            and np.isfinite(threshold)
            and threshold > 0.0
            and final <= threshold * (1.0 + 1.0e-9) + 1.0e-14,
            "{0} diagnostic {1} failed independent residual acceptance".format(filename, index),
        )
        nit = _integer(record.get("nonlinear_iterations"), "nonlinear iterations", 0)
        lit = _integer(record.get("linear_iterations"), "linear iterations", 0)
        rejected = _integer(record.get("function_domain_rejections", 0), "domain rejections", 0)
        if ranks is not None:
            require(
                _integer(record.get("ranks"), "diagnostic MPI ranks", 1) == ranks,
                "{0} diagnostic {1} has the wrong rank count".format(filename, index),
            )
        nonlinear.append(nit)
        linear.append(lit)
        rejections.append(rejected)
        ratios.append(final / threshold)
        snes_reasons[str(snes)] = snes_reasons.get(str(snes), 0) + 1
        ksp_reasons[str(ksp)] = ksp_reasons.get(str(ksp), 0) + 1
    return {
        "diagnostic_steps": len(diagnostics),
        "nonlinear_iterations": {
            "minimum": min(nonlinear),
            "maximum": max(nonlinear),
            "total": sum(nonlinear),
        },
        "linear_iterations": {
            "minimum": min(linear),
            "maximum": max(linear),
            "total": sum(linear),
        },
        "function_domain_rejections_total": sum(rejections),
        "maximum_final_residual_fraction": max(ratios),
        "snes_converged_reason_counts": snes_reasons,
        "ksp_converged_reason_counts": ksp_reasons,
    }


def _solver_configuration(archive, filename: str, *, ranks: Optional[int], ndof: int, n_element: int):
    configuration = _embedded_json(
        _scalar(archive, "solver_configuration_json", filename),
        filename + " solver configuration",
    )
    require(isinstance(configuration, dict), "{0} solver configuration is not an object".format(filename))
    common = {
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
    }
    for key, expected in common.items():
        require(configuration.get(key) == expected, "{0} solver has unexpected {1}".format(filename, key))

    evaluation = str(_scalar(archive, "element_evaluation_mode", filename))
    residual_only = _scalar(archive, "compiled_material_residual_only_available", filename)
    require(
        evaluation in {"joint", "split"}
        and configuration.get("element_evaluation_mode") == evaluation,
        "{0} element-evaluation metadata disagrees".format(filename),
    )
    require(
        isinstance(residual_only, (bool, np.bool_))
        and configuration.get("compiled_material_residual_only_available") is bool(residual_only),
        "{0} residual-only capability metadata disagrees".format(filename),
    )

    if ranks is None:
        require(
            str(_scalar(archive, "nonlinear_solver", filename)) == "petsc-snes"
            and str(_scalar(archive, "driver", filename)) == "examples/cardiac_benchmark/run.py"
            and configuration.get("name") == "petsc-snes",
            "{0} is not the required serial PETSc-SNES path".format(filename),
        )
        mpi_fields = [key for key in archive.files if key.startswith("mpi_")]
        require(not mpi_fields, "{0} serial result contains MPI-only metadata".format(filename))
        return configuration, None

    require(
        str(_scalar(archive, "nonlinear_solver", filename)) == "petsc-snes-mpi"
        and str(_scalar(archive, "driver", filename)) == "examples/cardiac_benchmark/run_mpi.py"
        and configuration.get("name") == "petsc-snes-mpi",
        "{0} is not the validated MPI companion".format(filename),
    )
    require(
        _scalar(archive, "mpi_enabled", filename) is True
        and _integer(_scalar(archive, "mpi_ranks", filename), "MPI ranks", 1) == ranks
        and _integer(_scalar(archive, "mpi_world_size", filename), "MPI world size", 1) == ranks
        and _integer(configuration.get("ranks"), "configured MPI ranks", 1) == ranks,
        "{0} rank metadata disagrees with its role".format(filename),
    )
    require(
        str(_scalar(archive, "mpi_implementation", filename)) == MPI_IMPLEMENTATION
        and configuration.get("implementation") == MPI_IMPLEMENTATION,
        "{0} has an unvalidated MPI implementation".format(filename),
    )
    expected_mpi = {
        "communicator": "PETSc.COMM_WORLD",
        "rank_local_element_assembly": True,
        "global_mesh_replicated": True,
        "collective_invalid_trial_policy": "all-owned-residual-entries-inf-for-bt",
        "commit_policy": "independent-global-residual-check-before-commit",
        "factor_solver_type": "superlu_dist",
        "configured_factor_solver_type": "superlu_dist",
        "mass_representation": "consistent_q1_hex8",
        "mass_partition": "owned-row-csr-all-touching-elements",
    }
    for key, expected in expected_mpi.items():
        require(configuration.get(key) == expected, "{0} MPI solver has unexpected {1}".format(filename, key))
    require(
        configuration.get("line_search_configuration_api")
        in {"SNES.getLineSearch", "namespaced PETSc option"},
        "{0} has unsupported line-search configuration".format(filename),
    )
    for key in ("petsc4py_version", "petsc_version"):
        require(isinstance(configuration.get(key), str) and configuration[key].strip(), "{0} lacks {1}".format(filename, key))

    counts = _array(archive, "mpi_local_element_counts", filename, kind="integer")
    quotient, remainder = divmod(n_element, ranks)
    expected_counts = np.full(ranks, quotient, dtype=np.int64)
    expected_counts[:remainder] += 1
    require(
        counts.shape == (ranks,) and np.array_equal(counts, expected_counts),
        "{0} element partition is incomplete".format(filename),
    )
    require(
        str(_scalar(archive, "mpi_partition", filename)) == "coupfe.partition_elements"
        and str(_scalar(archive, "mpi_build_layout", filename)) == "isolated-rank-directories"
        and str(_scalar(archive, "mpi_factor_solver_type", filename)) == "superlu_dist",
        "{0} MPI execution provenance is unsupported".format(filename),
    )
    require(
        str(_scalar(archive, "mpi_mass_partition", filename))
        == "owned-row-csr-all-touching-elements",
        "{0} lacks complete consistent-mass provenance".format(filename),
    )
    ranges = _array(archive, "mpi_mass_owned_row_ranges", filename, kind="integer")
    nnz = _array(archive, "mpi_mass_local_nnz", filename, kind="integer")
    touching = _array(archive, "mpi_mass_touching_element_counts", filename, kind="integer")
    require(
        ranges.shape == (ranks, 2)
        and nnz.shape == (ranks,)
        and touching.shape == (ranks,)
        and np.all(nnz > 0)
        and np.all(touching > 0),
        "{0} has malformed consistent-mass ownership arrays".format(filename),
    )
    expected_start = 0
    for start, stop in ranges:
        require(
            int(start) == expected_start and int(stop) > int(start),
            "{0} mass row ranges are not contiguous".format(filename),
        )
        expected_start = int(stop)
    require(expected_start == ndof, "{0} mass rows do not cover the system".format(filename))
    require(
        configuration.get("mass_owned_row_range") == ranges[0].astype(int).tolist()
        and _integer(configuration.get("mass_local_nnz"), "configured local mass nnz", 1)
        == int(nnz[0]),
        "{0} mass archive/configuration provenance disagrees".format(filename),
    )
    mpi = {
        "implementation": MPI_IMPLEMENTATION,
        "ranks": ranks,
        "partition": "coupfe.partition_elements",
        "local_element_counts": counts.astype(int).tolist(),
        "mass": {
            "representation": "consistent_q1_hex8",
            "partition": "owned-row-csr-all-touching-elements",
            "owned_row_ranges": ranges.astype(int).tolist(),
            "local_nnz": nnz.astype(int).tolist(),
            "touching_element_counts": touching.astype(int).tolist(),
        },
    }
    return configuration, mpi


def _expected_mesh_record():
    """Return the non-configurable retained rank-gate mesh."""
    plane_elements = (
        EXPECTED_N_CORE * EXPECTED_N_CORE
        + 4 * EXPECTED_N_CORE * EXPECTED_N_RADIAL
    )
    plane_nodes = (
        (EXPECTED_N_CORE + 1) ** 2
        + 4 * EXPECTED_N_CORE * EXPECTED_N_RADIAL
    )
    elements = EXPECTED_N_T * plane_elements
    nodes = (EXPECTED_N_T + 1) * plane_nodes
    return {
        "n_t": EXPECTED_N_T,
        "n_core": EXPECTED_N_CORE,
        "n_radial": EXPECTED_N_RADIAL,
        "core_half_width": EXPECTED_CORE_HALF_WIDTH,
        "expected_nodes": int(nodes),
        "expected_elements": int(elements),
        "expected_degrees_of_freedom": int(3 * nodes),
    }


def _load_run(path: Path, identity: Mapping[str, object], role: str, mesh: Mapping[str, object]):
    filename = str(identity["filename"])
    try:
        context = np.load(path.expanduser().resolve(), allow_pickle=False)
    except (OSError, ValueError) as error:
        raise RankGateInputError(
            "cannot load {0} without pickle: {1}".format(filename, error)
        ) from error
    with context as archive:
        completed = _integer(_scalar(archive, "completed_steps", filename), "completed steps", 0)
        expected = _integer(_scalar(archive, "expected_steps", filename), "expected steps", 0)
        converged = _scalar(archive, "converged", filename)
        require(
            isinstance(converged, (bool, np.bool_))
            and bool(converged)
            and completed == expected == EXPECTED_STEPS,
            "{0} is not a completed {1}-step solve".format(filename, EXPECTED_STEPS),
        )
        fixed = _validate_fixed_configuration(archive, filename, mesh)
        source = _source_identity(archive, filename)

        arrays = {}
        for field in EXACT_FLOAT_ARRAY_FIELDS + TOLERANT_FIELDS:
            arrays[field] = _array(archive, field, filename, kind="float")
        for field in EXACT_INTEGER_ARRAY_FIELDS:
            arrays[field] = _array(archive, field, filename, kind="integer")
        arrays["element_pressure_peak_pa"] = _array(
            archive, "element_pressure_peak_pa", filename, kind="float"
        )
        require(
            arrays["element_pressure_peak_pa"].shape == (0,),
            "{0} std-kappa result must have an empty element-pressure array".format(filename),
        )

        require(
            arrays["times"].shape == EXPECTED_TIMES.shape
            and np.array_equal(arrays["times"], EXPECTED_TIMES),
            "{0} time grid is not the exact 0--0.32 s prefix".format(filename),
        )
        require(
            arrays["tau"].shape == EXPECTED_TIMES.shape
            and np.array_equal(arrays["tau"], np.zeros_like(EXPECTED_TIMES)),
            "{0} activation history is not exactly zero".format(filename),
        )
        require(
            arrays["pres"].shape == EXPECTED_TIMES.shape
            and np.any(np.abs(arrays["pres"]) > 0.0),
            "{0} pressure history is invalid".format(filename),
        )
        expected_nodes = int(mesh["expected_nodes"])
        expected_elements = int(mesh["expected_elements"])
        ndof = int(mesh["expected_degrees_of_freedom"])
        require(arrays["nodes"].shape == (expected_nodes, 3), "{0} node count disagrees with mesh".format(filename))
        require(arrays["elems"].shape == (expected_elements, 8), "{0} element count disagrees with mesh".format(filename))
        require(
            np.all(arrays["elems"] >= 0) and np.all(arrays["elems"] < expected_nodes),
            "{0} connectivity is outside the mesh".format(filename),
        )
        require(arrays["u0"].shape == (EXPECTED_STEPS + 1, 3), "{0} u0 shape is invalid".format(filename))
        require(arrays["u1"].shape == (EXPECTED_STEPS + 1, 3), "{0} u1 shape is invalid".format(filename))
        require(arrays["U_peak"].shape == (ndof,), "{0} U_peak shape is invalid".format(filename))
        require(
            arrays["det_f_gauss_peak"].shape == (expected_elements, 8)
            and np.all(arrays["det_f_gauss_peak"] > 0.0),
            "{0} det(F) peak field is invalid".format(filename),
        )
        require(np.array_equal(arrays["p0"], np.array([0.025, 0.030, 0.0])), "{0} p0 is not the benchmark landmark".format(filename))
        require(np.array_equal(arrays["p1"], np.array([0.000, 0.030, 0.0])), "{0} p1 is not the benchmark landmark".format(filename))

        audit, audit_sha = _validate_audit(
            archive,
            filename,
            expected_nodes=expected_nodes,
            expected_elements=expected_elements,
        )
        ranks = EXPECTED_MPI_RANKS.get(role)
        configuration, mpi = _solver_configuration(
            archive,
            filename,
            ranks=ranks,
            ndof=ndof,
            n_element=expected_elements,
        )
        diagnostics = _diagnostic_summary(archive, filename, arrays["times"], ranks)

        invariant_dtype_fields = {}
        for key in archive.files:
            if key in RUN_DEPENDENT_COMMON_FIELDS or key.startswith("mpi_"):
                continue
            value = _array(archive, key, filename)
            if value.dtype.kind in "iubUS":
                invariant_dtype_fields[key] = value

        tbar = {
            "field_filename": str(_scalar(archive, "tbar_source_filename", filename)),
            "field_sha256": str(_scalar(archive, "tbar_source_sha256", filename)).lower(),
            "metadata_filename": str(_scalar(archive, "tbar_metadata_filename", filename)),
            "metadata_sha256": str(_scalar(archive, "tbar_metadata_sha256", filename)).lower(),
            "metadata_schema": TBAR_SCHEMA,
        }

    return {
        "arrays": arrays,
        "invariant_dtype_fields": invariant_dtype_fields,
        "audit": audit,
        "audit_sha256": audit_sha,
        "source": source,
        "fixed": fixed,
        "tbar": tbar,
        "completion": {
            "converged": True,
            "completed_steps": completed,
            "expected_steps": expected,
        },
        "solver": {
            "name": "petsc-snes" if role == "serial" else "petsc-snes-mpi",
            "diagnostics": diagnostics,
            "mpi": mpi,
        },
    }


def _assert_exact(reference: np.ndarray, candidate: np.ndarray, field: str, role: str) -> None:
    require(
        reference.shape == candidate.shape
        and reference.dtype.kind == candidate.dtype.kind
        and np.array_equal(reference, candidate),
        "{0} differs exactly in invariant field {1}".format(role, field),
    )


def _numerical_difference(reference: np.ndarray, candidate: np.ndarray, field: str):
    require(reference.shape == candidate.shape, "shape mismatch for {0}".format(field))
    difference = np.abs(candidate - reference)
    max_absolute = float(np.max(difference)) if difference.size else 0.0
    reference_scale = float(np.max(np.abs(reference))) if reference.size else 0.0
    if reference_scale == 0.0:
        max_relative = 0.0 if max_absolute == 0.0 else max_absolute / NUMERICAL_ATOL
        relative_definition = "max_abs/atol because the reference maximum is zero"
    else:
        max_relative = max_absolute / reference_scale
        relative_definition = "max_abs/max(abs(serial_reference))"
    allowed = NUMERICAL_ATOL + NUMERICAL_RTOL * np.abs(reference)
    fractions = np.divide(difference, allowed, out=np.zeros_like(difference), where=allowed > 0.0)
    maximum_fraction = float(np.max(fractions)) if fractions.size else 0.0
    passed = bool(np.all(difference <= allowed))
    return {
        "shape": list(reference.shape),
        "max_absolute_difference": max_absolute,
        "max_relative_difference": max_relative,
        "relative_difference_definition": relative_definition,
        "maximum_tolerance_fraction": maximum_fraction,
        "passed": passed,
    }


def build_report(paths: Mapping[str, Path], mesh: Mapping[str, object]):
    identities = {role: _identity(paths[role]) for role in INPUT_ROLES}
    resolved = [paths[role].expanduser().resolve() for role in INPUT_ROLES]
    require(len(set(resolved)) == len(resolved), "rank-gate inputs must be four distinct files")
    runs = {
        role: _load_run(paths[role], identities[role], role, mesh)
        for role in INPUT_ROLES
    }
    serial = runs["serial"]
    for role in ("mpi1", "mpi2", "mpi4"):
        candidate = runs[role]
        require(candidate["source"] == serial["source"], "{0} source revisions differ from serial".format(role))
        require(candidate["fixed"] == serial["fixed"], "{0} fixed setup differs from serial".format(role))
        require(candidate["tbar"] == serial["tbar"], "{0} Laplace-tbar identity differs from serial".format(role))
        require(candidate["audit"] == serial["audit"], "{0} pre-solve audits differ from serial".format(role))
        for field in EXACT_FLOAT_ARRAY_FIELDS + EXACT_INTEGER_ARRAY_FIELDS:
            _assert_exact(serial["arrays"][field], candidate["arrays"][field], field, role)
        require(
            set(candidate["invariant_dtype_fields"]) == set(serial["invariant_dtype_fields"]),
            "{0} invariant integer/text field inventory differs from serial".format(role),
        )
        for field, reference in serial["invariant_dtype_fields"].items():
            _assert_exact(reference, candidate["invariant_dtype_fields"][field], field, role)

    comparisons = {}
    for role in ("mpi1", "mpi2", "mpi4"):
        fields = {
            field: _numerical_difference(
                serial["arrays"][field], runs[role]["arrays"][field], field
            )
            for field in TOLERANT_FIELDS
        }
        passed = all(record["passed"] for record in fields.values())
        require(passed, "{0} exceeds the fixed numerical rank-gate tolerance".format(role))
        comparisons[role] = {"reference": "serial", "fields": fields, "passed": True}

    report_runs = {}
    for role in INPUT_ROLES:
        report_runs[role] = {
            "input": identities[role],
            "completion": runs[role]["completion"],
            "solver": runs[role]["solver"],
        }

    report = {
        "schema": REPORT_SCHEMA,
        "created_utc": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "status": "passed",
        "purpose": (
            "Closed Case B serial-versus-1/2/4-rank numerical equivalence "
            "through the snap interval."
        ),
        "expected_setup": {
            **dict(mesh),
            **serial["fixed"],
            "case": "B",
            "activation": "zero",
            "perturbation_m": 0.0,
            "serial_solver": "petsc-snes",
            "mpi_solver": "petsc-snes-mpi",
            "mpi_implementation": MPI_IMPLEMENTATION,
        },
        "source_identity": serial["source"],
        "laplace_tbar_identity": serial["tbar"],
        "pre_solve_audit": {
            "canonical_sha256": serial["audit_sha256"],
            "identical_across_all_runs": True,
            "geometry_passed": True,
            "pressure_passed": True,
            "robin_passed": True,
        },
        "comparison_definition": {
            "reference": "serial",
            "candidate_roles": ["mpi1", "mpi2", "mpi4"],
            "exact_float_arrays": list(EXACT_FLOAT_ARRAY_FIELDS),
            "exact_integer_arrays": list(EXACT_INTEGER_ARRAY_FIELDS),
            "all_other_common_integer_boolean_text_fields": "exact",
            "tolerant_fields": list(TOLERANT_FIELDS),
            "numerical_rule": "abs(candidate-serial) <= atol + rtol*abs(serial), elementwise",
            "rtol": NUMERICAL_RTOL,
            "atol": NUMERICAL_ATOL,
            "tolerance_origin": "historical CoupFE serial/MPI rank-equivalence gate precedent",
            "tolerances_are_cli_configurable": False,
            "element_pressure_peak_pa": "required empty for std-kappa",
        },
        "runs": report_runs,
        "comparisons": comparisons,
        "claim_boundary": (
            "Passing establishes numerical equivalence of these four retained "
            "archives for the declared setup and tolerance. It is not a mesh-"
            "convergence, performance-scaling, or benchmark-accuracy claim."
        ),
    }
    _finite_json(report, "rank-gate report")
    return report


def _write_json_atomic(path: Path, report: Mapping[str, object]) -> None:
    destination = path.expanduser().resolve()
    require(destination.suffix.lower() == ".json", "rank-gate output must be JSON")
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix="." + destination.name + ".",
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
            json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, destination)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def main(argv: Optional[Sequence[str]] = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--serial", type=Path, required=True)
    parser.add_argument("--mpi1", type=Path, required=True)
    parser.add_argument("--mpi2", type=Path, required=True)
    parser.add_argument("--mpi4", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args(argv)

    paths = {role: getattr(args, role) for role in INPUT_ROLES}
    output = args.report.expanduser().resolve()
    require(
        output not in {path.expanduser().resolve() for path in paths.values()},
        "report path must differ from every input",
    )
    mesh = _expected_mesh_record()
    report = build_report(paths, mesh)
    _write_json_atomic(args.report, report)
    print("rank gate passed: serial == MPI 1/2/4 within fixed tolerance")
    print("saved report -> {0}".format(args.report.name))
    return report


if __name__ == "__main__":
    main()
