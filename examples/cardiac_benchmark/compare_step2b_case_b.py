#!/usr/bin/env python3
"""Fail-closed comparison with the published Benchmark 1 Step-2 Case-B data.

The publisher archive contains eleven similarly named Case-B pickle files, but
``figures.py::TEAMS_DATASETS_B`` selects exactly ten of them.  This utility
loads those ten files in the publisher's order and deliberately excludes the
duplicate generic SimVascular alias.  Pickles are accepted only after a hash
manifest has matched and are decoded by an unpickler that permits the small
set of NumPy constructors used by the archive.

"Plumbing" in the JSON report means that identities, schemas, units, and the
CoupFE Step-2/Case-B physical configuration passed.  "Reproduction" is kept
separate: it reports numerical errors and event/branch metrics but does not
invent a pass tolerance that the paper did not define.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import pickle
import re
import sys
import tempfile
from typing import Dict, Iterable, Mapping, Optional, Sequence, Tuple

import numpy as np
import numpy.core.multiarray as _numpy_multiarray

try:  # package import
    from .activation import p_of_t, tau_of_t
    from .benchmark_parameters import (
        RUNTIME_SOURCE_FILES,
        benchmark_configuration,
        runtime_source_metadata,
    )
except ImportError:  # direct script import
    from activation import p_of_t, tau_of_t
    from benchmark_parameters import (
        RUNTIME_SOURCE_FILES,
        benchmark_configuration,
        runtime_source_metadata,
    )


REPORT_SCHEMA = "coupfe-cardiac-step2b-publisher-comparison-v1"
RESULT_SCHEMA = "coupfe-cardiac-result-v1"
HASH_MANIFEST_SCHEMA = "coupfe-cardiac-step2b-publisher-hashes-v1"
BENCHMARK_CONFIGURATION_ID = (
    "benchmark-1-step-2-case-B-active-stress-plus-pressure"
)
BENCHMARK_LOAD_CONTRACT = "active-stress-plus-pressure"
BENCHMARK_PEAK_LOAD_DEFINITION = "argmax(abs(active_tension_pa+pressure_pa))"
GENERALIZED_ALPHA_PARAMETERS = {
    "generalized_alpha_alpha_m": 0.2,
    "generalized_alpha_alpha_f": 0.4,
    "generalized_alpha_gamma": 0.7,
    "generalized_alpha_beta": 0.36,
}
GENERALIZED_ALPHA_STAGE_CONTRACT = "simula-source-matched-v1"

PUBLISHED_TIME_S = np.arange(101, dtype=float) * 0.01
COUPFE_TIME_S = np.arange(1001, dtype=float) * 0.001
POINTS = ("p0", "p1")
COMPONENTS = ("ux", "uy", "uz")
LANDMARKS_M = {
    "p0": np.array([0.025, 0.030, 0.0]),
    "p1": np.array([0.000, 0.030, 0.0]),
}

# Order and spelling are copied from figures.py::TEAMS_DATASETS_B.
PUBLISHER_SELECTION = (
    ("carpentry", "CARPentry", "monoventricular_blinded_B_group_carpentry.pkl"),
    ("ambit", "Ambit", "monoventricular_blinded_B_group_ambit.pkl"),
    ("4c", "4C", "monoventricular_blinded_B_group_4c.pkl"),
    ("simula", "Simula", "monoventricular_blinded_B_group_simula.pkl"),
    ("chimera", "CHimeRA", "monoventricular_blinded_B_group_chimera.pkl"),
    ("cheart", "CHeart", "monoventricular_blinded_B_group_cheart.pkl"),
    ("lifex", "lifeX", "monoventricular_blinded_B_group_lifex.pkl"),
    (
        "simvascular_p1",
        "SimVascular P1",
        "monoventricular_blinded_B_group_simvascular_p1p1.pkl",
    ),
    (
        "simvascular_p2",
        "SimVascular P2",
        "monoventricular_blinded_B_group_simvascular_p2.pkl",
    ),
    ("comsol", "COMSOL", "monoventricular_blinded_B_group_comsol.pkl"),
)
EXCLUDED_GENERIC_SIMVASCULAR = "monoventricular_blinded_B_group_simvascular.pkl"
DEFAULT_HASH_MANIFEST = Path(__file__).with_name(
    "step2b_case_b_reference_hashes.json"
)
OFFICIAL_HASH_MANIFEST_SHA256 = (
    "8392e28a1e6971b0c5c26e79aa2d0e86a48935cb7174c2454920b08a486631ac"
)
OFFICIAL_HASH_MANIFEST_SIZE_BYTES = 2658
OFFICIAL_SOURCE_DOI = "https://doi.org/10.5281/zenodo.14260459"
OFFICIAL_SOURCE_TITLE = "A software benchmark for cardiac elastodynamics"
OFFICIAL_SOURCE_CREATORS = (
    "Arostica Barrera, R.A.",
    "Bertoglio, Cristobal",
)
OFFICIAL_SOURCE_LICENSE = "CC-BY-4.0"
OFFICIAL_SOURCE_LICENSE_URL = "https://creativecommons.org/licenses/by/4.0/"
OFFICIAL_SELECTION_SOURCE = (
    "benchmark_article_data/results_time_curves/figures.py::TEAMS_DATASETS_B"
)
RUNTIME_SOURCE_MANIFEST_SCHEMA = (
    "coupfe-cardiac-step2b-runtime-source-hashes-v1"
)
DEFAULT_RUNTIME_SOURCE_MANIFEST = Path(__file__).with_name(
    "step2b_case_b_runtime_source_hashes.json"
)
REVIEWED_RUNTIME_SOURCE_MANIFEST_SHA256 = (
    "d39ccbd6e67d7517b31a536f0c34472afb770bfd76916a3572e20d52acf39a41"
)
REVIEWED_RUNTIME_SOURCE_MANIFEST_SIZE_BYTES = 2628
REVIEWED_RESULT_APP_REVISION = "e9b7d9084b24f7098170a221061eb159d0b090c1"
REVIEWED_RUNTIME_SOURCE_SNAPSHOT_REVISION = (
    "d06c3e9e827cdd7fc3208d38c6c2fda2fcb6c626"
)
REVIEWED_RUNTIME_SOURCE_SHA256 = (
    "6b96395761dd3203f0e9ffab90a77d6389dca13cdad43490a1deac95073480f1"
)

WINDOWS_S = (
    ("contraction_event", 0.16, 0.32),
    ("post_snap_plateau", 0.32, 0.48),
    ("relaxation", 0.48, 0.58),
    ("full_history", 0.00, 1.00),
)
FIXED_DOWNWARD_THRESHOLDS_M = (-5.0e-3, -15.0e-3)
NORMALIZED_DROP_FRACTIONS = (0.10, 0.50, 0.90)
FULL_SHA256 = re.compile(r"[0-9a-f]{64}")
FULL_SHA1 = re.compile(r"[0-9a-f]{40}")
SIMULA_TEAM_INDEX = next(
    index
    for index, selected in enumerate(PUBLISHER_SELECTION)
    if selected[0] == "simula"
)

EXPECTED_MATERIAL = {
    "a": 295.0,
    "b": 8.023,
    "a_f": 92360.0,
    "b_f": 16.026,
    "a_s": 12405.0,
    "b_s": 11.12,
    "a_fs": 1080.0,
    "b_fs": 11.436,
    "kappa": 1.0e6,
    "k_sw": 100.0,
    "eta": 100.0,
    "Ta": 0.0,
}
EXPECTED_ACTIVATION = {
    "t_sys": 0.16,
    "t_dias": 0.484,
    "gamma": 0.005,
    "a_max": 5.0,
    "a_min": -30.0,
    "sigma_0": 1.0e5,
}
EXPECTED_PRESSURE = {
    "t_sys_pre": 0.17,
    "t_dias_pre": 0.484,
    "gamma": 0.005,
    "a_max": 5.0,
    "a_min": -30.0,
    "alpha_pre": 5.0,
    "alpha_mid": 1.0,
    "sigma_pre": 7000.0,
    "sigma_mid": 16000.0,
}
EXPECTED_CORE_REVISION = "e2f42ed5772850a0a23a2ce434f430c287eae5c8"
EXPECTED_CORE_URL = "https://github.com/tengzhang48/CoupFE.git"
EXPECTED_MATERIAL_MODEL_ID = (
    "holzapfel-ogden-smooth-switch-complete-energy-derivative-v1"
)
FORMULATION_CONTRACTS = {
    "hex8_standard_pointwise_kappa": {
        "material_kappa_pa": 1.0e6,
        "local_pressure_bulk_modulus_pa": 0.0,
        "local_pressure_volume_law": "not-applicable",
        "mpi_implementation": (
            "cardiac-owned-distributed-closed-std-kappa-"
            "generalized-alpha-step0"
        ),
    },
    "hex8_local_pressure_p0_condensed_logj": {
        "material_kappa_pa": 0.0,
        "local_pressure_bulk_modulus_pa": 1.0e6,
        "local_pressure_volume_law": "linear-reference-volume-mean-log-j-v1",
        "mpi_implementation": (
            "cardiac-owned-distributed-closed-local-pressure-"
            "generalized-alpha-step0"
        ),
    },
    "hex8_local_pressure_p0_condensed_mean_logj_paper_j2": {
        "material_kappa_pa": 0.0,
        "local_pressure_bulk_modulus_pa": 1.0e6,
        "local_pressure_volume_law": (
            "paper-j2-of-reference-volume-weighted-geometric-mean-j-v1"
        ),
        "mpi_implementation": (
            "cardiac-owned-distributed-closed-local-pressure-mean-logj-"
            "paper-j2-generalized-alpha-step0"
        ),
    },
}


class ComparisonInputError(RuntimeError):
    """An input is unsafe, incomplete, or not the declared benchmark case."""


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
            "cannot read {0}: {1}".format(path.name, error)
        ) from error
    return digest.hexdigest(), size


def _identity(path: Path) -> Dict[str, object]:
    _payload, identity = _read_identity_payload(path)
    return identity


def _read_identity_payload(path: Path) -> Tuple[bytes, Dict[str, object]]:
    """Read once, then identify and consume the same immutable byte snapshot."""
    require(path.is_file(), "input file does not exist: {0}".format(path))
    try:
        payload = path.read_bytes()
    except OSError as error:
        raise ComparisonInputError(
            "cannot read {0}: {1}".format(path.name, error)
        ) from error
    return payload, {
        "filename": path.name,
        "sha256": hashlib.sha256(payload).hexdigest(),
        "size_bytes": len(payload),
    }


def _load_json_payload(payload: bytes, name: str, description: str):
    def reject_constant(value):
        raise ValueError("non-finite JSON constant {0}".format(value))

    try:
        text = payload.decode("utf-8")
        value = json.loads(text, parse_constant=reject_constant)
    except (UnicodeError, ValueError, json.JSONDecodeError) as error:
        raise ComparisonInputError(
            "cannot parse {0} {1}: {2}".format(description, name, error)
        ) from error
    return value


def _load_json(path: Path, description: str):
    payload, _identity_record = _read_identity_payload(path)
    return _load_json_payload(payload, path.name, description)


def load_hash_manifest(path: Path) -> Tuple[Mapping[str, object], Dict[str, object]]:
    """Load and structurally validate one exact ten-file manifest."""
    path = Path(path).expanduser().resolve()
    payload, identity = _read_identity_payload(path)
    require(
        identity["sha256"] == OFFICIAL_HASH_MANIFEST_SHA256
        and identity["size_bytes"] == OFFICIAL_HASH_MANIFEST_SIZE_BYTES,
        "publisher hash manifest is not the reviewed official manifest",
    )
    manifest = _load_json_payload(payload, path.name, "hash manifest")
    require(isinstance(manifest, dict), "hash manifest must be a JSON object")
    require(
        set(manifest) == {"schema", "source_doi", "selection_source", "files"},
        "hash manifest has missing or unexpected top-level fields",
    )
    require(
        manifest["schema"] == HASH_MANIFEST_SCHEMA,
        "unsupported hash manifest schema",
    )
    require(
        manifest["source_doi"] == OFFICIAL_SOURCE_DOI
        and manifest["selection_source"] == OFFICIAL_SELECTION_SOURCE,
        "hash manifest does not name the reviewed official source selection",
    )
    entries = manifest["files"]
    require(isinstance(entries, list), "hash manifest files must be a list")
    require(
        len(entries) == len(PUBLISHER_SELECTION),
        "hash manifest must contain exactly ten selected files",
    )
    by_name = {}
    for index, (entry, selected) in enumerate(zip(entries, PUBLISHER_SELECTION)):
        team_id, label, filename = selected
        require(isinstance(entry, dict), "hash manifest file entry must be an object")
        require(
            set(entry) == {"team_id", "label", "filename", "sha256", "size_bytes"},
            "hash manifest file entry has missing or unexpected fields",
        )
        require(
            entry["team_id"] == team_id
            and entry["label"] == label
            and entry["filename"] == filename,
            "hash manifest disagrees with publisher selection at index {0}".format(
                index
            ),
        )
        digest = entry["sha256"]
        size = entry["size_bytes"]
        require(
            isinstance(digest, str) and FULL_SHA256.fullmatch(digest) is not None,
            "hash manifest contains an invalid SHA-256",
        )
        require(
            isinstance(size, int) and not isinstance(size, bool) and size > 0,
            "hash manifest contains an invalid byte count",
        )
        require(filename not in by_name, "hash manifest contains a duplicate filename")
        by_name[filename] = entry
    return manifest, identity


def load_reviewed_runtime_source_manifest(
    path: Path = DEFAULT_RUNTIME_SOURCE_MANIFEST,
) -> Tuple[Dict[str, str], Dict[str, object]]:
    """Load the exact dirty-tree source identity used by the retained run."""
    path = Path(path).expanduser().resolve()
    payload, identity = _read_identity_payload(path)
    require(
        identity["sha256"] == REVIEWED_RUNTIME_SOURCE_MANIFEST_SHA256
        and identity["size_bytes"]
        == REVIEWED_RUNTIME_SOURCE_MANIFEST_SIZE_BYTES,
        "runtime source manifest is not the reviewed retained-run manifest",
    )
    manifest = _load_json_payload(payload, path.name, "runtime source manifest")
    require(
        isinstance(manifest, dict),
        "runtime source manifest must be a JSON object",
    )
    require(
        set(manifest)
        == {
            "schema",
            "result_application_revision",
            "result_application_tree_state",
            "repository_snapshot_revision",
            "runtime_source_sha256",
            "files",
        },
        "runtime source manifest has missing or unexpected top-level fields",
    )
    require(
        manifest["schema"] == RUNTIME_SOURCE_MANIFEST_SCHEMA
        and manifest["result_application_revision"]
        == REVIEWED_RESULT_APP_REVISION
        and manifest["result_application_tree_state"] == "dirty"
        and manifest["repository_snapshot_revision"]
        == REVIEWED_RUNTIME_SOURCE_SNAPSHOT_REVISION
        and manifest["runtime_source_sha256"]
        == REVIEWED_RUNTIME_SOURCE_SHA256,
        "runtime source manifest has altered retained-run provenance",
    )
    files = manifest["files"]
    # The hard-coded whole-file hash above already pins this historical file
    # inventory. Do not compare it with the evolving current runtime inventory.
    require(
        isinstance(files, dict) and bool(files),
        "runtime source manifest has no result-producing sources",
    )
    require(
        all(
            isinstance(digest, str)
            and FULL_SHA256.fullmatch(digest) is not None
            for digest in files.values()
        ),
        "runtime source manifest contains an invalid source SHA-256",
    )
    encoded = json.dumps(
        files, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    aggregate = hashlib.sha256(encoded.encode("utf-8")).hexdigest()
    require(
        aggregate == manifest["runtime_source_sha256"],
        "runtime source manifest aggregate does not match its source hashes",
    )
    return {
        "benchmark_runtime_source_manifest_json": encoded,
        "benchmark_runtime_source_sha256": aggregate,
    }, identity


class _RestrictedNumpyUnpickler(pickle.Unpickler):
    """Permit only the NumPy constructors present in the publisher pickles."""

    _ALLOWED = {
        ("numpy.core.multiarray", "_reconstruct"): _numpy_multiarray._reconstruct,
        ("numpy._core.multiarray", "_reconstruct"): _numpy_multiarray._reconstruct,
        ("numpy.core.multiarray", "scalar"): _numpy_multiarray.scalar,
        ("numpy._core.multiarray", "scalar"): _numpy_multiarray.scalar,
        ("numpy", "ndarray"): np.ndarray,
        ("numpy", "dtype"): np.dtype,
    }

    def find_class(self, module: str, name: str):
        try:
            return self._ALLOWED[(module, name)]
        except KeyError as error:
            raise pickle.UnpicklingError(
                "disallowed pickle global {0}.{1}".format(module, name)
            ) from error


def restricted_numpy_loads(payload: bytes, description: str = "pickle"):
    try:
        stream = io.BytesIO(payload)
        value = _RestrictedNumpyUnpickler(stream).load()
        require(stream.read(1) == b"", "{0} has trailing bytes".format(description))
        return value
    except ComparisonInputError:
        raise
    except (pickle.UnpicklingError, EOFError, ValueError, TypeError) as error:
        raise ComparisonInputError(
            "cannot safely decode {0}: {1}".format(description, error)
        ) from error


def _numeric_vector(value, *, length: int, description: str) -> np.ndarray:
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as error:
        raise ComparisonInputError(
            "{0} must be numeric".format(description)
        ) from error
    require(
        array.shape == (length,),
        "{0} must have shape ({1},)".format(description, length),
    )
    require(np.all(np.isfinite(array)), "{0} must be finite".format(description))
    return array.copy()


def _validate_publisher_dataset(
    value, *, team_id: str, filename: str
) -> Tuple[np.ndarray, np.ndarray]:
    require(isinstance(value, dict), "{0} must contain a dictionary".format(filename))
    require(
        set(value) == {"time", "displacement", "stress", "volume"},
        "{0} has an unexpected top-level schema".format(filename),
    )
    time = _numeric_vector(value["time"], length=101, description=filename + " time")
    expected_time = PUBLISHED_TIME_S.copy()
    if team_id == "ambit":
        # This is the one timestamp quirk in the hash-pinned publisher files.
        # Upstream statistics combine samples by index using CARPentry's grid.
        expected_time[0] = 0.001
    require(
        np.allclose(time, expected_time, rtol=0.0, atol=2.0e-12),
        "{0} time is not the publisher's 101-sample seconds grid".format(filename),
    )
    displacement = value["displacement"]
    require(
        isinstance(displacement, dict) and set(displacement) == set(POINTS),
        "{0} displacement must contain exactly p0 and p1".format(filename),
    )
    result = np.empty((2, 3, 101), dtype=float)
    for point_index, point in enumerate(POINTS):
        components = displacement[point]
        require(
            isinstance(components, dict)
            and set(components) == {"ux", "uy", "uz", "magnitude"},
            "{0} {1} displacement schema is invalid".format(filename, point),
        )
        for component_index, component in enumerate(COMPONENTS):
            result[point_index, component_index] = _numeric_vector(
                components[component],
                length=101,
                description="{0} {1} {2}".format(filename, point, component),
            )
        magnitude = _numeric_vector(
            components["magnitude"],
            length=101,
            description="{0} {1} magnitude".format(filename, point),
        )
        calculated = np.linalg.norm(result[point_index], axis=0)
        require(
            np.allclose(magnitude, calculated, rtol=2.0e-12, atol=2.0e-14),
            "{0} {1} magnitude disagrees with its SI components".format(
                filename, point
            ),
        )

    # The pickle has no unit tag.  Hash identity plus the documented one-second
    # grid and metre-scale physiological bounds make a mm-as-m error fail closed.
    require(
        float(np.max(np.abs(result))) < 0.10,
        "{0} displacement is not plausible in metres".format(filename),
    )
    require(
        -0.050 < float(np.min(result[0, 0])) < -0.020
        and -0.040 < float(np.min(result[1, 0])) < -0.015,
        "{0} x displacement is inconsistent with published SI metres".format(
            filename
        ),
    )
    return result, time


def load_publisher_reference(
    data_directory: Path, hash_manifest: Path
) -> Dict[str, object]:
    """Load the exact upstream Step-2 Case-B team selection in its fixed order."""
    data_directory = Path(data_directory).expanduser().resolve()
    require(data_directory.is_dir(), "publisher data directory does not exist")
    manifest, manifest_identity = load_hash_manifest(hash_manifest)
    entries = {entry["filename"]: entry for entry in manifest["files"]}
    expected_names = {selected[2] for selected in PUBLISHER_SELECTION}
    actual_names = {
        path.name
        for path in data_directory.glob("monoventricular_blinded_B_group_*.pkl")
        if path.is_file()
    }
    missing = sorted(expected_names - actual_names)
    unexpected = sorted(
        actual_names - expected_names - {EXCLUDED_GENERIC_SIMVASCULAR}
    )
    require(
        not missing,
        "missing selected publisher files: {0}".format(", ".join(missing)),
    )
    require(
        not unexpected,
        "unexpected publisher Case-B files: {0}".format(", ".join(unexpected)),
    )

    team_arrays = []
    team_times = []
    identities = []
    for team_id, label, filename in PUBLISHER_SELECTION:
        path = data_directory / filename
        payload, identity = _read_identity_payload(path)
        expected = entries[filename]
        require(
            identity["sha256"] == expected["sha256"]
            and identity["size_bytes"] == expected["size_bytes"],
            "selected publisher input changed: {0}".format(filename),
        )
        value = restricted_numpy_loads(payload, filename)
        displacement, time = _validate_publisher_dataset(
            value, team_id=team_id, filename=filename
        )
        team_arrays.append(displacement)
        team_times.append(time)
        identities.append(
            {
                "team_id": team_id,
                "label": label,
                **identity,
            }
        )

    generic = data_directory / EXCLUDED_GENERIC_SIMVASCULAR
    generic_identity = None
    if generic.is_file():
        _generic_payload, generic_identity = _read_identity_payload(generic)
        p2_identity = identities[8]
        require(
            generic_identity["sha256"] == p2_identity["sha256"]
            and generic_identity["size_bytes"] == p2_identity["size_bytes"],
            "excluded generic SimVascular file is not a byte-for-byte P2 duplicate",
        )

    teams = np.stack(team_arrays, axis=0)
    require(teams.shape == (10, 2, 3, 101), "publisher team stack is incomplete")
    native_times = np.stack(team_times, axis=0)
    require(
        native_times.shape == (10, 101),
        "publisher team time-grid stack is incomplete",
    )
    return {
        "time_s": PUBLISHED_TIME_S.copy(),
        "team_times_s": native_times,
        "teams_m": teams,
        "mean_m": np.mean(teams, axis=0),
        "std_m": np.std(teams, axis=0),
        "identities": identities,
        "manifest_identity": manifest_identity,
        "manifest_source_doi": manifest["source_doi"],
        "manifest_selection_source": manifest["selection_source"],
        "excluded_generic_identity": generic_identity,
    }


def _archive_scalar(archive, key: str, path: Path):
    require(key in archive, "{0} is missing {1!r}".format(path.name, key))
    try:
        value = np.asarray(archive[key])
    except (TypeError, ValueError) as error:
        raise ComparisonInputError(
            "cannot read {0!r} from {1}".format(key, path.name)
        ) from error
    require(value.shape == (), "{0!r} must be scalar".format(key))
    return value.item()


def _archive_array(archive, key: str, path: Path, shape) -> np.ndarray:
    require(key in archive, "{0} is missing {1!r}".format(path.name, key))
    try:
        value = np.asarray(archive[key], dtype=float)
    except (TypeError, ValueError) as error:
        raise ComparisonInputError(
            "{0!r} in {1} must be numeric".format(key, path.name)
        ) from error
    require(value.shape == shape, "{0!r} has the wrong shape".format(key))
    require(np.all(np.isfinite(value)), "{0!r} must be finite".format(key))
    return value.copy()


def _exact_parameter_json(archive, key: str, expected, path: Path) -> None:
    raw = _archive_scalar(archive, key, path)
    require(isinstance(raw, str), "{0!r} must be a JSON string".format(key))
    def reject_constant(value):
        raise ValueError("non-finite JSON constant {0}".format(value))

    try:
        observed = json.loads(raw, parse_constant=reject_constant)
    except (TypeError, ValueError, json.JSONDecodeError) as error:
        raise ComparisonInputError("{0!r} is invalid JSON".format(key)) from error
    require(
        isinstance(observed, dict) and set(observed) == set(expected),
        "{0!r} does not contain the exact Step-2 Case-B parameters".format(key),
    )
    for name, expected_value in expected.items():
        value = observed[name]
        require(
            isinstance(value, (int, float))
            and not isinstance(value, bool)
            and np.isfinite(float(value))
            and float(value) == expected_value,
            "{0!r} parameter {1!r} disagrees with Step-2 Case-B".format(
                key, name
            ),
        )


def _archive_json(archive, key: str, path: Path, expected_type):
    raw = _archive_scalar(archive, key, path)
    require(isinstance(raw, str), "{0!r} must be a JSON string".format(key))
    try:
        value = json.loads(
            raw,
            parse_constant=lambda constant: (_ for _ in ()).throw(
                ValueError("non-finite JSON constant {0}".format(constant))
            ),
        )
    except (TypeError, ValueError, json.JSONDecodeError) as error:
        raise ComparisonInputError("{0!r} is invalid JSON".format(key)) from error
    require(isinstance(value, expected_type), "{0!r} has the wrong JSON type".format(key))
    return value


def _validate_canonical_run_metadata(archive, path: Path) -> Dict[str, object]:
    """Reject a physically correct label on a noncanonical numerical run."""
    exact_scalars = {
        "driver": "examples/cardiac_benchmark/run_mpi.py",
        "benchmark_reproduction_profile": "paper-source-matched-full-cycle",
        "density": 1000.0,
        "material_kernel_formulation": "standard",
        "material_model_id": EXPECTED_MATERIAL_MODEL_ID,
        "mass_representation": "consistent_q1_hex8",
        "nonlinear_solver": "petsc-snes-mpi",
        "mesh_topology": "closed_multiblock_disk",
        "fiber_sampling": "gp_direct_rule",
        "fiber_sampling_option": "gp-direct",
        "tbar_definition": "laplace_presolved",
        "point_sampling": "hex8_reference_isoparametric",
        "viscous_rate": "velocity_consistent_green_lagrange_at_alpha_f_stage",
        "material_eta_pa_s": 100.0,
        "parameter_variant": "benchmark_eta",
        "a_top": 1.0e5,
        "b_top": 5.0e3,
        "a_epi": 1.0e8,
        "b_epi": 5.0e3,
        "perturb": 0.0,
        "apex_offset": 0.0,
    }
    for key, expected in exact_scalars.items():
        require(
            _archive_scalar(archive, key, path) == expected,
            "CoupFE archive has noncanonical {0!r}".format(key),
        )
    for key in ("mpi_enabled", "viscous_term_active", "flip_helix"):
        require(
            _archive_scalar(archive, key, path) is True,
            "CoupFE archive must record {0}=true".format(key),
        )
    require(
        _archive_scalar(archive, "isotropic", path) is False,
        "CoupFE archive must use anisotropic benchmark fibers",
    )
    ranks = int(_archive_scalar(archive, "mpi_ranks", path))
    world_size = int(_archive_scalar(archive, "mpi_world_size", path))
    require(
        ranks == world_size and ranks in {1, 2, 4, 8},
        "CoupFE archive has an unsupported MPI rank count",
    )
    formulation = _archive_scalar(archive, "formulation", path)
    require(
        formulation in FORMULATION_CONTRACTS,
        "CoupFE archive has an unsupported Step-2 formulation",
    )
    contract = FORMULATION_CONTRACTS[formulation]
    for key, expected in contract.items():
        require(
            _archive_scalar(archive, key, path) == expected,
            "CoupFE archive formulation metadata is inconsistent for {0!r}".format(
                key
            ),
        )
    for key in ("tbar_source_sha256", "tbar_metadata_sha256"):
        value = _archive_scalar(archive, key, path)
        require(
            isinstance(value, str) and FULL_SHA256.fullmatch(value) is not None,
            "CoupFE archive has invalid Laplace-field provenance",
        )
    require(
        _archive_scalar(archive, "tbar_metadata_schema", path)
        == "coupfe-cardiac-laplace-tbar-v1",
        "CoupFE archive has unsupported Laplace-field metadata",
    )

    archive_source = {
        key: _archive_scalar(archive, key, path)
        for key in (
            "benchmark_runtime_source_manifest_json",
            "benchmark_runtime_source_sha256",
        )
    }
    reviewed_source, _ = load_reviewed_runtime_source_manifest()
    current_source = runtime_source_metadata()
    require(
        archive_source in (reviewed_source, current_source),
        "CoupFE archive does not match the reviewed runtime source manifest "
        "or the current source tree",
    )
    uses_reviewed_source = archive_source == reviewed_source
    app_revision = _archive_scalar(archive, "app_revision", path)
    app_tree_state = _archive_scalar(archive, "app_tree_state", path)
    require(
        isinstance(app_revision, str)
        and FULL_SHA1.fullmatch(app_revision) is not None,
        "CoupFE archive has invalid application revision provenance",
    )
    require(
        _archive_scalar(archive, "app_source_kind", path) == "git-checkout"
        and app_tree_state in {"clean", "dirty"},
        "CoupFE archive has unsupported application source provenance",
    )
    if uses_reviewed_source:
        require(
            app_revision == REVIEWED_RESULT_APP_REVISION
            and app_tree_state == "dirty",
            "CoupFE archive retained runtime source does not match its "
            "reviewed application revision and tree state",
        )
    require(
        _archive_scalar(archive, "core_revision", path) == EXPECTED_CORE_REVISION
        and _archive_scalar(archive, "core_source_url", path) == EXPECTED_CORE_URL
        and _archive_scalar(archive, "core_source_kind", path)
        in {"git-checkout", "pep610-vcs"}
        and _archive_scalar(archive, "core_tree_state", path)
        in {"clean", "installed"},
        "CoupFE archive does not use the reviewed Core revision",
    )

    audit = _archive_json(archive, "pre_solve_audit_json", path, dict)
    require(
        set(audit) == {"geometry", "pressure", "robin"},
        "CoupFE archive has incomplete pre-solve audits",
    )
    expected_audit_schemas = {
        "geometry": "coupfe-cardiac-pre-solve-geometry-v1",
        "pressure": "coupfe-cardiac-pre-solve-pressure-v1",
        "robin": "coupfe-cardiac-pre-solve-robin-v1",
    }
    for name, schema in expected_audit_schemas.items():
        record = audit[name]
        require(
            isinstance(record, dict)
            and record.get("schema") == schema
            and record.get("passed") is True,
            "CoupFE archive has a failed or malformed {0} audit".format(name),
        )
    return {
        "profile": "paper-source-matched-full-cycle",
        "formulation": formulation,
        "mpi_ranks": ranks,
        "app_revision": app_revision,
        "app_tree_state": app_tree_state,
        "runtime_source_sha256": archive_source[
            "benchmark_runtime_source_sha256"
        ],
        "core_revision": EXPECTED_CORE_REVISION,
    }


def _validate_step_diagnostics(
    archive, path: Path, tau: np.ndarray, pressure: np.ndarray
) -> None:
    diagnostics = _archive_json(
        archive, "nonlinear_step_diagnostics_json", path, list
    )
    require(
        len(diagnostics) == 1000,
        "CoupFE archive must retain one diagnostic for every accepted step",
    )
    for index, record in enumerate(diagnostics, start=1):
        require(isinstance(record, dict), "CoupFE step diagnostic is malformed")
        required = {
            "time",
            "ranks",
            "final_residual_norm",
            "residual_acceptance_threshold",
            "snes_converged_reason",
            "active_tension_pa",
            "pressure_pa",
        }
        require(
            required.issubset(record),
            "CoupFE step diagnostic is missing acceptance/load evidence",
        )
        require(
            np.isclose(float(record["time"]), index * 0.001, rtol=0.0, atol=5.0e-13)
            and int(record["ranks"]) in {1, 2, 4, 8}
            and float(record["final_residual_norm"])
            <= float(record["residual_acceptance_threshold"])
            and int(record["snes_converged_reason"]) > 0
            and np.isclose(
                float(record["active_tension_pa"]), tau[index], rtol=0.0, atol=1.0e-9
            )
            and np.isclose(
                float(record["pressure_pa"]), pressure[index], rtol=0.0, atol=1.0e-9
            ),
            "CoupFE step diagnostic disagrees with accepted-state evidence",
        )


def _expected_step2b_load_histories() -> Tuple[np.ndarray, np.ndarray]:
    """Recreate the exact source-matched load histories on the 1 ms grid."""
    configuration = benchmark_configuration(2, "B")
    evaluation_times = COUPFE_TIME_S.copy()
    evaluation_times[1:] -= 0.4 * 0.001
    tau = np.zeros_like(COUPFE_TIME_S)
    pressure = np.zeros_like(COUPFE_TIME_S)
    tau[1:] = tau_of_t(
        evaluation_times[1:],
        p=configuration.activation_parameters,
        t_span=(0.0, 1.0),
    )
    pressure[1:] = p_of_t(
        evaluation_times[1:],
        p=configuration.pressure_parameters,
        t_span=(0.0, 1.0),
    )
    return tau, pressure


def load_coupfe_step2b(path: Path) -> Dict[str, object]:
    """Load one completed archive and reject Step-0/Case-B identity reuse."""
    path = Path(path).expanduser().resolve()
    payload, identity = _read_identity_payload(path)
    try:
        context = np.load(io.BytesIO(payload), allow_pickle=False)
    except (OSError, ValueError) as error:
        raise ComparisonInputError(
            "cannot load CoupFE archive {0}: {1}".format(path.name, error)
        ) from error
    try:
        with context as archive:
            require(
                _archive_scalar(archive, "result_schema", path) == RESULT_SCHEMA,
                "unsupported CoupFE result schema",
            )
            require(
                _archive_scalar(archive, "converged", path) is True,
                "CoupFE archive is not marked converged",
            )
            require(
                _archive_scalar(archive, "completed_steps", path) == 1000
                and _archive_scalar(archive, "expected_steps", path) == 1000,
                "CoupFE Step-2 Case-B archive must be complete at 1000/1000",
            )
            require(
                _archive_scalar(archive, "case", path) == "B"
                and _archive_scalar(archive, "benchmark_step", path) == 2,
                "CoupFE archive is not Benchmark 1 Step 2 Case B",
            )
            require(
                _archive_scalar(archive, "benchmark_configuration_id", path)
                == BENCHMARK_CONFIGURATION_ID,
                "CoupFE archive has the wrong benchmark configuration identity",
            )
            require(
                _archive_scalar(archive, "benchmark_load_contract", path)
                == BENCHMARK_LOAD_CONTRACT,
                "CoupFE archive does not apply active stress plus pressure",
            )
            require(
                _archive_scalar(archive, "benchmark_peak_load_definition", path)
                == BENCHMARK_PEAK_LOAD_DEFINITION,
                "CoupFE archive has the wrong combined-load peak definition",
            )
            require(
                _archive_scalar(archive, "benchmark_active_stress_enabled", path)
                is True
                and _archive_scalar(archive, "benchmark_pressure_enabled", path)
                is True,
                "CoupFE archive must enable both active stress and pressure",
            )
            require(
                float(_archive_scalar(archive, "dt", path)) == 0.001
                and float(_archive_scalar(archive, "t_end", path)) == 1.0,
                "CoupFE archive must use the published SI time setup",
            )
            require(
                _archive_scalar(archive, "integrator", path)
                == "generalized-alpha",
                "Step-2 Case-B reproduction requires generalized-alpha",
            )
            for key, expected_value in GENERALIZED_ALPHA_PARAMETERS.items():
                require(
                    float(_archive_scalar(archive, key, path)) == expected_value,
                    "CoupFE archive has the wrong source-matched {0}".format(key),
                )
            require(
                _archive_scalar(archive, "generalized_alpha_stage_contract", path)
                == GENERALIZED_ALPHA_STAGE_CONTRACT,
                "CoupFE archive has the wrong generalized-alpha stage contract",
            )
            _exact_parameter_json(
                archive,
                "benchmark_material_parameters_json",
                EXPECTED_MATERIAL,
                path,
            )
            _exact_parameter_json(
                archive,
                "benchmark_activation_parameters_json",
                EXPECTED_ACTIVATION,
                path,
            )
            _exact_parameter_json(
                archive,
                "benchmark_pressure_parameters_json",
                EXPECTED_PRESSURE,
                path,
            )
            run_contract = _validate_canonical_run_metadata(archive, path)
            times = _archive_array(archive, "times", path, (1001,))
            u0 = _archive_array(archive, "u0", path, (1001, 3))
            u1 = _archive_array(archive, "u1", path, (1001, 3))
            tau = _archive_array(archive, "tau", path, (1001,))
            pressure = _archive_array(archive, "pres", path, (1001,))
            load_evaluation_times = _archive_array(
                archive, "load_evaluation_times_s", path, (1001,)
            )
            _validate_step_diagnostics(archive, path, tau, pressure)
            for point in POINTS:
                landmark = _archive_array(archive, point, path, (3,))
                require(
                    np.allclose(landmark, LANDMARKS_M[point], rtol=0.0, atol=1.0e-14),
                    "CoupFE archive has the wrong {0} SI landmark".format(point),
                )
    except ComparisonInputError:
        raise
    except (KeyError, TypeError, ValueError) as error:
        raise ComparisonInputError(
            "malformed CoupFE archive {0}: {1}".format(path.name, error)
        ) from error

    require(
        np.allclose(times, COUPFE_TIME_S, rtol=0.0, atol=5.0e-13),
        "CoupFE archive must contain the published 0.001 s grid",
    )
    expected_load_times = COUPFE_TIME_S.copy()
    expected_load_times[1:] -= 0.4 * 0.001
    require(
        np.allclose(
            load_evaluation_times,
            expected_load_times,
            rtol=0.0,
            atol=5.0e-13,
        ),
        "CoupFE generalized-alpha loads are not staged at t_np1-alpha_f*dt",
    )
    displacement = np.stack((u0.T, u1.T), axis=0)
    require(
        float(np.max(np.abs(displacement))) < 0.20,
        "CoupFE displacement is not plausible in SI metres",
    )
    require(
        float(np.max(np.abs(displacement[:, :, 0]))) <= 1.0e-10,
        "CoupFE initial displacement must be zero in metres",
    )
    require(
        tau[0] == 0.0
        and pressure[0] == 0.0
        and np.all(tau >= 0.0)
        and np.all(pressure >= 0.0),
        "CoupFE load histories must start at zero and remain nonnegative",
    )
    expected_tau, expected_pressure = _expected_step2b_load_histories()
    require(
        np.allclose(tau, expected_tau, rtol=2.0e-6, atol=1.0e-3),
        "CoupFE active-stress history disagrees with the staged Step-2 load",
    )
    require(
        np.allclose(pressure, expected_pressure, rtol=2.0e-6, atol=1.0e-3),
        "CoupFE pressure history disagrees with the staged Step-2 load",
    )
    common = np.empty((2, 3, 101), dtype=float)
    for point_index in range(2):
        for component_index in range(3):
            common[point_index, component_index] = np.interp(
                PUBLISHED_TIME_S,
                times,
                displacement[point_index, component_index],
            )
    return {
        "identity": identity,
        "time_s": times,
        "displacement_m": displacement,
        "common_displacement_m": common,
        "run_contract": run_contract,
    }


def _window_mask(time_s: np.ndarray, start_s: float, end_s: float) -> np.ndarray:
    return (time_s >= start_s - 1.0e-12) & (time_s <= end_s + 1.0e-12)


def _relative_l2(actual: np.ndarray, reference: np.ndarray) -> Optional[float]:
    denominator = float(np.linalg.norm(reference.ravel()))
    if denominator <= np.finfo(float).tiny:
        return None
    return float(np.linalg.norm((actual - reference).ravel()) / denominator)


def paper_relative_discrepancy(
    actual: np.ndarray, reference: np.ndarray
) -> Dict[str, object]:
    """Return the benchmark paper Eq. 21 RED at p0 and p1."""
    actual = np.asarray(actual, dtype=float)
    reference = np.asarray(reference, dtype=float)
    require(
        actual.shape == (2, 3, 101) and reference.shape == (2, 3, 101),
        "paper RED inputs must be two 3-component 101-sample histories",
    )
    require(
        np.all(np.isfinite(actual)) and np.all(np.isfinite(reference)),
        "paper RED inputs must be finite",
    )
    points = {}
    for point_index, point in enumerate(POINTS):
        actual_history = actual[point_index].T
        reference_history = reference[point_index].T
        numerator = np.linalg.norm(actual_history - reference_history, axis=1)
        denominator = np.linalg.norm(reference_history, axis=1) + 1.0e-30
        points[point] = float(np.mean(numerator / denominator))
    return {
        "metric": "relative discrepancy (benchmark paper Eq. 21)",
        "definition": (
            "mean over the 101 common-grid times of the pointwise vector-error "
            "norm divided by the all-team-mean displacement-vector norm"
        ),
        "near_zero_denominator_offset": 1.0e-30,
        "points": points,
    }


def _error_metrics(actual: np.ndarray, reference: np.ndarray) -> Dict[str, object]:
    difference = actual - reference
    return {
        "rmse_m": float(np.sqrt(np.mean(difference * difference))),
        "rmse_mm": float(1.0e3 * np.sqrt(np.mean(difference * difference))),
        "relative_l2": _relative_l2(actual, reference),
        "bias_m": float(np.mean(difference)),
        "max_abs_error_m": float(np.max(np.abs(difference))),
    }


def trajectory_error_metrics(
    actual: np.ndarray, reference: np.ndarray, time_s: np.ndarray
) -> Dict[str, object]:
    """Return component, vector, and aggregate errors for named time windows."""
    result = {}
    for name, start_s, end_s in WINDOWS_S:
        mask = _window_mask(time_s, start_s, end_s)
        window = {
            "start_s": start_s,
            "end_s": end_s,
            "sample_count": int(np.count_nonzero(mask)),
            "points": {},
            "all_components": _error_metrics(actual[:, :, mask], reference[:, :, mask]),
        }
        for point_index, point in enumerate(POINTS):
            point_metrics = {"components": {}}
            for component_index, component in enumerate(COMPONENTS):
                point_metrics["components"][component] = _error_metrics(
                    actual[point_index, component_index, mask],
                    reference[point_index, component_index, mask],
                )
            # Index the point first so boolean indexing cannot move the sample
            # axis ahead of the component axis in NumPy's advanced rules.
            actual_vector = actual[point_index][:, mask].T
            reference_vector = reference[point_index][:, mask].T
            vector_difference = actual_vector - reference_vector
            point_metrics["vector"] = {
                "rmse_vector_m": float(
                    np.sqrt(np.mean(np.sum(vector_difference * vector_difference, axis=1)))
                ),
                "relative_l2_vector": _relative_l2(
                    actual_vector, reference_vector
                ),
                "max_vector_error_m": float(
                    np.max(np.linalg.norm(vector_difference, axis=1))
                ),
            }
            window["points"][point] = point_metrics
        result[name] = window
    return result


def _downward_crossings(
    time_s: np.ndarray,
    values: np.ndarray,
    threshold_m: float,
    *,
    start_s: float,
    end_s: float,
) -> Sequence[float]:
    crossings = []
    for index in range(len(time_s) - 1):
        left_time = float(time_s[index])
        right_time = float(time_s[index + 1])
        if left_time < start_s - 1.0e-12 or right_time > end_s + 1.0e-12:
            continue
        left = float(values[index])
        right = float(values[index + 1])
        if left > threshold_m and right <= threshold_m:
            fraction = (threshold_m - left) / (right - left)
            crossings.append(left_time + fraction * (right_time - left_time))
    return crossings


def _single_reference_crossing(
    time_s: np.ndarray,
    values: np.ndarray,
    threshold_m: float,
    description: str,
) -> float:
    crossings = _downward_crossings(
        time_s, values, threshold_m, start_s=0.16, end_s=0.32
    )
    require(
        len(crossings) == 1,
        "{0} must have exactly one contraction crossing".format(description),
    )
    return float(crossings[0])


def _candidate_crossing(
    time_s: np.ndarray, values: np.ndarray, threshold_m: float
) -> Dict[str, object]:
    crossings = _downward_crossings(
        time_s, values, threshold_m, start_s=0.16, end_s=0.32
    )
    return {
        "crossing_count": len(crossings),
        "crossing_time_s": float(crossings[0]) if len(crossings) == 1 else None,
    }


def _reference_summary(values: Iterable[float]) -> Dict[str, float]:
    array = np.asarray(tuple(values), dtype=float)
    return {
        "team_mean_s": float(np.mean(array)),
        "team_std_s": float(np.std(array)),
        "team_min_s": float(np.min(array)),
        "team_max_s": float(np.max(array)),
    }


def contraction_event_metrics(
    reference: Mapping[str, object], coupfe: Mapping[str, object]
) -> Dict[str, object]:
    """Measure fixed and normalized x-drop events on common and native grids."""
    teams = reference["teams_m"]
    common = coupfe["common_displacement_m"]
    native = coupfe["displacement_m"]
    native_time = coupfe["time_s"]
    report = {
        "definition": (
            "linear-interpolated downward ux crossings during 0.16--0.32 s; "
            "publisher summaries are statistics of ten per-team crossing times"
        ),
        "fixed_threshold_crossings": {},
        "normalized_drop_crossings": {},
    }
    for point_index, point in enumerate(POINTS):
        fixed = {}
        fixed_times = {}
        for threshold in FIXED_DOWNWARD_THRESHOLDS_M:
            team_times = [
                _single_reference_crossing(
                    PUBLISHED_TIME_S,
                    teams[team, point_index, 0],
                    threshold,
                    "{0} {1} m team {2}".format(point, threshold, team),
                )
                for team in range(10)
            ]
            summary = _reference_summary(team_times)
            common_event = _candidate_crossing(
                PUBLISHED_TIME_S, common[point_index, 0], threshold
            )
            native_event = _candidate_crossing(
                native_time, native[point_index, 0], threshold
            )
            common_time = common_event["crossing_time_s"]
            native_crossing_time = native_event["crossing_time_s"]
            fixed[str(int(round(threshold * 1000.0)))] = {
                "threshold_m": threshold,
                "threshold_mm": threshold * 1000.0,
                "reference": summary,
                "coupfe_common_grid": common_event,
                "coupfe_native_grid": native_event,
                "common_grid_delta_from_team_mean_s": (
                    None if common_time is None else common_time - summary["team_mean_s"]
                ),
                "simula_crossing_time_s": float(team_times[SIMULA_TEAM_INDEX]),
                "common_grid_delta_from_simula_s": (
                    None
                    if common_time is None
                    else common_time - team_times[SIMULA_TEAM_INDEX]
                ),
                "native_grid_delta_from_team_mean_s": (
                    None
                    if native_crossing_time is None
                    else native_crossing_time - summary["team_mean_s"]
                ),
                "native_grid_delta_from_simula_s": (
                    None
                    if native_crossing_time is None
                    else native_crossing_time - team_times[SIMULA_TEAM_INDEX]
                ),
            }
            fixed_times[threshold] = (summary, common_time)
        low = fixed_times[-5.0e-3]
        high = fixed_times[-15.0e-3]
        fixed["minus5_to_minus15_transition"] = {
            "reference_team_mean_duration_s": (
                high[0]["team_mean_s"] - low[0]["team_mean_s"]
            ),
            "coupfe_common_grid_duration_s": (
                None
                if low[1] is None or high[1] is None
                else high[1] - low[1]
            ),
        }
        report["fixed_threshold_crossings"][point] = fixed

        normalized = {}
        # Normalization is against zero displacement and each trajectory's most
        # negative ux through the end of the post-snap branch (0.48 s).
        team_minima = np.min(teams[:, point_index, 0, 16:49], axis=1)
        common_minimum = float(np.min(common[point_index, 0, 16:49]))
        native_branch = _window_mask(native_time, 0.16, 0.48)
        native_minimum = float(np.min(native[point_index, 0, native_branch]))
        for fraction in NORMALIZED_DROP_FRACTIONS:
            team_times = [
                _single_reference_crossing(
                    PUBLISHED_TIME_S,
                    teams[team, point_index, 0],
                    fraction * team_minima[team],
                    "{0} normalized {1:.0%} team {2}".format(
                        point, fraction, team
                    ),
                )
                for team in range(10)
            ]
            summary = _reference_summary(team_times)
            common_event = _candidate_crossing(
                PUBLISHED_TIME_S,
                common[point_index, 0],
                fraction * common_minimum,
            )
            native_event = _candidate_crossing(
                native_time,
                native[point_index, 0],
                fraction * native_minimum,
            )
            common_time = common_event["crossing_time_s"]
            native_crossing_time = native_event["crossing_time_s"]
            normalized[str(int(round(100.0 * fraction)))] = {
                "fraction_of_zero_to_post_snap_minimum": fraction,
                "reference": summary,
                "coupfe_common_grid": common_event,
                "coupfe_native_grid": native_event,
                "common_grid_delta_from_team_mean_s": (
                    None if common_time is None else common_time - summary["team_mean_s"]
                ),
                "simula_crossing_time_s": float(team_times[SIMULA_TEAM_INDEX]),
                "common_grid_delta_from_simula_s": (
                    None
                    if common_time is None
                    else common_time - team_times[SIMULA_TEAM_INDEX]
                ),
                "native_grid_delta_from_team_mean_s": (
                    None
                    if native_crossing_time is None
                    else native_crossing_time - summary["team_mean_s"]
                ),
                "native_grid_delta_from_simula_s": (
                    None
                    if native_crossing_time is None
                    else native_crossing_time - team_times[SIMULA_TEAM_INDEX]
                ),
            }
        report["normalized_drop_crossings"][point] = normalized
    return report


def _linear_slope(time_s: np.ndarray, values: np.ndarray) -> float:
    centered = time_s - float(np.mean(time_s))
    return float(np.dot(centered, values) / np.dot(centered, centered))


def _curve_features(time_s: np.ndarray, values: np.ndarray) -> Dict[str, float]:
    minimum_index = int(np.argmin(values))
    return {
        "mean_m": float(np.mean(values)),
        "range_m": float(np.ptp(values)),
        "total_variation_m": float(np.sum(np.abs(np.diff(values)))),
        "linear_slope_m_per_s": _linear_slope(time_s, values),
        "minimum_m": float(values[minimum_index]),
        "minimum_time_s": float(time_s[minimum_index]),
    }


def _branch_sign(value: float, tolerance: float = 1.0e-12) -> int:
    if value > tolerance:
        return 1
    if value < -tolerance:
        return -1
    return 0


def post_snap_branch_metrics(
    reference: Mapping[str, object], coupfe: Mapping[str, object]
) -> Dict[str, object]:
    """Describe the settled plateau without mixing in the snap transition."""
    teams = reference["teams_m"]
    mean = reference["mean_m"]
    actual = coupfe["common_displacement_m"]
    plateau = _window_mask(PUBLISHED_TIME_S, 0.32, 0.48)
    branch = plateau
    result = {
        "definition": (
            "settled post-snap branch and plateau metrics on 0.32--0.48 s; "
            "the 0.16--0.32 s transition is reported separately"
        ),
        "points": {},
    }
    for point_index, point in enumerate(POINTS):
        reference_branch = mean[point_index, 0, branch]
        actual_branch = actual[point_index, 0, branch]
        branch_time = PUBLISHED_TIME_S[branch]
        aligned_difference = (
            actual_branch - actual_branch[0]
        ) - (reference_branch - reference_branch[0])
        team_features = [
            _curve_features(branch_time, teams[team, point_index, 0, branch])
            for team in range(10)
        ]
        reference_features = _curve_features(branch_time, reference_branch)
        actual_features = _curve_features(branch_time, actual_branch)
        feature_comparison = {}
        for name in reference_features:
            team_values = np.array([entry[name] for entry in team_features])
            feature_comparison[name] = {
                "reference_team_mean": float(np.mean(team_values)),
                "reference_team_std": float(np.std(team_values)),
                "reference_mean_trajectory": reference_features[name],
                "reference_simula": team_features[SIMULA_TEAM_INDEX][name],
                "coupfe": actual_features[name],
                "delta_from_reference_mean_trajectory": (
                    actual_features[name] - reference_features[name]
                ),
                "delta_from_simula": (
                    actual_features[name]
                    - team_features[SIMULA_TEAM_INDEX][name]
                ),
            }
        levels = {}
        for sample_time in (0.32, 0.40, 0.48):
            index = int(round(sample_time / 0.01))
            values = teams[:, point_index, 0, index]
            levels["{0:.2f}".format(sample_time)] = {
                "time_s": sample_time,
                "reference_team_mean_m": float(np.mean(values)),
                "reference_team_std_m": float(np.std(values)),
                "reference_simula_m": float(values[SIMULA_TEAM_INDEX]),
                "coupfe_m": float(actual[point_index, 0, index]),
                "delta_m": float(actual[point_index, 0, index] - np.mean(values)),
            }
        plateau_time = PUBLISHED_TIME_S[plateau]
        reference_plateau = mean[point_index, 0, plateau]
        actual_plateau = actual[point_index, 0, plateau]
        component_envelopes = {}
        for component_index, component in enumerate(COMPONENTS):
            team_component = teams[:, point_index, component_index][:, plateau]
            actual_component = actual[point_index, component_index, plateau]
            lower = np.min(team_component, axis=0)
            upper = np.max(team_component, axis=0)
            outside_distance = np.maximum(lower - actual_component, 0.0) + np.maximum(
                actual_component - upper, 0.0
            )
            team_net_change = team_component[:, -1] - team_component[:, 0]
            reference_component = mean[point_index, component_index, plateau]
            reference_net_change = float(
                reference_component[-1] - reference_component[0]
            )
            coupfe_net_change = float(actual_component[-1] - actual_component[0])
            reference_sign = _branch_sign(reference_net_change)
            coupfe_sign = _branch_sign(coupfe_net_change)
            official_level_signs = sorted(
                {
                    _branch_sign(float(value))
                    for value in team_component.ravel()
                }
            )
            coupfe_level_signs = sorted(
                {_branch_sign(float(value)) for value in actual_component}
            )
            expected_level_sign = (
                official_level_signs[0]
                if len(official_level_signs) == 1
                else None
            )
            reference_increments = np.diff(reference_component)
            coupfe_increments = np.diff(actual_component)
            informative = np.abs(reference_increments) > 1.0e-12
            increment_sign_agreement = (
                float(
                    np.mean(
                        np.sign(coupfe_increments[informative])
                        == np.sign(reference_increments[informative])
                    )
                )
                if np.any(informative)
                else None
            )
            component_envelopes[component] = {
                "official_overall_min_m": float(np.min(team_component)),
                "official_overall_max_m": float(np.max(team_component)),
                "coupfe_min_m": float(np.min(actual_component)),
                "coupfe_max_m": float(np.max(actual_component)),
                "samplewise_official_envelope_coverage_fraction": float(
                    np.mean(outside_distance == 0.0)
                ),
                "outside_envelope_sample_count": int(
                    np.count_nonzero(outside_distance)
                ),
                "max_outside_envelope_distance_m": float(
                    np.max(outside_distance)
                ),
                "plateau_level_sign": {
                    "official_sample_signs": official_level_signs,
                    "official_unanimous_sign": expected_level_sign,
                    "coupfe_sample_signs": coupfe_level_signs,
                    "coupfe_mean_level_sign": _branch_sign(
                        float(np.mean(actual_component))
                    ),
                    "all_coupfe_samples_match_official_unanimous_sign": (
                        None
                        if expected_level_sign is None
                        else coupfe_level_signs == [expected_level_sign]
                    ),
                },
                "plateau_drift": {
                    "official_team_net_change_mean_m": float(
                        np.mean(team_net_change)
                    ),
                    "official_team_net_change_min_m": float(
                        np.min(team_net_change)
                    ),
                    "official_team_net_change_max_m": float(
                        np.max(team_net_change)
                    ),
                    "official_mean_trajectory_net_change_m": reference_net_change,
                    "coupfe_net_change_m": coupfe_net_change,
                    "official_mean_sign": reference_sign,
                    "coupfe_sign": coupfe_sign,
                    "net_change_sign_agrees": (
                        None
                        if reference_sign == 0
                        else coupfe_sign == reference_sign
                    ),
                    "increment_sign_agreement_fraction": increment_sign_agreement,
                },
            }
        official_team_octants = [
            [
                _branch_sign(
                    float(np.mean(teams[team, point_index, component, plateau]))
                )
                for component in range(3)
            ]
            for team in range(10)
        ]
        unique_official_octants = sorted(
            {tuple(octant) for octant in official_team_octants}
        )
        coupfe_octant = [
            _branch_sign(float(np.mean(actual[point_index, component, plateau])))
            for component in range(3)
        ]
        result["points"][point] = {
            "levels": levels,
            "branch_features": feature_comparison,
            "shape_rmse_after_0p32_level_alignment_m": float(
                np.sqrt(np.mean(aligned_difference * aligned_difference))
            ),
            "plateau_vector_octant": {
                "official_team_octants": official_team_octants,
                "official_unique_octants": [
                    list(octant) for octant in unique_official_octants
                ],
                "coupfe_octant": coupfe_octant,
                "matches_official_unanimous_octant": (
                    None
                    if len(unique_official_octants) != 1
                    else coupfe_octant == list(unique_official_octants[0])
                ),
            },
            "plateau": {
                "reference": _curve_features(plateau_time, reference_plateau),
                "coupfe": _curve_features(plateau_time, actual_plateau),
                "mean_level_bias_m": float(
                    np.mean(actual_plateau - reference_plateau)
                ),
            },
            "component_envelopes_0p32_to_0p48": component_envelopes,
        }
    return result


def relaxation_rebound_metrics(
    reference: Mapping[str, object], coupfe: Mapping[str, object]
) -> Dict[str, object]:
    """Report the post-0.484 s rebound without folding it into snap onset."""
    teams = reference["teams_m"]
    actual = coupfe["common_displacement_m"]
    native_time = coupfe["time_s"]
    native = coupfe["displacement_m"]
    result = {
        "definition": (
            "centered rebound speed at 0.51 s is [ux(0.52)-ux(0.50)]/0.02; "
            "recovery is measured from the 0.40--0.50 s ux minimum toward zero"
        ),
        "points": {},
    }
    for point_index, point in enumerate(POINTS):
        team_speeds = (
            teams[:, point_index, 0, 52] - teams[:, point_index, 0, 50]
        ) / 0.02
        common_speed = float(
            (actual[point_index, 0, 52] - actual[point_index, 0, 50]) / 0.02
        )
        native_050 = float(np.interp(0.50, native_time, native[point_index, 0]))
        native_052 = float(np.interp(0.52, native_time, native[point_index, 0]))
        team_recovery = []
        for team in range(10):
            curve = teams[team, point_index, 0]
            minimum = float(np.min(curve[40:51]))
            team_recovery.append(float((curve[58] - minimum) / (-minimum)))
        actual_minimum = float(np.min(actual[point_index, 0, 40:51]))
        common_recovery = float(
            (actual[point_index, 0, 58] - actual_minimum) / (-actual_minimum)
        ) if actual_minimum < 0.0 else None
        result["points"][point] = {
            "rebound_speed_centered_0p51_m_per_s": {
                "reference_team_mean": float(np.mean(team_speeds)),
                "reference_team_std": float(np.std(team_speeds)),
                "reference_team_min": float(np.min(team_speeds)),
                "reference_team_max": float(np.max(team_speeds)),
                "coupfe_common_grid": common_speed,
                "coupfe_native_grid": (native_052 - native_050) / 0.02,
                "common_grid_delta_from_team_mean": common_speed
                - float(np.mean(team_speeds)),
            },
            "recovery_fraction_by_0p58": {
                "reference_team_mean": float(np.mean(team_recovery)),
                "reference_team_std": float(np.std(team_recovery)),
                "coupfe_common_grid": common_recovery,
            },
        }
    return result


def _statistics_json(reference: Mapping[str, object]) -> Dict[str, object]:
    mean = reference["mean_m"]
    std = reference["std_m"]
    points = {}
    for point_index, point in enumerate(POINTS):
        points[point] = {}
        for component_index, component in enumerate(COMPONENTS):
            points[point][component] = {
                "mean_m": mean[point_index, component_index].tolist(),
                "std_m": std[point_index, component_index].tolist(),
            }
    return {
        "sample_count": 101,
        "team_count": 10,
        "time_s": PUBLISHED_TIME_S.tolist(),
        "points": points,
    }


def _comparison_curves_json(
    reference: Mapping[str, object], coupfe: Mapping[str, object]
) -> Dict[str, object]:
    """Return self-contained, reviewable curves for the public Step-2B plot."""
    teams = np.asarray(reference["teams_m"], dtype=float)
    mean = np.asarray(reference["mean_m"], dtype=float)
    candidate = np.asarray(coupfe["common_displacement_m"], dtype=float)
    require(
        teams.shape == (10, 2, 3, 101),
        "publisher comparison curves have the wrong shape",
    )
    require(
        mean.shape == (2, 3, 101) and candidate.shape == (2, 3, 101),
        "mean or CoupFE comparison curves have the wrong shape",
    )
    points = {}
    for point_index, point in enumerate(POINTS):
        points[point] = {}
        for component_index, component in enumerate(COMPONENTS):
            team_curves = teams[:, point_index, component_index]
            points[point][component] = {
                "coupfe_m": candidate[point_index, component_index].tolist(),
                "publisher_mean_m": mean[point_index, component_index].tolist(),
                "publisher_min_m": np.min(team_curves, axis=0).tolist(),
                "publisher_max_m": np.max(team_curves, axis=0).tolist(),
                "publisher_simula_m": team_curves[SIMULA_TEAM_INDEX].tolist(),
            }
    return {
        "sample_count": 101,
        "team_count": 10,
        "time_s": PUBLISHED_TIME_S.tolist(),
        "points": points,
        "meaning": (
            "derived SI-metre curves retained so the public figure can be "
            "regenerated without the CoupFE NPZ or publisher pickle files"
        ),
    }


def build_report(
    publisher_data_directory: Path,
    hash_manifest: Path,
    coupfe_run: Path,
) -> Dict[str, object]:
    reference = load_publisher_reference(publisher_data_directory, hash_manifest)
    coupfe = load_coupfe_step2b(coupfe_run)
    errors = trajectory_error_metrics(
        coupfe["common_displacement_m"], reference["mean_m"], PUBLISHED_TIME_S
    )
    simula_errors = trajectory_error_metrics(
        coupfe["common_displacement_m"],
        reference["teams_m"][SIMULA_TEAM_INDEX],
        PUBLISHED_TIME_S,
    )
    paper_red = paper_relative_discrepancy(
        coupfe["common_displacement_m"], reference["mean_m"]
    )
    return {
        "schema": REPORT_SCHEMA,
        "scope": {
            "benchmark": 1,
            "step": 2,
            "case": "B",
            "configuration_id": BENCHMARK_CONFIGURATION_ID,
            "loads": "active stress plus ventricular pressure",
            "comparison_grid": "101 samples, 0.00--1.00 s at 0.01 s",
            "displacement_unit": "m",
        },
        "inputs": {
            "publisher_hash_manifest": reference["manifest_identity"],
            "publisher_selected_files_in_upstream_order": reference["identities"],
            "excluded_generic_simvascular_alias": (
                {
                    "status": "excluded-byte-identical-P2-duplicate",
                    **reference["excluded_generic_identity"],
                }
                if reference["excluded_generic_identity"] is not None
                else {"status": "not-present"}
            ),
            "coupfe_run": coupfe["identity"],
            "coupfe_run_contract": coupfe["run_contract"],
        },
        "plumbing": {
            "status": "passed",
            "meaning": (
                "hashes, restricted decoding, schemas, SI units, completion, "
                "and exact Step-2 Case-B physical identity passed"
            ),
            "publisher_selection_source": reference["manifest_selection_source"],
            "publisher_source_doi": reference["manifest_source_doi"],
            "ambit_first_timestamp_handling": (
                "hash-pinned 0.001 s first timestamp mapped by upstream sample index "
                "to the common 0.000 s grid"
            ),
        },
        "publisher_material_attribution": {
            "source_title": OFFICIAL_SOURCE_TITLE,
            "creators": list(OFFICIAL_SOURCE_CREATORS),
            "source_doi": OFFICIAL_SOURCE_DOI,
            "license": OFFICIAL_SOURCE_LICENSE,
            "license_url": OFFICIAL_SOURCE_LICENSE_URL,
            "transformation_notice": (
                "CoupFE-Cardiac selected the ten hash-pinned Step 2 Case B "
                "publisher files, mapped the Ambit first sample to the "
                "publisher common grid by upstream sample index, and computed "
                "team mean, standard deviation, minimum/maximum envelope, "
                "selected Simula curves, event summaries, and comparison "
                "metrics. These are modified/derived data; raw publisher "
                "pickle files are not redistributed."
            ),
        },
        "reproduction": {
            "status": "quantified-no-paper-acceptance-threshold",
            "meaning": (
                "these are numerical agreement diagnostics, not a plumbing pass"
            ),
            "trajectory_errors": errors,
            "trajectory_errors_vs_named_simula": simula_errors,
            "paper_relative_discrepancy": paper_red,
            "contraction_events": contraction_event_metrics(reference, coupfe),
            "post_snap_branch": post_snap_branch_metrics(reference, coupfe),
            "relaxation_rebound": relaxation_rebound_metrics(reference, coupfe),
        },
        "publisher_mean_std": _statistics_json(reference),
        "comparison_curves": _comparison_curves_json(reference, coupfe),
    }


def _write_report(path: Path, report: Mapping[str, object]) -> None:
    destination = Path(path).expanduser()
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=".{0}.".format(destination.name),
            suffix=".tmp",
            dir=str(destination.parent),
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
            json.dump(report, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        os.replace(str(temporary), str(destination))
    except BaseException:
        if temporary is not None:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
        raise


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--publisher-data-dir",
        type=Path,
        required=True,
        help="directory containing the publisher results_time_curves/data files",
    )
    parser.add_argument(
        "--hash-manifest",
        type=Path,
        default=DEFAULT_HASH_MANIFEST,
        help="exact ten-file SHA-256 manifest",
    )
    parser.add_argument(
        "--coupfe-run", type=Path, required=True, help="completed Step-2 Case-B NPZ"
    )
    parser.add_argument(
        "--output", type=Path, help="write JSON atomically instead of stdout"
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = _parser().parse_args(argv)
    try:
        report = build_report(
            arguments.publisher_data_dir,
            arguments.hash_manifest,
            arguments.coupfe_run,
        )
        if arguments.output is None:
            json.dump(report, sys.stdout, indent=2, sort_keys=True, allow_nan=False)
            sys.stdout.write("\n")
        else:
            _write_report(arguments.output, report)
    except ComparisonInputError as error:
        print("comparison input error: {0}".format(error), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
