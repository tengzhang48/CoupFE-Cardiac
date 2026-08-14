"""Immutable physical/load definitions for Cardiac Benchmark 1.

Step 0 retains the historical public-driver split: Case A applies active
stress only and Case B applies ventricular pressure only.  The published
Benchmark 1 Step 2 Case B configuration applies *both* histories and uses the
Table 5 passive coefficients together with ``sigma_0 = 100000 Pa``.

The tuples stored by :class:`BenchmarkConfiguration` are the source of truth.
Compatibility dictionaries exported by ``material.py`` and ``activation.py``
are copies, so mutating one of those legacy names cannot alter a selected
benchmark configuration.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from types import MappingProxyType

import numpy as np


MATERIAL_PARAMETER_NAMES = (
    "a",
    "b",
    "a_f",
    "b_f",
    "a_s",
    "b_s",
    "a_fs",
    "b_fs",
    "kappa",
    "k_sw",
    "eta",
    "Ta",
)
ACTIVATION_PARAMETER_NAMES = (
    "t_sys",
    "t_dias",
    "gamma",
    "a_max",
    "a_min",
    "sigma_0",
)
PRESSURE_PARAMETER_NAMES = (
    "t_sys_pre",
    "t_dias_pre",
    "gamma",
    "a_max",
    "a_min",
    "alpha_pre",
    "alpha_mid",
    "sigma_pre",
    "sigma_mid",
)

# Exact application-owned inputs that can change a cardiac result. CoupFE Core
# is versioned separately by the driver; this manifest closes the provenance
# gap for an intentionally uncommitted application worktree on a remote VM.
RUNTIME_SOURCE_FILES = (
    "examples/cardiac_benchmark/activation.py",
    "examples/cardiac_benchmark/benchmark_parameters.py",
    "examples/cardiac_benchmark/boundary_audit.py",
    "examples/cardiac_benchmark/consistent_mass.py",
    "examples/cardiac_benchmark/distributed_local_pressure.py",
    "examples/cardiac_benchmark/distributed_mass.py",
    "examples/cardiac_benchmark/distributed_material.py",
    "examples/cardiac_benchmark/distributed_solver.py",
    "examples/cardiac_benchmark/generalized_alpha.py",
    "examples/cardiac_benchmark/geometry.py",
    "examples/cardiac_benchmark/local_pressure.py",
    "examples/cardiac_benchmark/material.py",
    "examples/cardiac_benchmark/pressure.py",
    "examples/cardiac_benchmark/result_io.py",
    "examples/cardiac_benchmark/robin.py",
    "examples/cardiac_benchmark/run.py",
    "examples/cardiac_benchmark/run_mpi.py",
    "examples/cardiac_benchmark/sampling.py",
    "examples/cardiac_benchmark/structural_directions.py",
    "examples/cardiac_benchmark/tbar_laplace.py",
)


STEP0_MATERIAL_PARAMETER_ITEMS = (
    ("a", 59.0),
    ("b", 8.023),
    ("a_f", 18472.0),
    ("b_f", 16.026),
    ("a_s", 2481.0),
    ("b_s", 11.12),
    ("a_fs", 216.0),
    ("b_fs", 11.436),
    ("kappa", 1.0e6),
    ("k_sw", 100.0),
    ("eta", 100.0),
    ("Ta", 0.0),
)
STEP2_CASE_B_MATERIAL_PARAMETER_ITEMS = (
    ("a", 295.0),
    ("b", 8.023),
    ("a_f", 92360.0),
    ("b_f", 16.026),
    ("a_s", 12405.0),
    ("b_s", 11.12),
    ("a_fs", 1080.0),
    ("b_fs", 11.436),
    ("kappa", 1.0e6),
    ("k_sw", 100.0),
    ("eta", 100.0),
    ("Ta", 0.0),
)

STEP0_ACTIVATION_PARAMETER_ITEMS = (
    ("t_sys", 0.16),
    ("t_dias", 0.484),
    ("gamma", 0.005),
    ("a_max", 5.0),
    ("a_min", -30.0),
    ("sigma_0", 1.5e5),
)
STEP2_CASE_B_ACTIVATION_PARAMETER_ITEMS = (
    ("t_sys", 0.16),
    ("t_dias", 0.484),
    ("gamma", 0.005),
    ("a_max", 5.0),
    ("a_min", -30.0),
    ("sigma_0", 1.0e5),
)

PRESSURE_PARAMETER_ITEMS = (
    ("t_sys_pre", 0.17),
    ("t_dias_pre", 0.484),
    ("gamma", 0.005),
    ("a_max", 5.0),
    ("a_min", -30.0),
    ("alpha_pre", 5.0),
    ("alpha_mid", 1.0),
    ("sigma_pre", 7000.0),
    ("sigma_mid", 16000.0),
)


@dataclass(frozen=True)
class BenchmarkConfiguration:
    """Complete immutable identity for one supported Benchmark 1 load case."""

    benchmark_step: int
    case: str
    identity: str
    material_parameter_items: tuple[tuple[str, float], ...]
    activation_parameter_items: tuple[tuple[str, float], ...]
    pressure_parameter_items: tuple[tuple[str, float], ...]
    active_stress_enabled: bool
    pressure_enabled: bool
    load_contract: str

    @property
    def material_parameters(self):
        return MappingProxyType(dict(self.material_parameter_items))

    @property
    def activation_parameters(self):
        return MappingProxyType(dict(self.activation_parameter_items))

    @property
    def pressure_parameters(self):
        return MappingProxyType(dict(self.pressure_parameter_items))


_STEP0_CASE_A = BenchmarkConfiguration(
    benchmark_step=0,
    case="A",
    identity="benchmark-1-step-0-case-A-active-stress-only",
    material_parameter_items=STEP0_MATERIAL_PARAMETER_ITEMS,
    activation_parameter_items=STEP0_ACTIVATION_PARAMETER_ITEMS,
    pressure_parameter_items=PRESSURE_PARAMETER_ITEMS,
    active_stress_enabled=True,
    pressure_enabled=False,
    load_contract="active-stress-only",
)
_STEP0_CASE_B = BenchmarkConfiguration(
    benchmark_step=0,
    case="B",
    identity="benchmark-1-step-0-case-B-pressure-only",
    material_parameter_items=STEP0_MATERIAL_PARAMETER_ITEMS,
    activation_parameter_items=STEP0_ACTIVATION_PARAMETER_ITEMS,
    pressure_parameter_items=PRESSURE_PARAMETER_ITEMS,
    active_stress_enabled=False,
    pressure_enabled=True,
    load_contract="pressure-only",
)
_STEP2_CASE_B = BenchmarkConfiguration(
    benchmark_step=2,
    case="B",
    identity="benchmark-1-step-2-case-B-active-stress-plus-pressure",
    material_parameter_items=STEP2_CASE_B_MATERIAL_PARAMETER_ITEMS,
    activation_parameter_items=STEP2_CASE_B_ACTIVATION_PARAMETER_ITEMS,
    pressure_parameter_items=PRESSURE_PARAMETER_ITEMS,
    active_stress_enabled=True,
    pressure_enabled=True,
    load_contract="active-stress-plus-pressure",
)

BENCHMARK_CONFIGURATIONS = MappingProxyType(
    {
        (0, "A"): _STEP0_CASE_A,
        (0, "B"): _STEP0_CASE_B,
        (2, "B"): _STEP2_CASE_B,
    }
)


def benchmark_configuration(benchmark_step, case):
    """Return a supported configuration or reject ambiguous combinations."""
    if isinstance(benchmark_step, bool) or not isinstance(benchmark_step, int):
        raise ValueError("benchmark step must be the integer 0 or 2")
    if case not in {"A", "B"}:
        raise ValueError("case must be 'A' or 'B'")
    if benchmark_step == 2 and case != "B":
        raise ValueError(
            "Benchmark 1 Step 2 is currently implemented only for Case B"
        )
    try:
        return BENCHMARK_CONFIGURATIONS[(benchmark_step, case)]
    except KeyError as error:
        raise ValueError(
            "supported Benchmark 1 modes are Step 0 Cases A/B and Step 2 Case B"
        ) from error


def _validated_parameter_dict(actual, expected_items, label):
    try:
        observed = dict(actual)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"malformed {label} parameter mapping") from error
    expected = dict(expected_items)
    if set(observed) != set(expected):
        raise RuntimeError(
            f"{label} parameter names disagree with the selected benchmark mode"
        )
    for name, expected_value in expected.items():
        value = observed[name]
        if isinstance(value, bool):
            raise RuntimeError(
                f"{label} parameter {name!r} disagrees with the selected benchmark mode"
            )
        try:
            numeric = float(value)
        except (TypeError, ValueError) as error:
            raise RuntimeError(
                f"{label} parameter {name!r} is not numeric"
            ) from error
        if numeric != expected_value:
            raise RuntimeError(
                f"{label} parameter {name!r} disagrees with the selected benchmark mode"
            )
        observed[name] = numeric
    return {name: observed[name] for name in expected}


def benchmark_metadata(
    configuration,
    *,
    material_parameters,
    activation_parameters,
    pressure_parameters,
):
    """Validate the applied mappings and return fail-closed archive fields."""
    selected = benchmark_configuration(
        configuration.benchmark_step, configuration.case
    )
    if configuration != selected:
        raise RuntimeError("benchmark configuration identity was modified")
    material = _validated_parameter_dict(
        material_parameters,
        configuration.material_parameter_items,
        "material",
    )
    activation = _validated_parameter_dict(
        activation_parameters,
        configuration.activation_parameter_items,
        "activation",
    )
    pressure = _validated_parameter_dict(
        pressure_parameters,
        configuration.pressure_parameter_items,
        "pressure",
    )
    compact = dict(sort_keys=True, separators=(",", ":"), allow_nan=False)
    return {
        "benchmark_step": configuration.benchmark_step,
        "benchmark_configuration_id": configuration.identity,
        "benchmark_identity_scope": (
            "physical-mode-defined-by-benchmark-step-and-configuration-id;"
            "mpi-implementation-is-numerical-provenance"
        ),
        "benchmark_load_contract": configuration.load_contract,
        "benchmark_peak_load_definition": (
            "argmax(abs(active_tension_pa+pressure_pa))"
        ),
        "benchmark_active_stress_enabled": configuration.active_stress_enabled,
        "benchmark_pressure_enabled": configuration.pressure_enabled,
        "benchmark_material_parameters_json": json.dumps(material, **compact),
        "benchmark_activation_parameters_json": json.dumps(activation, **compact),
        "benchmark_pressure_parameters_json": json.dumps(pressure, **compact),
        **runtime_source_metadata(),
    }


def runtime_source_metadata(repository_root=None):
    """Return a canonical content identity for result-producing app sources."""
    root = (
        Path(repository_root).expanduser().resolve()
        if repository_root is not None
        else Path(__file__).resolve().parents[2]
    )
    manifest = {}
    for relative in RUNTIME_SOURCE_FILES:
        path = root / relative
        try:
            payload = path.read_bytes()
        except OSError as error:
            raise RuntimeError(
                f"cannot identify benchmark runtime source {relative!r}"
            ) from error
        manifest[relative] = hashlib.sha256(payload).hexdigest()
    encoded = json.dumps(
        manifest, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    return {
        "benchmark_runtime_source_manifest_json": encoded,
        "benchmark_runtime_source_sha256": hashlib.sha256(
            encoded.encode("utf-8")
        ).hexdigest(),
    }


def validate_runtime_material_properties(
    configuration,
    property_names,
    property_values,
    *,
    condensed_local_pressure=False,
    active_tension_pa=0.0,
):
    """Reject a compiled-element property vector that drifts from its mode."""
    names = tuple(property_names)
    values = np.asarray(property_values, dtype=float)
    if values.shape != (len(names),) or len(set(names)) != len(names):
        raise RuntimeError("malformed compiled material property vector")
    # Shape and uniqueness were checked above, so plain zip is equally strict
    # here and keeps the declared Python 3.9 support.
    observed = dict(zip(names, values))
    expected = dict(configuration.material_parameters)
    expected["Ta"] = float(active_tension_pa)
    if condensed_local_pressure:
        expected["kappa"] = 0.0
    if set(observed) != set(expected):
        raise RuntimeError(
            "compiled material property names disagree with the benchmark mode"
        )
    for name, expected_value in expected.items():
        if not np.isfinite(observed[name]) or observed[name] != expected_value:
            raise RuntimeError(
                f"compiled material property {name!r} disagrees with the benchmark mode"
            )


def validate_load_histories(configuration, times, active_tension, pressure):
    """Validate load presence/absence and return the scalar peak-load proxy."""
    times = np.asarray(times, dtype=float)
    tau = np.asarray(active_tension, dtype=float)
    pres = np.asarray(pressure, dtype=float)
    if (
        times.ndim != 1
        or len(times) < 2
        or tau.shape != times.shape
        or pres.shape != times.shape
        or not np.all(np.isfinite(times))
        or not np.all(np.isfinite(tau))
        or not np.all(np.isfinite(pres))
    ):
        raise RuntimeError("malformed benchmark load histories")
    if tau[0] != 0.0 or pres[0] != 0.0:
        raise RuntimeError("benchmark load histories must start from zero")
    if not configuration.active_stress_enabled and np.any(tau != 0.0):
        raise RuntimeError("active stress is nonzero for a pressure-only mode")
    if not configuration.pressure_enabled and np.any(pres != 0.0):
        raise RuntimeError("pressure is nonzero for an active-stress-only mode")
    if np.any(tau < 0.0) or np.any(pres < 0.0):
        raise RuntimeError("benchmark load histories must be nonnegative")
    return tau + pres
