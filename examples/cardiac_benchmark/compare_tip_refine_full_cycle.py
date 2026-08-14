"""Step 0B two-/four-layer full-cycle comparison: JSON + accessible SVG.

Regenerates the retained figure
``docs/figures/step0b_tip_refine_full_cycle.svg`` and its bound metrics
report ``results/step0b_tip6p0_full_cycle_comparison.report.json`` from
the retained archives.  All inputs are required CLI arguments; no
machine-local path is written into any output.

The FEniCS record spans t = 0.001...0.999 s while both CoupFE archives
retain t = 0.000...1.000 s, so every FEniCS comparison uses matched
timestamps.  A separately retained four-layer 0.32 s prefix is required
as a continuity gate for the four-layer full cycle.  The figure is a
comparison rendering, not a validation or pass claim.
"""

from __future__ import annotations

import argparse
from collections.abc import Mapping
import hashlib
import io
import json
import math
import os
from pathlib import Path
import pickle
import tempfile

import numpy as np

try:  # package import
    from . import post
except ImportError:  # direct script execution
    import post

SCHEMA = "coupfe-cardiac-step0b-tip6p0-full-cycle-comparison-v2"
REPORT_NAME = "step0b_tip6p0_full_cycle_comparison.report.json"
FIGURE_TITLE = "Benchmark 1, Step 0 Case B - full-cycle displacement comparison"
FIGURE_SUBTITLE = (
    "Two and four wall layers, displacement in mm, t=0-1 s; published "
    "ten-team min-max envelope; FEniCS samples span t=0.001-0.999 s"
)
FIGURE_DESCRIPTION = (
    "Six line charts compare CoupFE-Cardiac 2x20x17 and 4x20x17 "
    "tip_refine=6.0 full-cycle trajectories with the retained local FEniCS "
    "Step 0B curves and the published ten-team envelope from zero to one "
    "second. Rows show p0 and p1; columns show x, y, and z displacement in "
    "millimetres. CoupFE-FEniCS comparisons use matched timestamps "
    "(FEniCS spans 0.001-0.999 s). The chart is not a validation or pass "
    "claim."
)
VISIBLE_LABELS = (
    "closed t/core/radial 2x20x17 and 4x20x17",
    "Hex8 Q1/P0 local pressure, tip_refine=6.0",
    "consistent-mass generalized-alpha",
    f"Report: {REPORT_NAME}",
    "benchmark DOI 10.5281/zenodo.14260459",
    "CC-BY-4.0",
)

EXPECTED_FENICS_SHA256 = {
    "parameters.json": "c1cd4c8d2521fd6c28774975843740a8af12568edd1240f5daa133d469e6fb76",
    "time_stamps.npy": "ddba330b1c8f8c1bb61282e187047f3aa99d0df37b2c4ed2139ea1b0e0ff0f0c",
    "componentwise_displacement_up0.npy": "4344a4f599a6eabb16159682339a735bff572eaa18eedd1fe2a97ebd3ee7f4a0",
    "componentwise_displacement_up1.npy": "88a679de2189bc137de5d64186c698f1702e9df20333849377cdfd01aac8bf1e",
}
EXPECTED_COUPFE_TWO_LAYER_SHA256 = (
    "5bb152c47b693af1dc2c0d650dde8b07ba28ef51d44e45db8c493cfe9a339375"
)
EXPECTED_COUPFE_FOUR_LAYER_SHA256 = (
    "1e333b29b05f01dedce9272b32b82ce6ccfda56036c1c8f57eb395b9b4494800"
)
EXPECTED_COUPFE_FOUR_LAYER_PREFIX_SHA256 = (
    "774a7dc5fc970bb744ff0188f0f428cff54c532b401a365643f4b626584d7acf"
)
EXPECTED_TBAR_IDENTITY_BY_NT = {
    2: {
        "tbar_source_filename": "tbar_nt2_core20_rad17_tip6p0.npy",
        "tbar_source_sha256": (
            "d848b6cafa0e74c6e8cf56ddd825e8b3d8c91fd8490237ed8b204539f7d3cbeb"
        ),
        "tbar_metadata_filename": "tbar_nt2_core20_rad17_tip6p0.meta.json",
        "tbar_metadata_sha256": (
            "9ae52c05c12c9c8d9c5ef5659e156c2eb573bc50dbb9895f417e785da52166b6"
        ),
        "tbar_metadata_schema": "coupfe-cardiac-laplace-tbar-v1",
    },
    4: {
        "tbar_source_filename": "tbar_nt4_core20_rad17_tip6p0.npy",
        "tbar_source_sha256": (
            "1578362593495b6fe48d6a2fd2e1332150121be4d6b361915d04f3980d78da8f"
        ),
        "tbar_metadata_filename": "tbar_nt4_core20_rad17_tip6p0.meta.json",
        "tbar_metadata_sha256": (
            "ecdadd335e41922eab459f4d0d6a17cf7fd4a3add496ddc0b179ca9f25daceeb"
        ),
        "tbar_metadata_schema": "coupfe-cardiac-laplace-tbar-v1",
    },
}
# Compatibility alias for code that imported the v1 two-layer constant.
EXPECTED_COUPFE_SHA256 = EXPECTED_COUPFE_TWO_LAYER_SHA256
TEAM_FILE_SHA256 = {
    "4C": "d066fdc316d92c5e76c2da08aeb8612b9c3f7d5fc2a6f364a970adaeed505a99",
    "ambit": "e9d287ffad46fd4a96f7b0f9f6187067d1ff7e392dbf188368ed9ce857e0895c",
    "carpentry": "5a223507bdd54daf0d775a8c14b07ab26933923063b171cd16daf61625b2a7cf",
    "cheart": "08460c6708673b57f07ad1475db88381da9608b4ef5f34248c2558a3bec38da6",
    "chimera": "9062c2204d21b6bc3146711f061cd82af5be19543f4fcce2ff2cfde351b1cbcd",
    "comsol": "80121d17056e46cf7a2419fe25f70c14afac7f1599e2de08db33124582191c88",
    "lifex": "f7e72fae62055f2b202196e18f5e8527fc4476383d17b376f6c105f13e7aa4f4",
    "simula": "4dd3aace9580dccd803e0e70bb74524bc061015b8dfc27274a00f747369f16c4",
    "simvascular_p1p1": "d82772ff96ce900cf45cc8cca50985f5bb9042476742637adc77765b363db9e3",
    "simvascular_p2": "1ca6f64f7eb63f1d7aef2023211c9f14d2d7fdfc4c9d232fc914d7a3c3ed6cf5",
}
TEAM_FILENAME_PREFIX = "monoventricular_nonblinded_step_0B_group_"
TEAM_ALIAS = "simvascular"
TEAM_ALIAS_TARGET = "simvascular_p2"
REFERENCE_DOI = "10.5281/zenodo.14260459"
REFERENCE_LICENSE = "CC-BY-4.0"
EXPECTED_TIMES = np.arange(0, 1001, dtype=float) * 1.0e-3
EXPECTED_PREFIX_TIMES = np.arange(0, 321, dtype=float) * 1.0e-3
EXPECTED_FENICS_TIMES = np.arange(1, 1000, dtype=float) * 1.0e-3
PREFIX_CONTINUITY_TOLERANCE_MM = 1.0e-9
PHASES = {
    "snap_window": (0.200, 0.320),
    "peak": (0.350, 0.484),
    "relaxation": (0.484, 1.000),
    "late_relaxation": (0.750, 1.000),
}
FENICS_WINDOWS = {
    "full_shared": (0.001, 0.999),
    "pre_snap": (0.001, 0.199),
    "snap_window": (0.200, 0.320),
    "post_snap_shared": (0.321, 0.999),
    "published_peak": (0.350, 0.484),
    "relaxation_shared": (0.484, 0.999),
    "late_relaxation": (0.750, 0.999),
}
CROSS_MESH_WINDOWS = {
    "full_cycle": (0.000, 1.000),
    "pre_snap": (0.001, 0.199),
    "snap_window": (0.200, 0.320),
    "post_snap_shared": (0.321, 0.999),
    "published_peak": (0.350, 0.484),
    "relaxation_shared": (0.484, 0.999),
    "late_relaxation": (0.750, 0.999),
}
TRANSIENT_WINDOWS = ((0.200, 0.320), (0.484, 0.750))
SEMANTIC_CONTRACT_FIELDS = (
    "benchmark_activation_parameters_json",
    "benchmark_identity_scope",
    "benchmark_load_contract",
    "benchmark_material_parameters_json",
    "benchmark_peak_load_definition",
    "benchmark_pressure_parameters_json",
    "benchmark_reproduction_profile",
    "fiber_direction_reconstruction",
    "generalized_alpha_stage_contract",
    "local_pressure_volume_law",
    "material_kernel_formulation",
    "material_model_id",
    "parameter_variant",
    "point_sampling",
    "tbar_definition",
    "viscous_rate",
)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _scalar(archive, key: str):
    _require(key in archive, f"CoupFE archive is missing {key!r}")
    value = np.asarray(archive[key])
    _require(value.shape == (), f"CoupFE archive field {key!r} is not scalar")
    return value.item()


def _same_scalar(actual, expected) -> bool:
    if isinstance(expected, float):
        try:
            return math.isclose(float(actual), expected, rel_tol=0.0, abs_tol=1.0e-12)
        except (TypeError, ValueError):
            return False
    return actual == expected


def _load_coupfe(
    coupfe_npz: Path,
    *,
    role: str,
    expected_hash: str,
    n_t: int,
    completed_steps: int,
    t_end: float,
):
    _require(
        coupfe_npz.is_file(),
        f"CoupFE {role} archive does not exist: {coupfe_npz}",
    )
    payload = coupfe_npz.read_bytes()
    observed_hash = _sha256_bytes(payload)
    _require(
        observed_hash == expected_hash,
        f"CoupFE {role} archive has SHA-256 {observed_hash}, expected "
        f"{expected_hash}",
    )
    with np.load(io.BytesIO(payload), allow_pickle=False) as archive:
        result = post._load_validated_result(
            archive, coupfe_npz, requested_case="step_0B"
        )
        expected_scalars = {
            "result_schema": "coupfe-cardiac-result-v1",
            "completed_steps": completed_steps,
            "expected_steps": completed_steps,
            "converged": True,
            "case": "B",
            "benchmark_step": 0,
            "benchmark_configuration_id": "benchmark-1-step-0-case-B-pressure-only",
            "benchmark_active_stress_enabled": False,
            "benchmark_pressure_enabled": True,
            "benchmark_load_contract": "pressure-only",
            "benchmark_reproduction_profile": "not-applicable",
            "mpi_enabled": True,
            "mpi_ranks": 8,
            "mpi_linear_solver_profile": "fgmres-gamg-rigid-rebuild",
            "integrator": "generalized-alpha",
            "generalized_alpha_alpha_m": 0.2,
            "generalized_alpha_alpha_f": 0.4,
            "generalized_alpha_gamma": 0.7,
            "generalized_alpha_beta": 0.36,
            "generalized_alpha_stage_contract": "simula-source-matched-v1",
            "formulation": "hex8_local_pressure_p0_condensed_logj",
            "material_kappa_pa": 0.0,
            "local_pressure_bulk_modulus_pa": 1.0e6,
            "local_pressure_volume_law": "linear-reference-volume-mean-log-j-v1",
            "mass_representation": "consistent_q1_hex8",
            "material_eta_pa_s": 100.0,
            "material_kernel_formulation": "standard",
            "material_model_id": (
                "holzapfel-ogden-smooth-switch-complete-energy-derivative-v1"
            ),
            "parameter_variant": "benchmark_eta",
            "density": 1000.0,
            "mesh_topology": "closed_multiblock_disk",
            "a_top": 1.0e5,
            "b_top": 5.0e3,
            "a_epi": 1.0e8,
            "b_epi": 5.0e3,
            "core_half_width": 0.36,
            "apex_offset": 0.0,
            "perturb": 0.0,
            "flip_helix": True,
            "isotropic": False,
            "n_t": n_t,
            "n_core": 20,
            "n_radial": 17,
            "tip_refine": 6.0,
            "dt": 1.0e-3,
            "t_end": t_end,
            "load_horizon": 1.0,
            "fiber_sampling_option": "gp-direct",
            "fiber_sampling": "gp_direct_rule",
            "fiber_direction_reconstruction": "toolkit-physical-coordinate-u-v-v1",
            "point_sampling": "hex8_reference_isoparametric",
            "tbar_definition": "laplace_presolved",
            "viscous_rate": (
                "velocity_consistent_green_lagrange_at_alpha_f_stage"
            ),
            "viscous_term_active": True,
            "element_evaluation_mode": "joint",
            "app_tree_state": "clean",
            "core_tree_state": "clean",
        }
        expected_scalars.update(EXPECTED_TBAR_IDENTITY_BY_NT[n_t])
        for key, expected in expected_scalars.items():
            actual = _scalar(archive, key)
            _require(
                _same_scalar(actual, expected),
                f"CoupFE {role} archive field {key!r} is {actual!r}, "
                f"expected {expected!r}",
            )

        expected_points = {
            "p0": np.array([0.025, 0.030, 0.0]),
            "p1": np.array([0.000, 0.030, 0.0]),
        }
        for key, expected in expected_points.items():
            _require(
                key in archive,
                f"CoupFE {role} archive is missing {key!r}",
            )
            actual = np.asarray(archive[key], dtype=float)
            _require(
                actual.shape == (3,)
                and np.allclose(actual, expected, rtol=0.0, atol=1.0e-15),
                f"CoupFE {role} archive has the wrong physical landmark {key}",
            )

        times = np.asarray(result["times"], dtype=float)
        expected_times = (
            EXPECTED_TIMES if completed_steps == 1000 else EXPECTED_PREFIX_TIMES
        )
        _require(
            times.shape == expected_times.shape
            and np.allclose(times, expected_times, rtol=0.0, atol=1.0e-15),
            f"CoupFE {role} archive time grid is not the exact 1 ms grid "
            f"through t={t_end:.2f} s",
        )
        ours = {}
        for label, key in (("p0", "u0"), ("p1", "u1")):
            values = np.asarray(result["histories"][key], dtype=float)
            _require(
                values.shape == (len(times), 3) and np.all(np.isfinite(values)),
                f"CoupFE {role} archive {key!r} is not a finite "
                f"({len(times)}, 3) history",
            )
            ours[label] = values.copy()
        det_f = np.asarray(result.get("det_f_gauss_peak"), dtype=float)
        _require(
            det_f.ndim == 2
            and det_f.shape[1] == 8
            and np.all(np.isfinite(det_f))
            and np.all(det_f > 0.0),
            f"CoupFE {role} archive lacks positive finite 8-GP det(F) evidence",
        )
        solver_diagnostics = result.get("solver_diagnostics")
        _require(
            isinstance(solver_diagnostics, list)
            and len(solver_diagnostics) == completed_steps,
            f"CoupFE {role} archive lacks one solver diagnostic per step",
        )
        _require(
            all(
                "function_domain_rejections" in record
                and "nonlinear_iterations" in record
                for record in solver_diagnostics
            ),
            f"CoupFE {role} archive has incomplete solver diagnostics",
        )
        rejection_counts = [
            int(record["function_domain_rejections"])
            for record in solver_diagnostics
        ]
        _require(
            min(rejection_counts) >= 0,
            f"CoupFE {role} archive has invalid domain-rejection diagnostics",
        )
        domain_rejections = sum(rejection_counts)
        _require(
            domain_rejections == 0,
            f"CoupFE {role} archive records {domain_rejections} domain rejections",
        )
        nonlinear_iterations = [
            int(record["nonlinear_iterations"])
            for record in solver_diagnostics
        ]
        _require(
            min(nonlinear_iterations) >= 0,
            f"CoupFE {role} archive has incomplete nonlinear diagnostics",
        )
        selected_state_index = int(result.get("n_peak", -1))
        _require(
            0 <= selected_state_index < len(times),
            f"CoupFE {role} archive has an invalid selected-state index",
        )
        validation_evidence = {
            "det_f_saved_selected_state_min": float(det_f.min()),
            "det_f_saved_selected_state_max": float(det_f.max()),
            "selected_state_index": selected_state_index,
            "selected_state_time_s": float(times[selected_state_index]),
            "selection_note": (
                "The archive stores det(F) at its selected peak-load state; for a "
                "truncated prefix this can be the final retained state."
            ),
            "function_domain_rejections": domain_rejections,
            "nonlinear_iterations_min": min(nonlinear_iterations),
            "nonlinear_iterations_max": max(nonlinear_iterations),
        }
        _require(
            "nodes" in archive and "elems" in archive,
            f"CoupFE {role} lacks mesh continuity arrays",
        )
        nodes = np.ascontiguousarray(np.asarray(archive["nodes"]))
        elements = np.ascontiguousarray(np.asarray(archive["elems"]))
        _require(
            nodes.ndim == 2
            and nodes.shape[1] == 3
            and np.issubdtype(nodes.dtype, np.floating)
            and elements.ndim == 2
            and np.issubdtype(elements.dtype, np.integer)
            and np.all(np.isfinite(nodes)),
            f"CoupFE {role} has invalid mesh continuity arrays",
        )
        sampling_metadata = result.get("sampling_metadata")
        _require(
            isinstance(sampling_metadata, Mapping),
            f"CoupFE {role} lacks sampling metadata",
        )
        sampling_payload = json.dumps(
            sampling_metadata,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ).encode("utf-8")
        mesh_sampling_identity = {
            "nodes_shape": list(nodes.shape),
            "nodes_dtype": nodes.dtype.str,
            "nodes_sha256": _sha256_bytes(nodes.tobytes()),
            "elements_shape": list(elements.shape),
            "elements_dtype": elements.dtype.str,
            "elements_sha256": _sha256_bytes(elements.tobytes()),
            "sampling_metadata_sha256": _sha256_bytes(sampling_payload),
        }
        runtime_versions = dict(result.get("runtime_versions") or {})
        solver_configuration = result.get("solver_configuration")
        _require(
            runtime_versions
            and isinstance(solver_configuration, Mapping)
            and "petsc_version" in solver_configuration
            and "petsc4py_version" in solver_configuration,
            f"CoupFE {role} lacks runtime-version evidence",
        )
        environment_identity = {
            "runtime_versions": runtime_versions,
            "petsc_version": str(solver_configuration["petsc_version"]),
            "petsc4py_version": str(solver_configuration["petsc4py_version"]),
        }
        identity = {
            "app_revision": str(_scalar(archive, "app_revision")),
            "app_tree_state": str(_scalar(archive, "app_tree_state")),
            "core_revision": str(_scalar(archive, "core_revision")),
            "core_tree_state": str(_scalar(archive, "core_tree_state")),
            "runtime_source_sha256": str(
                _scalar(archive, "benchmark_runtime_source_sha256")
            ),
        }
        semantic_contract = {
            key: _scalar(archive, key) for key in SEMANTIC_CONTRACT_FIELDS
        }
    archive_identity = {
        "role": role,
        "name": coupfe_npz.name,
        "sha256": observed_hash,
        "size_bytes": len(payload),
        "source_identity": identity,
        "semantic_contract": semantic_contract,
        "mesh": {
            "n_t": n_t,
            "n_core": 20,
            "n_radial": 17,
            "tip_refine": 6.0,
        },
        "tbar_identity": {
            key: expected
            for key, expected in EXPECTED_TBAR_IDENTITY_BY_NT[n_t].items()
        },
        "completed_steps": completed_steps,
        "t_end_s": t_end,
        "validation_evidence": validation_evidence,
        "mesh_sampling_identity": mesh_sampling_identity,
        "environment_identity": environment_identity,
    }
    return times, ours, archive_identity


def _load_fenics(fenics_dir: Path):
    identities = {}
    payloads = {}
    for name, expected_hash in EXPECTED_FENICS_SHA256.items():
        path = fenics_dir / name
        _require(path.is_file(), f"FEniCS input does not exist: {path}")
        payload = path.read_bytes()
        observed_hash = _sha256_bytes(payload)
        _require(
            observed_hash == expected_hash,
            f"FEniCS input {name} has SHA-256 {observed_hash}, expected "
            f"{expected_hash}",
        )
        identities[name] = {
            "name": name,
            "sha256": observed_hash,
            "size_bytes": len(payload),
        }
        payloads[name] = payload
    times = np.asarray(
        np.load(io.BytesIO(payloads["time_stamps.npy"]), allow_pickle=False),
        dtype=float,
    )
    p0 = np.asarray(
        np.load(
            io.BytesIO(payloads["componentwise_displacement_up0.npy"]),
            allow_pickle=False,
        ),
        dtype=float,
    )
    p1 = np.asarray(
        np.load(
            io.BytesIO(payloads["componentwise_displacement_up1.npy"]),
            allow_pickle=False,
        ),
        dtype=float,
    )
    _require(
        times.shape == EXPECTED_FENICS_TIMES.shape
        and np.allclose(times, EXPECTED_FENICS_TIMES, rtol=0.0, atol=1.0e-15),
        "FEniCS time grid is not t=0.001...0.999 s at dt=0.001 s",
    )
    for label, values in (("p0", p0), ("p1", p1)):
        _require(
            values.shape == (len(times), 3) and np.all(np.isfinite(values)),
            f"FEniCS {label} is not a finite (999, 3) history",
        )
    return times, p0, p1, identities


def _load_teams(teams_dir: Path, times: np.ndarray):
    expected_names = {
        f"{TEAM_FILENAME_PREFIX}{team}.pickle" for team in TEAM_FILE_SHA256
    }
    alias_name = f"{TEAM_FILENAME_PREFIX}{TEAM_ALIAS}.pickle"
    matching = {
        path.name: path
        for path in teams_dir.glob(f"{TEAM_FILENAME_PREFIX}*.pickle")
    }
    _require(
        set(matching) == expected_names | {alias_name},
        "ten-team directory does not contain the exact reviewed Step 0B manifest",
    )
    alias_payload = matching[alias_name].read_bytes()
    target_name = f"{TEAM_FILENAME_PREFIX}{TEAM_ALIAS_TARGET}.pickle"
    _require(
        alias_payload == matching[target_name].read_bytes(),
        "excluded SimVascular alias is not byte-identical to selected SimVascular P2",
    )

    teams = {"p0": [], "p1": []}
    provenance = []
    for team, expected_hash in TEAM_FILE_SHA256.items():
        name = f"{TEAM_FILENAME_PREFIX}{team}.pickle"
        path = matching[name]
        payload = path.read_bytes()
        observed_hash = _sha256_bytes(payload)
        _require(
            observed_hash == expected_hash,
            f"team input {name} has SHA-256 {observed_hash}, expected {expected_hash}",
        )
        data = pickle.loads(payload)  # trusted, hash-pinned Zenodo input
        _require(isinstance(data, Mapping), f"team input {name} is not a mapping")
        _require(
            "time" in data and "displacement" in data,
            f"team input {name} lacks time or displacement data",
        )
        team_times = np.asarray(data["time"], dtype=float)
        _require(
            team_times.ndim == 1
            and len(team_times) >= 2
            and np.all(np.isfinite(team_times))
            and np.all(np.diff(team_times) > 0.0)
            and team_times[0] <= times[0] + 1.000001e-3
            and team_times[-1] >= times[-1] - 1.000001e-3,
            f"team input {name} has an invalid or incomplete time grid",
        )
        displacement = data["displacement"]
        _require(
            isinstance(displacement, Mapping),
            f"team input {name} has malformed displacement data",
        )
        for lab in teams:
            _require(lab in displacement, f"team input {name} lacks {lab}")
            disp = displacement[lab]
            if isinstance(disp, np.ndarray) and disp.shape == ():
                disp = disp.item()
            _require(
                isinstance(disp, Mapping),
                f"team input {name} has malformed {lab}",
            )
            components = []
            for component in ("ux", "uy", "uz"):
                _require(
                    component in disp,
                    f"team input {name} lacks {lab}.{component}",
                )
                values = np.asarray(disp[component], dtype=float)
                _require(
                    values.shape == team_times.shape and np.all(np.isfinite(values)),
                    f"team input {name} has malformed {lab}.{component}",
                )
                components.append(np.interp(times, team_times, values))
            column = np.column_stack(
                components
            )
            teams[lab].append(column)
        provenance.append(
            {
                "team": team,
                "name": name,
                "sha256": observed_hash,
                "size_bytes": len(payload),
                "source_sample_count": int(len(team_times)),
                "source_time_start_s": float(team_times[0]),
                "source_time_end_s": float(team_times[-1]),
            }
        )
    selection = {
        "selected_count": len(provenance),
        "selected_files": [item["name"] for item in provenance],
        "excluded_alias": {
            "name": alias_name,
            "sha256": _sha256_bytes(alias_payload),
            "identical_to": target_name,
        },
        "interpolation": (
            "linear interpolation to the CoupFE 1 ms grid; endpoint hold is used "
            "only where a published team curve begins at 0.001 s or ends at 0.999 s"
        ),
    }
    return {lab: np.asarray(cur) for lab, cur in teams.items()}, provenance, selection


def _require_comparable_coupfe_runs(runs: Mapping[str, dict]) -> dict:
    reference = runs["two_layer"]["identity"]
    runtime_identity = {
        "runtime_source_sha256": reference["source_identity"][
            "runtime_source_sha256"
        ],
        "core_revision": reference["source_identity"]["core_revision"],
    }
    semantic_contract = reference["semantic_contract"]
    environment_identity = reference["environment_identity"]
    for role, run in runs.items():
        identity = run["identity"]
        observed_runtime = {
            "runtime_source_sha256": identity["source_identity"][
                "runtime_source_sha256"
            ],
            "core_revision": identity["source_identity"]["core_revision"],
        }
        _require(
            observed_runtime == runtime_identity,
            f"CoupFE runtime source identity mismatch for {role}",
        )
        _require(
            identity["semantic_contract"] == semantic_contract,
            f"CoupFE semantic contract mismatch for {role}",
        )
        _require(
            identity["environment_identity"] == environment_identity,
            f"CoupFE runtime environment mismatch for {role}",
        )
    return {
        "runtime_source_identity": runtime_identity,
        "semantic_contract": semantic_contract,
        "environment_identity": environment_identity,
        "app_revisions_may_differ": True,
        "comparison": (
            "runtime-source SHA-256 and Core revision are exact matches; "
            "all declared semantic-contract fields are exact matches"
        ),
    }


def _prefix_continuity(runs: Mapping[str, dict]) -> dict:
    full = runs["four_layer"]
    prefix = runs["four_layer_prefix"]
    _require(
        np.array_equal(full["times"][: len(EXPECTED_PREFIX_TIMES)], prefix["times"]),
        "four-layer full cycle does not share the exact retained prefix time grid",
    )
    _require(
        full["identity"]["mesh_sampling_identity"]
        == prefix["identity"]["mesh_sampling_identity"],
        "four-layer full cycle and retained prefix differ in mesh or point sampling",
    )
    point_metrics = {}
    maximum_vector_mm = 0.0
    maximum_component_mm = 0.0
    for label in ("p0", "p1"):
        delta_mm = (
            full["histories"][label][: len(EXPECTED_PREFIX_TIMES)]
            - prefix["histories"][label]
        ) * 1.0e3
        vector_mm = np.linalg.norm(delta_mm, axis=1)
        local_component = float(np.max(np.abs(delta_mm)))
        local_vector = float(np.max(vector_mm))
        maximum_component_mm = max(maximum_component_mm, local_component)
        maximum_vector_mm = max(maximum_vector_mm, local_vector)
        point_metrics[label] = {
            "max_abs_component_difference_mm": local_component,
            "max_vector_difference_mm": local_vector,
        }
    _require(
        maximum_vector_mm <= PREFIX_CONTINUITY_TOLERANCE_MM,
        "four-layer full cycle fails retained-prefix continuity: maximum vector "
        f"difference {maximum_vector_mm:.17g} mm exceeds "
        f"{PREFIX_CONTINUITY_TOLERANCE_MM:.17g} mm",
    )
    return {
        "sample_count": len(EXPECTED_PREFIX_TIMES),
        "time_window_s": [0.0, 0.32],
        "tolerance_mm": PREFIX_CONTINUITY_TOLERANCE_MM,
        "max_abs_component_difference_mm": maximum_component_mm,
        "max_vector_difference_mm": maximum_vector_mm,
        "points": point_metrics,
        "passed": True,
    }


def _load_inputs(
    coupfe_npz: Path,
    coupfe_four_layer: Path,
    coupfe_four_layer_prefix: Path,
    fenics_dir: Path,
    teams_dir: Path,
) -> dict:
    two_times, two_histories, two_identity = _load_coupfe(
        coupfe_npz,
        role="two_layer",
        expected_hash=EXPECTED_COUPFE_SHA256,
        n_t=2,
        completed_steps=1000,
        t_end=1.0,
    )
    four_times, four_histories, four_identity = _load_coupfe(
        coupfe_four_layer,
        role="four_layer",
        expected_hash=EXPECTED_COUPFE_FOUR_LAYER_SHA256,
        n_t=4,
        completed_steps=1000,
        t_end=1.0,
    )
    prefix_times, prefix_histories, prefix_identity = _load_coupfe(
        coupfe_four_layer_prefix,
        role="four_layer_prefix",
        expected_hash=EXPECTED_COUPFE_FOUR_LAYER_PREFIX_SHA256,
        n_t=4,
        completed_steps=320,
        t_end=0.32,
    )
    runs = {
        "two_layer": {
            "times": two_times,
            "histories": two_histories,
            "identity": two_identity,
        },
        "four_layer": {
            "times": four_times,
            "histories": four_histories,
            "identity": four_identity,
        },
        "four_layer_prefix": {
            "times": prefix_times,
            "histories": prefix_histories,
            "identity": prefix_identity,
        },
    }
    _require(
        np.array_equal(two_times, four_times),
        "two- and four-layer CoupFE full-cycle time grids differ",
    )
    comparison_contract = _require_comparable_coupfe_runs(runs)
    prefix_continuity = _prefix_continuity(runs)
    f_times, f_p0, f_p1, fenics_identities = _load_fenics(fenics_dir)
    teams, team_provenance, team_selection = _load_teams(teams_dir, two_times)
    return {
        "times": two_times,
        "runs": runs,
        "comparison_contract": comparison_contract,
        "prefix_continuity": prefix_continuity,
        "fenics_times": f_times,
        "fenics": {"p0": f_p0, "p1": f_p1},
        "fenics_identities": fenics_identities,
        "teams": teams,
        "team_provenance": team_provenance,
        "team_selection": team_selection,
    }


def _first_threshold_crossing(
    times: np.ndarray,
    u: np.ndarray,
    *,
    direction: str,
    start_index: int = 1,
) -> tuple[float | None, int | None]:
    z = u[:, 2]
    _require(direction in {"downward", "upward"}, "invalid crossing direction")
    for i in range(max(1, start_index), len(times)):
        crossed = (
            z[i - 1] > -0.005 and z[i] <= -0.005
            if direction == "downward"
            else z[i - 1] < -0.005 and z[i] >= -0.005
        )
        if crossed:
            frac = (-0.005 - z[i - 1]) / (z[i] - z[i - 1])
            crossing = float(times[i - 1] + frac * (times[i] - times[i - 1]))
            return crossing, i
    return None, None


def _uz_minus_5mm_events(times: np.ndarray, u: np.ndarray) -> dict:
    downward, downward_index = _first_threshold_crossing(
        times, u, direction="downward"
    )
    if downward_index is None:
        upward = None
    else:
        upward, _ = _first_threshold_crossing(
            times,
            u,
            direction="upward",
            start_index=downward_index + 1,
        )
    return {"downward_s": downward, "upward_s": upward}


def _run_reference_metrics(
    times: np.ndarray,
    histories: Mapping[str, np.ndarray],
    f_times: np.ndarray,
    fenics: Mapping[str, np.ndarray],
    teams: Mapping[str, np.ndarray],
) -> dict:
    index = np.searchsorted(times, f_times)
    _require(
        np.allclose(times[index], f_times, rtol=0.0, atol=1.0e-12),
        "CoupFE and FEniCS time grids do not align",
    )
    metrics: dict = {"fenics": {}, "ten_team_envelope": {}}
    for label in ("p0", "p1"):
        diff = (histories[label][index] - fenics[label]) * 1.0e3
        norms = np.linalg.norm(diff, axis=1)
        metrics["fenics"][label] = {
            "alignment": "matched timestamps (fenics 0.001-0.999 s)",
            "difference_direction": "CoupFE minus FEniCS",
            "last_shared_time_s": float(f_times[-1]),
            "last_shared_time_gap_mm": float(norms[-1]),
            "max_gap_mm": float(norms.max()),
            "max_gap_time_s": float(f_times[int(np.argmax(norms))]),
            "vector_rmse_mm": float(np.sqrt(np.mean(norms**2))),
            "component_rmse_mm": {
                component: float(np.sqrt(np.mean(diff[:, ci] ** 2)))
                for ci, component in enumerate(("x", "y", "z"))
            },
            "windows": {
                name: {
                    "time_window_s": [start, end],
                    **_difference_metrics(
                        f_times,
                        diff,
                        (f_times >= start - 1.0e-15)
                        & (f_times <= end + 1.0e-15),
                    ),
                }
                for name, (start, end) in FENICS_WINDOWS.items()
            },
        }
        band_lo = teams[label].min(axis=0) * 1.0e3
        band_hi = teams[label].max(axis=0) * 1.0e3
        own = histories[label] * 1.0e3
        inside_components = np.column_stack(
            [
                (own[:, ci] >= band_lo[:, ci])
                & (own[:, ci] <= band_hi[:, ci])
                for ci in range(3)
            ]
        )
        phases = {}
        for phase, (start, end) in PHASES.items():
            mask = (times >= start - 1.0e-15) & (times <= end + 1.0e-15)
            phases[phase] = {
                "time_window_s": [start, end],
                "sample_count": int(mask.sum()),
                "component_inside_fraction": {
                    component: float(inside_components[mask, ci].mean())
                    for ci, component in enumerate(("x", "y", "z"))
                },
                "all_components_inside_fraction": float(
                    inside_components[mask].all(axis=1).mean()
                ),
            }
        metrics["ten_team_envelope"][label] = {
            "component_inside_fraction": {
                component: {
                    "inside_fraction": float(inside_components[:, ci].mean())
                }
                for ci, component in enumerate(("x", "y", "z"))
            },
            "all_components_inside_fraction": float(
                inside_components.all(axis=1).mean()
            ),
            "phases": phases,
        }
    coupfe_events = {
        label: _uz_minus_5mm_events(times, histories[label])
        for label in ("p0", "p1")
    }
    fenics_events = {
        label: _uz_minus_5mm_events(f_times, fenics[label])
        for label in ("p0", "p1")
    }
    event_deltas = {}
    for label in ("p0", "p1"):
        event_deltas[label] = {}
        for direction in ("downward_s", "upward_s"):
            coupfe_event = coupfe_events[label][direction]
            fenics_event = fenics_events[label][direction]
            event_deltas[label][direction] = (
                None
                if coupfe_event is None or fenics_event is None
                else float(coupfe_event - fenics_event)
            )
    metrics["uz_minus_5mm_events"] = {
        "definition": (
            "first linearly interpolated downward u_z=-5 mm crossing and first "
            "later upward crossing"
        ),
        "coupfe": coupfe_events,
        "fenics": fenics_events,
        "coupfe_minus_fenics_s": event_deltas,
    }
    return metrics


def _difference_metrics(
    times: np.ndarray,
    difference_mm: np.ndarray,
    mask: np.ndarray,
) -> dict:
    selected = difference_mm[mask]
    selected_times = times[mask]
    norms = np.linalg.norm(selected, axis=1)
    maximum = int(np.argmax(norms))
    return {
        "sample_count": int(mask.sum()),
        "vector_rmse_mm": float(np.sqrt(np.mean(norms**2))),
        "max_vector_difference_mm": float(norms[maximum]),
        "max_vector_difference_time_s": float(selected_times[maximum]),
        "last_sample_time_s": float(selected_times[-1]),
        "last_sample_delta_mm": {
            component: float(selected[-1, ci])
            for ci, component in enumerate(("x", "y", "z"))
        },
        "last_sample_vector_difference_mm": float(norms[-1]),
        "component_rmse_mm": {
            component: float(np.sqrt(np.mean(selected[:, ci] ** 2)))
            for ci, component in enumerate(("x", "y", "z"))
        },
    }


def _cross_mesh_metrics(loaded: dict) -> dict:
    times = loaded["times"]
    two = loaded["runs"]["two_layer"]["histories"]
    four = loaded["runs"]["four_layer"]["histories"]
    windows = {}
    for name, (start, end) in CROSS_MESH_WINDOWS.items():
        mask = (times >= start - 1.0e-15) & (times <= end + 1.0e-15)
        windows[name] = {
            "time_window_s": [start, end],
            "points": {
                label: _difference_metrics(
                    times,
                    (four[label] - two[label]) * 1.0e3,
                    mask,
                )
                for label in ("p0", "p1")
            },
        }
    events = {
        role: {
            label: _uz_minus_5mm_events(
                times, loaded["runs"][role]["histories"][label]
            )
            for label in ("p0", "p1")
        }
        for role in ("two_layer", "four_layer")
    }
    event_delta = {}
    sensitivity = {}
    transient_mask = np.zeros(times.shape, dtype=bool)
    for start, end in TRANSIENT_WINDOWS:
        transient_mask |= (
            (times >= start - 1.0e-15) & (times <= end + 1.0e-15)
        )
    for label in ("p0", "p1"):
        event_delta[label] = {}
        matching_events = True
        event_shifts = []
        for direction in ("downward_s", "upward_s"):
            two_event = events["two_layer"][label][direction]
            four_event = events["four_layer"][label][direction]
            if two_event is None and four_event is None:
                event_delta[label][direction] = None
            elif two_event is None or four_event is None:
                matching_events = False
                event_delta[label][direction] = None
            else:
                shift = float(four_event - two_event)
                event_delta[label][direction] = shift
                event_shifts.append(abs(shift))

        difference_mm = (four[label] - two[label]) * 1.0e3
        norms = np.linalg.norm(difference_mm, axis=1)
        transient_maximum_mm = float(norms[transient_mask].max())
        pre_rmse_mm = windows["pre_snap"]["points"][label]["vector_rmse_mm"]
        late_rmse_mm = windows["late_relaxation"]["points"][label][
            "vector_rmse_mm"
        ]
        endpoint_mm = float(norms[-1])
        if transient_maximum_mm == 0.0:
            _require(
                pre_rmse_mm == 0.0 and late_rmse_mm == 0.0 and endpoint_mm == 0.0,
                "zero transient maximum is inconsistent with nonzero differences",
            )
            ratios = {
                "pre_snap_rmse_over_transient_max": 0.0,
                "late_relaxation_rmse_over_transient_max": 0.0,
                "endpoint_over_transient_max": 0.0,
            }
        else:
            ratios = {
                "pre_snap_rmse_over_transient_max": (
                    pre_rmse_mm / transient_maximum_mm
                ),
                "late_relaxation_rmse_over_transient_max": (
                    late_rmse_mm / transient_maximum_mm
                ),
                "endpoint_over_transient_max": endpoint_mm / transient_maximum_mm,
            }
        complete_events = matching_events and len(event_shifts) == 2
        sensitivity[label] = {
            "transient_max_difference_mm": transient_maximum_mm,
            "transient_windows_s": [list(window) for window in TRANSIENT_WINDOWS],
            "ratios": ratios,
            "matching_event_presence": matching_events,
            "complete_downward_and_upward_events": complete_events,
            "maximum_absolute_event_shift_s": (
                float(max(event_shifts)) if event_shifts else None
            ),
        }
    return {
        "direction": "four_layer minus two_layer",
        "windows": windows,
        "uz_minus_5mm_events": {
            "definition": (
                "first linearly interpolated downward crossing and first later "
                "upward crossing"
            ),
            "event_time_s": events,
            "four_minus_two_s": event_delta,
        },
        "cross_mesh_decay_diagnostic": {
            "points": sensitivity,
            "claim_boundary": (
                "These descriptive ratios quantify whether the observed cross-mesh "
                "difference decays. They are not a preregistered decision rule and "
                "cannot establish a physical solution branch or mesh convergence."
            ),
        },
        "four_layer_prefix_continuity": loaded["prefix_continuity"],
    }


def compute_metrics(
    coupfe_npz: Path,
    coupfe_four_layer: Path,
    coupfe_four_layer_prefix: Path,
    fenics_dir: Path,
    teams_dir: Path,
    *,
    _loaded: dict | None = None,
) -> dict:
    loaded = _loaded or _load_inputs(
        coupfe_npz,
        coupfe_four_layer,
        coupfe_four_layer_prefix,
        fenics_dir,
        teams_dir,
    )
    runs = {}
    for role in ("two_layer", "four_layer"):
        run = loaded["runs"][role]
        runs[role] = {
            "mesh": run["identity"]["mesh"],
            **_run_reference_metrics(
                loaded["times"],
                run["histories"],
                loaded["fenics_times"],
                loaded["fenics"],
                loaded["teams"],
            ),
        }
    fenics_change = {}
    for label in ("p0", "p1"):
        two_metrics = runs["two_layer"]["fenics"][label]
        four_metrics = runs["four_layer"]["fenics"][label]
        fenics_change[label] = {
            "vector_rmse_reduction_percent": float(
                100.0
                * (two_metrics["vector_rmse_mm"] - four_metrics["vector_rmse_mm"])
                / two_metrics["vector_rmse_mm"]
            ),
            "maximum_gap_reduction_percent": float(
                100.0
                * (two_metrics["max_gap_mm"] - four_metrics["max_gap_mm"])
                / two_metrics["max_gap_mm"]
            ),
            "full_rmse_and_maximum_both_improve": bool(
                four_metrics["vector_rmse_mm"] < two_metrics["vector_rmse_mm"]
                and four_metrics["max_gap_mm"] < two_metrics["max_gap_mm"]
            ),
        }
    return {
        "schema": SCHEMA,
        "runs": runs,
        "cross_mesh": _cross_mesh_metrics(loaded),
        "four_vs_two_fenics_change": {
            "positive_percent_means_the_four_layer_error_is_smaller": True,
            "points": fenics_change,
            "claim_boundary": (
                "An error reduction against this one FEniCS reference is an "
                "accuracy diagnostic, not a convergence proof. Window-level "
                "changes can be mixed."
            ),
        },
        "comparison_contract": loaded["comparison_contract"],
        "inputs": {
            "coupfe": {
                role: loaded["runs"][role]["identity"]
                for role in ("two_layer", "four_layer", "four_layer_prefix")
            },
            "fenics_files": loaded["fenics_identities"],
            "ten_team_dataset": {
                "doi": REFERENCE_DOI,
                "license": REFERENCE_LICENSE,
                "selection": loaded["team_selection"],
                "files": loaded["team_provenance"],
            },
        },
        "method": {
            "coupfe_fenics": (
                "Euclidean displacement-vector differences at exact shared "
                "timestamps t=0.001...0.999 s"
            ),
            "ten_team_envelope": (
                "Per-component minimum and maximum of the exact ten "
                "upstream-selected Step 0B team curves after interpolation to "
                "the CoupFE 1 ms grid"
            ),
            "cross_mesh": (
                "Four-layer minus two-layer displacement on their exact common "
                "1 ms grid; named windows include both endpoints"
            ),
            "units": "displacement metrics are millimetres; time is seconds",
        },
        "claim_boundary": (
            "Same-time comparison metrics and ten-team envelope fractions for "
            "two full-cycle CoupFE trajectories plus a cross-mesh sensitivity "
            "record. This is not a validation, pass claim, convergence proof, "
            "or evidence of distinct physical solution branches."
        ),
    }


def _accessible_svg(fig, labels: tuple[str, ...]) -> bytes:
    import io
    import re

    buffer = io.StringIO()
    fig.savefig(buffer, format="svg")
    text = buffer.getvalue()
    # Strip the XML declaration and DOCTYPE (the release contract is inert).
    text = re.sub(r"^\s*<\?xml[^?]*\?>\s*", "", text, count=1)
    text = re.sub(
        r"<!DOCTYPE[^>]*(\[[^]]*\])?>", "", text, count=1, flags=re.DOTALL
    )
    # Drop the metadata block (RDF/foreign namespaces; also holds a
    # nondeterministic date).
    text = re.sub(r"<metadata>.*?</metadata>\s*", "", text, flags=re.DOTALL)
    title = (
        f'<title id="figure-title">{FIGURE_TITLE}</title>\n'
        f'<desc id="figure-description">{FIGURE_DESCRIPTION}</desc>'
    )
    match = re.search(r"<svg[^>]*>", text)
    if match is None:
        raise RuntimeError("matplotlib SVG root not found")
    opening = match.group(0)
    if "role=" not in opening:
        opening = (
            opening[:-1]
            + ' role="img" aria-labelledby="figure-title figure-description">'
        )
    text = text[: match.start()] + opening + "\n" + title + "\n" + text[match.end():]
    text = "\n".join(line.rstrip() for line in text.splitlines()) + "\n"
    return text.encode("utf-8")


def render_figure(
    coupfe_npz: Path,
    coupfe_four_layer: Path,
    coupfe_four_layer_prefix: Path,
    fenics_dir: Path,
    teams_dir: Path,
    *,
    _loaded: dict | None = None,
) -> bytes:
    import matplotlib

    matplotlib.use("Agg")
    matplotlib.rcParams["svg.fonttype"] = "none"
    matplotlib.rcParams["svg.hashsalt"] = SCHEMA
    from matplotlib import pyplot as plt

    loaded = _loaded or _load_inputs(
        coupfe_npz,
        coupfe_four_layer,
        coupfe_four_layer_prefix,
        fenics_dir,
        teams_dir,
    )
    times = loaded["times"]
    two_m = loaded["runs"]["two_layer"]["histories"]
    four_m = loaded["runs"]["four_layer"]["histories"]
    two = {label: values * 1.0e3 for label, values in two_m.items()}
    four = {label: values * 1.0e3 for label, values in four_m.items()}
    f_times = loaded["fenics_times"]
    fenics = {
        label: values * 1.0e3 for label, values in loaded["fenics"].items()
    }
    teams = loaded["teams"]

    fig, axes = plt.subplots(2, 3, figsize=(15, 7.5), sharex=True)
    for row, lab in enumerate(("p0", "p1")):
        array = teams[lab] * 1.0e3
        band_lo, band_hi = array.min(axis=0), array.max(axis=0)
        for col, comp in enumerate(("x", "y", "z")):
            axis = axes[row, col]
            axis.fill_between(
                times, band_lo[:, col], band_hi[:, col],
                color="#D0D5DD", alpha=0.65,
                label="ten-team envelope" if (row, col) == (0, 0) else None,
            )
            axis.plot(
                f_times,
                fenics[lab][:, col],
                color="#1F2937",
                lw=1.4,
                ls="-",
                label="local FEniCS" if (row, col) == (0, 0) else None,
            )
            axis.plot(
                times,
                two[lab][:, col],
                color="#D97706",
                lw=1.5,
                ls="--",
                label=(
                    "CoupFE 2x20x17 tip_refine=6.0"
                    if (row, col) == (0, 0)
                    else None
                ),
            )
            axis.plot(
                times,
                four[lab][:, col],
                color="#2563EB",
                lw=1.5,
                ls="-.",
                label=(
                    "CoupFE 4x20x17 tip_refine=6.0"
                    if (row, col) == (0, 0)
                    else None
                ),
            )
            axis.set_ylabel(f"{lab}  u_{comp} (mm)")
            axis.grid(alpha=0.3)
            axis.set_xlim(0.0, 1.0)
            if row == 0:
                axis.set_title(f"u_{comp}")
            if row == 1:
                axis.set_xlabel("t (s)")
    axes[0, 0].legend(loc="upper right", fontsize=9)
    two_source = loaded["runs"]["two_layer"]["identity"]["source_identity"]
    four_source = loaded["runs"]["four_layer"]["identity"]["source_identity"]
    two_app_revision = two_source["app_revision"][:7]
    four_app_revision = four_source["app_revision"][:7]
    core_revision = two_source["core_revision"][:7]
    fig.suptitle(FIGURE_TITLE, y=0.995)
    fig.text(0.5, 0.958, FIGURE_SUBTITLE, ha="center", fontsize=9, color="#475467")
    fig.text(
        0.01,
        0.005,
        "  |  ".join(VISIBLE_LABELS)
        + f"  |  app 2L {two_app_revision}, 4L {four_app_revision}"
        + f"  |  Core {core_revision}",
        fontsize=7,
    )
    fig.tight_layout(rect=(0.0, 0.025, 1.0, 0.94))
    payload = _accessible_svg(fig, VISIBLE_LABELS)
    plt.close(fig)
    return payload


def _write_bytes_atomic(path: Path, payload: bytes) -> None:
    path = path.expanduser()
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=path.parent,
            prefix=f".{path.name}.",
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.chmod(temporary, 0o644)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def _write_artifact_pair_with_rollback(
    first: tuple[Path, bytes], second: tuple[Path, bytes]
) -> None:
    """Replace both outputs, restoring the old pair after a write failure."""

    items = (first, second)
    snapshots = {}
    for path, _ in items:
        if path.exists():
            snapshots[path] = {
                "existed": True,
                "payload": path.read_bytes(),
                "mode": path.stat().st_mode & 0o777,
            }
        else:
            snapshots[path] = {"existed": False, "payload": None, "mode": None}
    try:
        for path, payload in items:
            _write_bytes_atomic(path, payload)
    except Exception as write_error:
        rollback_errors = []
        for path, _ in reversed(items):
            snapshot = snapshots[path]
            try:
                if snapshot["existed"]:
                    _write_bytes_atomic(path, snapshot["payload"])
                    os.chmod(path, snapshot["mode"])
                elif path.exists():
                    path.unlink()
            except Exception as rollback_error:  # pragma: no cover - OS failure
                rollback_errors.append(f"{path}: {rollback_error}")
        if rollback_errors:  # pragma: no cover - requires repeated OS failure
            raise RuntimeError(
                "artifact-pair write failed and rollback was incomplete: "
                + "; ".join(rollback_errors)
            ) from write_error
        raise


def _validate_cli_paths(args: argparse.Namespace) -> None:
    coupfe_npz = args.coupfe_npz.expanduser().resolve()
    coupfe_four_layer = args.coupfe_four_layer.expanduser().resolve()
    coupfe_four_layer_prefix = args.coupfe_four_layer_prefix.expanduser().resolve()
    fenics_dir = args.fenics_dir.expanduser().resolve()
    teams_dir = args.teams_dir.expanduser().resolve()
    out_metrics = args.out_metrics.expanduser().resolve()
    out_figure = args.out_figure.expanduser().resolve()

    args.coupfe_npz = coupfe_npz
    args.coupfe_four_layer = coupfe_four_layer
    args.coupfe_four_layer_prefix = coupfe_four_layer_prefix
    args.fenics_dir = fenics_dir
    args.teams_dir = teams_dir
    args.out_metrics = out_metrics
    args.out_figure = out_figure

    _require(out_metrics != out_figure, "output paths must be distinct")
    coupfe_inputs = {
        coupfe_npz,
        coupfe_four_layer,
        coupfe_four_layer_prefix,
    }
    _require(
        len(coupfe_inputs) == 3,
        "two-layer, four-layer, and four-layer-prefix inputs must be distinct",
    )
    _require(
        out_metrics.suffix.lower() == ".json",
        "--out-metrics must name a .json file",
    )
    _require(
        out_figure.suffix.lower() == ".svg",
        "--out-figure must name an .svg file",
    )
    _require(
        out_metrics not in coupfe_inputs and out_figure not in coupfe_inputs,
        "an output path must not overwrite any CoupFE input archive",
    )
    for output in (out_metrics, out_figure):
        _require(
            not output.is_relative_to(fenics_dir),
            "outputs must not be written inside the FEniCS input directory",
        )
        _require(
            not output.is_relative_to(teams_dir),
            "outputs must not be written inside the ten-team input directory",
        )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--coupfe-npz",
        "--coupfe-two-layer",
        dest="coupfe_npz",
        type=Path,
        required=True,
        help="retained 2x20x17 full-cycle archive",
    )
    parser.add_argument("--coupfe-four-layer", type=Path, required=True)
    parser.add_argument("--coupfe-four-layer-prefix", type=Path, required=True)
    parser.add_argument("--fenics-dir", type=Path, required=True)
    parser.add_argument("--teams-dir", type=Path, required=True)
    parser.add_argument("--out-metrics", type=Path, required=True)
    parser.add_argument("--out-figure", type=Path, required=True)
    args = parser.parse_args(argv)
    _validate_cli_paths(args)

    loaded = _load_inputs(
        args.coupfe_npz,
        args.coupfe_four_layer,
        args.coupfe_four_layer_prefix,
        args.fenics_dir,
        args.teams_dir,
    )
    metrics = compute_metrics(
        args.coupfe_npz,
        args.coupfe_four_layer,
        args.coupfe_four_layer_prefix,
        args.fenics_dir,
        args.teams_dir,
        _loaded=loaded,
    )
    metrics_bytes = (
        json.dumps(metrics, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    figure = render_figure(
        args.coupfe_npz,
        args.coupfe_four_layer,
        args.coupfe_four_layer_prefix,
        args.fenics_dir,
        args.teams_dir,
        _loaded=loaded,
    )
    # Build and validate both artifacts before replacing either destination.
    _write_artifact_pair_with_rollback(
        (args.out_metrics, metrics_bytes),
        (args.out_figure, figure),
    )
    print(
        json.dumps(
            {
                "metrics": {
                    "path": str(args.out_metrics),
                    "sha256": _sha256_bytes(metrics_bytes),
                },
                "figure": {
                    "path": str(args.out_figure),
                    "sha256": _sha256_bytes(figure),
                    "size_bytes": len(figure),
                },
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
