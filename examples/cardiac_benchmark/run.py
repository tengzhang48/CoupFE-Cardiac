"""Benchmark 1, Step 0 Cases A/B and Step 2 Case B.

Backward Euler is the default coherent time-discretization path.
``--integrator newmark`` selects an experimental Newmark inertia/Robin path;
the material viscous strain rate remains a backward difference. The driver composes:
  * Holzapfel-Ogden + active-stress + viscous element (Hex8 F-bar, or
    standard Hex8 plus an element-local condensed pressure)
  * consistent Q1 Hex8 mass by default (ρ = 1000 kg/m³), with the
    historical row-summed lumped path available explicitly
  * pericardial Robin spring-dashpot operator (epi normal-only + base full)

Active tension τ(t) is updated as a per-step material property. Step 2 Case B
combines that active stress with endocardial pressure and its published Table
5 parameters. No Dirichlet
BCs -- the Robin springs hold the body.  Displacement history at the reference
points p0, p1 is recorded for comparison with the benchmark reference data.
"""
from __future__ import annotations

import argparse
import atexit
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

import coupfe

sys.path.insert(0, os.path.dirname(__file__))
import geometry as geom
import tbar_laplace
from activation import tau_of_t, p_of_t
from benchmark_parameters import (
    benchmark_configuration,
    benchmark_metadata,
    validate_load_histories,
    validate_runtime_material_properties,
)
from material import (
    CardiacHex8,
    HO_PARAMS,
    MATERIAL_MODEL_ID,
    _PER_GP,
    struct_tensors,
)
from local_pressure import LocalPressureHex8Operator
from boundary_audit import audit_geometry, audit_pressure, audit_robin
from consistent_mass import (
    ConsistentMassInertia,
    ConsistentNewmarkInertia,
    consistent_mass_coo,
)
from robin import RobinOperator
from pressure import FollowerPressureOperator
from newmark import NewmarkInertia
from result_io import save_completed
from sampling import interpolate_displacement, locate_hex8_point
from solver import PetscSnesSolver, SnesSolveError, checked_newton_solve

from coupfe.operators.element_group import ElementGroup
from coupfe.operators.inertia import InertiaOperator
from coupfe.runtime.compiled_element import CompiledElement

DENSITY = 1.0e3
A_TOP, B_TOP = 1.0e5, 5.0e3       # base Robin (Pa/m, Pa·s/m)
A_EPI, B_EPI = 1.0e8, 5.0e3       # epi  Robin (normal-only)

_KERNEL = {}
PUBLIC_CORE_URL = "https://github.com/tengzhang48/CoupFE.git"


def _sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _package_version(name):
    try:
        return importlib.metadata.version(name)
    except importlib.metadata.PackageNotFoundError:
        return "unknown"


def _checkout_root(start, *, require_core_package=False):
    start = Path(start).resolve()
    for path in (start, *start.parents):
        if not (path / ".git").exists():
            continue
        if require_core_package:
            checkout_package = (path / "coupfe").resolve()
            if start != checkout_package or not (
                checkout_package / "__init__.py"
            ).is_file():
                continue
        return path
    return None


def _source_revision(start, override, *, require_core_package=False):
    """Return an explicit or nearest-checkout Git revision without recording a path."""
    value = os.environ.get(override)
    if value:
        return value
    checkout = _checkout_root(start, require_core_package=require_core_package)
    if checkout is None:
        return "unknown"
    try:
        completed = subprocess.run(
            ["git", "-C", str(checkout), "rev-parse", "HEAD"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=2,
        )
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    return completed.stdout.strip() or "unknown"


def _source_tree_state(start, override, *, require_core_package=False):
    """Return clean/dirty/unknown without recording the checkout path."""
    value = os.environ.get(override)
    if value:
        return value
    checkout = _checkout_root(start, require_core_package=require_core_package)
    if checkout is None:
        return "unknown"
    try:
        completed = subprocess.run(
            ["git", "-C", str(checkout), "status", "--porcelain", "--untracked-files=normal"],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=2,
        )
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    return "dirty" if completed.stdout else "clean"


def _installed_vcs_identity(distribution_name="coupfe", *, module_file=None):
    """Return a PEP 610 VCS identity for a non-worktree Core installation."""
    try:
        distribution = importlib.metadata.distribution(distribution_name)
        raw = distribution.read_text("direct_url.json")
    except importlib.metadata.PackageNotFoundError:
        return None
    if raw is None:
        return None
    resolved_module = module_file or getattr(coupfe, "__file__", None)
    if resolved_module is None:
        return None
    try:
        imported_path = Path(resolved_module).resolve(strict=True)
        distributed_path = Path(
            distribution.locate_file("coupfe/__init__.py")
        ).resolve(strict=True)
    except (AttributeError, OSError, RuntimeError):
        return None
    if imported_path != distributed_path:
        return None
    try:
        record = json.loads(raw)
    except json.JSONDecodeError:
        return None
    vcs = record.get("vcs_info", {})
    revision = vcs.get("commit_id")
    if (
        record.get("url") != PUBLIC_CORE_URL
        or vcs.get("vcs") != "git"
        or not isinstance(revision, str)
        or len(revision) != 40
        or any(character not in "0123456789abcdefABCDEF" for character in revision)
    ):
        return None
    return {
        "revision": revision.lower(),
        "tree_state": "installed",
        "source_kind": "pep610-vcs",
        "source_url": PUBLIC_CORE_URL,
    }


def _runtime_metadata():
    app_source = Path(__file__).resolve().parents[2]
    core_file = getattr(coupfe, "__file__", None)
    core_source = Path(core_file).resolve().parent if core_file else None
    app_revision = _source_revision(app_source, "COUPFE_CARDIAC_APP_REVISION")
    app_tree_state = _source_tree_state(app_source, "COUPFE_CARDIAC_TREE_STATE")
    app_source_kind = (
        os.environ.get("COUPFE_CARDIAC_SOURCE_KIND", "asserted")
        if os.environ.get("COUPFE_CARDIAC_APP_REVISION")
        else ("git-checkout" if app_revision != "unknown" else "unknown")
    )
    core_revision = (
        _source_revision(
            core_source,
            "COUPFE_CORE_REVISION",
            require_core_package=True,
        )
        if core_source else os.environ.get("COUPFE_CORE_REVISION", "unknown")
    )
    core_tree_state = (
        _source_tree_state(
            core_source,
            "COUPFE_CORE_TREE_STATE",
            require_core_package=True,
        )
        if core_source else os.environ.get("COUPFE_CORE_TREE_STATE", "unknown")
    )
    core_source_kind = (
        os.environ.get("COUPFE_CORE_SOURCE_KIND", "asserted")
        if os.environ.get("COUPFE_CORE_REVISION")
        else ("git-checkout" if core_revision != "unknown" else "unknown")
    )
    core_source_url = PUBLIC_CORE_URL if core_revision != "unknown" else "unknown"
    if core_revision == "unknown":
        installed = _installed_vcs_identity()
        if installed is not None:
            core_revision = installed["revision"]
            core_tree_state = installed["tree_state"]
            core_source_kind = installed["source_kind"]
            core_source_url = installed["source_url"]

    return {
        "result_schema": "coupfe-cardiac-result-v1",
        "driver": "examples/cardiac_benchmark/run.py",
        "app_revision": app_revision,
        "app_tree_state": app_tree_state,
        "app_source_kind": app_source_kind,
        "core_revision": core_revision,
        "core_tree_state": core_tree_state,
        "core_source_kind": core_source_kind,
        "core_source_url": core_source_url,
        "coupfe_version": _package_version("coupfe"),
        "numpy_version": np.__version__,
        "scipy_version": _package_version("scipy"),
        "python_version": platform.python_version(),
    }


def _hex_lumped_mass(nodes, elems, density, ndof, dof_per_node=3):
    nodal = np.zeros(len(nodes))
    for e in elems:
        vol = geom._hex_volume(nodes[e])
        np.add.at(nodal, e, density * vol / 8.0)
    M = np.zeros(ndof)
    for c in range(3):
        M[c::dof_per_node] = nodal
    return M


# Hex8 node natural coords (Abaqus order; matches shape_hex8.for + our elems)
_NAT = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                 [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)


def _hex8_gp_shape():
    """Shape functions N[gp, node] at the 8 Gauss points, in gauss_hex8 order
    (ii fastest, then jj, then kk; sign +g1 for the 1st of each pair)."""
    g1 = 1.0 / np.sqrt(3.0)
    GP = np.empty((8, 3))
    for g in range(8):
        GP[g, 0] = g1 if g % 2 == 0 else -g1
        GP[g, 1] = g1 if (g // 2) % 2 == 0 else -g1
        GP[g, 2] = g1 if g // 4 == 0 else -g1
    N = np.empty((8, 8))
    for g in range(8):
        for a in range(8):
            N[g, a] = 0.125 * np.prod(1.0 + _NAT[a] * GP[g])
    return N


def _init_fiber_svars(
    elem,
    mesh,
    schema,
    per_gp,
    *,
    sampling="cg1",
    tbar=None,
    asign=-1.0,
):
    """Initialize the structural tensors at the eight Hex8 Gauss points.

    ``cg1`` interpolates nodal directions and applies the historical local
    Gram--Schmidt projection. ``gp-direct`` evaluates the topology's declared
    rule at each Gauss point: physical-coordinate toolkit reconstruction for
    the closed mesh, and the retained parametric rule for the historical mesh.
    """
    if sampling not in {"cg1", "gp-direct"}:
        raise ValueError(f"unsupported fiber sampling {sampling!r}")
    n_gp = elem.svars.shape[1] // per_gp
    o = {k: schema[k]["offset"] for k in schema}
    N = _hex8_gp_shape()
    if n_gp != len(N):
        raise ValueError(f"expected 8 Hex8 Gauss points, found {n_gp}")
    if tbar is None:
        tbar = mesh.param[:, 0]
    tbar = np.asarray(tbar, dtype=float)
    if tbar.shape != (mesh.n_node,) or not np.all(np.isfinite(tbar)):
        raise ValueError("fiber transmural coordinates must be finite per-node data")
    for e in range(mesh.n_elem):
        element_nodes = mesh.elems[e]
        if sampling == "cg1":
            fn = mesh.fiber_node[element_nodes]  # (8,3) nodal fibers
            sn = mesh.sheet_node[element_nodes]
        for g in range(n_gp):
            if sampling == "cg1":
                f0 = N[g] @ fn
                s0 = N[g] @ sn
            else:
                f0, s0, _ = geom.sample_structural_frame(
                    mesh,
                    element_nodes,
                    N[g],
                    tbar,
                    asign=asign,
                )
            ff, ss, fssym = struct_tensors(f0, s0)
            base = g * per_gp
            elem.svars[e, base + o["ff"]:base + o["ff"] + 9] = ff.ravel()
            elem.svars[e, base + o["ss"]:base + o["ss"] + 9] = ss.ravel()
            elem.svars[e, base + o["fssym"]:base + o["fssym"] + 9] = fssym.ravel()
    elem.svars_trial = elem.svars.copy()


def build_group(
    mesh,
    dt,
    build_dir=None,
    *,
    formulation="fbar",
    evaluation_mode="joint",
    fiber_sampling="cg1",
    fiber_tbar=None,
    fiber_asign=-1.0,
    material_parameters=None,
):
    """Build the stateful material group and optional condensed pressure term."""
    from material import build_kernel
    kernel_formulations = {
        "fbar": "fbar_mechanics",
        "local-pressure": "standard",
        "std-kappa": "standard",
    }
    if formulation not in kernel_formulations:
        raise ValueError(f"unsupported cardiac formulation {formulation!r}")
    kernel_formulation = kernel_formulations[formulation]
    if kernel_formulation not in _KERNEL:
        module_tag = {
            "fbar": "fbar",
            "local-pressure": "local_pressure_standard",
            "std-kappa": "std_kappa",
        }[formulation]
        _, _KERNEL[kernel_formulation] = build_kernel(
            tmpdir=build_dir,
            formulation=kernel_formulation,
            module_name=f"cardiac_hex8_run_{module_tag}",
        )
    mod = _KERNEL[kernel_formulation]
    problem = CardiacHex8()
    dpn, comps = 3, (0, 1, 2)
    schema = problem._mat.state_schema
    per_gp = sum(v["size"] for v in schema.values())
    props = np.asarray(problem._mat.props_array, float).copy()
    prop_names = list(problem._mat.props.keys())
    if material_parameters is not None:
        selected = dict(material_parameters)
        if set(selected) != set(prop_names):
            raise RuntimeError(
                "selected material parameters do not match the generated kernel"
            )
        props[:] = [float(selected[name]) for name in prop_names]
    eta_idx = prop_names.index("eta")
    if props[eta_idx] != 100.0:
        raise RuntimeError("the benchmark material must retain eta=100 Pa s")
    kappa_idx = prop_names.index("kappa")
    bulk_modulus = float(props[kappa_idx])
    if formulation == "local-pressure":
        # The standard material kernel supplies only the deviatoric/anisotropic,
        # active, and viscous terms. The one-pressure-per-element operator below
        # supplies the complete volumetric contribution exactly once.
        props[kappa_idx] = 0.0
    elem = CompiledElement(
        mod, props=props,
        dof_per_node=dpn, n_svars=per_gp * 8, mcrd=3,
        n_elem=mesh.n_elem, dt=dt, backend="native",
        state_schema=schema,
    )
    _init_fiber_svars(
        elem,
        mesh,
        schema,
        per_gp,
        sampling=fiber_sampling,
        tbar=fiber_tbar,
        asign=fiber_asign,
    )
    grp = ElementGroup(
        elem,
        mesh.nodes,
        mesh.elems,
        dof_per_node=dpn,
        comps=comps,
        evaluation_mode=evaluation_mode,
    )
    ta_idx = prop_names.index("Ta")
    local_pressure = None
    if formulation == "local-pressure":
        local_pressure = LocalPressureHex8Operator(
            mesh.nodes,
            mesh.elems,
            dpn * mesh.n_node,
            bulk_modulus=bulk_modulus,
            dof_per_node=dpn,
        )
    return grp, elem, ta_idx, local_pressure


def _apply_laplace_tbar(mesh, path, *, asign):
    """Replace the analytic parametric t̃ with a precomputed Laplace t̄ field.

    Recomputes per-node and per-element fiber/sheet/normal frames using the
    topology's declared reconstruction rule.  Closed meshes use physical
    coordinates exactly as the pinned toolkit does; historical polar meshes
    retain their parametric convention. The field must cover every mesh node
    in mesh order and satisfy t̄=0 on endocardial and t̄=1 on epicardial nodes
    (which also validates the node ordering).
    """
    tbar = np.load(path, allow_pickle=False)
    if tbar.shape != (mesh.n_node,):
        raise ValueError(
            f"Laplace tbar field shape {tbar.shape} does not match mesh "
            f"node count {mesh.n_node}"
        )
    if (
        not np.all(np.isfinite(tbar))
        or tbar.min() < -1e-9
        or tbar.max() > 1.0 + 1e-9
    ):
        raise ValueError("Laplace tbar field must be finite and within [0, 1]")
    endo = mesh.param[:, 0] == 0.0
    epi = mesh.param[:, 0] == 1.0
    if (np.abs(tbar[endo]).max() > 1e-6) or (np.abs(tbar[epi] - 1.0).max() > 1e-6):
        raise ValueError(
            "Laplace tbar field does not satisfy endo=0/epi=1 on this mesh "
            "(wrong mesh or node ordering?)"
        )
    delta = tbar - mesh.param[:, 0]
    print(
        f"tbar: injected Laplace field {path} "
        f"(max |delta_tbar| vs analytic = {np.abs(delta).max():.4f}, implied "
        f"max fiber-angle change {120.0 * np.abs(delta).max():.2f} deg)",
        flush=True,
    )
    geom.update_structural_frames(mesh, tbar, asign=asign)
    return np.asarray(tbar, dtype=float).copy()


def _load_laplace_tbar_metadata(mesh, field_path):
    """Validate the portable sidecar emitted with a Laplace tbar field."""
    field_path = Path(field_path)
    metadata_path = field_path.with_suffix(".meta.json")
    if not metadata_path.is_file():
        raise ValueError(
            f"Laplace tbar field requires sibling metadata {metadata_path.name}"
        )
    def reject_nonfinite_constant(value):
        raise ValueError(f"non-finite JSON constant {value!r}")

    try:
        metadata = json.loads(
            metadata_path.read_text(encoding="utf-8"),
            parse_constant=reject_nonfinite_constant,
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError(
            f"cannot read Laplace tbar metadata {metadata_path.name}: {error}"
        ) from error
    if not isinstance(metadata, dict):
        raise ValueError("Laplace tbar metadata must be a JSON object")

    def finite_tree(value, location="metadata"):
        if value is None or isinstance(value, (str, bool)):
            return
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            if not np.isfinite(float(value)):
                raise ValueError(f"Laplace tbar {location} is non-finite")
            return
        if isinstance(value, list):
            for index, child in enumerate(value):
                finite_tree(child, f"{location}[{index}]")
            return
        if isinstance(value, dict) and all(
            isinstance(key, str) for key in value
        ):
            for key, child in value.items():
                finite_tree(child, f"{location}.{key}")
            return
        raise ValueError(f"Laplace tbar {location} has unsupported JSON data")

    finite_tree(metadata)
    schema = tbar_laplace.SCHEMA
    if metadata.get("schema") != schema:
        raise ValueError("Laplace tbar metadata has an unsupported schema")
    if metadata.get("mesh_topology") != mesh.topology:
        raise ValueError("Laplace tbar metadata does not match mesh topology")
    if metadata.get("output_npy") != field_path.name:
        raise ValueError("Laplace tbar metadata does not name this field")
    field_sha256 = _sha256_file(field_path)
    if metadata.get("sha256") != field_sha256:
        raise ValueError("Laplace tbar field SHA-256 disagrees with metadata")

    mesh_parameters = metadata.get("mesh_parameters")
    if not isinstance(mesh_parameters, dict):
        raise ValueError("Laplace tbar metadata is missing mesh parameters")
    expected_mesh = {
        "n_t": int(mesh.n_t),
        "n_core": int(mesh.n_core),
        "n_radial": int(mesh.n_radial),
        "nodes": int(mesh.n_node),
        "elements": int(mesh.n_elem),
    }
    for name, expected in expected_mesh.items():
        observed = mesh_parameters.get(name)
        if (
            isinstance(observed, bool)
            or not isinstance(observed, (int, float))
            or not float(observed).is_integer()
            or int(observed) != expected
        ):
            raise ValueError(
                f"Laplace tbar metadata mesh parameter {name!r} disagrees"
            )
    core_half_width = mesh_parameters.get("core_half_width")
    if (
        isinstance(core_half_width, bool)
        or not isinstance(core_half_width, (int, float))
        or not np.isclose(
            float(core_half_width),
            float(mesh.core_half_width),
            rtol=0.0,
            atol=1.0e-15,
        )
    ):
        raise ValueError(
            "Laplace tbar metadata core-half-width disagrees with the mesh"
        )
    # Sidecars predate the tip-refine control omit the key; 1.0 is the
    # uniform-mesh default, so absence remains fail-closed for graded meshes.
    tip_refine = mesh_parameters.get("tip_refine", 1.0)
    if (
        isinstance(tip_refine, bool)
        or not isinstance(tip_refine, (int, float))
        or not np.isclose(
            float(tip_refine),
            float(mesh.tip_refine),
            rtol=0.0,
            atol=1.0e-15,
        )
    ):
        raise ValueError(
            "Laplace tbar metadata tip-refine disagrees with the mesh"
        )

    mesh_identity = metadata.get("mesh_identity")
    if not isinstance(mesh_identity, dict):
        raise ValueError(
            "Laplace tbar metadata is not bound to mesh geometry; regenerate "
            "the field and sidecar"
        )
    expected_identity = tbar_laplace.closed_mesh_identity(
        mesh.nodes, mesh.elems
    )
    if set(mesh_identity) != set(expected_identity):
        raise ValueError("Laplace tbar mesh identity has invalid fields")
    if mesh_identity.get("schema") != expected_identity["schema"]:
        raise ValueError("Laplace tbar mesh identity has an unsupported schema")
    for name in (
        "node_coordinates_sha256",
        "element_connectivity_sha256",
    ):
        if mesh_identity.get(name) != expected_identity[name]:
            raise ValueError(
                f"Laplace tbar mesh identity {name!r} disagrees with the mesh"
            )

    expected_boundaries = {
        "endocardium": "Dirichlet tbar=0",
        "epicardium": "Dirichlet tbar=1",
        "base": "natural homogeneous Neumann",
    }
    if metadata.get("boundary_conditions") != expected_boundaries:
        raise ValueError("Laplace tbar metadata has different boundary conditions")

    solver_record = metadata.get("solver")
    if not isinstance(solver_record, dict):
        raise ValueError("Laplace tbar metadata is missing solver diagnostics")
    expected_solver_text = {
        "method": "Q1 Hex8 Galerkin; SciPy sparse direct solve",
        "linear_solver": "scipy.sparse.linalg.spsolve",
        "natural_boundary": "homogeneous Neumann on base",
    }
    if any(
        solver_record.get(name) != expected
        for name, expected in expected_solver_text.items()
    ):
        raise ValueError("Laplace tbar metadata has different solver semantics")
    for name, expected in {
        "n_node": mesh.n_node,
        "n_element": mesh.n_elem,
    }.items():
        observed = solver_record.get(name)
        if (
            isinstance(observed, bool)
            or not isinstance(observed, (int, float))
            or not float(observed).is_integer()
            or int(observed) != expected
        ):
            raise ValueError(f"Laplace tbar solver {name} disagrees with mesh")
    matrix_nnz = solver_record.get("matrix_nnz")
    if (
        isinstance(matrix_nnz, bool)
        or not isinstance(matrix_nnz, (int, float))
        or not float(matrix_nnz).is_integer()
        or int(matrix_nnz) <= 0
    ):
        raise ValueError("Laplace tbar metadata has invalid matrix_nnz")
    numeric = {}
    for name in (
        "linear_residual_l2",
        "linear_residual_inf",
        "linear_residual_limit",
        "linear_rhs_l2",
        "max_abs_boundary_endo",
        "max_abs_boundary_epi_minus_one",
        "minimum",
        "maximum",
        "minimum_gauss_det_j_m3",
        "maximum_gauss_det_j_m3",
    ):
        value = solver_record.get(name)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"Laplace tbar metadata is missing numeric {name}")
        numeric[name] = float(value)
    if (
        numeric["linear_residual_l2"] < 0.0
        or numeric["linear_residual_inf"] < 0.0
        or numeric["linear_residual_limit"] <= 0.0
        or numeric["linear_residual_l2"] > numeric["linear_residual_limit"]
        or numeric["max_abs_boundary_endo"] > 1.0e-12
        or numeric["max_abs_boundary_epi_minus_one"] > 1.0e-12
        or numeric["minimum"] < -1.0e-9
        or numeric["maximum"] > 1.0 + 1.0e-9
        or numeric["minimum_gauss_det_j_m3"] <= 0.0
        or numeric["maximum_gauss_det_j_m3"]
        < numeric["minimum_gauss_det_j_m3"]
    ):
        raise ValueError("Laplace tbar metadata records a failed solver check")
    return {
        "field_sha256": field_sha256,
        "metadata_filename": metadata_path.name,
        "metadata_sha256": _sha256_file(metadata_path),
        "metadata_schema": schema,
    }


def _time_grid(dt, t_end):
    """Return an exact, increasing grid whose final point is ``t_end``."""
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("--dt must be finite and positive")
    if not np.isfinite(t_end) or t_end <= 0.0:
        raise ValueError("--tend must be finite and positive")
    intervals = int(round(t_end / dt))
    if intervals < 1 or not np.isclose(
        intervals * dt, t_end, rtol=1.0e-12, atol=1.0e-15
    ):
        raise ValueError("--tend must be an integer multiple of --dt")
    return np.linspace(0.0, t_end, intervals + 1)


def _load_schedule_grid(dt, t_end, load_horizon=None):
    """Return the solve prefix and the grid used to integrate its load.

    ``solve_ivp`` is adaptive, so integrating the same ODE over different final
    intervals can change roundoff-level values in their common prefix.  An
    explicit horizon makes truncated serial/MPI rank gates use the exact prefix
    of the production load history.
    """
    solve_grid = _time_grid(dt, t_end)
    horizon = t_end if load_horizon is None else float(load_horizon)
    if not np.isfinite(horizon) or horizon <= 0.0:
        raise ValueError("--load-horizon must be finite and positive")
    schedule_grid = _time_grid(dt, horizon)
    if len(schedule_grid) < len(solve_grid):
        raise ValueError("--load-horizon must be greater than or equal to --tend")
    solve_grid = schedule_grid[: len(solve_grid)].copy()
    if not np.isclose(solve_grid[-1], t_end, rtol=1.0e-12, atol=1.0e-15):
        raise ValueError("--load-horizon and --tend must share the --dt grid")
    return solve_grid, schedule_grid, horizon


def _benchmark_load_histories(
    case, dt, t_end, load_horizon=None, *, benchmark_step=0
):
    """Integrate the selected benchmark loads and return their solve prefix."""
    configuration = benchmark_configuration(benchmark_step, case)
    times, schedule_times, horizon = _load_schedule_grid(
        dt, t_end, load_horizon
    )
    count = len(times)
    tau = np.zeros(count, dtype=float)
    pressure = np.zeros(count, dtype=float)
    if configuration.active_stress_enabled:
        tau = tau_of_t(
            schedule_times, p=configuration.activation_parameters
        )[:count].copy()
    if configuration.pressure_enabled:
        pressure = p_of_t(
            schedule_times, p=configuration.pressure_parameters
        )[:count].copy()
    validate_load_histories(configuration, times, tau, pressure)
    return times, tau, pressure, horizon


def _parser():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nt", type=int, default=2)
    ap.add_argument("--nmu", type=int, default=12)
    ap.add_argument("--ntheta", type=int, default=16)
    ap.add_argument(
        "--mesh-topology",
        choices=["polar-ring", "closed-multiblock"],
        default="closed-multiblock",
        help="polar-ring = historical open/collapsed latitude-longitude mesh; "
             "closed-multiblock = noncollapsed five-block benchmark domain "
             "(default)",
    )
    ap.add_argument("--ncore", type=int, default=20,
                    help="even central-square element count per side for the closed mesh")
    ap.add_argument("--nradial", type=int, default=17,
                    help="radial element layers in each closed-mesh outer block")
    ap.add_argument("--core-half-width", type=float, default=0.36,
                    help="central-square half-width in unit-disk coordinates")
    ap.add_argument("--tip-refine", type=float, default=1.0,
                    help="closed-mesh local refinement near the apex: the meridian "
                         "spacing adjacent to the apex is scaled by 1/tip-refine "
                         "while element counts, connectivity, boundary labels, and "
                         "the base rim are unchanged; 1.0 keeps the uniform "
                         "benchmark mesh")
    ap.add_argument("--dt", type=float, default=2e-3)
    ap.add_argument("--tend", type=float, default=1.0)
    ap.add_argument(
        "--load-horizon",
        type=float,
        default=None,
        metavar="SECONDS",
        help=(
            "integrate the benchmark load through this fixed horizon, then use "
            "the exact prefix through --tend; use 1.0 s for shortened rank gates"
        ),
    )
    ap.add_argument("--raw-helix", action="store_true",
                    help="use asign=+1 instead of the pinned-formula-matched "
                         "default; this changes prescribed helix handedness")
    ap.add_argument(
        "--case",
        choices=["A", "B"],
        default="A",
        help=(
            "case identity within --benchmark-step; Step 0 A is active-only, "
            "Step 0 B is pressure-only, and Step 2 B combines both"
        ),
    )
    ap.add_argument(
        "--benchmark-step",
        type=int,
        choices=[0, 2],
        default=0,
        help=(
            "Benchmark 1 protocol step; Step 0 preserves the historical A/B "
            "split, while Step 2 currently supports Case B with joint active "
            "stress and pressure"
        ),
    )
    ap.add_argument("--perturb", type=float, default=0.0,
                    help="break mesh symmetry with N(0,σ) node noise [m] (diagnostic)")
    ap.add_argument("--integrator", choices=["newmark", "be"], default="be",
                    help="be = coherent backward Euler; newmark = experimental inertia/Robin path")
    ap.add_argument(
        "--nonlinear-solver",
        choices=["core-newton", "petsc-snes"],
        default="core-newton",
        help="core-newton = guarded Core serial Newton; petsc-snes = persistent "
             "application-owned SNES with the recovered 2026-06-27 settings",
    )
    ap.add_argument(
        "--formulation",
        choices=["fbar", "local-pressure", "std-kappa"],
        default="fbar",
        help="fbar = historical single-field F-bar path; local-pressure = "
             "standard Hex8 material with one P0 pressure condensed per element; "
             "std-kappa = standard Hex8 material with the pointwise kappa "
             "penalty (paper Eq. 3 volumetric form)",
    )
    ap.add_argument(
        "--isotropic",
        action="store_true",
        help="retired forensic control; the public benchmark driver rejects "
             "this physical-parameter change",
    )
    ap.add_argument(
        "--mass",
        choices=["consistent", "lumped"],
        default="consistent",
        help="consistent = assembled Q1-Hex8 mass matrix (default; the "
             "local Simula reference convention); lumped = row-summed diagonal "
             "mass (historical CoupFE path)",
    )
    ap.add_argument(
        "--material-eta",
        type=float,
        default=100.0,
        metavar="PA_S",
        help="passive material viscosity eta in Pa s; the public benchmark "
             "driver requires the paper value 100",
    )
    ap.add_argument(
        "--tbar-laplace",
        default=None,
        metavar="PATH.npy",
        help="precomputed per-node Laplace transmural field (paper Eq. 13) "
             "replacing the analytic parametric t~ for fiber generation",
    )
    ap.add_argument(
        "--fiber-sampling",
        choices=["cg1", "gp-direct"],
        default="cg1",
        help="cg1 = interpolate nodal fiber/sheet directions then apply "
             "Gram--Schmidt (historical); gp-direct = evaluate the analytic "
             "fiber rule at every Hex8 Gauss point",
    )
    ap.add_argument(
        "--element-evaluation",
        choices=["joint", "split"],
        default="joint",
        help=(
            "joint = cached paired material R/K evaluation (default); "
            "split = residual-only material evaluation until a tangent is needed"
        ),
    )
    ap.add_argument(
        "--viscous-evidence-out",
        default=None,
        metavar="PATH.npz",
        help=(
            "retired forensic eta-split spelling; the public benchmark driver "
            "rejects this physical-parameter sweep"
        ),
    )
    ap.add_argument(
        "--viscous-evidence-start",
        type=float,
        default=None,
        metavar="SECONDS",
        help=argparse.SUPPRESS,
    )
    ap.add_argument(
        "--viscous-evidence-end",
        type=float,
        default=None,
        metavar="SECONDS",
        help=argparse.SUPPRESS,
    )
    ap.add_argument("--apex-offset", type=float, default=0.2, dest="apex_offset",
                    help="non-degenerate apex: last ring at mu=-pi+offset (rad), small open hole; "
                         "0 = collapsed-apex research variant (degenerate)")
    ap.add_argument("--out", default=None)
    ap.add_argument(
        "--build-dir",
        help="persistent generated-kernel build directory (default: process temporary directory)",
    )
    ap.add_argument(
        "--epicardial-normal-mode",
        choices=["facet", "nodal-smoothed"],
        default="facet",
        dest="epicardial_normal_mode",
        help="epicardial Robin normal field: 'facet' (default, historical per-Gauss facet normal, "
             "bit-identical) or 'nodal-smoothed' (area-weighted Q1-interpolated nodal normal; "
             "a mechanism diagnostic, not a benchmark-faithful operator)",
    )
    return ap


def main(argv=None):
    ap = _parser()
    args = ap.parse_args(argv)
    try:
        configuration = benchmark_configuration(args.benchmark_step, args.case)
    except ValueError as error:
        ap.error(str(error))
    if not np.isfinite(args.apex_offset) or args.apex_offset < 0.0:
        ap.error("--apex-offset must be finite and nonnegative")
    if args.ncore < 4 or args.ncore % 2:
        ap.error("--ncore must be an even integer >= 4")
    if args.nradial < 1:
        ap.error("--nradial must be positive")
    if not 0.1 < args.core_half_width < 1.0 / np.sqrt(2.0):
        ap.error("--core-half-width must lie in (0.1, 1/sqrt(2))")
    if not np.isfinite(args.tip_refine) or not 1.0 <= args.tip_refine <= 8.0:
        ap.error("--tip-refine must lie in [1, 8]")
    if not np.isfinite(args.material_eta):
        ap.error("--material-eta must be finite")
    if args.material_eta != 100.0:
        ap.error(
            "the benchmark physical parameters are fixed: --material-eta "
            "must equal the paper value 100 Pa s; parameter sweeps are not run"
        )
    if args.isotropic:
        ap.error(
            "the benchmark physical parameters are fixed: --isotropic is a "
            "historical forensic control and is not run by the public driver"
        )
    if (
        args.mesh_topology == "polar-ring"
        and args.formulation == "local-pressure"
        and args.apex_offset <= 0.0
    ):
        ap.error(
            "--formulation local-pressure requires a nondegenerate open apex "
            "(--apex-offset > 0)"
        )
    if (
        args.viscous_evidence_out is not None
        or args.viscous_evidence_start is not None
        or args.viscous_evidence_end is not None
    ):
        ap.error(
            "the accepted-state eta-split diagnostic is a historical forensic "
            "control and is not run by the public benchmark driver; the paper "
            "physical parameters are fixed and parameter sweeps are not run"
        )
    try:
        _load_schedule_grid(args.dt, args.tend, args.load_horizon)
    except ValueError as error:
        ap.error(str(error))
    if args.out is None:
        args.out = (
            f"benchmark1_step{args.benchmark_step}_case{args.case}_full.npz"
            if args.benchmark_step != 0
            else f"case{args.case}_full.npz"
        )
    if args.integrator == "newmark":
        print(
            "WARNING: Newmark is an experimental discretization: inertia and Robin "
            "damping use Newmark kinematics, while material viscosity uses a backward "
            "strain difference.",
            flush=True,
        )
    if args.mesh_topology == "polar-ring" and args.apex_offset == 0.0:
        print(
            "WARNING: apex_offset=0 selects collapsed, degenerate apex elements; "
            "this is a research geometry variant.",
            flush=True,
        )
    if args.mesh_topology == "polar-ring" and args.apex_offset > 0.0:
        print(
            "WARNING: polar-ring with a positive apex offset selects the "
            "historical truncated open-tip domain, which is not the Benchmark 1 "
            "geometry; retained results from it are archive records, not "
            "current benchmark evidence.",
            flush=True,
        )

    if args.build_dir is None:
        temporary_build = tempfile.TemporaryDirectory(prefix="coupfe-cardiac-build-")
        atexit.register(temporary_build.cleanup)
        build_dir = temporary_build.name
    else:
        build_dir = str(Path(args.build_dir).expanduser())

    if args.mesh_topology == "closed-multiblock":
        mesh = geom.build_closed_mesh(
            args.nt,
            args.ncore,
            args.nradial,
            core_half_width=args.core_half_width,
            tip_refine=args.tip_refine,
            flip_helix=not args.raw_helix,
        )
        print(
            "geometry: closed square-to-disk topology; no apex cut surface",
            flush=True,
        )
    else:
        mesh = geom.build_mesh(
            args.nt,
            args.nmu,
            args.ntheta,
            flip_helix=not args.raw_helix,
            apex_offset=args.apex_offset,
        )
    fiber_tbar = np.asarray(mesh.param[:, 0], dtype=float).copy()
    tbar_definition = "analytic_parametric"
    tbar_source_filename = ""
    tbar_source_sha256 = ""
    tbar_metadata_filename = ""
    tbar_metadata_sha256 = ""
    tbar_metadata_schema = ""
    if args.tbar_laplace:
        tbar_path = Path(args.tbar_laplace).expanduser().resolve()
        if not tbar_path.is_file():
            ap.error(f"--tbar-laplace does not name a file: {args.tbar_laplace}")
        try:
            tbar_provenance = _load_laplace_tbar_metadata(mesh, tbar_path)
            fiber_tbar = _apply_laplace_tbar(
                mesh, tbar_path, asign=-1.0 if not args.raw_helix else 1.0
            )
        except (OSError, ValueError) as error:
            ap.error(str(error))
        tbar_definition = "laplace_presolved"
        tbar_source_filename = tbar_path.name
        tbar_source_sha256 = tbar_provenance["field_sha256"]
        tbar_metadata_filename = tbar_provenance["metadata_filename"]
        tbar_metadata_sha256 = tbar_provenance["metadata_sha256"]
        tbar_metadata_schema = tbar_provenance["metadata_schema"]
    if args.perturb > 0:
        rng = np.random.default_rng(0)
        mesh.nodes[:] = mesh.nodes + rng.normal(0.0, args.perturb, mesh.nodes.shape)
        if mesh.topology == "closed_multiblock_disk":
            geom.update_structural_frames(
                mesh,
                fiber_tbar,
                asign=-1.0 if not args.raw_helix else 1.0,
            )
        print(f"perturbed mesh nodes with σ={args.perturb:.1e} m", flush=True)
    dpn = 3
    ndof = dpn * mesh.n_node
    print(
        f"mesh: {mesh.n_node} nodes, {mesh.n_elem} hexes, ndof={ndof}  "
        "(displacement system)"
    )
    pre_solve_audit = {
        "geometry": audit_geometry(
            mesh,
            require_closed=args.mesh_topology == "closed-multiblock",
        )
    }
    print(
        "pre-solve geometry audit: "
        + json.dumps(pre_solve_audit["geometry"], sort_keys=True),
        flush=True,
    )

    grp, elem, ta_idx, local_pressure = build_group(
        mesh,
        args.dt,
        build_dir=build_dir,
        formulation=args.formulation,
        evaluation_mode=args.element_evaluation,
        fiber_sampling=args.fiber_sampling,
        fiber_tbar=fiber_tbar,
        fiber_asign=-1.0 if not args.raw_helix else 1.0,
        material_parameters=configuration.material_parameters,
    )
    material_property_names = tuple(CardiacHex8()._mat.props.keys())
    fiber_sampling_label = {
        "cg1": "cg1_gram_schmidt",
        "gp-direct": "gp_direct_rule",
    }[args.fiber_sampling]
    print(f"fiber sampling: {fiber_sampling_label}", flush=True)
    parameter_variant = "benchmark_eta"
    print(
        f"material viscosity: eta={args.material_eta:.9g} Pa s "
        f"({parameter_variant}; paper parameter fixed)",
        flush=True,
    )
    if args.mass == "consistent":
        Mcoo = consistent_mass_coo(
            mesh.nodes, mesh.elems, DENSITY, dof_per_node=dpn
        )
        if args.integrator == "newmark":
            inertia = ConsistentNewmarkInertia(*Mcoo, ndof, beta=0.25, gamma=0.5)
        else:
            inertia = ConsistentMassInertia(*Mcoo, ndof, damping=0.0)
        print("mass: consistent Q1-Hex8 (assembled, default)", flush=True)
    else:
        M = _hex_lumped_mass(mesh.nodes, mesh.elems, DENSITY, ndof, dof_per_node=dpn)
        if args.integrator == "newmark":
            inertia = NewmarkInertia(M, ndof, beta=0.25, gamma=0.5)
        else:
            inertia = InertiaOperator(M, ndof, damping=0.0)
        print("mass: lumped (row-summed)", flush=True)
    _epi_mode = {
        "facet": "normal",
        "nodal-smoothed": "normal-smoothed",
    }[args.epicardial_normal_mode]
    robin = RobinOperator(
        mesh.nodes, ndof,
        [(mesh.facets_base, A_TOP, B_TOP, "full"),
         (mesh.facets_epi, A_EPI, B_EPI, _epi_mode)],
        dof_per_node=dpn,
        kinematics=inertia if args.integrator == "newmark" else None,
    )
    pre_solve_audit["robin"] = audit_robin(robin)
    ops = [grp]
    if local_pressure is not None:
        ops.append(local_pressure)
        print(
            "formulation: standard Hex8 material (kappa=0) + one exactly "
            "condensed element pressure, K=1.000e+06 Pa",
            flush=True,
        )
    elif args.formulation == "std-kappa":
        print(
            "formulation: standard Hex8 material with pointwise kappa "
            f"penalty, kappa={HO_PARAMS['kappa']:.3e} Pa "
            "(paper Eq. 3 volumetric form)",
            flush=True,
        )
    else:
        print("formulation: historical Hex8 F-bar penalty", flush=True)
    print(
        "compiled material element evaluation: "
        f"{args.element_evaluation} "
        f"(residual-only available="
        f"{bool(getattr(elem, 'has_element_r_batch', False))})",
        flush=True,
    )
    ops.extend([inertia, robin])

    # Step 0 retains the historical split; Step 2 Case B applies both loads.
    times, tau, pres, load_horizon = _benchmark_load_histories(
        args.case,
        args.dt,
        args.tend,
        args.load_horizon,
        benchmark_step=args.benchmark_step,
    )
    pressure_op = None
    if configuration.pressure_enabled:
        # parent-element centroid per endo facet -> robust out-of-solid normal
        interior = mesh.nodes[mesh.elems[mesh.facets_endo_elem]].mean(axis=1)
        pressure_op = FollowerPressureOperator(mesh.nodes, ndof, mesh.facets_endo,
                                               dof_per_node=dpn, interior=interior)
        ops.append(pressure_op)
        pre_solve_audit["pressure"] = audit_pressure(
            mesh, dof_per_node=dpn
        )
    print(
        "pre-solve boundary audit: "
        + json.dumps(pre_solve_audit, sort_keys=True),
        flush=True,
    )
    print(
        f"BENCHMARK 1 STEP {args.benchmark_step} CASE {args.case} "
        f"({configuration.load_contract}): "
        f"tau_peak={tau.max():.1f}Pa  p_peak={pres.max():.1f}Pa",
        flush=True,
    )

    petsc_solver = None
    if args.nonlinear_solver == "petsc-snes":
        petsc_solver = PetscSnesSolver()
        atexit.register(petsc_solver.close)
        print(
            "nonlinear solver: PETSc SNES newtonls + bt, preonly + LU "
            "(recovered 2026-06-27 settings)",
            flush=True,
        )
    else:
        print("nonlinear solver: guarded Core Newton", flush=True)
    p0_location = locate_hex8_point(mesh.nodes, mesh.elems, geom.P0)
    p1_location = locate_hex8_point(mesh.nodes, mesh.elems, geom.P1)
    print(
        "Hex8 output sampling: "
        f"p0 elem={p0_location.element_index} "
        f"|recon|={p0_location.reconstruction_error * 1e3:.3e}mm; "
        f"p1 elem={p1_location.element_index} "
        f"|recon|={p1_location.reconstruction_error * 1e3:.3e}mm",
        flush=True,
    )

    U = np.zeros(ndof)
    hist_u0 = np.zeros((len(times), 3))
    hist_u1 = np.zeros((len(times), 3))
    load = validate_load_histories(configuration, times, tau, pres)
    n_peak = int(np.argmax(np.abs(load)))
    U_peak = np.zeros(ndof)
    step_diagnostics = []
    t0 = time.time()
    completed_steps = 0
    for n in range(1, len(times)):
        elem.props[ta_idx] = tau[n]
        if pressure_op is not None:
            pressure_op.p = float(pres[n])
        U_pred = inertia.predictor(args.dt)
        # Solver failures propagate to the CLI before a completed archive can
        # be written.
        if petsc_solver is None:
            U, _state, nit = checked_newton_solve(
                ops,
                U_pred,
                None,
                ndof,
                {},
                t=times[n],
                dt=args.dt,
                rtol=1e-8,
                maxit=40,
            )
            step_diagnostics.append(
                {
                    "time": float(times[n]),
                    "dt": float(args.dt),
                    "nonlinear_iterations": int(nit),
                    "active_tension_pa": float(tau[n]),
                    "pressure_pa": float(pres[n]),
                }
            )
        else:
            try:
                U, _state, diagnostics = petsc_solver.solve(
                    ops,
                    U_pred,
                    None,
                    ndof,
                    {},
                    t=times[n],
                    dt=args.dt,
                )
            except SnesSolveError as error:
                print(
                    "FAILED before state commit: "
                    f"step={n}/{len(times) - 1}, t={times[n]:.9g}s, "
                    f"material_eta={args.material_eta:.9g} Pa s, error={error}",
                    file=sys.stderr,
                    flush=True,
                )
                if error.diagnostics is not None:
                    print(
                        "failure diagnostics: "
                        + json.dumps(
                            error.diagnostics.as_dict(),
                            sort_keys=True,
                            separators=(",", ":"),
                            allow_nan=False,
                        ),
                        file=sys.stderr,
                        flush=True,
                    )
                raise
            nit = diagnostics.nonlinear_iterations
            diagnostic_record = diagnostics.as_dict()
            diagnostic_record.update(
                {
                    "active_tension_pa": float(tau[n]),
                    "pressure_pa": float(pres[n]),
                }
            )
            # Serialize each accepted step immediately. This catches any
            # nonfinite diagnostic at the step that produced it instead of
            # discovering it only after a long simulation has completed.
            json.dumps(diagnostic_record, allow_nan=False)
            step_diagnostics.append(diagnostic_record)
        completed_steps = n
        hist_u0[n] = interpolate_displacement(
            U, p0_location, dof_per_node=dpn
        )
        hist_u1[n] = interpolate_displacement(
            U, p1_location, dof_per_node=dpn
        )
        if n == n_peak:
            U_peak = U.copy()
        if n % max(1, (len(times) // 25)) == 0 or n < 3:
            umax = np.abs(U.reshape(-1, dpn)[:, :3]).max()   # over u comps only
            if configuration.load_contract == "active-stress-plus-pressure":
                load_label = f"tau={tau[n]:9.1f}Pa p={pres[n]:8.1f}Pa"
            else:
                load_label = f"load={load[n]:9.1f}Pa"
            print(
                f"  step {n:4d} t={times[n]:.3f}s  {load_label}  "
                f"nit={nit:2d}  |u|max={umax*1e3:7.3f}mm  "
                f"u0=({hist_u0[n,0]*1e3:+.2f},{hist_u0[n,1]*1e3:+.2f},"
                f"{hist_u0[n,2]*1e3:+.2f})mm",
                flush=True,
            )
    print(f"elapsed {time.time()-t0:.1f}s", flush=True)

    if petsc_solver is None:
        solver_configuration = {
            "name": "core-newton",
            "rtol": 1.0e-8,
            "max_it": 40,
            "independent_acceptance_atol": 1.0e-14,
        }
    else:
        solver_configuration = petsc_solver.configuration()
    solver_configuration.update(
        {
            "element_evaluation_mode": args.element_evaluation,
            "compiled_material_residual_only_available": bool(
                getattr(elem, "has_element_r_batch", False)
            ),
            "epicardial_normal_mode": args.epicardial_normal_mode,
        }
    )
    solver_configuration_json = json.dumps(
        solver_configuration, sort_keys=True, separators=(",", ":"), allow_nan=False
    )
    step_diagnostics_json = json.dumps(
        step_diagnostics, sort_keys=True, separators=(",", ":"), allow_nan=False
    )

    # The same checked Q1 kinematics are retained for both formulations. The
    # projected pressure is retained only when it actually participates in the
    # solve; it must not be mistaken for the F-bar penalty pressure.
    kinematic_probe = local_pressure or LocalPressureHex8Operator(
        mesh.nodes,
        mesh.elems,
        ndof,
        bulk_modulus=HO_PARAMS["kappa"],
        dof_per_node=dpn,
    )
    det_f_peak = kinematic_probe.deformation_jacobians(U_peak)
    element_pressure_peak = (
        local_pressure.element_pressure(U_peak)
        if local_pressure is not None
        else np.empty(0, dtype=float)
    )
    formulation_label = {
        "fbar": "hex8_fbar",
        "local-pressure": "hex8_local_pressure_p0_condensed_logj",
        "std-kappa": "hex8_standard_pointwise_kappa",
    }[args.formulation]
    validate_runtime_material_properties(
        configuration,
        material_property_names,
        elem.props,
        condensed_local_pressure=args.formulation == "local-pressure",
        active_tension_pa=tau[-1],
    )
    benchmark_archive_metadata = benchmark_metadata(
        configuration,
        material_parameters=configuration.material_parameters,
        activation_parameters=configuration.activation_parameters,
        pressure_parameters=configuration.pressure_parameters,
    )

    out_path = save_completed(
        args.out,
        completed_steps=completed_steps,
        expected_steps=len(times) - 1,
        **_runtime_metadata(),
        **benchmark_archive_metadata,
        benchmark_reproduction_profile=(
            "diagnostic-serial-backward-euler"
            if args.benchmark_step == 2
            else "not-applicable"
        ),
        times=times, tau=tau, pres=pres, u0=hist_u0, u1=hist_u1,
        p0=geom.P0, p1=geom.P1, U_peak=U_peak, n_peak=n_peak,
        nodes=mesh.nodes, elems=mesh.elems, fiber=mesh.fiber,
        facets_endo=mesh.facets_endo,
        case=args.case, dt=args.dt, integrator=args.integrator,
        formulation=formulation_label, density=DENSITY,
        material_kernel_formulation=(
            "fbar_mechanics" if args.formulation == "fbar" else "standard"
        ),
        material_model_id=MATERIAL_MODEL_ID,
        material_kappa_pa=(
            0.0
            if args.formulation == "local-pressure"
            else HO_PARAMS["kappa"]
        ),
        local_pressure_bulk_modulus_pa=(
            HO_PARAMS["kappa"] if args.formulation == "local-pressure" else 0.0
        ),
        isotropic=bool(args.isotropic),
        mass_representation=(
            "consistent_q1_hex8" if args.mass == "consistent" else "lumped_row_sum"
        ),
        det_f_gauss_peak=det_f_peak,
        element_pressure_peak_pa=element_pressure_peak,
        element_evaluation_mode=args.element_evaluation,
        compiled_material_residual_only_available=bool(
            getattr(elem, "has_element_r_batch", False)
        ),
        epicardial_normal_mode=args.epicardial_normal_mode,
        nonlinear_solver=args.nonlinear_solver,
        solver_configuration_json=solver_configuration_json,
        nonlinear_step_diagnostics_json=step_diagnostics_json,
        pre_solve_audit_json=json.dumps(
            pre_solve_audit,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        ),
        fiber_sampling=fiber_sampling_label,
        fiber_sampling_option=args.fiber_sampling,
        fiber_direction_reconstruction=geom.structural_direction_reconstruction(
            mesh.topology
        ),
        tbar_definition=tbar_definition,
        tbar_source_filename=tbar_source_filename,
        tbar_source_sha256=tbar_source_sha256,
        tbar_metadata_filename=tbar_metadata_filename,
        tbar_metadata_sha256=tbar_metadata_sha256,
        tbar_metadata_schema=tbar_metadata_schema,
        point_sampling="hex8_reference_isoparametric",
        p0_sampling_element=p0_location.element_index,
        p0_sampling_natural=np.asarray(p0_location.natural_coordinates),
        p0_sampling_weights=np.asarray(p0_location.weights),
        p0_sampling_reconstruction_error_m=p0_location.reconstruction_error,
        p1_sampling_element=p1_location.element_index,
        p1_sampling_natural=np.asarray(p1_location.natural_coordinates),
        p1_sampling_weights=np.asarray(p1_location.weights),
        p1_sampling_reconstruction_error_m=p1_location.reconstruction_error,
        viscous_rate="backward_difference",
        material_eta_pa_s=args.material_eta,
        viscous_term_active=bool(args.material_eta > 0.0),
        parameter_variant=parameter_variant,
        mesh_topology=mesh.topology,
        n_side=mesh.n_side,
        n_core=mesh.n_core,
        n_radial=mesh.n_radial,
        core_half_width=mesh.core_half_width,
        tip_refine=mesh.tip_refine,
        apex_offset=(args.apex_offset if args.mesh_topology == "polar-ring" else 0.0),
        flip_helix=not args.raw_helix,
        perturb=args.perturb, t_end=args.tend, load_horizon=load_horizon,
        a_top=A_TOP, b_top=B_TOP, a_epi=A_EPI, b_epi=B_EPI,
        n_t=mesh.n_t, n_mu=mesh.n_mu, n_theta=mesh.n_theta,
    )
    print(f"saved -> {out_path}", flush=True)


if __name__ == "__main__":
    main()
