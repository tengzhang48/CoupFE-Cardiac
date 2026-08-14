"""Generate the Q1-Hex Laplace transmural field for a closed LV mesh.

The cardiac benchmark defines the scalar transmural coordinate ``tbar`` by

    div(grad(tbar)) = 0          in the ventricular wall,
    tbar = 0                     on the endocardium,
    tbar = 1                     on the epicardium,
    grad(tbar) . n = 0           on the base.

This application-owned utility assembles that problem with the same Q1 Hex8
mesh used by the serial cardiac driver.  The zero-flux base condition is the
natural boundary condition of the weak form, so no base-boundary contribution
is added.  It deliberately uses SciPy rather than requiring FEniCSx/PETSc just
to prepare this driver input.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile

import numpy as np

import geometry


# Keep the established result/archive schema label.  Newly generated sidecars
# add the required nested mesh identity below; the execution loader rejects
# legacy v1 sidecars that lack it, while historical result records stay
# readable without rewriting their provenance.
SCHEMA = "coupfe-cardiac-laplace-tbar-v1"
MESH_IDENTITY_SCHEMA = "coupfe-cardiac-closed-mesh-identity-v1"
_GAUSS = 1.0 / np.sqrt(3.0)
_GAUSS_POINTS = tuple(
    (a, b, c)
    for a in (-_GAUSS, _GAUSS)
    for b in (-_GAUSS, _GAUSS)
    for c in (-_GAUSS, _GAUSS)
)
_XI = np.array([-1, 1, 1, -1, -1, 1, 1, -1], dtype=float)
_ETA = np.array([-1, -1, 1, 1, -1, -1, 1, 1], dtype=float)
_ZETA = np.array([-1, -1, -1, -1, 1, 1, 1, 1], dtype=float)


def _as_mesh_arrays(nodes, elems):
    nodes = np.asarray(nodes, dtype=float)
    elems = np.asarray(elems)
    if nodes.ndim != 2 or nodes.shape[1] != 3 or not np.all(np.isfinite(nodes)):
        raise ValueError("nodes must be a finite array with shape (n_node, 3)")
    if elems.ndim != 2 or elems.shape[1] != 8:
        raise ValueError("elems must have shape (n_elem, 8)")
    if not np.issubdtype(elems.dtype, np.integer):
        if not np.all(np.isfinite(elems)) or not np.all(elems == np.floor(elems)):
            raise ValueError("elems must contain integer node indices")
    elems = elems.astype(np.int64, copy=False)
    if len(nodes) == 0 or len(elems) == 0:
        raise ValueError("nodes and elems must be nonempty")
    if np.min(elems) < 0 or np.max(elems) >= len(nodes):
        raise ValueError("element connectivity contains an out-of-range node index")
    return nodes, elems


def _canonical_array_sha256(values, dtype):
    """Hash shape, dtype, and C-order bytes in a platform-stable encoding."""
    canonical_dtype = np.dtype(dtype)
    canonical = np.ascontiguousarray(
        np.asarray(values).astype(canonical_dtype, copy=False)
    )
    header = json.dumps(
        {
            "dtype": canonical_dtype.str,
            "shape": list(canonical.shape),
        },
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("ascii")
    digest = hashlib.sha256()
    digest.update(b"coupfe-cardiac-canonical-array-v1\0")
    digest.update(header)
    digest.update(b"\0")
    digest.update(canonical.tobytes(order="C"))
    return digest.hexdigest()


def closed_mesh_identity(nodes, elems):
    """Return an exact, portable identity for closed-mesh geometry/topology."""
    nodes, elems = _as_mesh_arrays(nodes, elems)
    return {
        "schema": MESH_IDENTITY_SCHEMA,
        "node_coordinates_sha256": _canonical_array_sha256(nodes, "<f8"),
        "element_connectivity_sha256": _canonical_array_sha256(elems, "<i8"),
    }


def _boundary_nodes(values, n_node, name):
    nodes = np.asarray(values)
    if nodes.ndim != 1:
        nodes = nodes.ravel()
    if not np.issubdtype(nodes.dtype, np.integer):
        if not np.all(np.isfinite(nodes)) or not np.all(nodes == np.floor(nodes)):
            raise ValueError(f"{name} must contain integer node indices")
    nodes = np.unique(nodes.astype(np.int64, copy=False))
    if len(nodes) == 0:
        raise ValueError(f"{name} must not be empty")
    if nodes[0] < 0 or nodes[-1] >= n_node:
        raise ValueError(f"{name} contains an out-of-range node index")
    return nodes


def solve_q1_hex_laplace(
    nodes,
    elems,
    endocardial_nodes,
    epicardial_nodes,
    *,
    bound_tolerance: float = 1.0e-9,
    boundary_tolerance: float = 1.0e-12,
    residual_tolerance: float = 1.0e-10,
):
    """Solve the nodal Q1 Hex8 Laplace problem and return field, diagnostics.

    Dirichlet values are zero on ``endocardial_nodes`` and one on
    ``epicardial_nodes``.  All otherwise unlabeled exterior faces, including
    the base of the closed benchmark mesh, receive the homogeneous natural
    Neumann condition.  The returned diagnostics are JSON-compatible.
    """
    from scipy import sparse
    from scipy.sparse.linalg import spsolve

    nodes, elems = _as_mesh_arrays(nodes, elems)
    endocardial_nodes = _boundary_nodes(
        endocardial_nodes, len(nodes), "endocardial_nodes"
    )
    epicardial_nodes = _boundary_nodes(
        epicardial_nodes, len(nodes), "epicardial_nodes"
    )
    if np.intersect1d(endocardial_nodes, epicardial_nodes).size:
        raise ValueError("endocardial and epicardial node sets must be disjoint")
    for name, value in (
        ("bound_tolerance", bound_tolerance),
        ("boundary_tolerance", boundary_tolerance),
        ("residual_tolerance", residual_tolerance),
    ):
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and positive")

    rows = []
    cols = []
    values = []
    minimum_det_j = np.inf
    maximum_det_j = -np.inf
    for element_index, connectivity in enumerate(elems):
        coordinates = nodes[connectivity]
        element_matrix = np.zeros((8, 8), dtype=float)
        for xi_gp, eta_gp, zeta_gp in _GAUSS_POINTS:
            derivatives = 0.125 * np.array(
                [
                    _XI * (1.0 + _ETA * eta_gp) * (1.0 + _ZETA * zeta_gp),
                    (1.0 + _XI * xi_gp) * _ETA * (1.0 + _ZETA * zeta_gp),
                    (1.0 + _XI * xi_gp) * (1.0 + _ETA * eta_gp) * _ZETA,
                ]
            )
            jacobian_transpose = derivatives @ coordinates
            det_j = float(np.linalg.det(jacobian_transpose))
            minimum_det_j = min(minimum_det_j, det_j)
            maximum_det_j = max(maximum_det_j, det_j)
            if not np.isfinite(det_j) or det_j <= 0.0:
                raise ValueError(
                    "nonpositive or nonfinite Gauss-point Jacobian in "
                    f"element {element_index}: detJ={det_j!r}"
                )
            gradients = np.linalg.solve(jacobian_transpose, derivatives)
            element_matrix += (gradients.T @ gradients) * det_j
        rows.append(np.repeat(connectivity, 8))
        cols.append(np.tile(connectivity, 8))
        values.append(element_matrix.ravel())

    row = np.concatenate(rows)
    col = np.concatenate(cols)
    value = np.concatenate(values)
    constrained = np.concatenate([endocardial_nodes, epicardial_nodes])
    prescribed = np.zeros(len(nodes), dtype=float)
    prescribed[epicardial_nodes] = 1.0
    rhs = np.zeros(len(nodes), dtype=float)

    constrained_row = np.isin(row, constrained)
    constrained_col = np.isin(col, constrained)
    contribution = (~constrained_row) & constrained_col
    np.add.at(
        rhs,
        row[contribution],
        -value[contribution] * prescribed[col[contribution]],
    )
    keep = (
        (~constrained_row & ~constrained_col)
        | (constrained_row & (row == col))
    )
    matrix = sparse.coo_matrix(
        (value[keep], (row[keep], col[keep])),
        shape=(len(nodes), len(nodes)),
    ).tolil()
    matrix[constrained, :] = 0.0
    matrix[constrained, constrained] = 1.0
    rhs[constrained] = prescribed[constrained]
    matrix = matrix.tocsr()

    field = np.asarray(spsolve(matrix, rhs), dtype=float)
    if field.shape != (len(nodes),) or not np.all(np.isfinite(field)):
        raise RuntimeError("Laplace solve returned a nonfinite or malformed field")

    residual = np.asarray(matrix @ field - rhs, dtype=float)
    residual_l2 = float(np.linalg.norm(residual))
    residual_inf = float(np.linalg.norm(residual, ord=np.inf))
    rhs_l2 = float(np.linalg.norm(rhs))
    residual_limit = residual_tolerance * max(1.0, rhs_l2)
    if residual_l2 > residual_limit:
        raise RuntimeError(
            "Laplace linear residual exceeds the declared tolerance: "
            f"{residual_l2:.6e} > {residual_limit:.6e}"
        )

    endocardial_error = float(np.max(np.abs(field[endocardial_nodes])))
    epicardial_error = float(
        np.max(np.abs(field[epicardial_nodes] - 1.0))
    )
    minimum = float(np.min(field))
    maximum = float(np.max(field))
    if endocardial_error > boundary_tolerance:
        raise RuntimeError("Laplace field does not satisfy tbar=0 on endocardium")
    if epicardial_error > boundary_tolerance:
        raise RuntimeError("Laplace field does not satisfy tbar=1 on epicardium")
    if minimum < -bound_tolerance or maximum > 1.0 + bound_tolerance:
        raise RuntimeError(
            f"Laplace field lies outside [0, 1]: [{minimum:.6e}, {maximum:.6e}]"
        )

    diagnostics = {
        "method": "Q1 Hex8 Galerkin; SciPy sparse direct solve",
        "linear_solver": "scipy.sparse.linalg.spsolve",
        "natural_boundary": "homogeneous Neumann on base",
        "n_node": int(len(nodes)),
        "n_element": int(len(elems)),
        "matrix_nnz": int(matrix.nnz),
        "minimum_gauss_det_j_m3": float(minimum_det_j),
        "maximum_gauss_det_j_m3": float(maximum_det_j),
        "linear_residual_l2": residual_l2,
        "linear_residual_inf": residual_inf,
        "linear_rhs_l2": rhs_l2,
        "linear_residual_limit": float(residual_limit),
        "minimum": minimum,
        "maximum": maximum,
        "max_abs_boundary_endo": endocardial_error,
        "max_abs_boundary_epi_minus_one": epicardial_error,
    }
    return field, diagnostics


def solve_closed_mesh_tbar(
    *,
    n_t: int = 2,
    n_core: int = 20,
    n_radial: int = 17,
    core_half_width: float = 0.36,
    tip_refine: float = 1.0,
):
    """Build the closed multiblock mesh and solve its Laplace field."""
    mesh = geometry.build_closed_mesh(
        n_t=n_t,
        n_core=n_core,
        n_radial=n_radial,
        core_half_width=core_half_width,
        tip_refine=tip_refine,
        flip_helix=True,
    )
    field, solver = solve_q1_hex_laplace(
        mesh.nodes,
        mesh.elems,
        np.unique(mesh.facets_endo),
        np.unique(mesh.facets_epi),
    )
    delta = field - mesh.param[:, 0]
    metadata = {
        "schema": SCHEMA,
        "mesh_topology": mesh.topology,
        "mesh_identity": closed_mesh_identity(mesh.nodes, mesh.elems),
        "mesh_parameters": {
            "n_t": int(mesh.n_t),
            "n_core": int(mesh.n_core),
            "n_radial": int(mesh.n_radial),
            "core_half_width": float(mesh.core_half_width),
            "tip_refine": float(mesh.tip_refine),
            "nodes": int(mesh.n_node),
            "elements": int(mesh.n_elem),
        },
        "boundary_conditions": {
            "endocardium": "Dirichlet tbar=0",
            "epicardium": "Dirichlet tbar=1",
            "base": "natural homogeneous Neumann",
        },
        "solver": solver,
        "max_abs_difference_from_layer_coordinate": float(
            np.max(np.abs(delta))
        ),
        "rms_difference_from_layer_coordinate": float(
            np.sqrt(np.mean(delta**2))
        ),
    }
    return field, metadata


def _sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_outputs_atomic(output_path, field, metadata):
    destination = Path(output_path).expanduser()
    if destination.suffix != ".npy":
        raise ValueError("output path must end in .npy")
    destination.parent.mkdir(parents=True, exist_ok=True)
    metadata_destination = destination.with_suffix(".meta.json")
    npy_temporary = None
    json_temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        ) as stream:
            npy_temporary = Path(stream.name)
            np.save(stream, np.asarray(field, dtype=float))
            stream.flush()
            os.fsync(stream.fileno())
        digest = _sha256_file(npy_temporary)
        record = dict(metadata)
        record["output_npy"] = destination.name
        record["sha256"] = digest
        payload = json.dumps(
            record, indent=2, sort_keys=True, allow_nan=False
        ) + "\n"
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=f".{metadata_destination.name}.",
            suffix=".tmp",
            dir=metadata_destination.parent,
            delete=False,
        ) as stream:
            json_temporary = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(npy_temporary, destination)
        npy_temporary = None
        os.replace(json_temporary, metadata_destination)
        json_temporary = None
    except BaseException:
        if npy_temporary is not None:
            npy_temporary.unlink(missing_ok=True)
        if json_temporary is not None:
            json_temporary.unlink(missing_ok=True)
        raise
    return destination, metadata_destination, record


def generate_closed_mesh_tbar(
    output_path,
    *,
    n_t: int = 2,
    n_core: int = 20,
    n_radial: int = 17,
    core_half_width: float = 0.36,
    tip_refine: float = 1.0,
):
    """Solve and atomically write the field plus its JSON metadata record."""
    field, metadata = solve_closed_mesh_tbar(
        n_t=n_t,
        n_core=n_core,
        n_radial=n_radial,
        core_half_width=core_half_width,
        tip_refine=tip_refine,
    )
    _, _, written_metadata = _write_outputs_atomic(
        output_path, field, metadata
    )
    return field, written_metadata


def _parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nt", type=int, default=2)
    parser.add_argument("--ncore", type=int, default=20)
    parser.add_argument("--nradial", type=int, default=17)
    parser.add_argument("--core-half-width", type=float, default=0.36)
    parser.add_argument("--tip-refine", type=float, default=1.0)
    parser.add_argument("--out", type=Path, required=True)
    return parser


def main(argv=None):
    args = _parser().parse_args(argv)
    _, metadata = generate_closed_mesh_tbar(
        args.out,
        n_t=args.nt,
        n_core=args.ncore,
        n_radial=args.nradial,
        core_half_width=args.core_half_width,
        tip_refine=args.tip_refine,
    )
    print(json.dumps(metadata, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
