"""Physics diagnostics at peak load for a completed cardiac result archive.

The report covers wall inflation, fiber stretch, and ``det(F)``. Current runs
store ``det(F)`` at all eight Hex8 Gauss points; that retained field is used in
preference to the older centroid-only reconstruction.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np


RESULT_SCHEMA = "coupfe-cardiac-result-v1"
_NAT = np.array(
    [
        [-1, -1, -1],
        [1, -1, -1],
        [1, 1, -1],
        [-1, 1, -1],
        [-1, -1, 1],
        [1, -1, 1],
        [1, 1, 1],
        [-1, 1, 1],
    ],
    dtype=float,
)


def F_centroid(X, u):
    """Return the deformation gradient at a Hex8 centroid.

    With ``J = dx/dxi``, a row-form natural gradient transforms as
    ``dN/dx = dN/dxi @ inv(J)``. The former transpose on ``inv(J)`` was only
    harmless for orthogonal mappings and is deliberately covered by a
    distorted affine-element regression test.
    """
    X = np.asarray(X, dtype=float)
    u = np.asarray(u, dtype=float)
    if X.shape != (8, 3) or u.shape != (8, 3):
        raise ValueError("Hex8 centroid kinematics require X and u with shape (8, 3)")
    if not np.all(np.isfinite(X)) or not np.all(np.isfinite(u)):
        raise ValueError("Hex8 centroid kinematics require finite inputs")
    dN_nat = _NAT / 8.0
    jacobian = X.T @ dN_nat
    determinant = float(np.linalg.det(jacobian))
    if not np.isfinite(determinant) or determinant <= 0.0:
        raise RuntimeError("Hex8 reference Jacobian is nonpositive or non-finite")
    dN_phys = dN_nat @ np.linalg.inv(jacobian)
    deformation_gradient = np.eye(3) + u.T @ dN_phys
    if not np.all(np.isfinite(deformation_gradient)):
        raise RuntimeError("Hex8 centroid deformation gradient is non-finite")
    return deformation_gradient


def _scalar(archive, key, path):
    if key not in archive:
        raise RuntimeError(f"{path} is missing required result field {key!r}")
    value = np.asarray(archive[key])
    if value.shape != ():
        raise RuntimeError(f"{path} has non-scalar result field {key!r}")
    return value.item()


def _finite(archive, key, path, *, shape=None):
    if key not in archive:
        raise RuntimeError(f"{path} is missing required result field {key!r}")
    value = np.asarray(archive[key], dtype=float)
    if shape is not None and value.shape != shape:
        raise RuntimeError(
            f"{path} has invalid {key!r} shape {value.shape}; expected {shape}"
        )
    if not np.all(np.isfinite(value)):
        raise RuntimeError(f"{path} has non-finite values in {key!r}")
    return value


def _nonnegative_integer_scalar(archive, key, path):
    value = _scalar(archive, key, path)
    if isinstance(value, (bool, np.bool_)):
        raise RuntimeError(f"{path} has invalid integer field {key!r}")
    try:
        numeric = float(value)
    except (TypeError, ValueError) as error:
        raise RuntimeError(f"{path} has invalid integer field {key!r}") from error
    if not np.isfinite(numeric) or not numeric.is_integer() or numeric < 0.0:
        raise RuntimeError(f"{path} has invalid integer field {key!r}")
    return int(numeric)


def _indices(archive, key, path, *, width=None):
    if key not in archive:
        raise RuntimeError(f"{path} is missing required result field {key!r}")
    raw = np.asarray(archive[key])
    if raw.ndim != 2 or (width is not None and raw.shape[1] != width):
        raise RuntimeError(f"{path} has invalid {key!r} connectivity")
    if raw.dtype.kind not in "iu":
        raise RuntimeError(f"{path} has non-integer {key!r} connectivity")
    return raw.astype(np.int64, copy=True)


def _validate_complete_history(archive, path):
    schema = str(_scalar(archive, "result_schema", path))
    if schema != RESULT_SCHEMA:
        raise RuntimeError(
            f"{path} uses unsupported result schema {schema!r}; expected {RESULT_SCHEMA!r}"
        )
    converged = _scalar(archive, "converged", path)
    if not isinstance(converged, (bool, np.bool_)):
        raise RuntimeError(f"{path} has non-boolean result field 'converged'")
    if not converged:
        raise RuntimeError(f"{path} is not marked as a completed solve")
    completed = _nonnegative_integer_scalar(archive, "completed_steps", path)
    expected = _nonnegative_integer_scalar(archive, "expected_steps", path)
    times = _finite(archive, "times", path)
    if (
        completed < 1
        or completed != expected
        or times.ndim != 1
        or len(times) != expected + 1
        or not np.all(np.diff(times) > 0.0)
    ):
        raise RuntimeError(
            f"{path} has an incomplete or inconsistent completed-step history"
        )
    n_peak = _nonnegative_integer_scalar(archive, "n_peak", path)
    if n_peak < 0 or n_peak >= len(times):
        raise RuntimeError(f"{path} has out-of-range n_peak")
    return times.copy(), n_peak


def analyze_result(npz):
    """Validate and return peak-load diagnostics without printing them."""
    path = Path(npz).expanduser()
    with np.load(path, allow_pickle=False) as archive:
        times, n_peak = _validate_complete_history(archive, path)
        nodes = _finite(archive, "nodes", path)
        if nodes.ndim != 2 or nodes.shape[1] != 3 or len(nodes) < 8:
            raise RuntimeError(f"{path} has invalid nodes")
        elems = _indices(archive, "elems", path, width=8)
        if len(elems) < 1 or np.any(elems < 0) or np.any(elems >= len(nodes)):
            raise RuntimeError(f"{path} has invalid Hex8 element connectivity")
        fiber = _finite(archive, "fiber", path, shape=(len(elems), 3))
        if np.any(np.linalg.norm(fiber, axis=1) <= 0.0):
            raise RuntimeError(f"{path} has zero-length element fibers")
        facets_endo = _indices(archive, "facets_endo", path)
        if (
            facets_endo.shape[1] < 3
            or len(facets_endo) < 1
            or np.any(facets_endo < 0)
            or np.any(facets_endo >= len(nodes))
        ):
            raise RuntimeError(f"{path} has invalid endocardial facets")
        displacement_flat = _finite(archive, "U_peak", path)
        if displacement_flat.ndim != 1 or len(displacement_flat) % len(nodes):
            raise RuntimeError(f"{path} has invalid U_peak")
        dof_per_node = len(displacement_flat) // len(nodes)
        if dof_per_node < 3:
            raise RuntimeError(f"{path} has fewer than three displacement components")
        displacement = displacement_flat.reshape(len(nodes), dof_per_node)[:, :3].copy()

        stored_det_f = None
        if "det_f_gauss_peak" in archive:
            stored_det_f = _finite(
                archive,
                "det_f_gauss_peak",
                path,
                shape=(len(elems), 8),
            ).copy()
            if np.any(stored_det_f <= 0.0):
                raise RuntimeError(f"{path} has nonpositive stored det_f_gauss_peak")
        pressure = None
        if "element_pressure_peak_pa" in archive:
            candidate = _finite(
                archive,
                "element_pressure_peak_pa",
                path,
            )
            if candidate.shape not in {(0,), (len(elems),)}:
                raise RuntimeError(f"{path} has invalid element_pressure_peak_pa shape")
            if candidate.shape == (len(elems),):
                pressure = candidate.copy()

    i4f = []
    centroid_det_f = []
    for element_index, connectivity in enumerate(elems):
        deformation_gradient = F_centroid(
            nodes[connectivity], displacement[connectivity]
        )
        right_cauchy_green = deformation_gradient.T @ deformation_gradient
        f0 = fiber[element_index] / np.linalg.norm(fiber[element_index])
        i4f.append(float(f0 @ right_cauchy_green @ f0))
        centroid_det_f.append(float(np.linalg.det(deformation_gradient)))
    i4f = np.asarray(i4f, dtype=float)
    centroid_det_f = np.asarray(centroid_det_f, dtype=float)
    if not np.all(np.isfinite(i4f)) or not np.all(np.isfinite(centroid_det_f)):
        raise RuntimeError(f"{path} produced non-finite peak diagnostics")

    if stored_det_f is None:
        det_f = centroid_det_f
        det_f_sampling = "reconstructed element centroid"
    else:
        det_f = stored_det_f
        det_f_sampling = "stored 8-Gauss-point field"

    endocardial_nodes = np.unique(facets_endo)
    reference = nodes[endocardial_nodes]
    endocardial_u = displacement[endocardial_nodes]
    radial_norm = np.linalg.norm(reference[:, 1:], axis=1, keepdims=True)
    if np.any(radial_norm <= 0.0):
        raise RuntimeError(f"{path} has an undefined endocardial radial direction")
    radial_direction = reference[:, 1:] / radial_norm
    radial_displacement = np.sum(endocardial_u[:, 1:] * radial_direction, axis=1)
    if not np.all(np.isfinite(radial_displacement)):
        raise RuntimeError(f"{path} produced non-finite radial displacement")

    return {
        "filename": path.name,
        "peak_index": n_peak,
        "peak_time_s": float(times[n_peak]),
        "fiber_stretch_i4f": i4f,
        "det_f": det_f,
        "centroid_det_f": centroid_det_f,
        "det_f_sampling": det_f_sampling,
        "endocardial_radial_displacement_m": radial_displacement,
        "element_pressure_peak_pa": pressure,
    }


def main(npz):
    diagnostics = analyze_result(npz)
    i4f = diagnostics["fiber_stretch_i4f"]
    det_f = diagnostics["det_f"]
    radial = diagnostics["endocardial_radial_displacement_m"]
    print(
        f"file: {diagnostics['filename']}   "
        f"(peak-load field, t={diagnostics['peak_time_s']:.6g} s)"
    )
    print(
        f"  fiber stretch I4f:  min={i4f.min():.4f}  "
        f"mean={i4f.mean():.4f}  max={i4f.max():.4f}"
    )
    print(
        "     fraction of elems with I4f>1 (fibers in tension): "
        f"{np.mean(i4f > 1.0):.2f}"
    )
    print(
        f"  J=det F ({diagnostics['det_f_sampling']}): "
        f"min={det_f.min():.4f}  mean={det_f.mean():.4f}  max={det_f.max():.4f}"
    )
    pressure = diagnostics["element_pressure_peak_pa"]
    if pressure is not None:
        print(
            "  condensed pressure:  "
            f"min={pressure.min():+.6e}  mean={pressure.mean():+.6e}  "
            f"max={pressure.max():+.6e} Pa"
        )
    print(
        f"  endo radial disp:   mean={radial.mean() * 1e3:+.3f} mm  "
        f"(>0 = inflation)  [{radial.min() * 1e3:+.2f},"
        f"{radial.max() * 1e3:+.2f}]"
    )
    return diagnostics


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "caseB_full.npz")
