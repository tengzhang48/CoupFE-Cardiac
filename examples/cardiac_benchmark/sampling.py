"""Hex8-aware point location and displacement interpolation.

The benchmark output points are material points in the reference
configuration.  This module therefore inverts each candidate element's
*reference* trilinear isoparametric map and uses the resulting Hex8 shape
functions to sample the nodal displacement.  Natural coordinates follow the
standard CoupFE/Abaqus node order::

    0 (-1,-1,-1)   1 (+1,-1,-1)   2 (+1,+1,-1)   3 (-1,+1,-1)
    4 (-1,-1,+1)   5 (+1,-1,+1)   6 (+1,+1,+1)   7 (-1,+1,+1)

Candidate elements are selected by reference-coordinate axis-aligned bounding
boxes.  A damped Newton solve then recovers ``(xi, eta, zeta)``.  A location is
accepted only when its inputs, Jacobian, natural coordinates, shape weights,
and reconstructed physical point are finite; the local Jacobian is
orientation-preserving and nonsingular; the natural point is inside the
closed reference cube (within the requested tolerance); and reconstruction
meets a scale-aware tolerance.  If a point lies on a shared face/edge/node,
all valid candidates are considered and the lowest element index is selected.
That rule makes boundary sampling deterministic without depending on search
or floating-point sort order.

This is intentionally application-owned code.  It does not tetrahedralize the
mesh and it does not change the finite-element interpolation space.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# Standard Hex8 natural coordinates in the connectivity order documented above.
HEX8_NATURAL_COORDINATES = np.array(
    [
        [-1.0, -1.0, -1.0],
        [+1.0, -1.0, -1.0],
        [+1.0, +1.0, -1.0],
        [-1.0, +1.0, -1.0],
        [-1.0, -1.0, +1.0],
        [+1.0, -1.0, +1.0],
        [+1.0, +1.0, +1.0],
        [-1.0, +1.0, +1.0],
    ],
    dtype=float,
)


@dataclass(frozen=True)
class Hex8PointLocation:
    """A checked inverse-isoparametric location in one Hex8 element.

    ``node_ids`` and ``weights`` use the element connectivity order above.
    ``reconstruction_error`` is the Euclidean reference-coordinate error in
    the caller's length units.  ``iterations`` counts accepted Newton updates.
    Tuple fields keep the record immutable and straightforward to serialize.
    """

    element_index: int
    node_ids: tuple[int, int, int, int, int, int, int, int]
    natural_coordinates: tuple[float, float, float]
    weights: tuple[float, float, float, float, float, float, float, float]
    reconstruction_error: float
    iterations: int


class _DegenerateHex8Error(RuntimeError):
    pass


class _InverseMapError(RuntimeError):
    pass


def hex8_shape(natural_coordinates):
    """Return the eight trilinear shape functions at ``(xi, eta, zeta)``."""
    natural = np.asarray(natural_coordinates, dtype=float)
    if natural.shape != (3,) or not np.all(np.isfinite(natural)):
        raise ValueError("natural coordinates must be a finite length-3 vector")
    return 0.125 * np.prod(
        1.0 + HEX8_NATURAL_COORDINATES * natural[None, :], axis=1
    )


def hex8_shape_derivatives(natural_coordinates):
    """Return ``dN_a/d(xi,eta,zeta)`` with shape ``(8, 3)``."""
    natural = np.asarray(natural_coordinates, dtype=float)
    if natural.shape != (3,) or not np.all(np.isfinite(natural)):
        raise ValueError("natural coordinates must be a finite length-3 vector")
    xi, eta, zeta = natural
    signs = HEX8_NATURAL_COORDINATES
    derivatives = np.empty((8, 3), dtype=float)
    derivatives[:, 0] = (
        0.125 * signs[:, 0] * (1.0 + signs[:, 1] * eta)
        * (1.0 + signs[:, 2] * zeta)
    )
    derivatives[:, 1] = (
        0.125 * (1.0 + signs[:, 0] * xi) * signs[:, 1]
        * (1.0 + signs[:, 2] * zeta)
    )
    derivatives[:, 2] = (
        0.125 * (1.0 + signs[:, 0] * xi)
        * (1.0 + signs[:, 1] * eta) * signs[:, 2]
    )
    return derivatives


def _checked_mesh(nodes, elements):
    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements)
    if nodes.ndim != 2 or nodes.shape[1] != 3 or len(nodes) == 0:
        raise ValueError("nodes must have finite shape (n_node, 3)")
    if not np.all(np.isfinite(nodes)):
        raise ValueError("nodes must have finite shape (n_node, 3)")
    if elements.ndim != 2 or elements.shape[1] != 8 or len(elements) == 0:
        raise ValueError("elements must have shape (n_element, 8)")
    if not np.issubdtype(elements.dtype, np.integer):
        raise ValueError("element connectivity must contain integer node ids")
    elements = elements.astype(np.int64, copy=False)
    if np.any(elements < 0) or np.any(elements >= len(nodes)):
        raise ValueError("element connectivity contains an out-of-range node id")
    return nodes, elements


def _checked_point(point):
    point = np.asarray(point, dtype=float)
    if point.shape != (3,) or not np.all(np.isfinite(point)):
        raise ValueError("point must be a finite length-3 vector")
    return point


def candidate_hex8_elements(nodes, elements, point, *, bbox_tolerance=None):
    """Return sorted Hex8 indices whose reference AABBs contain ``point``.

    ``bbox_tolerance`` is an absolute length in the coordinate units.  The
    default is a roundoff allowance based on the magnitudes of the supplied
    mesh and point; it is only a candidate-search expansion, not the final
    inside or reconstruction criterion.
    """
    nodes, elements = _checked_mesh(nodes, elements)
    point = _checked_point(point)
    coordinate_scale = max(
        float(np.max(np.abs(nodes))),
        float(np.max(np.abs(point))),
        np.finfo(float).tiny,
    )
    if bbox_tolerance is None:
        bbox_tolerance = 128.0 * np.finfo(float).eps * coordinate_scale
    if not np.isfinite(bbox_tolerance) or bbox_tolerance < 0.0:
        raise ValueError("bbox_tolerance must be finite and nonnegative")
    element_nodes = nodes[elements]
    lower = np.min(element_nodes, axis=1) - bbox_tolerance
    upper = np.max(element_nodes, axis=1) + bbox_tolerance
    inside = np.all((point >= lower) & (point <= upper), axis=1)
    return tuple(int(index) for index in np.flatnonzero(inside))


def _map_and_jacobian(element_nodes, natural):
    weights = hex8_shape(natural)
    derivatives = hex8_shape_derivatives(natural)
    mapped = weights @ element_nodes
    jacobian = element_nodes.T @ derivatives
    return weights, mapped, jacobian


def _check_jacobian(jacobian, *, element_index, jacobian_rtol):
    if jacobian.shape != (3, 3) or not np.all(np.isfinite(jacobian)):
        raise _DegenerateHex8Error(
            f"candidate Hex8 element {element_index} has a non-finite Jacobian"
        )
    determinant = float(np.linalg.det(jacobian))
    singular_values = np.linalg.svd(jacobian, compute_uv=False)
    if (
        not np.isfinite(determinant)
        or not np.all(np.isfinite(singular_values))
        or determinant <= 0.0
        or singular_values[0] <= 0.0
        or singular_values[-1] <= jacobian_rtol * singular_values[0]
    ):
        raise _DegenerateHex8Error(
            f"candidate Hex8 element {element_index} has a degenerate or "
            "orientation-reversing Jacobian"
        )


def _inverse_map(
    element_nodes,
    point,
    *,
    element_index,
    max_iterations,
    reconstruction_atol,
    reconstruction_rtol,
    jacobian_rtol,
):
    diameter = float(np.linalg.norm(np.ptp(element_nodes, axis=0)))
    coordinate_magnitude = max(
        float(np.max(np.abs(element_nodes))),
        float(np.max(np.abs(point))),
        np.finfo(float).tiny,
    )
    reconstruction_tolerance = (
        reconstruction_atol
        + reconstruction_rtol * diameter
        + 256.0 * np.finfo(float).eps * coordinate_magnitude
    )
    natural = np.zeros(3, dtype=float)
    updates = 0
    for _ in range(max_iterations + 1):
        weights, mapped, jacobian = _map_and_jacobian(element_nodes, natural)
        _check_jacobian(
            jacobian, element_index=element_index, jacobian_rtol=jacobian_rtol
        )
        residual = mapped - point
        residual_norm = float(np.linalg.norm(residual))
        if not np.isfinite(residual_norm):
            raise _InverseMapError("inverse map produced a non-finite residual")
        if residual_norm <= reconstruction_tolerance:
            return natural, weights, residual_norm, updates, reconstruction_tolerance
        if updates >= max_iterations:
            break
        try:
            step = np.linalg.solve(jacobian, residual)
        except np.linalg.LinAlgError as error:
            raise _DegenerateHex8Error(
                f"candidate Hex8 element {element_index} has a singular Jacobian"
            ) from error
        if not np.all(np.isfinite(step)):
            raise _InverseMapError("inverse map produced a non-finite Newton step")

        # Damping keeps a distorted trilinear map from accepting an explosive
        # full update.  The first trial is the ordinary Newton step.
        accepted = False
        best_natural = None
        best_norm = np.inf
        damping = 1.0
        for _line_search in range(16):
            trial = natural - damping * step
            if np.all(np.isfinite(trial)) and np.linalg.norm(trial, ord=np.inf) < 1.0e8:
                trial_weights = hex8_shape(trial)
                trial_mapped = trial_weights @ element_nodes
                trial_norm = float(np.linalg.norm(trial_mapped - point))
                if np.isfinite(trial_norm) and trial_norm < best_norm:
                    best_natural, best_norm = trial, trial_norm
                if np.isfinite(trial_norm) and (
                    trial_norm <= reconstruction_tolerance
                    or trial_norm < residual_norm * (1.0 - 1.0e-4 * damping)
                ):
                    natural = trial
                    accepted = True
                    break
            damping *= 0.5
        if not accepted:
            if best_natural is None or best_norm >= residual_norm:
                raise _InverseMapError("damped Newton inverse map did not make progress")
            natural = best_natural
        updates += 1
    raise _InverseMapError(
        f"inverse map did not converge within {max_iterations} Newton updates"
    )


def locate_hex8_point(
    nodes,
    elements,
    point,
    *,
    bbox_tolerance=None,
    inside_tolerance=1.0e-9,
    reconstruction_atol=0.0,
    reconstruction_rtol=1.0e-11,
    jacobian_rtol=1.0e-12,
    max_iterations=30,
):
    """Locate a reference point and return checked Hex8 interpolation data.

    Natural coordinates may exceed ``[-1, 1]`` only by ``inside_tolerance``.
    ``reconstruction_atol`` uses physical coordinate units, while
    ``reconstruction_rtol`` scales with the candidate element diameter.
    ``jacobian_rtol`` bounds the smallest-to-largest singular-value ratio.
    Degenerate/inverted candidates, points outside the mesh, and failed inverse
    maps raise ``RuntimeError`` rather than silently extrapolating.
    """
    nodes, elements = _checked_mesh(nodes, elements)
    point = _checked_point(point)
    for name, value, strictly_positive in (
        ("inside_tolerance", inside_tolerance, False),
        ("reconstruction_atol", reconstruction_atol, False),
        ("reconstruction_rtol", reconstruction_rtol, True),
        ("jacobian_rtol", jacobian_rtol, True),
    ):
        if not np.isfinite(value) or value < 0.0 or (strictly_positive and value == 0.0):
            qualifier = "positive" if strictly_positive else "nonnegative"
            raise ValueError(f"{name} must be finite and {qualifier}")
    if (
        isinstance(max_iterations, (bool, np.bool_))
        or not isinstance(max_iterations, (int, np.integer))
        or max_iterations < 1
    ):
        raise ValueError("max_iterations must be a positive integer")

    candidates = candidate_hex8_elements(
        nodes, elements, point, bbox_tolerance=bbox_tolerance
    )
    if not candidates:
        raise RuntimeError("point is outside the Hex8 mesh bounding boxes")

    valid = []
    failures = []
    for element_index in candidates:
        node_ids = elements[element_index]
        element_nodes = nodes[node_ids]
        try:
            natural, weights, _error, iterations, tolerance = _inverse_map(
                element_nodes,
                point,
                element_index=element_index,
                max_iterations=int(max_iterations),
                reconstruction_atol=float(reconstruction_atol),
                reconstruction_rtol=float(reconstruction_rtol),
                jacobian_rtol=float(jacobian_rtol),
            )
            if not np.all(np.isfinite(natural)):
                raise _InverseMapError("inverse map produced non-finite natural coordinates")
            if np.any(natural < -1.0 - inside_tolerance) or np.any(
                natural > 1.0 + inside_tolerance
            ):
                raise _InverseMapError("inverse-map solution is outside the reference cube")
            if (
                not np.all(np.isfinite(weights))
                or not np.isclose(np.sum(weights), 1.0, rtol=0.0, atol=2.0e-13)
            ):
                raise _InverseMapError("inverse map produced invalid Hex8 weights")
            reconstructed = weights @ element_nodes
            checked_error = float(np.linalg.norm(reconstructed - point))
            if not np.isfinite(checked_error) or checked_error > tolerance:
                raise _InverseMapError("inverse-map reconstruction check failed")
            _, _, jacobian = _map_and_jacobian(element_nodes, natural)
            _check_jacobian(
                jacobian,
                element_index=element_index,
                jacobian_rtol=float(jacobian_rtol),
            )
            valid.append(
                Hex8PointLocation(
                    element_index=element_index,
                    node_ids=tuple(int(value) for value in node_ids),
                    natural_coordinates=tuple(float(value) for value in natural),
                    weights=tuple(float(value) for value in weights),
                    reconstruction_error=checked_error,
                    iterations=iterations,
                )
            )
        except (_DegenerateHex8Error, _InverseMapError) as error:
            failures.append(str(error))

    if not valid:
        detail = "; ".join(failures[:3])
        suffix = f": {detail}" if detail else ""
        raise RuntimeError(
            "point is outside all candidate Hex8 cells or the candidates are "
            f"degenerate{suffix}"
        )
    return min(valid, key=lambda location: location.element_index)


def interpolate_displacement(displacement, location, *, dof_per_node=3):
    """Interpolate the first three nodal DOFs at a checked Hex8 location.

    ``displacement`` may be a flat global vector with ``dof_per_node`` entries
    per node, or an ``(n_node, dof_per_node)`` array.  The first three entries
    are interpreted as Cartesian displacement components.  Only the eight
    selected nodes are read, and non-finite selected values fail closed.
    """
    if not isinstance(location, Hex8PointLocation):
        raise TypeError("location must be a Hex8PointLocation")
    if (
        isinstance(dof_per_node, (bool, np.bool_))
        or not isinstance(dof_per_node, (int, np.integer))
        or dof_per_node < 3
    ):
        raise ValueError("dof_per_node must be an integer of at least 3")
    values = np.asarray(displacement, dtype=float)
    node_ids = np.asarray(location.node_ids, dtype=np.int64)
    if values.ndim == 1:
        required = int(dof_per_node) * (int(np.max(node_ids)) + 1)
        if len(values) < required:
            raise ValueError("flat displacement vector is too short for the location")
        component_ids = (
            int(dof_per_node) * node_ids[:, None] + np.arange(3, dtype=np.int64)
        )
        selected = values[component_ids]
    elif values.ndim == 2:
        if values.shape[1] < 3 or values.shape[0] <= int(np.max(node_ids)):
            raise ValueError("nodal displacement array is too small for the location")
        selected = values[node_ids, :3]
    else:
        raise ValueError("displacement must be a flat vector or a 2-D nodal array")
    if not np.all(np.isfinite(selected)):
        raise RuntimeError("selected nodal displacement contains non-finite values")
    result = np.asarray(location.weights, dtype=float) @ selected
    if result.shape != (3,) or not np.all(np.isfinite(result)):
        raise RuntimeError("Hex8 displacement interpolation produced a non-finite result")
    return result
