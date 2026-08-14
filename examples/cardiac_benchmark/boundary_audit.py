"""Fail-fast geometry and boundary checks for Cardiac Benchmark 1.

These checks run before nonlinear dynamics.  They distinguish setup identity
from solver convergence: every exterior face must have exactly one intended
boundary role, except for the explicitly identified traction-free tip of the
historical open polar-ring mesh.  The closed benchmark topology must have no
apex cut, reference geometry measures must agree within a declared meshing
tolerance, and the undeformed pressure resultant must match the signed analytic
base projection for the closed topology.  The historical open polar-ring
topology is instead checked against the projected boundary of its declared
truncated pressure surface; it is never treated as the closed benchmark.
"""

from __future__ import annotations

import numpy as np

from pressure import FollowerPressureOperator


_HEX_FACES = np.asarray(
    [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [0, 1, 5, 4],
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 0, 4, 7],
    ],
    dtype=int,
)
_NAT = np.asarray(
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
_REFERENCE = {
    "wall_volume_cm3": 177.694656,
    "endocardial_area_cm2": 155.263738,
    "epicardial_area_cm2": 228.505027,
    "base_area_cm2": 18.577290,
}


def _shape_gradients(xi, eta, zeta):
    gradients = np.empty((8, 3))
    for a, (sx, sy, sz) in enumerate(_NAT):
        gradients[a, 0] = sx * (1 + sy * eta) * (1 + sz * zeta) / 8
        gradients[a, 1] = sy * (1 + sx * xi) * (1 + sz * zeta) / 8
        gradients[a, 2] = sz * (1 + sx * xi) * (1 + sy * eta) / 8
    return gradients


def _gauss_jacobians(nodes, elements):
    gauss = 1.0 / np.sqrt(3.0)
    points = [
        (xi, eta, zeta)
        for zeta in (gauss, -gauss)
        for eta in (gauss, -gauss)
        for xi in (gauss, -gauss)
    ]
    determinant = np.empty((len(elements), 8))
    condition = np.empty_like(determinant)
    for e, connectivity in enumerate(elements):
        coordinates = nodes[connectivity]
        for g, point in enumerate(points):
            jacobian = coordinates.T @ _shape_gradients(*point)
            determinant[e, g] = np.linalg.det(jacobian)
            condition[e, g] = np.linalg.cond(jacobian)
    return determinant, condition


def _extended_jacobians(nodes, elements):
    """Sample vertices, edge/face midpoints, center, and Gauss locations."""
    gauss = 1.0 / np.sqrt(3.0)
    levels = (-1.0, -gauss, 0.0, gauss, 1.0)
    determinant = []
    condition = []
    scaled = []
    worst_locations = []
    for e, connectivity in enumerate(elements):
        coordinates = nodes[connectivity]
        for zeta in levels:
            for eta in levels:
                for xi in levels:
                    jacobian = coordinates.T @ _shape_gradients(xi, eta, zeta)
                    det = np.linalg.det(jacobian)
                    column_norms = np.linalg.norm(jacobian, axis=0)
                    determinant.append(det)
                    condition.append(np.linalg.cond(jacobian))
                    scaled.append(det / np.prod(column_norms))
                    worst_locations.append((e, xi, eta, zeta))
    return (
        np.asarray(determinant),
        np.asarray(condition),
        np.asarray(scaled),
        worst_locations,
    )


def _quad_area(nodes, facets):
    gauss = 1.0 / np.sqrt(3.0)
    area = 0.0
    for facet in facets:
        coordinates = nodes[facet]
        for xi, eta in (
            (-gauss, -gauss),
            (gauss, -gauss),
            (gauss, gauss),
            (-gauss, gauss),
        ):
            dxi = np.asarray(
                [-(1 - eta), 1 - eta, 1 + eta, -(1 + eta)]
            ) / 4
            deta = np.asarray(
                [-(1 - xi), -(1 + xi), 1 + xi, 1 - xi]
            ) / 4
            area += np.linalg.norm(np.cross(dxi @ coordinates, deta @ coordinates))
    return float(area)


def _facet_set(facets):
    return {tuple(sorted(np.asarray(face, dtype=int).tolist())) for face in facets}


def _historical_polar_tip(mesh, exterior):
    """Identify the named polar-ring topology's intentional open tip.

    The release policy is deliberately narrower than "allow an unlabeled
    boundary": the face set must be the complete terminal annulus of a
    noncollapsed ``polar_ring`` mesh, with the expected face/node counts and
    parametric structure.  The caller decides whether that historical policy
    is enabled via ``require_closed=False``.
    """
    errors = []
    if getattr(mesh, "topology", None) != "polar_ring":
        return set(), None, ["mesh topology is not polar_ring"]

    param = np.asarray(getattr(mesh, "param", np.empty((0, 3))), dtype=float)
    if param.shape != (mesh.n_node, 3) or not np.all(np.isfinite(param[:, :2])):
        return set(), None, ["polar-ring parametric coordinates are invalid"]
    if mesh.n_t < 1 or mesh.n_theta < 3:
        return set(), None, ["polar-ring topology counts are invalid"]

    tip_mu = float(param[:, 1].min())
    tolerance = 64.0 * np.finfo(float).eps * max(1.0, abs(tip_mu))
    tip_faces = {
        face
        for face in exterior
        if np.all(np.abs(param[np.asarray(face, dtype=int), 1] - tip_mu) <= tolerance)
    }
    expected_faces = int(mesh.n_t * mesh.n_theta)
    expected_nodes = int((mesh.n_t + 1) * mesh.n_theta)
    tip_nodes = set().union(*(set(face) for face in tip_faces)) if tip_faces else set()

    if tip_mu <= -np.pi + tolerance:
        errors.append("terminal ring is collapsed at the apex")
    if len(tip_faces) != expected_faces:
        errors.append(
            f"terminal annulus has {len(tip_faces)} faces; expected {expected_faces}"
        )
    if len(tip_nodes) != expected_nodes:
        errors.append(
            f"terminal annulus has {len(tip_nodes)} nodes; expected {expected_nodes}"
        )
    for face in tip_faces:
        face_param = param[np.asarray(face, dtype=int)]
        if len(set(face)) != 4:
            errors.append("terminal annulus contains a collapsed quadrilateral")
            break
        if len(np.unique(np.round(face_param[:, 0], decimals=14))) != 2:
            errors.append("terminal quadrilateral does not span one transmural layer")
            break
        theta = face_param[:, 2]
        if not np.all(np.isfinite(theta)) or len(np.unique(np.round(theta, 14))) != 2:
            errors.append("terminal quadrilateral has invalid circumferential nodes")
            break

    return tip_faces, tip_mu, errors


def audit_geometry(mesh, *, require_closed, reference_tolerance=0.01):
    """Return a JSON-safe audit or fail before element compilation."""
    face_counts = {}
    for element in mesh.elems:
        for local_face in _HEX_FACES:
            face = tuple(sorted(element[local_face].tolist()))
            face_counts[face] = face_counts.get(face, 0) + 1
    exterior = {face for face, count in face_counts.items() if count == 1}
    nonmanifold = sum(count > 2 for count in face_counts.values())
    by_label = {
        "endocardium": _facet_set(mesh.facets_endo),
        "epicardium": _facet_set(mesh.facets_epi),
        "base": _facet_set(mesh.facets_base),
    }
    labeled = set().union(*by_label.values())
    unclassified = exterior - labeled
    overlap = sum(
        len(by_label[left] & by_label[right])
        for left, right in (
            ("endocardium", "epicardium"),
            ("endocardium", "base"),
            ("epicardium", "base"),
        )
    )
    determinant, condition = _gauss_jacobians(mesh.nodes, mesh.elems)
    extended_det, extended_condition, scaled_jacobian, locations = (
        _extended_jacobians(mesh.nodes, mesh.elems)
    )
    worst_condition = int(np.argmax(extended_condition))
    worst_scaled = int(np.argmin(scaled_jacobian))
    measures = {
        "wall_volume_cm3": float(determinant.sum() * 1.0e6),
        "endocardial_area_cm2": float(
            _quad_area(mesh.nodes, mesh.facets_endo) * 1.0e4
        ),
        "epicardial_area_cm2": float(
            _quad_area(mesh.nodes, mesh.facets_epi) * 1.0e4
        ),
        "base_area_cm2": float(_quad_area(mesh.nodes, mesh.facets_base) * 1.0e4),
    }
    relative = {
        key: float(measures[key] / reference - 1.0)
        for key, reference in _REFERENCE.items()
    }
    tip_faces = set()
    tip_mu = None
    tip_policy_errors = []
    boundary_policy = "closed_three_boundaries" if require_closed else "three_boundaries"
    if not require_closed and mesh.topology == "polar_ring":
        boundary_policy = "historical_polar_ring_traction_free_tip"
        tip_faces, tip_mu, tip_policy_errors = _historical_polar_tip(mesh, exterior)
    permitted_tip_faces = tip_faces if not tip_policy_errors else set()
    unexpected_unclassified = unclassified - permitted_tip_faces
    mislabeled_tip = tip_faces & labeled

    failures = []
    if nonmanifold:
        failures.append(f"{nonmanifold} faces have more than two owners")
    if require_closed and unclassified:
        failures.append(
            f"closed benchmark topology has {len(unclassified)} unclassified "
            "exterior faces"
        )
    elif unexpected_unclassified:
        failures.append(
            f"{len(unexpected_unclassified)} exterior faces are unclassified "
            "outside the declared boundary policy"
        )
    if tip_policy_errors:
        failures.append(
            "historical polar-ring traction-free tip is invalid: "
            + "; ".join(tip_policy_errors)
        )
    if mislabeled_tip:
        failures.append(
            f"{len(mislabeled_tip)} historical tip faces have a boundary label; "
            "the tip must remain explicitly traction-free"
        )
    if labeled - exterior:
        failures.append(f"{len(labeled - exterior)} labeled faces are not exterior")
    if overlap:
        failures.append(f"{overlap} faces have multiple boundary labels")
    if np.count_nonzero(determinant <= 0.0):
        failures.append("one or more Gauss-point reference Jacobians are nonpositive")
    if not np.all(np.isfinite(condition)):
        failures.append("one or more Gauss-point Jacobian conditions are nonfinite")
    if np.count_nonzero(extended_det <= 0.0):
        failures.append(
            "one or more vertex/edge/face/center Jacobians are nonpositive"
        )
    if not np.all(np.isfinite(extended_condition)):
        failures.append("one or more extended Jacobian conditions are nonfinite")
    if extended_condition.max() > 50.0:
        failures.append(
            f"extended Jacobian condition {extended_condition.max():.3f} exceeds 50"
        )
    if scaled_jacobian.min() < 0.05:
        failures.append(
            f"extended scaled Jacobian {scaled_jacobian.min():.3f} is below 0.05"
        )
    if require_closed:
        for key, error in relative.items():
            if abs(error) > reference_tolerance:
                failures.append(
                    f"{key} differs from the reference mesh by {100*error:.3f}%"
                )
    record = {
        "schema": "coupfe-cardiac-pre-solve-geometry-v1",
        "mesh_topology": mesh.topology,
        "require_closed": bool(require_closed),
        "boundary_policy": boundary_policy,
        "nodes": int(mesh.n_node),
        "elements": int(mesh.n_elem),
        "exterior_faces": len(exterior),
        "labeled_exterior_faces": len(labeled),
        "unclassified_exterior_faces": len(unclassified),
        "intentional_traction_free_tip_faces": len(permitted_tip_faces),
        "unexpected_unclassified_exterior_faces": len(unexpected_unclassified),
        "traction_free_tip_mu_rad": tip_mu,
        "nonexterior_labeled_faces": len(labeled - exterior),
        "multiply_labeled_faces": overlap,
        "nonmanifold_faces": nonmanifold,
        "gauss_jacobian_min_m3": float(determinant.min()),
        "gauss_jacobian_max_m3": float(determinant.max()),
        "gauss_jacobian_condition_max": float(condition.max()),
        "nonpositive_gauss_jacobians": int(np.count_nonzero(determinant <= 0.0)),
        "extended_jacobian_min_m3": float(extended_det.min()),
        "extended_jacobian_condition_max": float(extended_condition.max()),
        "extended_condition_worst_location": locations[worst_condition],
        "extended_scaled_jacobian_min": float(scaled_jacobian.min()),
        "extended_scaled_jacobian_worst_location": locations[worst_scaled],
        "nonpositive_extended_jacobians": int(
            np.count_nonzero(extended_det <= 0.0)
        ),
        "measures": measures,
        "reference_relative_differences": relative,
        "x_extent_mm": [
            float(mesh.nodes[:, 0].min() * 1.0e3),
            float(mesh.nodes[:, 0].max() * 1.0e3),
        ],
        "passed": not failures,
        "failures": failures,
    }
    if failures:
        raise RuntimeError("pre-solve geometry audit failed: " + "; ".join(failures))
    return record


def _polar_ring_area_vector_and_moment(mesh, node_ids):
    """Integrate an oriented virtual cap bounded by one polar ring.

    The historical pressure surface is open at both its basal and truncated
    terminal rings.  A constant-pressure surface resultant can therefore be
    checked independently by closing it with those two polygonal rings.  The
    ordering comes from the mesh's declared circumferential coordinate, not
    from the pressure operator or its facet normals.
    """
    node_ids = np.asarray(node_ids, dtype=int)
    theta = np.asarray(mesh.param[node_ids, 2], dtype=float)
    if len(node_ids) != mesh.n_theta or not np.all(np.isfinite(theta)):
        raise RuntimeError("polar pressure boundary ring has invalid nodes")
    order = np.argsort(theta)
    coordinates = np.asarray(mesh.nodes[node_ids[order]], dtype=float)
    center = coordinates.mean(axis=0)
    area_vector = np.zeros(3)
    moment = np.zeros(3)
    for current, following in zip(coordinates, np.roll(coordinates, -1, axis=0)):
        triangle_area_vector = 0.5 * np.cross(
            current - center,
            following - center,
        )
        triangle_centroid = (center + current + following) / 3.0
        area_vector += triangle_area_vector
        moment += np.cross(triangle_centroid, triangle_area_vector)
    if area_vector[0] < 0.0:
        area_vector *= -1.0
        moment *= -1.0
    if not np.all(np.isfinite(area_vector)) or area_vector[0] <= 0.0:
        raise RuntimeError("polar pressure boundary ring has invalid projection")
    return area_vector, moment


def _historical_pressure_reference(mesh):
    """Return the actual truncated polar surface's pressure reference.

    The endocardial side, a positive-x basal cap, and a negative-x terminal
    cap form a closed surface.  Their vector-area and moment identities give
    an independent reference for the side load on the *discrete* historical
    mesh, including the polygonal circumferential approximation.
    """
    if getattr(mesh, "topology", None) != "polar_ring":
        raise RuntimeError("historical pressure reference requires polar_ring")
    param = np.asarray(getattr(mesh, "param", np.empty((0, 3))), dtype=float)
    if param.shape != (mesh.n_node, 3):
        raise RuntimeError("polar pressure reference has invalid coordinates")
    endocardial_nodes = np.unique(np.asarray(mesh.facets_endo, dtype=int))
    endocardial_mu = param[endocardial_nodes, 1]
    if not len(endocardial_nodes) or not np.all(np.isfinite(endocardial_mu)):
        raise RuntimeError("polar pressure surface has invalid endocardial nodes")
    base_mu = float(endocardial_mu.max())
    tip_mu = float(endocardial_mu.min())
    tolerance = 64.0 * np.finfo(float).eps * max(
        1.0,
        abs(base_mu),
        abs(tip_mu),
    )
    base_nodes = endocardial_nodes[
        np.abs(param[endocardial_nodes, 1] - base_mu) <= tolerance
    ]
    tip_nodes = endocardial_nodes[
        np.abs(param[endocardial_nodes, 1] - tip_mu) <= tolerance
    ]
    if base_mu <= tip_mu + tolerance:
        raise RuntimeError("polar pressure surface has no axial extent")
    base_area, base_moment = _polar_ring_area_vector_and_moment(mesh, base_nodes)
    tip_area, tip_moment = _polar_ring_area_vector_and_moment(mesh, tip_nodes)
    resultant = base_area - tip_area
    moment = base_moment - tip_moment
    if resultant[0] <= 0.0 or not np.all(np.isfinite(resultant)):
        raise RuntimeError("truncated polar pressure surface has invalid projection")
    return {
        "resultant": resultant,
        "moment": moment,
        "base_area_vector": base_area,
        "tip_area_vector": tip_area,
        "base_mu": base_mu,
        "tip_mu": tip_mu,
    }


def audit_pressure(mesh, *, dof_per_node=3):
    """Check the signed zero-state unit-pressure resultant and moment.

    Closed meshes retain the analytic projected-base gate.  Historical
    polar-ring meshes use their independently reconstructed polygonal base and
    terminal rings, so a deliberately truncated pressure surface is checked as
    declared rather than compared to a different, closed domain.
    """
    ndof = dof_per_node * mesh.n_node
    interior = mesh.nodes[mesh.elems[mesh.facets_endo_elem]].mean(axis=1)
    operator = FollowerPressureOperator(
        mesh.nodes,
        ndof,
        mesh.facets_endo,
        p=1.0,
        dof_per_node=dof_per_node,
        interior=interior,
    )
    residual = operator.residual(np.zeros(ndof), None, t=0.0, dt=1.0)
    nodal = np.zeros((mesh.n_node, 3))
    np.add.at(
        nodal.reshape(-1),
        residual.gdofs,
        residual.values,
    )
    resultant = nodal.sum(axis=0)
    moment = np.cross(mesh.nodes, nodal).sum(axis=0)
    # The closed benchmark's endocardium plus a flat base disk is a closed
    # cavity.  Its side-surface unit-pressure resultant is therefore +x times
    # the projected base area under this residual convention.  The historical
    # polar mesh has an additional terminal ring, so it requires the separate
    # discrete-surface closure above.  Checking only magnitude in either case
    # would allow a completely reversed pressure load to pass.
    from geometry import MU_BASE_ENDO, RS_ENDO

    base_radius = RS_ENDO * abs(np.sin(MU_BASE_ENDO))
    analytic_closed_area = np.pi * base_radius**2
    historical_reference = None
    if mesh.topology == "closed_multiblock_disk":
        pressure_surface_policy = "closed_analytic_base_projection"
        expected_resultant = np.asarray([analytic_closed_area, 0.0, 0.0])
        expected_moment = np.zeros(3)
    elif mesh.topology == "polar_ring":
        pressure_surface_policy = "historical_truncated_polar_surface"
        historical_reference = _historical_pressure_reference(mesh)
        expected_resultant = historical_reference["resultant"]
        expected_moment = historical_reference["moment"]
    else:
        raise RuntimeError(
            "pre-solve pressure audit does not support mesh topology "
            f"{mesh.topology!r}"
        )
    expected_scale = np.linalg.norm(expected_resultant)
    expected_axial = expected_resultant[0]
    relative_magnitude_error = abs(
        np.linalg.norm(resultant) / expected_scale - 1.0
    )
    signed_axial_ratio = resultant[0] / expected_axial
    relative_signed_axial_error = abs(signed_axial_ratio - 1.0)
    relative_resultant_error = (
        np.linalg.norm(resultant - expected_resultant) / expected_scale
    )
    transverse = np.linalg.norm(
        resultant[1:] - expected_resultant[1:]
    ) / expected_scale
    moment_scale = expected_scale * 0.1
    normalized_moment = np.linalg.norm(moment - expected_moment) / moment_scale
    failures = []
    if not np.all(np.isfinite(resultant)) or not np.all(np.isfinite(moment)):
        failures.append("pressure resultant or moment is nonfinite")
    if relative_signed_axial_error > 5.0e-3:
        failures.append(
            "pressure resultant signed axial error is "
            f"{relative_signed_axial_error:.3e}"
        )
    if transverse > 5.0e-10:
        failures.append(f"pressure resultant transverse fraction is {transverse:.3e}")
    if normalized_moment > 5.0e-10:
        failures.append(f"pressure resultant normalized moment is {normalized_moment:.3e}")
    record = {
        "schema": "coupfe-cardiac-pre-solve-pressure-v1",
        "pressure_surface_policy": pressure_surface_policy,
        "unit_pressure_resultant_N": resultant.tolist(),
        "expected_unit_pressure_resultant_N": expected_resultant.tolist(),
        "unit_pressure_moment_Nm": moment.tolist(),
        "expected_unit_pressure_moment_Nm": expected_moment.tolist(),
        # Retained as the closed-benchmark comparator for both topologies.  It
        # is the audit target only when pressure_surface_policy is closed.
        "analytic_projected_base_area_m2": float(analytic_closed_area),
        "declared_surface_expected_axial_resultant_N": float(expected_axial),
        "expected_to_closed_analytic_axial_ratio": float(
            expected_axial / analytic_closed_area
        ),
        "relative_magnitude_error": float(relative_magnitude_error),
        "signed_axial_ratio": float(signed_axial_ratio),
        "relative_signed_axial_error": float(relative_signed_axial_error),
        "relative_resultant_error": float(relative_resultant_error),
        "transverse_fraction": float(transverse),
        "normalized_moment": float(normalized_moment),
        "passed": not failures,
        "failures": failures,
    }
    if historical_reference is not None:
        record.update(
            {
                "historical_base_mu_rad": historical_reference["base_mu"],
                "historical_tip_mu_rad": historical_reference["tip_mu"],
                "historical_base_ring_area_vector_m2": historical_reference[
                    "base_area_vector"
                ].tolist(),
                "historical_traction_free_tip_ring_area_vector_m2": (
                    historical_reference["tip_area_vector"].tolist()
                ),
            }
        )
    if failures:
        raise RuntimeError("pre-solve pressure audit failed: " + "; ".join(failures))
    return record


def audit_robin(operator):
    """Check assembled reference-surface spring/dashpot matrices."""
    k_error = float(abs(operator.Kmat - operator.Kmat.T).max())
    c_error = float(abs(operator.Cmat - operator.Cmat.T).max())
    failures = []
    if not np.isfinite(k_error) or k_error > 1.0e-6:
        failures.append(f"Robin spring symmetry error is {k_error:.3e}")
    if not np.isfinite(c_error) or c_error > 1.0e-9:
        failures.append(f"Robin dashpot symmetry error is {c_error:.3e}")
    spring_groups = operator.rigid_spring_groups
    translation = sum(
        (group["translation_matrix"] for group in spring_groups),
        np.zeros((3, 3)),
    )
    base_rotation = sum(
        group["long_axis_rotation"]
        for group in spring_groups
        if group["mode"] == "full"
    )
    epicardial_rotation = sum(
        group["long_axis_rotation"]
        for group in spring_groups
        if group["mode"] == "normal"
    )
    rigid_body_stiffness = {
        "schema": "coupfe-cardiac-robin-rigid-stiffness-v1",
        "diagnostic_only": True,
        "interpretation": (
            "discrete boundary geometry diagnostic; not a validation or "
            "pass/fail target"
        ),
        "coordinate_system": "global Cartesian; long axis is +x",
        "translation": {
            "units": "N/m",
            "global_diagonal": [float(value) for value in np.diag(translation)],
        },
        "long_axis_rotation": {
            "units": "N m/rad",
            "kinematic_mode": "u = e_x cross X per unit radian",
            "base_full_vector": float(base_rotation),
            "epicardial_normal_only": float(epicardial_rotation),
            "total": float(base_rotation + epicardial_rotation),
        },
    }
    record = {
        "schema": "coupfe-cardiac-pre-solve-robin-v1",
        "active_dofs": int(len(operator.dofs)),
        "spring_symmetry_error": k_error,
        "dashpot_symmetry_error": c_error,
        "rigid_body_stiffness": rigid_body_stiffness,
        "passed": not failures,
        "failures": failures,
    }
    if failures:
        raise RuntimeError("pre-solve Robin audit failed: " + "; ".join(failures))
    return record
