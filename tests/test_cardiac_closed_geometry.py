"""Pre-solve geometry and boundary-audit regression checks."""

from __future__ import annotations

import copy

import numpy as np
import pytest

import boundary_audit
import geometry
from robin import RobinOperator


_A_TOP = 1.0e5
_B_TOP = 5.0e3
_A_EPI = 1.0e8
_B_EPI = 5.0e3


def _robin_operator(mesh, groups=None):
    if groups is None:
        groups = [
            (mesh.facets_base, _A_TOP, _B_TOP, "full"),
            (mesh.facets_epi, _A_EPI, _B_EPI, "normal"),
        ]
    return RobinOperator(mesh.nodes, 3 * mesh.n_node, groups)


@pytest.fixture(scope="module")
def closed_mesh():
    return geometry.build_closed_mesh(
        n_t=1,
        n_core=12,
        n_radial=6,
        core_half_width=0.36,
    )


@pytest.fixture(scope="module")
def geometry_record(closed_mesh):
    # One percent is an explicit discretization tolerance against the retained
    # closed reference mesh, not a nonlinear-solver acceptance tolerance.
    return boundary_audit.audit_geometry(
        closed_mesh,
        require_closed=True,
        reference_tolerance=0.01,
    )


@pytest.fixture(scope="module")
def closed_mesh_with_interior_wall_layer():
    return geometry.build_closed_mesh(
        n_t=2,
        n_core=12,
        n_radial=6,
        core_half_width=0.36,
    )


def test_closed_wall_mapping_preserves_exact_endo_and_epi_ellipses():
    for rho in (0.0, 0.17, 0.63, 1.0):
        for theta in (-2.1, 0.0, 1.4):
            mu_endo = -np.pi + rho * (np.pi + geometry.MU_BASE_ENDO)
            mu_epi = -np.pi + rho * (np.pi + geometry.MU_BASE_EPI)
            np.testing.assert_allclose(
                geometry.closed_wall_point(0.0, rho, theta),
                geometry.point(0.0, mu_endo, theta),
                rtol=0.0,
                atol=0.0,
            )
            np.testing.assert_allclose(
                geometry.closed_wall_point(1.0, rho, theta),
                geometry.point(1.0, mu_epi, theta),
                rtol=0.0,
                atol=0.0,
            )


@pytest.mark.parametrize("rho", [0.0, 1.0], ids=["apex", "base"])
def test_closed_wall_apex_and_base_are_straight_cartesian_segments(rho):
    theta = 0.73
    endo = geometry.closed_wall_point(0.0, rho, theta)
    epi = geometry.closed_wall_point(1.0, rho, theta)
    for t in (0.125, 0.5, 0.875):
        np.testing.assert_allclose(
            geometry.closed_wall_point(t, rho, theta),
            (1.0 - t) * endo + t * epi,
            rtol=0.0,
            atol=2.0e-17,
        )


def test_four_layer_base_has_the_toolkit_p0_positive_x_reach():
    mesh = geometry.build_closed_mesh(
        n_t=4,
        n_core=12,
        n_radial=6,
        core_half_width=0.36,
    )
    base_nodes = np.unique(mesh.facets_base)
    positive_y_meridian = base_nodes[
        (mesh.nodes[base_nodes, 1] > 0.0)
        & (np.abs(mesh.nodes[base_nodes, 2]) < 1.0e-14)
    ]
    positive_y_meridian = positive_y_meridian[
        np.argsort(mesh.nodes[positive_y_meridian, 1])
    ]
    assert len(positive_y_meridian) == mesh.n_t + 1

    meridian = mesh.nodes[positive_y_meridian]
    segment = np.flatnonzero(
        (meridian[:-1, 1] <= geometry.P0[1])
        & (geometry.P0[1] <= meridian[1:, 1])
    )
    assert segment.tolist() == [2]
    lower, upper = meridian[segment[0]:segment[0] + 2]
    fraction = (geometry.P0[1] - lower[1]) / (upper[1] - lower[1])
    ray_intersection = lower + fraction * (upper - lower)

    endo = geometry.point(0.0, geometry.MU_BASE_ENDO, np.pi)
    epi = geometry.point(1.0, geometry.MU_BASE_EPI, np.pi)
    reference_fraction = (geometry.P0[1] - endo[1]) / (epi[1] - endo[1])
    reference_intersection = endo + reference_fraction * (epi - endo)
    np.testing.assert_allclose(
        ray_intersection,
        reference_intersection,
        rtol=0.0,
        atol=2.0e-17,
    )
    assert 1.0e6 * (ray_intersection[0] - geometry.P0[0]) == pytest.approx(
        113.98499726388012,
        abs=2.0e-9,
    )


def test_closed_mapping_with_interior_wall_layer_retains_positive_quality(
    closed_mesh_with_interior_wall_layer,
):
    gauss_determinants, _ = boundary_audit._gauss_jacobians(
        closed_mesh_with_interior_wall_layer.nodes,
        closed_mesh_with_interior_wall_layer.elems,
    )
    extended_determinants, conditions, scaled, _ = (
        boundary_audit._extended_jacobians(
            closed_mesh_with_interior_wall_layer.nodes,
            closed_mesh_with_interior_wall_layer.elems,
        )
    )
    assert np.min(gauss_determinants) > 0.0
    assert np.min(extended_determinants) > 0.0
    assert np.max(conditions) < 10.0
    assert np.min(scaled) > 0.25


def test_closed_multiblock_topology_has_only_the_three_benchmark_boundaries(
    closed_mesh, geometry_record
):
    assert closed_mesh.topology == "closed_multiblock_disk"
    assert (closed_mesh.n_node, closed_mesh.n_elem) == (914, 432)
    assert len(closed_mesh.facets_endo) == 432
    assert len(closed_mesh.facets_epi) == 432
    assert len(closed_mesh.facets_base) == 48

    assert geometry_record["passed"] is True
    assert geometry_record["exterior_faces"] == 912
    assert geometry_record["labeled_exterior_faces"] == 912
    assert geometry_record["unclassified_exterior_faces"] == 0
    assert geometry_record["nonexterior_labeled_faces"] == 0
    assert geometry_record["multiply_labeled_faces"] == 0
    assert geometry_record["nonmanifold_faces"] == 0
    np.testing.assert_allclose(
        geometry_record["x_extent_mm"],
        [-97.0, 26.47058823529412],
        rtol=0.0,
        atol=2.0e-12,
    )


def test_closed_multiblock_extended_jacobians_and_measures_pass_declared_gate(
    geometry_record,
):
    assert geometry_record["nonpositive_gauss_jacobians"] == 0
    assert geometry_record["nonpositive_extended_jacobians"] == 0
    assert geometry_record["gauss_jacobian_min_m3"] > 0.0
    assert geometry_record["extended_jacobian_min_m3"] > 0.0
    assert geometry_record["extended_jacobian_condition_max"] < 10.0
    assert geometry_record["extended_scaled_jacobian_min"] > 0.25

    relative = geometry_record["reference_relative_differences"]
    assert set(relative) == {
        "wall_volume_cm3",
        "endocardial_area_cm2",
        "epicardial_area_cm2",
        "base_area_cm2",
    }
    assert max(abs(value) for value in relative.values()) < 0.01


def test_closed_multiblock_unit_pressure_matches_projected_base_resultant(
    closed_mesh,
):
    record = boundary_audit.audit_pressure(closed_mesh)

    assert record["passed"] is True
    assert record["pressure_surface_policy"] == "closed_analytic_base_projection"
    assert record["expected_to_closed_analytic_axial_ratio"] == pytest.approx(1.0)
    assert record["relative_magnitude_error"] < 5.0e-3
    assert record["relative_signed_axial_error"] < 5.0e-3
    assert record["signed_axial_ratio"] > 0.0
    assert record["transverse_fraction"] < 5.0e-10
    assert record["normalized_moment"] < 5.0e-10
    np.testing.assert_allclose(
        record["expected_unit_pressure_resultant_N"],
        [record["analytic_projected_base_area_m2"], 0.0, 0.0],
        rtol=0.0,
        atol=0.0,
    )
    np.testing.assert_array_equal(
        record["expected_unit_pressure_moment_Nm"],
        [0.0, 0.0, 0.0],
    )
    np.testing.assert_allclose(
        np.linalg.norm(record["unit_pressure_resultant_N"]),
        record["analytic_projected_base_area_m2"],
        rtol=5.0e-3,
        atol=0.0,
    )


def test_historical_pressure_audit_uses_the_declared_truncated_surface():
    # This deliberately coarse mesh is the focused serial/MPI smoke setup.  Its
    # actual polygonal, open-tip pressure surface has only about 61% of the
    # closed analytic projection, so comparing it to the closed domain is not a
    # meaningful pressure-operator gate.
    mesh = geometry.build_mesh(
        n_t=1,
        n_mu=2,
        n_theta=4,
        apex_offset=0.2,
    )

    record = boundary_audit.audit_pressure(mesh)

    assert record["passed"] is True
    assert record["pressure_surface_policy"] == (
        "historical_truncated_polar_surface"
    )
    assert record["expected_to_closed_analytic_axial_ratio"] == pytest.approx(
        0.6091132493442383,
        rel=2.0e-14,
    )
    assert record["relative_signed_axial_error"] < 1.0e-12
    assert record["relative_resultant_error"] < 1.0e-12
    base = np.asarray(record["historical_base_ring_area_vector_m2"])
    tip = np.asarray(
        record["historical_traction_free_tip_ring_area_vector_m2"]
    )
    np.testing.assert_allclose(
        record["expected_unit_pressure_resultant_N"],
        base - tip,
        rtol=0.0,
        atol=1.0e-18,
    )
    np.testing.assert_allclose(
        record["unit_pressure_resultant_N"],
        record["expected_unit_pressure_resultant_N"],
        rtol=1.0e-12,
        atol=1.0e-18,
    )
    assert record["declared_surface_expected_axial_resultant_N"] < (
        0.7 * record["analytic_projected_base_area_m2"]
    )


def test_historical_polar_ring_records_only_its_traction_free_tip():
    mesh = geometry.build_mesh(
        n_t=1,
        n_mu=4,
        n_theta=8,
        apex_offset=0.2,
    )

    record = boundary_audit.audit_geometry(mesh, require_closed=False)

    assert record["passed"] is True
    assert record["boundary_policy"] == "historical_polar_ring_traction_free_tip"
    assert record["unclassified_exterior_faces"] == 8
    assert record["intentional_traction_free_tip_faces"] == 8
    assert record["unexpected_unclassified_exterior_faces"] == 0
    assert record["traction_free_tip_mu_rad"] == pytest.approx(-np.pi + 0.2)


def test_historical_tip_policy_does_not_allow_an_unlabeled_base_face():
    mesh = geometry.build_mesh(
        n_t=1,
        n_mu=4,
        n_theta=8,
        apex_offset=0.2,
    )
    broken = copy.deepcopy(mesh)
    broken.facets_base = broken.facets_base[1:]

    with pytest.raises(RuntimeError, match="outside the declared boundary policy"):
        boundary_audit.audit_geometry(broken, require_closed=False)


def test_historical_tip_policy_rejects_a_labeled_tip_face():
    mesh = geometry.build_mesh(
        n_t=1,
        n_mu=4,
        n_theta=8,
        apex_offset=0.2,
    )
    broken = copy.deepcopy(mesh)
    terminal_mu = broken.param[:, 1].min()
    tip_face = next(
        face
        for element in broken.elems
        for face in element[boundary_audit._HEX_FACES]
        if np.allclose(broken.param[face, 1], terminal_mu)
    )
    broken.facets_base = np.vstack([broken.facets_base, tip_face])

    with pytest.raises(RuntimeError, match="tip faces have a boundary label"):
        boundary_audit.audit_geometry(broken, require_closed=False)


def test_closed_policy_rejects_an_unclassified_exterior_face(closed_mesh):
    broken = copy.deepcopy(closed_mesh)
    broken.facets_base = broken.facets_base[1:]

    with pytest.raises(RuntimeError, match="closed benchmark topology has 1 unclassified"):
        boundary_audit.audit_geometry(broken, require_closed=True)


def test_geometry_audit_retains_nonmanifold_broken_control(closed_mesh):
    broken = copy.deepcopy(closed_mesh)
    broken.elems = np.vstack([broken.elems, broken.elems[0]])

    with pytest.raises(RuntimeError, match="faces have more than two owners"):
        boundary_audit.audit_geometry(broken, require_closed=True)


def test_geometry_audit_retains_multiply_labeled_broken_control(closed_mesh):
    broken = copy.deepcopy(closed_mesh)
    broken.facets_epi = np.vstack([broken.facets_epi, broken.facets_endo[0]])

    with pytest.raises(RuntimeError, match="faces have multiple boundary labels"):
        boundary_audit.audit_geometry(broken, require_closed=True)


def test_geometry_audit_retains_nonexterior_label_broken_control(closed_mesh):
    face_counts = {}
    for element in closed_mesh.elems:
        for local_face in boundary_audit._HEX_FACES:
            face = tuple(sorted(element[local_face].tolist()))
            face_counts[face] = face_counts.get(face, 0) + 1
    interior_face = next(face for face, count in face_counts.items() if count == 2)
    broken = copy.deepcopy(closed_mesh)
    broken.facets_base = np.vstack([broken.facets_base, interior_face])

    with pytest.raises(RuntimeError, match="labeled faces are not exterior"):
        boundary_audit.audit_geometry(broken, require_closed=True)


def test_geometry_audit_retains_invalid_jacobian_broken_control(closed_mesh):
    broken = copy.deepcopy(closed_mesh)
    broken.elems[0, [0, 1]] = broken.elems[0, [1, 0]]

    with pytest.raises(RuntimeError, match="Jacobians are nonpositive"):
        boundary_audit.audit_geometry(broken, require_closed=True)


def test_pressure_audit_rejects_reversed_resultant_with_correct_magnitude(
    closed_mesh, monkeypatch
):
    original = boundary_audit.FollowerPressureOperator

    class ReversedPressureOperator(original):
        def residual(self, U, state, t, dt):
            result = super().residual(U, state, t, dt)
            return type(result)(result.gdofs, -result.values)

    monkeypatch.setattr(
        boundary_audit,
        "FollowerPressureOperator",
        ReversedPressureOperator,
    )

    with pytest.raises(RuntimeError, match="signed axial error"):
        boundary_audit.audit_pressure(closed_mesh)


def test_historical_pressure_audit_rejects_reversed_surface_orientation(
    monkeypatch,
):
    mesh = geometry.build_mesh(
        n_t=1,
        n_mu=2,
        n_theta=4,
        apex_offset=0.2,
    )
    original = boundary_audit.FollowerPressureOperator

    class ReversedPressureOperator(original):
        def residual(self, U, state, t, dt):
            result = super().residual(U, state, t, dt)
            return type(result)(result.gdofs, -result.values)

    monkeypatch.setattr(
        boundary_audit,
        "FollowerPressureOperator",
        ReversedPressureOperator,
    )

    with pytest.raises(RuntimeError, match="signed axial error"):
        boundary_audit.audit_pressure(mesh)


def test_robin_audit_reports_rigid_stiffness_from_assembled_operator():
    mesh = geometry.build_closed_mesh(
        n_t=1,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )
    operator = _robin_operator(mesh)
    rigid = boundary_audit.audit_robin(operator)["rigid_body_stiffness"]

    assert rigid["schema"] == "coupfe-cardiac-robin-rigid-stiffness-v1"
    assert rigid["diagnostic_only"] is True
    assert rigid["translation"]["units"] == "N/m"
    assert rigid["long_axis_rotation"]["units"] == "N m/rad"

    translation = []
    for component in range(3):
        displacement = np.zeros(operator.ndof)
        displacement[component::3] = 1.0
        translation.append(float(displacement @ operator.Kmat @ displacement))
    np.testing.assert_allclose(
        rigid["translation"]["global_diagonal"],
        translation,
        rtol=2.0e-15,
        atol=1.0e-9,
    )

    rotation = np.zeros((mesh.n_node, 3))
    rotation[:, 1] = -mesh.nodes[:, 2]
    rotation[:, 2] = mesh.nodes[:, 1]
    rotation = rotation.ravel()
    expected_total = float(rotation @ operator.Kmat @ rotation)
    assert rigid["long_axis_rotation"]["total"] == pytest.approx(
        expected_total, rel=2.0e-15, abs=1.0e-13
    )

    base = _robin_operator(
        mesh,
        [(mesh.facets_base, _A_TOP, _B_TOP, "full")],
    )
    epicardium = _robin_operator(
        mesh,
        [(mesh.facets_epi, _A_EPI, _B_EPI, "normal")],
    )
    assert rigid["long_axis_rotation"]["base_full_vector"] == pytest.approx(
        float(rotation @ base.Kmat @ rotation), rel=2.0e-15, abs=1.0e-13
    )
    assert rigid["long_axis_rotation"][
        "epicardial_normal_only"
    ] == pytest.approx(
        float(rotation @ epicardium.Kmat @ rotation), rel=2.0e-15, abs=1.0e-13
    )


def test_robin_rigid_stiffness_is_invariant_to_through_wall_refinement():
    records = []
    for n_t in (1, 3):
        mesh = geometry.build_closed_mesh(
            n_t=n_t,
            n_core=4,
            n_radial=1,
            core_half_width=0.36,
        )
        records.append(
            boundary_audit.audit_robin(_robin_operator(mesh))[
                "rigid_body_stiffness"
            ]
        )

    np.testing.assert_allclose(
        records[0]["translation"]["global_diagonal"],
        records[1]["translation"]["global_diagonal"],
        rtol=2.0e-15,
        atol=1.0e-9,
    )
    for component in (
        "base_full_vector",
        "epicardial_normal_only",
        "total",
    ):
        assert records[0]["long_axis_rotation"][component] == pytest.approx(
            records[1]["long_axis_rotation"][component],
            rel=2.0e-15,
            abs=1.0e-13,
        )


def test_surface_refinement_reduces_epicardial_faceting_stiffness():
    rotations = []
    for n_core, n_radial in ((4, 1), (8, 2)):
        mesh = geometry.build_closed_mesh(
            n_t=1,
            n_core=n_core,
            n_radial=n_radial,
            core_half_width=0.36,
        )
        rigid = boundary_audit.audit_robin(_robin_operator(mesh))[
            "rigid_body_stiffness"
        ]
        rotations.append(
            rigid["long_axis_rotation"]["epicardial_normal_only"]
        )

    assert rotations[0] > 0.0
    assert 0.0 < rotations[1] < rotations[0]


def _epi_meridian_edges(mesh):
    """Epi facet edges that connect meridian rings, sorted apex-first."""
    mu = mesh.param[:, 1]
    edges = []
    for facet in mesh.facets_epi:
        for a, b in zip(facet, np.roll(facet, -1)):
            if abs(mu[a] - mu[b]) > 1.0e-12:
                edges.append(
                    (
                        0.5 * (mu[a] + mu[b]),
                        float(np.linalg.norm(mesh.nodes[a] - mesh.nodes[b])),
                    )
                )
    edges.sort(key=lambda item: item[0])
    return edges


def test_tip_refine_default_reproduces_the_uniform_mesh_bit_for_bit():
    reference = geometry.build_closed_mesh(
        n_t=2, n_core=12, n_radial=6, core_half_width=0.36,
    )
    explicit = geometry.build_closed_mesh(
        n_t=2, n_core=12, n_radial=6, core_half_width=0.36, tip_refine=1.0,
    )

    np.testing.assert_array_equal(reference.nodes, explicit.nodes)
    np.testing.assert_array_equal(reference.elems, explicit.elems)
    np.testing.assert_array_equal(reference.param, explicit.param)
    np.testing.assert_array_equal(reference.facets_endo, explicit.facets_endo)
    np.testing.assert_array_equal(reference.facets_epi, explicit.facets_epi)
    np.testing.assert_array_equal(reference.facets_base, explicit.facets_base)
    assert reference.tip_refine == 1.0


def test_tip_refine_shrinks_only_apex_adjacent_meridian_elements():
    uniform = geometry.build_closed_mesh(
        n_t=2, n_core=12, n_radial=6, core_half_width=0.36,
    )
    graded = geometry.build_closed_mesh(
        n_t=2, n_core=12, n_radial=6, core_half_width=0.36, tip_refine=2.5,
    )

    # Identical topology and labels; only node coordinates may move.
    np.testing.assert_array_equal(uniform.elems, graded.elems)
    np.testing.assert_array_equal(uniform.facets_endo, graded.facets_endo)
    np.testing.assert_array_equal(uniform.facets_epi, graded.facets_epi)
    np.testing.assert_array_equal(uniform.facets_base, graded.facets_base)
    assert uniform.nodes.shape == graded.nodes.shape
    assert not np.array_equal(uniform.nodes, graded.nodes)
    assert graded.tip_refine == 2.5

    uniform_edges = _epi_meridian_edges(uniform)
    graded_edges = _epi_meridian_edges(graded)
    assert len(uniform_edges) == len(graded_edges)
    uniform_tip = np.mean([length for _, length in uniform_edges[:8]])
    graded_tip = np.mean([length for _, length in graded_edges[:8]])
    assert 2.0 < uniform_tip / graded_tip < 2.6

    # The base rim stays fixed to round-off.
    uniform_rim = np.sort(uniform.nodes[np.unique(uniform.facets_base)], axis=0)
    graded_rim = np.sort(graded.nodes[np.unique(graded.facets_base)], axis=0)
    np.testing.assert_allclose(uniform_rim, graded_rim, rtol=0.0, atol=1.0e-15)


def test_tip_refine_graded_mesh_passes_closed_geometry_audit():
    graded = geometry.build_closed_mesh(
        n_t=2, n_core=12, n_radial=6, core_half_width=0.36, tip_refine=2.5,
    )
    # The structural gates stay at production strictness.  The reference-
    # measures comparison uses a declared coarse-mesh bound: remapping the
    # meridian coordinate legitimately moves the faceted volume of this
    # coarse 12x6 test mesh by ~1.3%, while the benchmark-resolution
    # 20x17 mesh passes the production 1% gate (measured <=0.06% up to
    # tip_refine=4).
    record = boundary_audit.audit_geometry(
        graded, require_closed=True, reference_tolerance=0.02,
    )

    assert record["nonpositive_gauss_jacobians"] == 0
    assert record["nonpositive_extended_jacobians"] == 0
    assert record["extended_scaled_jacobian_min"] > 0.2
    assert record["extended_jacobian_condition_max"] < 15.0
    assert record["unclassified_exterior_faces"] == 0
    assert record["labeled_exterior_faces"] == record["exterior_faces"]
    relative = record["reference_relative_differences"]
    assert max(abs(value) for value in relative.values()) < 0.02


def test_tip_refine_rejects_invalid_strength():
    for bad in (0.5, 0.0, np.nan, np.inf, 8.5):
        with pytest.raises(ValueError, match="tip_refine"):
            geometry.build_closed_mesh(
                n_t=1, n_core=4, n_radial=1, tip_refine=bad,
            )


def test_closed_disk_seams_deduplicate_at_fine_resolution():
    # Regression: at some resolutions adjacent blocks computed shared seam
    # columns with last-bit-different floats, the bare rounded dedup key
    # split them, and the resulting overlapping rim quads inflated the base
    # area by 3.62%.  The builder now fails closed on nonconforming seams.
    mesh = geometry.build_closed_mesh(
        n_t=2, n_core=64, n_radial=72, core_half_width=0.36,
    )

    assert len(mesh.facets_base) == mesh.n_t * 4 * mesh.n_core
    base_area = boundary_audit._quad_area(mesh.nodes, mesh.facets_base) * 1.0e4
    reference = boundary_audit._REFERENCE["base_area_cm2"]
    assert abs(base_area / reference - 1.0) < 1.0e-3
