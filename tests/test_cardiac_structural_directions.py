"""Regression gates for physical-coordinate structural directions."""
from __future__ import annotations

import numpy as np

import fiber_crosscheck
import geometry
import material
import run as cardiac_run
import structural_directions


def _angle_degrees(left, right):
    return float(
        np.rad2deg(
            np.arccos(np.clip(np.dot(left, right), -1.0, 1.0))
        )
    )


def test_toolkit_apex_rule_is_applied_during_coordinate_reconstruction():
    tbar = 0.4
    positive_axis = np.array([geometry.r_long(tbar), 0.0, 0.0])

    u, v = structural_directions.ellipsoid_parameters(tbar, positive_axis)

    assert u == 0.0
    assert v == 0.0


def test_topologies_have_distinct_reconstruction_provenance():
    assert geometry.structural_direction_reconstruction(
        "closed_multiblock_disk"
    ) == "toolkit-physical-coordinate-u-v-v1"
    assert geometry.structural_direction_reconstruction(
        "polar_ring"
    ) == "historical-parametric-mu-theta-v1"


def test_physical_frame_matches_pinned_oracle_at_straight_wall_point():
    tbar, rho, theta = 0.37, 0.64, 0.81
    physical_point = geometry.closed_wall_point(tbar, rho, theta)
    stored_mu = -np.pi + rho * (np.pi + geometry.mu_base(tbar))
    interpolated_ellipsoid_point = geometry.point(tbar, stored_mu, theta)
    assert np.linalg.norm(physical_point - interpolated_ellipsoid_point) > 1.0e-8

    expected = fiber_crosscheck.authors_frame(physical_point, tbar)
    observed = geometry.toolkit_fiber_frame(
        tbar, physical_point, asign=-1.0
    )

    for actual, reference in zip(observed, expected):
        np.testing.assert_allclose(actual, reference, rtol=0.0, atol=3.0e-15)


def test_closed_gp_direct_reconstructs_at_q1_physical_point():
    mesh = geometry.build_closed_mesh(
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )
    shapes = cardiac_run._hex8_gp_shape()
    tbar = np.asarray(mesh.param[:, 0], dtype=float)
    largest_old_angle = -1.0
    selected = None
    for connectivity in mesh.elems:
        parametric = mesh.param[connectivity]
        for weights in shapes:
            observed = geometry.sample_structural_frame(
                mesh, connectivity, weights, tbar, asign=-1.0
            )
            t_gp = float(weights @ tbar[connectivity])
            physical_point = weights @ mesh.nodes[connectivity]
            expected = fiber_crosscheck.authors_frame(physical_point, t_gp)
            for actual, reference in zip(observed, expected):
                np.testing.assert_allclose(
                    actual, reference, rtol=0.0, atol=4.0e-15
                )

            mu_gp = float(weights @ parametric[:, 1])
            theta_gp = geometry._circular_shape_average(
                parametric[:, 2], weights
            )
            historical = geometry.fiber_frame(
                t_gp, mu_gp, theta_gp, asign=-1.0
            )[0]
            angle = _angle_degrees(historical, observed[0])
            if angle > largest_old_angle:
                largest_old_angle = angle
                selected = physical_point

    assert selected is not None
    assert largest_old_angle > 4.0


def test_injected_tbar_rebuilds_closed_nodal_and_element_frames(tmp_path):
    mesh = geometry.build_closed_mesh(
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )
    tbar = np.asarray(mesh.param[:, 0], dtype=float).copy()
    interior = np.flatnonzero(tbar == 0.5)
    transverse = mesh.nodes[interior, 1]
    tbar[interior] += 0.08 * transverse / np.max(np.abs(transverse))
    path = tmp_path / "tbar.npy"
    np.save(path, tbar)

    returned = cardiac_run._apply_laplace_tbar(mesh, path, asign=-1.0)

    np.testing.assert_array_equal(returned, tbar)
    node = int(interior[np.argmax(np.abs(transverse))])
    expected_node = fiber_crosscheck.authors_frame(mesh.nodes[node], tbar[node])
    np.testing.assert_allclose(
        mesh.fiber_node[node], expected_node[0], rtol=0.0, atol=3.0e-15
    )
    np.testing.assert_allclose(
        mesh.sheet_node[node], expected_node[1], rtol=0.0, atol=3.0e-15
    )

    connectivity = mesh.elems[0]
    weights = np.full(8, 1.0 / 8.0)
    t_center = float(weights @ tbar[connectivity])
    point_center = weights @ mesh.nodes[connectivity]
    expected_element = fiber_crosscheck.authors_frame(point_center, t_center)
    np.testing.assert_allclose(
        mesh.fiber[0], expected_element[0], rtol=0.0, atol=3.0e-15
    )
    np.testing.assert_allclose(
        mesh.sheet[0], expected_element[1], rtol=0.0, atol=3.0e-15
    )
    np.testing.assert_allclose(
        mesh.normal[0], expected_element[2], rtol=0.0, atol=3.0e-15
    )


def test_historical_polar_frames_retain_the_parametric_rule():
    mesh = geometry.build_mesh(
        n_t=1,
        n_mu=4,
        n_theta=8,
        flip_helix=True,
        apex_offset=0.2,
    )
    before = tuple(
        field.copy()
        for field in (
            mesh.fiber,
            mesh.sheet,
            mesh.normal,
            mesh.fiber_node,
            mesh.sheet_node,
        )
    )

    geometry.update_structural_frames(
        mesh, np.asarray(mesh.param[:, 0]), asign=-1.0
    )

    for actual, expected in zip(
        (
            mesh.fiber,
            mesh.sheet,
            mesh.normal,
            mesh.fiber_node,
            mesh.sheet_node,
        ),
        before,
    ):
        np.testing.assert_array_equal(actual, expected)


def test_sheet_sign_reversal_is_invariant_in_current_material_law():
    model = material.HolzapfelOgdenActive()
    state = material.verification_state()

    positive, _ = model.stress_PK1(
        state["F"],
        state["ff"],
        state["ss"],
        state["fssym"],
        state["E_prev"],
        state["dt"],
    )
    negative, _ = model.stress_PK1(
        state["F"],
        state["ff"],
        state["ss"],
        -state["fssym"],
        state["E_prev"],
        state["dt"],
    )

    np.testing.assert_allclose(positive, negative, rtol=0.0, atol=2.0e-13)
