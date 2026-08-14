from __future__ import annotations

import numpy as np
import pytest

import geometry
import sampling


def _unit_hex(x0=0.0, x1=1.0):
    natural = sampling.HEX8_NATURAL_COORDINATES
    nodes = 0.5 * (natural + 1.0)
    nodes[:, 0] = x0 + (x1 - x0) * nodes[:, 0]
    return nodes


def test_affine_hex_location_and_vector_interpolation():
    nodes = _unit_hex()
    elements = np.arange(8, dtype=int)[None, :]
    point = np.array([0.2, 0.3, 0.4])

    location = sampling.locate_hex8_point(nodes, elements, point)

    assert location.element_index == 0
    np.testing.assert_allclose(
        location.natural_coordinates, [-0.6, -0.4, -0.2], atol=2.0e-14
    )
    np.testing.assert_allclose(
        np.asarray(location.weights) @ nodes, point, rtol=0.0, atol=2.0e-14
    )
    # Every component is affine and must therefore be reproduced exactly by Hex8.
    nodal_u = np.column_stack(
        (nodes[:, 0] + 2.0 * nodes[:, 1], -3.0 * nodes[:, 2], nodes.sum(axis=1))
    )
    expected = np.array([0.8, -1.2, 0.9])
    np.testing.assert_allclose(
        sampling.interpolate_displacement(nodal_u, location), expected, atol=2.0e-14
    )
    flat_with_extra_dof = np.column_stack((nodal_u, np.full(8, 99.0))).ravel()
    np.testing.assert_allclose(
        sampling.interpolate_displacement(
            flat_with_extra_dof, location, dof_per_node=4
        ),
        expected,
        atol=2.0e-14,
    )


def test_warped_hex_uses_inverse_isoparametric_map():
    nodes = _unit_hex()
    nodes[6] += np.array([0.15, -0.05, 0.2])
    nodes[7] += np.array([-0.04, 0.03, 0.08])
    target_natural = np.array([0.23, -0.31, 0.44])
    point = sampling.hex8_shape(target_natural) @ nodes

    location = sampling.locate_hex8_point(nodes, np.arange(8)[None, :], point)

    np.testing.assert_allclose(
        location.natural_coordinates, target_natural, rtol=0.0, atol=2.0e-11
    )
    assert location.iterations >= 1
    assert location.reconstruction_error < 1.0e-11


def test_shared_face_choice_is_lowest_element_index():
    nodes = np.vstack((_unit_hex(0.0, 1.0), _unit_hex(1.0, 2.0)))
    elements = np.vstack((np.arange(8), np.arange(8, 16)))
    point = np.array([1.0, 0.37, 0.61])

    assert sampling.candidate_hex8_elements(nodes, elements, point) == (0, 1)
    location = sampling.locate_hex8_point(nodes, elements, point)

    assert location.element_index == 0
    np.testing.assert_allclose(
        np.asarray(location.weights) @ nodes[np.asarray(location.node_ids)],
        point,
        atol=2.0e-14,
    )


def test_cardiac_probe_points_are_located_in_hex8_elements():
    mesh = geometry.build_mesh(
        n_t=2, n_mu=12, n_theta=16, flip_helix=True, apex_offset=0.2
    )

    for point in (geometry.P0, geometry.P1):
        location = sampling.locate_hex8_point(mesh.nodes, mesh.elems, point)
        assert 0 <= location.element_index < mesh.n_elem
        natural = np.asarray(location.natural_coordinates)
        assert np.all(natural >= -1.0 - 1.0e-9)
        assert np.all(natural <= +1.0 + 1.0e-9)
        np.testing.assert_allclose(
            np.asarray(location.weights) @ mesh.nodes[np.asarray(location.node_ids)],
            point,
            rtol=0.0,
            atol=2.0e-11,
        )


def test_outside_point_and_degenerate_cell_fail_closed():
    nodes = _unit_hex()
    elements = np.arange(8)[None, :]
    with pytest.raises(RuntimeError, match="outside"):
        sampling.locate_hex8_point(nodes, elements, [2.0, 2.0, 2.0])

    collapsed = nodes.copy()
    collapsed[:, 2] = 0.0
    with pytest.raises(RuntimeError, match="degenerate"):
        sampling.locate_hex8_point(collapsed, elements, [0.5, 0.5, 0.0])


def test_interpolator_rejects_nonfinite_selected_displacement():
    nodes = _unit_hex()
    location = sampling.locate_hex8_point(
        nodes, np.arange(8)[None, :], [0.5, 0.5, 0.5]
    )
    displacement = np.zeros((8, 3))
    displacement[0, 0] = np.nan
    with pytest.raises(RuntimeError, match="non-finite"):
        sampling.interpolate_displacement(displacement, location)
