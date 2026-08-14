from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

import geometry
from local_pressure import (
    InvalidDeformationError,
    LocalPressureHex8Operator,
    PAPER_MEAN_DILATATION_LAW,
    _NATURAL_NODES,
)


def _unit_cube():
    nodes = 0.5 * (_NATURAL_NODES + 1.0)
    return nodes, np.arange(8, dtype=int)[None, :]


def _affine_displacement(nodes, deformation_gradient):
    return (nodes @ (deformation_gradient - np.eye(3)).T).ravel()


def _residual(operator, displacement):
    contribution = operator.residual(displacement, None, 0.0, 1.0)
    result = np.zeros(operator.ndof)
    np.add.at(result, contribution.gdofs, contribution.values)
    return result


def _tangent(operator, displacement):
    contribution = operator.tangent(displacement, None, 0.0, 1.0)
    return sp.coo_matrix(
        (contribution.values, (contribution.rows, contribution.cols)),
        shape=(operator.ndof, operator.ndof),
    ).toarray()


def test_uniform_affine_state_recovers_k_log_j_and_isochoric_zero():
    nodes, elements = _unit_cube()
    bulk = 1.0e6
    operator = LocalPressureHex8Operator(
        nodes, elements, 24, bulk_modulus=bulk
    )
    deformation = np.array(
        [[1.08, 0.03, 0.0], [0.0, 0.97, 0.02], [0.0, 0.0, 1.04]]
    )
    displacement = _affine_displacement(nodes, deformation)
    np.testing.assert_allclose(
        operator.element_pressure(displacement),
        [bulk * np.log(np.linalg.det(deformation))],
        rtol=2.0e-14,
        atol=1.0e-10,
    )
    np.testing.assert_allclose(
        operator.deformation_jacobians(displacement),
        np.linalg.det(deformation),
        rtol=2.0e-14,
        atol=2.0e-14,
    )

    isochoric = np.diag([1.1, 1.0 / 1.1, 1.0])
    isochoric_u = _affine_displacement(nodes, isochoric)
    np.testing.assert_allclose(operator.element_pressure(isochoric_u), 0.0, atol=2e-10)
    np.testing.assert_allclose(_residual(operator, isochoric_u), 0.0, atol=2e-10)


def test_uniform_affine_state_recovers_paper_mean_dilatation_pressure():
    nodes, elements = _unit_cube()
    bulk = 1.0e6
    operator = LocalPressureHex8Operator(
        nodes,
        elements,
        24,
        bulk_modulus=bulk,
        pressure_law=PAPER_MEAN_DILATATION_LAW,
    )
    deformation = np.array(
        [[1.08, 0.03, 0.0], [0.0, 0.97, 0.02], [0.0, 0.0, 1.04]]
    )
    determinant = np.linalg.det(deformation)
    displacement = _affine_displacement(nodes, deformation)
    np.testing.assert_allclose(
        operator.element_pressure(displacement),
        [0.5 * bulk * (determinant**2 - 1.0)],
        rtol=2.0e-14,
        atol=1.0e-10,
    )
    compression = np.diag([0.95, 0.96, 0.97])
    compression_u = _affine_displacement(nodes, compression)
    assert operator.element_pressure(compression_u)[0] < 0.0


def test_paper_law_is_applied_after_mean_log_j_not_pointwise():
    nodes, elements = _unit_cube()
    bulk = 1.0e6
    operator = LocalPressureHex8Operator(
        nodes,
        elements,
        24,
        bulk_modulus=bulk,
        pressure_law=PAPER_MEAN_DILATATION_LAW,
    )
    rng = np.random.default_rng(20260804)
    displacement = rng.normal(scale=2.5e-2, size=24)
    determinant = operator.deformation_jacobians(displacement)[0]
    mean_log_j = np.sum(operator._weight[0] * np.log(determinant)) / (
        operator._volume[0]
    )
    expected = 0.5 * bulk * np.expm1(2.0 * mean_log_j)
    pointwise_average = np.sum(
        operator._weight[0] * 0.5 * bulk * (determinant**2 - 1.0)
    ) / operator._volume[0]
    np.testing.assert_allclose(
        operator.element_pressure(displacement), expected, rtol=2.0e-14
    )
    assert abs(expected - pointwise_average) > 100.0


def test_unknown_local_pressure_law_is_rejected():
    nodes, elements = _unit_cube()
    with pytest.raises(ValueError, match="pressure_law must be one of"):
        LocalPressureHex8Operator(
            nodes,
            elements,
            24,
            bulk_modulus=1.0e6,
            pressure_law="unknown",
        )


def test_default_local_pressure_law_is_bit_identical_to_explicit_log():
    nodes, elements = _unit_cube()
    default = LocalPressureHex8Operator(
        nodes, elements, 24, bulk_modulus=1.0e6
    )
    explicit = LocalPressureHex8Operator(
        nodes, elements, 24, bulk_modulus=1.0e6, pressure_law="log"
    )
    displacement = np.random.default_rng(913).normal(scale=7.0e-3, size=24)
    np.testing.assert_array_equal(
        default.element_pressure(displacement),
        explicit.element_pressure(displacement),
    )
    np.testing.assert_array_equal(
        _residual(default, displacement), _residual(explicit, displacement)
    )
    np.testing.assert_array_equal(
        _tangent(default, displacement), _tangent(explicit, displacement)
    )


def test_paper_law_rejects_nonfinite_tangent_slope_even_if_pressure_is_finite():
    nodes, elements = _unit_cube()
    bulk = 1.0e6
    operator = LocalPressureHex8Operator(
        nodes,
        elements,
        24,
        bulk_modulus=bulk,
        pressure_law=PAPER_MEAN_DILATATION_LAW,
    )
    twice_mean_log_j = np.log(np.finfo(float).max / bulk) + np.log(1.5)
    stretch = np.exp(twice_mean_log_j / 6.0)
    displacement = _affine_displacement(nodes, stretch * np.eye(3))
    with pytest.raises(InvalidDeformationError, match="pressure"):
        operator.element_pressure(displacement)


@pytest.mark.parametrize("pressure_law", ["log", "paper"])
def test_condensed_tangent_matches_finite_difference_and_is_symmetric(
    pressure_law,
):
    nodes, elements = _unit_cube()
    operator = LocalPressureHex8Operator(
        nodes,
        elements,
        24,
        bulk_modulus=1.0e4,
        pressure_law=pressure_law,
    )
    rng = np.random.default_rng(18)
    displacement = rng.normal(scale=8.0e-3, size=24)
    direction = rng.normal(size=24)
    direction /= np.linalg.norm(direction)
    tangent = _tangent(operator, displacement)
    step = 2.0e-7
    finite_difference = (
        _residual(operator, displacement + step * direction)
        - _residual(operator, displacement - step * direction)
    ) / (2.0 * step)
    np.testing.assert_allclose(
        tangent @ direction, finite_difference, rtol=4.0e-8, atol=3.0e-7
    )
    np.testing.assert_allclose(tangent, tangent.T, rtol=2.0e-13, atol=2.0e-10)


@pytest.mark.parametrize("pressure_law", ["log", "paper"])
def test_condensed_tangent_on_curved_cardiac_hex8_matches_finite_difference(
    pressure_law,
):
    mesh = geometry.build_mesh(
        n_t=1, n_mu=2, n_theta=4, flip_helix=True, apex_offset=0.2
    )
    nodes = mesh.nodes[mesh.elems[-1]]

    # Every face cross-difference vanishes for an affine Hex8 map.  Protect the
    # purpose of this regression by requiring a materially non-affine cell from
    # the actual curved cardiac mesh rather than merely trusting its origin.
    face_cross_differences = np.array(
        [
            nodes[0] - nodes[1] + nodes[2] - nodes[3],
            nodes[4] - nodes[5] + nodes[6] - nodes[7],
            nodes[0] - nodes[1] - nodes[4] + nodes[5],
            nodes[3] - nodes[2] - nodes[7] + nodes[6],
            nodes[0] - nodes[3] - nodes[4] + nodes[7],
            nodes[1] - nodes[2] - nodes[5] + nodes[6],
        ]
    )
    diameter = max(
        np.linalg.norm(nodes[a] - nodes[b])
        for a in range(8)
        for b in range(a)
    )
    assert np.max(np.linalg.norm(face_cross_differences, axis=1)) / diameter > 1.0e-2

    operator = LocalPressureHex8Operator(
        nodes,
        np.arange(8, dtype=int)[None, :],
        24,
        bulk_modulus=1.0e6,
        pressure_law=pressure_law,
    )
    rng = np.random.default_rng(314159)
    displacement = rng.normal(scale=1.0e-5, size=24)
    direction = rng.normal(size=24)
    direction /= np.linalg.norm(direction)
    tangent = _tangent(operator, displacement)
    step = 3.0e-8
    finite_difference = (
        _residual(operator, displacement + step * direction)
        - _residual(operator, displacement - step * direction)
    ) / (2.0 * step)

    tangent_error = np.linalg.norm(tangent @ direction - finite_difference) / max(
        np.linalg.norm(finite_difference), 1.0
    )
    symmetry_error = np.linalg.norm(tangent - tangent.T) / max(
        np.linalg.norm(tangent), 1.0
    )
    # The central difference is limited by truncation/roundoff; these bounds are
    # deliberately looser than the observed O(1e-10) and O(1e-16) errors while
    # remaining tight enough to catch a missing condensed or geometric term.
    assert tangent_error < 2.0e-9
    assert symmetry_error < 5.0e-13


def test_local_pressure_fails_closed_for_bad_reference_or_deformation():
    nodes, elements = _unit_cube()
    degenerate = nodes.copy()
    degenerate[4:] = degenerate[:4]
    with pytest.raises(ValueError, match="non-positive Gauss-point"):
        LocalPressureHex8Operator(degenerate, elements, 24, bulk_modulus=1.0e6)

    operator = LocalPressureHex8Operator(
        nodes, elements, 24, bulk_modulus=1.0e6
    )
    inverted = _affine_displacement(nodes, np.diag([-1.0, 1.0, 1.0]))
    with pytest.raises(InvalidDeformationError, match="inverted or non-finite"):
        operator.residual(inverted, None, 0.0, 1.0)

    nonfinite = np.zeros(24)
    nonfinite[0] = np.nan
    with pytest.raises(RuntimeError, match="non-finite"):
        operator.tangent(nonfinite, None, 0.0, 1.0)

    overflow_reference = nodes * 1.0e200
    with pytest.raises(ValueError, match="non-finite Jacobian determinant"):
        LocalPressureHex8Operator(
            overflow_reference, elements, 24, bulk_modulus=1.0e6
        )

    overflow_u = np.full(24, 1.0e200)
    overflow_u[0::3] *= np.asarray(_NATURAL_NODES[:, 0])
    overflow_u[1::3] *= np.asarray(_NATURAL_NODES[:, 1])
    overflow_u[2::3] *= np.asarray(_NATURAL_NODES[:, 2])
    with pytest.raises(RuntimeError, match="inverted or non-finite"):
        operator.residual(overflow_u, None, 0.0, 1.0)


def test_local_pressure_max_step_keeps_core_newton_trial_in_domain():
    nodes, elements = _unit_cube()
    operator = LocalPressureHex8Operator(
        nodes, elements, 24, bulk_modulus=1.0e6
    )
    displacement = np.zeros(24)
    # The full correction reflects the cube and alpha=0.5 collapses it. The
    # operator must shorten the correction to the first valid halving, while a
    # direct evaluation of the invalid full trial remains fail-closed.
    increment = _affine_displacement(nodes, np.diag([-1.0, 1.0, 1.0]))
    with pytest.raises(InvalidDeformationError):
        operator.residual(displacement + increment, None, 0.0, 1.0)

    alpha = operator.max_step(displacement, increment, 0.0)
    assert alpha == pytest.approx(0.25)
    assert np.min(
        operator.deformation_jacobians(displacement + alpha * increment)
    ) > 0.0


def test_distorted_affine_patch_and_shared_node_scatter():
    nodes, elements = _unit_cube()
    affine_map = np.array(
        [[1.2, 0.13, -0.04], [0.08, 0.91, 0.06], [-0.03, 0.11, 1.07]]
    )
    distorted = nodes @ affine_map.T + np.array([0.3, -0.2, 0.5])
    deformation = np.array(
        [[1.04, 0.02, -0.01], [0.01, 0.98, 0.03], [0.0, -0.02, 1.03]]
    )
    operator = LocalPressureHex8Operator(
        distorted, elements, 24, bulk_modulus=8.0e5
    )
    displacement = _affine_displacement(distorted, deformation)
    np.testing.assert_allclose(
        operator.deformation_jacobians(displacement),
        np.linalg.det(deformation),
        rtol=3.0e-14,
        atol=3.0e-14,
    )

    # Two elements share a face. Assembly must sum both element contributions
    # at the shared global nodes rather than overwrite repeated COO entries.
    right = distorted.copy()
    shift = distorted[1] - distorted[0]
    right += shift
    combined_nodes = np.vstack((distorted, right[[1, 2, 5, 6]]))
    second = np.array([1, 8, 9, 2, 5, 10, 11, 6], dtype=int)
    combined_elements = np.vstack((np.arange(8), second))
    combined = LocalPressureHex8Operator(
        combined_nodes,
        combined_elements,
        3 * len(combined_nodes),
        bulk_modulus=8.0e5,
    )
    combined_u = _affine_displacement(combined_nodes, deformation)
    residual = _residual(combined, combined_u).reshape(-1, 3)
    first_only = LocalPressureHex8Operator(
        combined_nodes, combined_elements[:1], 3 * len(combined_nodes),
        bulk_modulus=8.0e5,
    )
    second_only = LocalPressureHex8Operator(
        combined_nodes, combined_elements[1:], 3 * len(combined_nodes),
        bulk_modulus=8.0e5,
    )
    expected = (
        _residual(first_only, combined_u) + _residual(second_only, combined_u)
    ).reshape(-1, 3)
    np.testing.assert_allclose(residual, expected, rtol=2.0e-14, atol=2.0e-10)
    assert np.linalg.norm(residual[[1, 2, 5, 6]]) > 0.0
