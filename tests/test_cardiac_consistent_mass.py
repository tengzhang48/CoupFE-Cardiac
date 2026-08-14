"""Checks for the application-owned consistent Hex8 inertia operators."""

from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from consistent_mass import (
    ConsistentMassInertia,
    ConsistentNewmarkInertia,
    consistent_mass_coo,
)


def _unit_cube():
    nodes = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [1.0, 1.0, 1.0],
            [0.0, 1.0, 1.0],
        ]
    )
    return nodes, np.arange(8, dtype=int).reshape(1, 8)


def _matrix(density=3.25):
    nodes, elements = _unit_cube()
    rows, cols, values = consistent_mass_coo(nodes, elements, density)
    matrix = sp.coo_matrix((values, (rows, cols)), shape=(24, 24)).tocsr()
    return rows, cols, values, matrix


def _tangent_matrix(tangent, ndof=24):
    return sp.coo_matrix(
        (tangent.values, (tangent.rows, tangent.cols)), shape=(ndof, ndof)
    ).tocsr()


def test_consistent_mass_is_symmetric_and_has_the_expected_row_sum_and_mass():
    density = 3.25
    _, _, _, matrix = _matrix(density)

    np.testing.assert_allclose(matrix.toarray(), matrix.toarray().T, atol=2e-16)
    np.testing.assert_allclose(
        np.asarray(matrix.sum(axis=1)).ravel(),
        np.full(24, density / 8.0),
        rtol=2e-15,
        atol=2e-16,
    )
    # Each displacement component carries rho * volume; the full block has
    # three identical component copies.
    assert float(matrix.sum()) == pytest.approx(3.0 * density, rel=2e-15)


def test_backward_euler_residual_tangent_predictor_and_commit_share_one_mass():
    rows, cols, values, matrix = _matrix()
    u0 = np.linspace(-0.02, 0.03, 24)
    v0 = np.linspace(0.04, -0.01, 24)
    trial = np.linspace(0.01, 0.07, 24)
    dt = 0.125
    damping = 0.3
    inertia = ConsistentMassInertia(
        rows,
        cols,
        values,
        24,
        u0=u0,
        v0=v0,
        damping=damping,
    )

    np.testing.assert_allclose(inertia.predictor(dt), u0 + dt * v0)
    expected = matrix @ ((trial - u0 - dt * v0) / dt**2)
    expected += damping * (matrix @ ((trial - u0) / dt))
    residual = inertia.residual(trial, None, 0.125, dt)
    np.testing.assert_array_equal(residual.gdofs, np.arange(24))
    np.testing.assert_allclose(residual.values, expected, rtol=2e-15, atol=2e-15)

    tangent = _tangent_matrix(inertia.tangent(trial, None, 0.125, dt))
    expected_tangent = matrix / dt**2 + damping * matrix / dt
    np.testing.assert_allclose(
        tangent.toarray(), expected_tangent.toarray(), rtol=2e-15, atol=2e-15
    )
    direction = np.linspace(-0.3, 0.2, 24)
    # The residual is affine in U. A moderate centered-difference step avoids
    # subtractive cancellation while checking the assembled COO tangent.
    step = 1.0e-4
    finite_difference = (
        inertia.residual(trial + step * direction, None, 0.125, dt).values
        - inertia.residual(trial - step * direction, None, 0.125, dt).values
    ) / (2.0 * step)
    np.testing.assert_allclose(
        finite_difference, tangent @ direction, rtol=2e-9, atol=2e-10
    )

    state = {"sentinel": 1}
    assert inertia.commit(trial, state, 0.125, dt) is state
    np.testing.assert_array_equal(inertia.u_prev, trial)
    np.testing.assert_allclose(inertia.v_prev, (trial - u0) / dt)


def test_newmark_residual_tangent_velocity_and_commit_are_consistent():
    rows, cols, values, matrix = _matrix()
    u0 = np.linspace(-0.01, 0.02, 24)
    v0 = np.linspace(0.03, -0.02, 24)
    a0 = np.linspace(-0.04, 0.01, 24)
    trial = np.linspace(0.02, 0.06, 24)
    dt = 0.05
    beta = 0.25
    gamma = 0.5
    inertia = ConsistentNewmarkInertia(
        rows,
        cols,
        values,
        24,
        beta=beta,
        gamma=gamma,
        u0=u0,
        v0=v0,
        a0=a0,
    )

    u_pred = u0 + dt * v0 + (0.5 - beta) * dt**2 * a0
    acceleration = (trial - u_pred) / (beta * dt**2)
    velocity = v0 + dt * ((1.0 - gamma) * a0 + gamma * acceleration)
    np.testing.assert_allclose(
        inertia.predictor(dt), u0 + dt * v0 + 0.5 * dt**2 * a0
    )
    np.testing.assert_allclose(inertia.velocity(trial, dt), velocity)
    assert inertia.velocity_tangent(dt) == pytest.approx(gamma / (beta * dt))
    np.testing.assert_allclose(
        inertia.residual(trial, None, dt, dt).values,
        matrix @ acceleration,
        rtol=2e-15,
        atol=2e-15,
    )
    tangent = _tangent_matrix(inertia.tangent(trial, None, dt, dt))
    np.testing.assert_allclose(
        tangent.toarray(), (matrix / (beta * dt**2)).toarray(), rtol=2e-15
    )

    state = np.arange(3)
    assert inertia.commit(trial, state, dt, dt) is state
    np.testing.assert_array_equal(inertia.u_prev, trial)
    np.testing.assert_allclose(inertia.v_prev, velocity)
    np.testing.assert_allclose(inertia.a_prev, acceleration)


def test_consistent_mass_rejects_invalid_density_and_element_orientation():
    nodes, elements = _unit_cube()
    with pytest.raises(ValueError, match="density"):
        consistent_mass_coo(nodes, elements, 0.0)
    inverted = elements.copy()
    inverted[0, [0, 1]] = inverted[0, [1, 0]]
    with pytest.raises(ValueError, match="nonpositive Jacobian"):
        consistent_mass_coo(nodes, inverted, 1.0)
