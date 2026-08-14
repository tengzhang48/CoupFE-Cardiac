"""Element-local condensed-pressure contribution for cardiac Hex8 meshes.

This application operator supplies only the volumetric part of the stress. It
is composed with the standard (non-F-bar) compiled cardiac element, whose
``kappa`` property is set to zero and which continues to own all Gauss-point
fiber, viscous-history, and active-stress state.

For each element, one pressure is eliminated algebraically from the selected
mean-dilatation energy.  The historical log law is

``p_e = K / V_e * integral(log(det(F))) dV``

and the paper-law variant is

``p_e = K/2 * (exp(2 m_e) - 1)``,
``m_e = 1/V_e * integral(log(det(F))) dV``.

In both cases ``P_vol = p_e F**-T``. The returned displacement tangent is the
exact derivative, including the dependence of ``p_e`` on every element
displacement degree of freedom. No pressure state is committed because the
local response is algebraic and is evaluated exactly on every call.

The formulation is intentionally application-owned until it has broader use
and validation. It does not change CoupFE Core's generic local-pressure code
generator or pretend to reproduce the older F-bar Case B campaign.
"""
from __future__ import annotations

import numpy as np

from coupfe.operators.base import Residual, Tangent


_NATURAL_NODES = np.array(
    [
        [-1.0, -1.0, -1.0],
        [1.0, -1.0, -1.0],
        [1.0, 1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
        [1.0, -1.0, 1.0],
        [1.0, 1.0, 1.0],
        [-1.0, 1.0, 1.0],
    ],
    dtype=float,
)


def _gauss_shape_gradients() -> np.ndarray:
    """Return ``dN/dxi`` at the eight 2x2x2 Gauss points."""
    coordinate = 1.0 / np.sqrt(3.0)
    points = np.array(
        [
            [sx * coordinate, sy * coordinate, sz * coordinate]
            for sz in (1.0, -1.0)
            for sy in (1.0, -1.0)
            for sx in (1.0, -1.0)
        ],
        dtype=float,
    )
    gradients = np.empty((8, 8, 3), dtype=float)
    for gp, xi in enumerate(points):
        for node, sign in enumerate(_NATURAL_NODES):
            for direction in range(3):
                other = [axis for axis in range(3) if axis != direction]
                gradients[gp, node, direction] = (
                    0.125
                    * sign[direction]
                    * (1.0 + sign[other[0]] * xi[other[0]])
                    * (1.0 + sign[other[1]] * xi[other[1]])
                )
    return gradients


_DN_DXI = _gauss_shape_gradients()

LOG_MEAN_DILATATION_LAW = "log"
PAPER_MEAN_DILATATION_LAW = "paper"
LOCAL_PRESSURE_LAWS = frozenset(
    {LOG_MEAN_DILATATION_LAW, PAPER_MEAN_DILATATION_LAW}
)


def _pressure_response(mean_log_j, bulk_modulus, pressure_law):
    """Return ``(p, dp/dm)`` for ``m = <log(J)>``.

    ``pressure_law="paper"`` applies the benchmark paper's scalar volumetric
    energy to the element mean dilatation.  ``expm1`` preserves accuracy close
    to the reference state, where both laws have tangent ``K``.
    """
    mean_log_j = np.asarray(mean_log_j, dtype=float)
    bulk_modulus = float(bulk_modulus)
    pressure_law = str(pressure_law)
    if pressure_law == LOG_MEAN_DILATATION_LAW:
        pressure = bulk_modulus * mean_log_j
        derivative = np.full_like(mean_log_j, bulk_modulus)
    elif pressure_law == PAPER_MEAN_DILATATION_LAW:
        with np.errstate(over="ignore", invalid="ignore"):
            exponential = np.exp(2.0 * mean_log_j)
            pressure = 0.5 * bulk_modulus * np.expm1(2.0 * mean_log_j)
            derivative = bulk_modulus * exponential
    else:
        raise ValueError(
            f"pressure_law must be one of {sorted(LOCAL_PRESSURE_LAWS)}"
        )
    return pressure, derivative


class InvalidDeformationError(RuntimeError):
    """A recoverable trial displacement lies outside the finite-strain domain.

    Direct callers still receive a fail-closed exception. Nonlinear solvers may
    catch this exact type to shorten or reject a trial step; other numerical or
    programming errors must continue to propagate.
    """


def _first_bad_element(values: np.ndarray) -> int:
    """Return the first element containing a nonfinite or nonpositive value."""
    bad = np.any(~np.isfinite(values) | (values <= 0.0), axis=1)
    indices = np.flatnonzero(bad)
    if len(indices) == 0:
        raise RuntimeError("internal local-pressure validity check is inconsistent")
    return int(indices[0])


class LocalPressureHex8Operator:
    """Vectorized one-pressure-per-element volumetric Hex8 operator."""

    def __init__(
        self,
        nodes,
        elements,
        ndof,
        *,
        bulk_modulus,
        dof_per_node=3,
        pressure_law=LOG_MEAN_DILATATION_LAW,
    ):
        self.nodes = np.asarray(nodes, dtype=float)
        self.elements = np.asarray(elements, dtype=int)
        self.ndof = int(ndof)
        self.bulk_modulus = float(bulk_modulus)
        self.dof_per_node = int(dof_per_node)
        self.pressure_law = str(pressure_law)
        if self.nodes.ndim != 2 or self.nodes.shape[1] != 3:
            raise ValueError("nodes must have shape (n_node, 3)")
        if self.elements.ndim != 2 or self.elements.shape[1] != 8:
            raise ValueError("elements must have shape (n_element, 8)")
        if not np.all(np.isfinite(self.nodes)):
            raise ValueError("nodes must be finite")
        if len(self.elements) < 1:
            raise ValueError("at least one Hex8 element is required")
        if np.min(self.elements) < 0 or np.max(self.elements) >= len(self.nodes):
            raise ValueError("element connectivity is outside the node array")
        if self.dof_per_node != 3:
            raise ValueError("the cardiac local-pressure operator requires 3 DOF/node")
        if self.ndof != self.dof_per_node * len(self.nodes):
            raise ValueError("ndof is inconsistent with nodes and dof_per_node")
        if not np.isfinite(self.bulk_modulus) or self.bulk_modulus <= 0.0:
            raise ValueError("bulk_modulus must be finite and positive")
        if self.pressure_law not in LOCAL_PRESSURE_LAWS:
            raise ValueError(
                f"pressure_law must be one of {sorted(LOCAL_PRESSURE_LAWS)}"
            )

        coordinates = self.nodes[self.elements]
        # J[i,j] = dx_i/dxi_j. Rows of the physical gradient are
        # dN/dx = dN/dxi @ inv(J).
        jacobian = np.einsum("eai,gaj->egij", coordinates, _DN_DXI)
        with np.errstate(over="ignore", invalid="ignore"):
            determinant = np.linalg.det(jacobian)
        if not np.all(np.isfinite(determinant)) or np.any(determinant <= 0.0):
            bad = _first_bad_element(determinant)
            raise ValueError(
                "local-pressure reference Hex8 has a non-positive Gauss-point "
                f"or non-finite Jacobian determinant in element {bad}"
            )
        try:
            inverse = np.linalg.inv(jacobian)
        except np.linalg.LinAlgError as error:
            raise ValueError("local-pressure reference Hex8 Jacobian is singular") from error
        self._gradient = np.einsum("gaj,egjk->egak", _DN_DXI, inverse)
        self._weight = determinant
        self._volume = determinant.sum(axis=1)
        if (
            not np.all(np.isfinite(inverse))
            or not np.all(np.isfinite(self._gradient))
            or not np.all(np.isfinite(self._volume))
            or np.any(self._volume <= 0.0)
        ):
            raise ValueError(
                "local-pressure reference Hex8 produced non-finite gradients or volume"
            )
        self._element_dofs = (
            self.elements[:, :, None] * self.dof_per_node
            + np.arange(self.dof_per_node)[None, None, :]
        ).reshape(len(self.elements), 24)

    def _kinematics(self, U):
        displacement = np.asarray(U, dtype=float)
        if displacement.shape != (self.ndof,):
            raise ValueError(
                f"displacement has shape {displacement.shape}; expected {(self.ndof,)}"
            )
        if not np.all(np.isfinite(displacement)):
            raise RuntimeError("local-pressure displacement is non-finite")
        element_u = displacement.reshape(-1, 3)[self.elements]
        deformation = np.eye(3)[None, None, :, :] + np.einsum(
            "eai,egaj->egij", element_u, self._gradient
        )
        with np.errstate(over="ignore", invalid="ignore"):
            determinant = np.linalg.det(deformation)
        if not np.all(np.isfinite(determinant)) or np.any(determinant <= 0.0):
            bad = _first_bad_element(determinant)
            finite = determinant[bad][np.isfinite(determinant[bad])]
            minimum = float(np.min(finite)) if len(finite) else float("nan")
            raise InvalidDeformationError(
                "local-pressure Hex8 deformation is inverted or non-finite: "
                f"element={bad}, min(det(F))={minimum:.6e}"
            )
        try:
            inverse_transpose = np.swapaxes(np.linalg.inv(deformation), -1, -2)
        except np.linalg.LinAlgError as error:
            raise InvalidDeformationError(
                "local-pressure Hex8 deformation is singular"
            ) from error
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            mean_log_j = (
                (self._weight * np.log(determinant)).sum(axis=1) / self._volume
            )
            pressure, pressure_derivative = _pressure_response(
                mean_log_j, self.bulk_modulus, self.pressure_law
            )
        if not np.all(np.isfinite(inverse_transpose)) or not np.all(
            np.isfinite(pressure)
        ) or not np.all(np.isfinite(pressure_derivative)):
            raise InvalidDeformationError(
                "local element pressure or inverse deformation is non-finite"
            )
        return inverse_transpose, determinant, pressure, pressure_derivative

    def element_pressure(self, U) -> np.ndarray:
        """Return the algebraically eliminated pressure for every element."""
        return self._kinematics(U)[2].copy()

    def deformation_jacobians(self, U) -> np.ndarray:
        """Return ``det(F)`` at every element/Gauss-point pair."""
        return self._kinematics(U)[1].copy()

    def max_step(self, U, dU, t, dt=None) -> float:
        """Bound a Core Newton correction to an orientation-preserving trial.

        The Cardiac PETSc adapter maps this exact exception to a nonfinite
        residual so PETSc's ``newtonls``/``bt`` line search rejects the trial.
        Core's serial Newton instead consumes this optional operator hook before
        its residual backtracking loop. Halving is intentionally conservative:
        the base iterate is already required to be valid, so continuity
        guarantees a valid positive step for a finite correction.
        """
        displacement = np.asarray(U, dtype=float)
        increment = np.asarray(dU, dtype=float)
        if displacement.shape != (self.ndof,) or increment.shape != (self.ndof,):
            raise ValueError(
                "local-pressure Newton iterate and increment must match ndof"
            )
        if not np.all(np.isfinite(increment)):
            raise RuntimeError("local-pressure Newton increment is non-finite")
        # Establish that the accepted/base iterate itself is in the domain.
        self._kinematics(displacement)
        alpha = 1.0
        for _ in range(60):
            try:
                self._kinematics(displacement + alpha * increment)
            except InvalidDeformationError:
                alpha *= 0.5
                continue
            return alpha
        raise InvalidDeformationError(
            "local-pressure Newton correction has no resolvable valid trial step"
        )

    def residual(self, U, state, t, dt):
        inverse_transpose, _determinant, pressure, _pressure_derivative = (
            self._kinematics(U)
        )
        stress = pressure[:, None, None, None] * inverse_transpose
        element_residual = np.einsum(
            "egij,egaj,eg->eai", stress, self._gradient, self._weight
        )
        if not np.all(np.isfinite(element_residual)):
            raise RuntimeError("local-pressure residual is non-finite")
        return Residual(self._element_dofs.ravel(), element_residual.ravel())

    def tangent(self, U, state, t, dt):
        inverse_transpose, _determinant, pressure, pressure_derivative = (
            self._kinematics(U)
        )

        # q_ai = integral(F^-T_ij N_a,j) dV and
        # dp/du_bk = dp/dm * q_bk/V for m=<log(J)>.
        q = np.einsum(
            "egij,egaj,eg->eai",
            inverse_transpose,
            self._gradient,
            self._weight,
        )
        dp_du = (
            pressure_derivative[:, None, None]
            * q
            / self._volume[:, None, None]
        )
        condensed = np.einsum("eai,ebk->eaibk", q, dp_du)

        # d(F^-T)_ij/dF_kl = -(F^-T)_kj (F^-T)_il.
        geometric = -pressure[:, None, None, None, None] * np.einsum(
            "egkj,egil,egaj,egbl,eg->eaibk",
            inverse_transpose,
            inverse_transpose,
            self._gradient,
            self._gradient,
            self._weight,
        )
        element_tangent = condensed + geometric
        if not np.all(np.isfinite(element_tangent)):
            raise RuntimeError("local-pressure tangent is non-finite")
        rows = np.broadcast_to(
            self._element_dofs[:, :, None], element_tangent.reshape(-1, 24, 24).shape
        )
        columns = np.broadcast_to(
            self._element_dofs[:, None, :], element_tangent.reshape(-1, 24, 24).shape
        )
        return Tangent(rows.ravel(), columns.ravel(), element_tangent.ravel())

    def commit(self, U, state, t, dt):
        # The eliminated pressure is an algebraic function of U, not history.
        return state
