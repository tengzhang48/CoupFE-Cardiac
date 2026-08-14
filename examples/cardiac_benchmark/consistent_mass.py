"""Consistent (non-lumped) Q1-Hex8 mass for implicit cardiac dynamics.

The local Simula/FEniCS reference uses a consistent mass matrix. The historical
CoupFE driver used a row-summed lumped mass, which defines a different spatial
discretization of inertia. This module assembles the standard
consistent mass

    M_e = rho * sum_gp w_gp * det(J_gp) * N_gp^T N_gp   (2x2x2 Gauss)

as COO triplets and provides backward-Euler and Newmark inertia operators
with the same ``(residual, tangent, commit)`` contract as the lumped
variants (``coupfe.operators.inertia.InertiaOperator`` and
``newmark.NewmarkInertia``), so they drop into the existing driver
unchanged.

Row sums of the assembled matrix equal the corresponding row-summed mass
vector (partition of unity), which is one of the correctness invariants checked
by the tests.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from coupfe.operators.base import Residual, Tangent

_NAT = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                 [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)


def _gp_table():
    """N and dN/dxi at the eight 2x2x2 Gauss points (run.py ordering)."""
    g1 = 1.0 / np.sqrt(3.0)
    GP = np.empty((8, 3))
    for g in range(8):
        GP[g, 0] = g1 if g % 2 == 0 else -g1
        GP[g, 1] = g1 if (g // 2) % 2 == 0 else -g1
        GP[g, 2] = g1 if g // 4 == 0 else -g1
    N = np.empty((8, 8))
    dN = np.empty((8, 8, 3))
    for g in range(8):
        xi, eta, zeta = GP[g]
        for a in range(8):
            sa, ta, ua = _NAT[a]
            N[g, a] = 0.125 * (1 + sa * xi) * (1 + ta * eta) * (1 + ua * zeta)
            dN[g, a, 0] = 0.125 * sa * (1 + ta * eta) * (1 + ua * zeta)
            dN[g, a, 1] = 0.125 * (1 + sa * xi) * ta * (1 + ua * zeta)
            dN[g, a, 2] = 0.125 * (1 + sa * xi) * (1 + ta * eta) * ua
    return N, dN


_N_GP, _DN_GP = _gp_table()


def consistent_mass_coo(nodes, elems, density, dof_per_node=3):
    """Assemble the consistent Hex8 mass as (rows, cols, vals) COO triplets.

    Per element, the 8x8 nodal block ``rho * sum_gp det(J) N^T N`` is applied
    to each of the ``dof_per_node`` components (node-major global DOF order,
    dof = node * dof_per_node + c, matching the driver layout).
    """
    nodes = np.asarray(nodes, float)
    elems = np.asarray(elems, int)
    if not np.isfinite(density) or density <= 0.0:
        raise ValueError("density must be finite and positive")
    rows, cols, vals = [], [], []
    eye = np.eye(dof_per_node)
    for e in elems:
        xe = nodes[e]                                   # (8, 3)
        Me = np.zeros((8, 8))
        for g in range(8):
            J = _DN_GP[g].T @ xe                        # d x / d xi
            detJ = float(np.linalg.det(J))
            if detJ <= 0.0:
                raise ValueError(f"nonpositive Jacobian in element {e.tolist()}")
            Me += detJ * np.outer(_N_GP[g], _N_GP[g])
        Me *= density
        block = np.kron(Me, eye)                        # node-major 24x24
        gdof = (e[:, None] * dof_per_node
                + np.arange(dof_per_node)).ravel()
        rr, cc = np.meshgrid(gdof, gdof, indexing="ij")
        rows.append(rr.ravel())
        cols.append(cc.ravel())
        vals.append(block.ravel())
    return (np.concatenate(rows), np.concatenate(cols),
            np.concatenate(vals))


class _ConsistentMassBase:
    def __init__(self, rows, cols, vals, ndof):
        self.ndof = int(ndof)
        rows = np.asarray(rows, dtype=np.int64)
        cols = np.asarray(cols, dtype=np.int64)
        vals = np.asarray(vals, dtype=float)
        if rows.shape != cols.shape or rows.shape != vals.shape:
            raise ValueError("COO arrays must have matching shapes")
        self.M = sp.coo_matrix((vals, (rows, cols)),
                               shape=(self.ndof, self.ndof)).tocsr()
        self._coo = (rows, cols, vals)
        self._idx = np.arange(self.ndof)

    def _check_state(self, name, value):
        if value.shape != (self.ndof,) or not np.all(np.isfinite(value)):
            raise ValueError(f"{name} must be a finite vector with shape (ndof,)")


class ConsistentMassInertia(_ConsistentMassBase):
    """Backward-Euler inertia with a consistent mass matrix.

    residual ``M/dt^2 (u - uhat)``, tangent ``M/dt^2``; ``commit`` updates
    ``v = (u - u_prev)/dt`` then ``u_prev = u``. Mirrors
    ``coupfe.operators.inertia.InertiaOperator`` semantics with the diagonal
    mass replaced by the consistent matrix. ``damping`` is the optional
    Rayleigh mass-proportional coefficient (C = alpha * M).
    """

    def __init__(self, rows, cols, vals, ndof, *, u0=None, v0=None,
                 damping=0.0):
        super().__init__(rows, cols, vals, ndof)
        self.u_prev = (np.zeros(ndof) if u0 is None
                       else np.asarray(u0, dtype=float).copy())
        self.v_prev = (np.zeros(ndof) if v0 is None
                       else np.asarray(v0, dtype=float).copy())
        self._check_state("u0", self.u_prev)
        self._check_state("v0", self.v_prev)
        self.damping = float(damping)

    def predictor(self, dt):
        """Inertial predictor ``u_prev + dt * v_prev`` (Newton warm start)."""
        return self.u_prev + dt * self.v_prev

    def residual(self, U, state, t, dt) -> Residual:
        U = np.asarray(U, dtype=float)
        uhat = self.u_prev + dt * self.v_prev
        r = self.M @ ((U - uhat) / dt ** 2)
        if self.damping:
            r = r + self.damping * (self.M @ ((U - self.u_prev) / dt))
        return Residual(self._idx, r)

    def tangent(self, U, state, t, dt) -> Tangent:
        rows, cols, vals = self._coo
        v = vals / dt ** 2
        if self.damping:
            v = v + self.damping * vals / dt
        return Tangent(rows, cols, v)

    def commit(self, U, state, t, dt):
        U = np.asarray(U, dtype=float)
        self.v_prev = (U - self.u_prev) / dt
        self.u_prev = U.copy()
        return state


class ConsistentNewmarkInertia(_ConsistentMassBase):
    """Newmark-beta inertia with a consistent mass matrix.

    Same kinematics as ``newmark.NewmarkInertia`` (beta=0.25, gamma=0.5
    default); the mass enters as the consistent matrix:
    residual ``M @ a_{n+1}``, tangent ``M / (beta dt^2)``.
    """

    def __init__(self, rows, cols, vals, ndof, *, beta=0.25, gamma=0.5,
                 u0=None, v0=None, a0=None):
        super().__init__(rows, cols, vals, ndof)
        self.beta = float(beta)
        self.gamma = float(gamma)
        if not np.isfinite(self.beta) or self.beta <= 0.0:
            raise ValueError("beta must be finite and positive")
        if not np.isfinite(self.gamma) or self.gamma <= 0.0:
            raise ValueError("gamma must be finite and positive")
        self.u_prev = (np.zeros(ndof) if u0 is None
                       else np.asarray(u0, dtype=float).copy())
        self.v_prev = (np.zeros(ndof) if v0 is None
                       else np.asarray(v0, dtype=float).copy())
        self.a_prev = (np.zeros(ndof) if a0 is None
                       else np.asarray(a0, dtype=float).copy())
        self._check_state("u0", self.u_prev)
        self._check_state("v0", self.v_prev)
        self._check_state("a0", self.a_prev)

    def _u_pred(self, dt):
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt must be finite and positive")
        return (self.u_prev + dt * self.v_prev
                + (0.5 - self.beta) * dt * dt * self.a_prev)

    def predictor(self, dt):
        """Warm-start displacement predictor (constant-average-acceleration)."""
        return (self.u_prev + dt * self.v_prev
                + 0.5 * dt * dt * self.a_prev)

    def _accel(self, U, dt):
        return (np.asarray(U, dtype=float) - self._u_pred(dt)) / (
            self.beta * dt * dt)

    def velocity(self, U, dt):
        """Velocity at the trial displacement under the Newmark update."""
        acceleration = self._accel(U, dt)
        return self.v_prev + dt * (
            (1.0 - self.gamma) * self.a_prev + self.gamma * acceleration)

    def velocity_tangent(self, dt):
        """Scalar ``dv_{n+1}/du_{n+1}`` for boundary dashpot tangents."""
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt must be finite and positive")
        return self.gamma / (self.beta * dt)

    def residual(self, U, state, t, dt) -> Residual:
        a = self._accel(U, dt)
        return Residual(self._idx, self.M @ a)

    def tangent(self, U, state, t, dt) -> Tangent:
        rows, cols, vals = self._coo
        return Tangent(rows, cols, vals / (self.beta * dt * dt))

    def commit(self, U, state, t, dt):
        U = np.asarray(U, dtype=float)
        a_new = self._accel(U, dt)
        self.v_prev = self.velocity(U, dt)
        self.a_prev = a_new
        self.u_prev = U.copy()
        return state
