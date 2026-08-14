"""Pericardial Robin (spring-dashpot) boundary operator.

``A_EPI=1e8 Pa/m`` is the published Benchmark 1 value. Local operator checks
do not establish whole-case benchmark agreement; see
``docs/BENCHMARK_COMPARISON.md``.

Paper Eq. 1 boundary terms, acting in the *reference* configuration (so the
spring/dashpot matrices are constant and assembled once):

    base Γ_top :   J·T·F⁻ᵀ·N + α_top·u   + β_top·u̇   = 0   (full vector)
    epi  Γ_epi :   (J·T·F⁻ᵀ·N)·N + α_epi·(u·N) + β_epi·(u̇·N) = 0   (normal only)

Discretised over boundary quad facets with 2×2 Gauss quadrature.  Yields
constant spring (K) and dashpot (C) sparse matrices.  With backward-Euler
(u̇ = (u - u_prev)/dt) the operator contributes

    residual = K·u + (C/dt)·(u - u_prev)
    tangent  = K + C/dt

α_top, α_epi in Pa/m ; β_top, β_epi in Pa·s/m.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from coupfe.operators.base import Residual, Tangent

_G = 1.0 / np.sqrt(3.0)
_GP = [(-_G, -_G), (_G, -_G), (_G, _G), (-_G, _G)]   # 2×2 Gauss, weight 1 each


def _shape(xi, eta):
    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                         (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
    dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
    deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    return N, dxi, deta


_ROBIN_MODES = ("full", "normal", "normal-smoothed")


def _nodal_smoothed_normals(nodes, facets):
    """Area-weighted unit nodal normals over a consistently-oriented surface.

    Each facet contributes its area-weighted unit normal (averaged over the
    2x2 Gauss points) to its four nodes; each nodal accumulator is then
    normalized.  Facet orientation is assumed consistent -- the same
    assumption the default facet-normal mode relies on, that
    ``cross(dX/dxi, dX/deta)`` points along one coherent side of the surface.
    Nodes that receive no contribution keep a zero vector and fall back to the
    facet normal at the quadrature point during assembly.  This is a mechanism
    diagnostic, not a benchmark-faithful operator.
    """
    nodes = np.asarray(nodes, float)
    acc = np.zeros((nodes.shape[0], 3))
    for f in facets:
        X = nodes[f]
        fvec = np.zeros(3)
        farea = 0.0
        for xi, eta in _GP:
            _, dxi, deta = _shape(xi, eta)
            cx = np.cross(dxi @ X, deta @ X)
            detA = np.linalg.norm(cx)
            if not np.isfinite(detA) or detA <= 0.0:
                raise ValueError("Robin boundary contains a degenerate facet")
            fvec += detA * (cx / detA)
            farea += detA
        nrm = np.linalg.norm(fvec)
        if nrm > 0.0:
            fvec = fvec / nrm
        for a in f:
            acc[int(a)] += farea * fvec
    norms = np.linalg.norm(acc, axis=1)
    unit = np.zeros_like(acc)
    nz = norms > 0.0
    unit[nz] = acc[nz] / norms[nz, None]
    return unit


def _assemble(nodes, facet_groups, ndof, dof_per_node=3):
    """Build constant spring K and dashpot C COO triplets.

    facet_groups: list of (facets (M,4) int, alpha, beta, mode) with
    mode in {"full", "normal", "normal-smoothed"}. ``dof_per_node`` is the
    global stride; the springs act on displacement components 0..2 only.

    "normal" (default) uses the per-facet-Gauss-point normal
    ``cross(dX/dxi, dX/deta)/detA``.  "normal-smoothed" replaces only that
    normal field with a Q1-interpolated, area-weighted nodal-smoothed normal;
    the facet area weights (``detA``), 2x2 quadrature, coefficients, and
    assembly structure are identical, so the default mode is bit-identical to
    the historical operator.
    """
    dpn = dof_per_node
    rows, cols, kv, cv = [], [], [], []
    for facets, alpha, beta, mode in facet_groups:
        if mode not in _ROBIN_MODES:
            raise ValueError(f"unsupported Robin projection mode {mode!r}")
        nodal_n = (_nodal_smoothed_normals(nodes, facets)
                   if mode == "normal-smoothed" else None)
        for f in facets:
            X = nodes[f]                      # (4,3) reference coords
            for xi, eta in _GP:
                N, dxi, deta = _shape(xi, eta)
                g_xi = dxi @ X                # (3,)
                g_eta = deta @ X
                cx = np.cross(g_xi, g_eta)
                detA = np.linalg.norm(cx)
                if not np.isfinite(detA) or detA <= 0.0:
                    raise ValueError("Robin boundary contains a degenerate facet")
                if mode == "full":
                    proj = np.eye(3)
                elif mode == "normal":
                    nrm = cx / detA
                    proj = np.outer(nrm, nrm)
                else:                          # normal-smoothed: Q1-interpolated nodal normal
                    nq = nodal_n[f].T @ N      # Σ_a N_a(ξ,η) n_a
                    nn_q = np.linalg.norm(nq)
                    nq = nq / nn_q if nn_q > 0.0 else cx / detA
                    proj = np.outer(nq, nq)
                for a in range(4):
                    for b in range(4):
                        w = N[a] * N[b] * detA
                        for ci in range(3):
                            for cj in range(3):
                                pij = proj[ci, cj]
                                if pij == 0.0:
                                    continue
                                rows.append(dpn * f[a] + ci)
                                cols.append(dpn * f[b] + cj)
                                kv.append(alpha * pij * w)
                                cv.append(beta * pij * w)
    rows = np.asarray(rows, int); cols = np.asarray(cols, int)
    kv = np.asarray(kv, float); cv = np.asarray(cv, float)
    Kmat = sp.csr_matrix((kv, (rows, cols)), shape=(ndof, ndof))
    Cmat = sp.csr_matrix((cv, (rows, cols)), shape=(ndof, ndof))
    return Kmat, Cmat, rows, cols, kv, cv


def _rigid_spring_reduction(nodes, facets, alpha, mode):
    """Reduce one spring group onto rigid translations and x-axis rotation.

    This is a diagnostic of the *discrete* quadrilateral boundary geometry.
    It reuses the operator quadrature but neither changes the assembled Robin
    condition nor defines a mesh-validation target.
    """
    nodes = np.asarray(nodes, float)
    nodal_n = (_nodal_smoothed_normals(nodes, facets)
               if mode == "normal-smoothed" else None)
    translation = np.zeros((3, 3))
    long_axis_rotation = 0.0
    for facet in facets:
        coordinates = nodes[facet]
        for xi, eta in _GP:
            shape, dxi, deta = _shape(xi, eta)
            cross = np.cross(dxi @ coordinates, deta @ coordinates)
            area_weight = np.linalg.norm(cross)
            if not np.isfinite(area_weight) or area_weight <= 0.0:
                raise ValueError("Robin boundary contains a degenerate facet")
            if mode == "full":
                projection = np.eye(3)
            elif mode == "normal":
                normal = cross / area_weight
                projection = np.outer(normal, normal)
            elif mode == "normal-smoothed":
                nq = nodal_n[facet].T @ shape
                nn_q = np.linalg.norm(nq)
                nq = nq / nn_q if nn_q > 0.0 else cross / area_weight
                projection = np.outer(nq, nq)
            else:
                raise ValueError(f"unsupported Robin projection mode {mode!r}")
            translation += alpha * projection * area_weight
            point = shape @ coordinates
            rotation = np.asarray([0.0, -point[2], point[1]])
            long_axis_rotation += (
                alpha * float(rotation @ projection @ rotation) * area_weight
            )
    return {
        "mode": mode,
        "translation_matrix": translation,
        "long_axis_rotation": float(long_axis_rotation),
    }


class RobinOperator:
    """Constant spring-dashpot surface operator (Operator contract)."""

    def __init__(
        self,
        nodes,
        ndof,
        facet_groups,
        u0=None,
        dof_per_node=3,
        kinematics=None,
    ):
        self.ndof = int(ndof)
        nodes = np.asarray(nodes, float)
        facet_groups = tuple(facet_groups)
        self.Kmat, self.Cmat, self._rows, self._cols, self._kv, self._cv = \
            _assemble(nodes, facet_groups, self.ndof, dof_per_node)
        self.rigid_spring_groups = tuple(
            _rigid_spring_reduction(nodes, facets, alpha, mode)
            for facets, alpha, _beta, mode in facet_groups
        )
        self.u_prev = np.zeros(ndof) if u0 is None else np.asarray(u0, float).copy()
        self.kinematics = kinematics
        self._dofs = np.unique(self._rows)
        self.dofs = self._dofs          # public name for the core solve_dynamics_distributed robin= contract

    def residual(self, U, state, t, dt) -> Residual:
        U = np.asarray(U, dtype=float)
        if self.kinematics is None:
            if not np.isfinite(dt) or dt <= 0.0:
                raise ValueError("dt must be finite and positive")
            velocity = (U - self.u_prev) / dt
        else:
            velocity = self.kinematics.velocity(U, dt)
        r = self.Kmat @ U + self.Cmat @ velocity
        return Residual(self._dofs, r[self._dofs])

    def tangent(self, U, state, t, dt) -> Tangent:
        velocity_factor = (
            1.0 / dt
            if self.kinematics is None
            else self.kinematics.velocity_tangent(dt)
        )
        return Tangent(self._rows, self._cols, self._kv + self._cv * velocity_factor)

    def commit(self, U, state, t, dt):
        if self.kinematics is None:
            self.u_prev = np.asarray(U, dtype=float).copy()
        return state


if __name__ == "__main__":
    # self-check: uniform displacement of the base produces a restoring force
    # equal to alpha_top * base_area * u (full-vector springs), and the spring
    # matrix is symmetric.
    import os, sys
    sys.path.insert(0, os.path.dirname(__file__))
    from geometry import build_mesh

    m = build_mesh(n_t=2, n_mu=10, n_theta=16, apex_offset=0.2)
    ndof = 3 * m.n_node
    A_TOP, B_TOP, A_EPI, B_EPI = 1e5, 5e3, 1e8, 5e3
    op = RobinOperator(
        m.nodes, ndof,
        [(m.facets_base, A_TOP, B_TOP, "full"),
         (m.facets_epi, A_EPI, B_EPI, "normal")],
    )

    # base area (reference) from facet quadrature
    base_area = 0.0
    for f in m.facets_base:
        X = m.nodes[f]
        for xi, eta in _GP:
            _, dxi, deta = _shape(xi, eta)
            base_area += np.linalg.norm(np.cross(dxi @ X, deta @ X))
    print(f"base area = {base_area*1e4:.3f} cm^2")

    # uniform u_z displacement, very large dt so dashpot negligible
    u = np.zeros(ndof)
    uz = 1e-3
    u[2::3] = uz
    R = op.residual(u, None, t=0.0, dt=1e12)
    fz = R.values[(R.gdofs % 3) == 2].sum()
    # only base contributes a z-reaction for uniform u_z if epi normals cancel?
    # epi also resists the normal component, so compute base-only separately:
    op_base = RobinOperator(m.nodes, ndof, [(m.facets_base, A_TOP, B_TOP, "full")])
    Rb = op_base.residual(u, None, t=0.0, dt=1e12)
    fz_base = Rb.values[(Rb.gdofs % 3) == 2].sum()
    print(f"base z-reaction      = {fz_base:.4f} N")
    print(f"alpha_top*area*uz    = {A_TOP*base_area*uz:.4f} N")
    sym_err = abs((op.Kmat - op.Kmat.T)).max()
    print(f"K symmetry error     = {sym_err:.2e}")
    reaction_error = abs(fz_base - A_TOP * base_area * uz) / (A_TOP * base_area * uz)
    if reaction_error >= 1e-10:
        raise RuntimeError(f"Robin reaction check failed: error={reaction_error:.3e}")
    if sym_err >= 1e-6:
        raise RuntimeError(f"Robin symmetry check failed: error={sym_err:.3e}")
    print("ROBIN SELF-CHECK PASSED")
