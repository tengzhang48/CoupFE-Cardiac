"""Endocardial follower-pressure boundary operator.

Paper Eq. 1 endocardial Neumann term: J·T·F⁻ᵀ·N = p·J·F⁻ᵀ·N on Γ_endo — a
*follower* (deformation-dependent) pressure load.  By Nanson's formula
J·F⁻ᵀ·N dA = n da, the external nodal force is

    f_a = p ∫ N_a (a_ξ × a_η) dξ dη ,   a_ξ = Σ_b dN_b/dξ (X_b+u_b)  (deformed)

integrated over the deformed boundary quad facets.  The tangent is
deformation-dependent (the classic follower-load / pressure stiffness)
→ complex-stepped per facet.

Pressure p [Pa] is a positive cavity pressure updated per step by the driver.
The operator uses `+p*(n da)` with each facet area vector oriented out of the
solid. Local tests cover orientation reversal and tangent consistency. Earlier
full-case cavity-volume observations were not retained and are not a release
validation claim; see `docs/CASE_B_STATUS.md`.
"""
from __future__ import annotations

import numpy as np

from coupfe.operators.base import Residual, Tangent, complex_step_tangent

_G = 1.0 / np.sqrt(3.0)
_GP = [(-_G, -_G), (_G, -_G), (_G, _G), (-_G, _G)]


def _shape(xi, eta):
    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                         (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
    dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
    deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
    return N, dxi, deta


def _facet_residual(Xf, ue, p, orient):
    """Residual (12,) for one quad facet; ue (12,) facet displacements.

    Implementation of the benchmark term in paper Eq. 1,
    p · J·F⁻ᵀ·N · δu (= p·n·da in the current configuration), with N the
    **out-of-solid** normal (into the cavity on the
    endocardium).  ``orient`` (±1) fixes the facet's a_ξ×a_η sense to that
    consistently-out-of-solid direction; the residual is **+p·n·da** under the
    repository's residual convention. Analytic in ue → complex-step safe.
    """
    xf = Xf + ue.reshape(4, 3)              # deformed facet coords
    r = np.zeros(12, dtype=ue.dtype)
    for xi, eta in _GP:
        N, dxi, deta = _shape(xi, eta)
        a_xi = dxi @ xf                      # (3,) deformed tangent
        a_eta = deta @ xf
        nda = orient * np.cross(a_xi, a_eta)  # out-of-solid (into cavity) n·da
        for a in range(4):
            # Reversing orient reverses the complete local pressure load; the
            # paired broken-control test protects this convention.
            r[3 * a:3 * a + 3] += p * N[a] * nda
    return r


def _facet_residual_batch(Xf, ue, p, orient):
    """Vectorised over ALL facets at once.  Xf,ue: (Nf,4,3); orient: (Nf,);
    returns r: (Nf,12).  Same `+p·n·da` term as `_facet_residual`, just batched
    (so the complex-step tangent costs 12 vector ops, not 12·Nf scalar ones).
    Analytic in ue (poly + cross) → complex-step safe (ue may be complex)."""
    xf = Xf + ue                                  # (Nf,4,3)
    r = np.zeros((xf.shape[0], 12), dtype=ue.dtype)
    for xi, eta in _GP:
        N, dxi, deta = _shape(xi, eta)
        a_xi = np.einsum('a,nad->nd', dxi, xf)    # (Nf,3) deformed tangents
        a_eta = np.einsum('a,nad->nd', deta, xf)
        nda = orient[:, None] * np.cross(a_xi, a_eta)   # (Nf,3) out-of-solid n·da
        for a in range(4):
            r[:, 3 * a:3 * a + 3] += p * N[a] * nda
    return r


class FollowerPressureOperator:
    """Deformation-dependent endocardial pressure load (Operator contract).

    Normals are oriented consistently **out of the solid** (into the cavity on the
    endocardium) per facet. The residual is
    **+p·(n da)** under the repository residual convention. Pass the driver's
    parent-element-centroid ``interior`` so each facet's orientation matches the
    solve. This local convention is not by itself a full Case B validation.
    """

    def __init__(self, nodes, ndof, facets, p=0.0, dof_per_node=3, interior=None):
        self.nodes = np.asarray(nodes, float)
        self.ndof = int(ndof)
        self.facets = np.asarray(facets, int)
        self.p = float(p)              # updated per step by the driver
        # Acts on displacement components 0..2; dof_per_node is the global stride.
        self.gm = (self.facets[:, :, None] * dof_per_node
                   + np.arange(3)[None, None, :]).reshape(len(self.facets), 12)
        # Orient each facet's reference normal OUT of the solid (into the cavity):
        # away from its parent-element centroid. A single cavity centroid can
        # mis-sign tapered facets, so this input is mandatory.
        if interior is None:
            raise ValueError(
                "interior parent-element centroids are required for facet orientation"
            )
        interior = np.asarray(interior, float)
        if interior.shape != (len(self.facets), 3) or not np.all(np.isfinite(interior)):
            raise ValueError(
                f"interior must have finite shape {(len(self.facets), 3)}"
            )
        self._orient = np.empty(len(self.facets))
        N0, dxi0, deta0 = _shape(0.0, 0.0)
        for i, f in enumerate(self.facets):
            X = self.nodes[f]
            n_ref = np.cross(dxi0 @ X, deta0 @ X)
            fc = X.mean(axis=0)
            # out-of-solid points from the interior (element centroid) toward the facet
            self._orient[i] = 1.0 if np.dot(n_ref, fc - interior[i]) > 0 else -1.0
        self._Xf = self.nodes[self.facets]          # (Nf,4,3) reference facet coords

    def residual(self, U, state, t, dt) -> Residual:
        U = np.asarray(U, float)
        ue = U[self.gm].reshape(-1, 4, 3)
        r = _facet_residual_batch(self._Xf, ue, self.p, self._orient)   # (Nf,12)
        return Residual(self.gm.ravel(), r.ravel())

    def tangent(self, U, state, t, dt) -> Tangent:
        # Vectorised complex-step over all facets: 12 batched passes (one per local
        # dof) instead of 12·Nf scalar facet evals — ~100x faster on a fine mesh.
        U = np.asarray(U, float)
        ue = U[self.gm].reshape(-1, 4, 3).astype(complex)
        Nf = ue.shape[0]
        h = 1e-30
        K = np.empty((Nf, 12, 12))
        for j in range(12):                          # perturb local dof j (node j//3, comp j%3)
            up = ue.copy()
            up[:, j // 3, j % 3] += 1j * h
            rp = _facet_residual_batch(self._Xf, up, self.p, self._orient)
            K[:, :, j] = rp.imag / h
        g = self.gm                                  # (Nf,12) global dofs
        rows = np.repeat(g, 12, axis=1)              # (Nf,144): g[a] repeated over b
        cols = np.tile(g, (1, 12))                   # (Nf,144): g[b] tiled over a
        data = K.reshape(Nf, 144)                    # K[n,a,b] row-major matches rows/cols
        return Tangent(rows.ravel(), cols.ravel(), data.ravel())

    def commit(self, U, state, t, dt):
        return state


if __name__ == "__main__":
    # Self-check: finiteness + tangent consistency, built with the DRIVER's
    # orientation (parent-element-centroid `interior`). This self-check covers
    # finiteness and tangent construction only; full-case cavity response remains
    # outside this local operator check.
    import os, sys
    sys.path.insert(0, os.path.dirname(__file__))
    from geometry import build_mesh

    m = build_mesh(n_t=2, n_mu=10, n_theta=16, apex_offset=0.2)
    ndof = 3 * m.n_node
    interior = m.nodes[m.elems[m.facets_endo_elem]].mean(axis=1)   # = the driver's
    op = FollowerPressureOperator(m.nodes, ndof, m.facets_endo, p=1000.0,
                                  interior=interior)
    R = op.residual(np.zeros(ndof), None, 0.0, 1.0)
    f = np.zeros(ndof)
    np.add.at(f, R.gdofs, R.values)
    net = f.reshape(-1, 3).sum(axis=0)
    print(f"apex x = {m.nodes[:, 0].min():+.4f} m  (apex at -x, base at +x)")
    print(f"net pressure residual (report only): fx={net[0]:+.3f} fy={net[1]:+.3f} fz={net[2]:+.3f}")
    print(f"residual finite: {np.all(np.isfinite(R.values))}")
    T = op.tangent(np.zeros(ndof), None, 0.0, 1.0)
    print(f"tangent entries: {len(T.values)}, finite: {np.all(np.isfinite(T.values))}")
    if not np.all(np.isfinite(R.values)) or not np.all(np.isfinite(T.values)):
        raise RuntimeError("pressure residual or tangent contains non-finite values")
    print("PRESSURE SELF-CHECK PASSED (finite residual + tangent construction)")
