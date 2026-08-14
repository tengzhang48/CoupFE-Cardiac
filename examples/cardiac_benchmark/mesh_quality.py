"""Research diagnostic for the structured Hex8 LV-wall mesh.

The ``apex_offset=0`` variant studied here collapses the apex ring to a single
node, and curved thin-wall elements can be distorted. The diagnostic reports,
per element, the Jacobian determinant at all eight Gauss points (min/max gives
the scaled Jacobian; <=0 is inverted/degenerate), plus volume and aspect ratio,
and localises the worst elements (apex versus bulk).
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import geometry as geom

_NAT = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                 [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)
_G = 1.0 / np.sqrt(3.0)


def _gp_coords():
    GP = np.empty((8, 3))
    for g in range(8):
        GP[g] = [_G if g % 2 == 0 else -_G,
                 _G if (g // 2) % 2 == 0 else -_G,
                 _G if g // 4 == 0 else -_G]
    return GP


def _dN(xi, eta, zeta):
    """dN_a/dξ_i at a natural point (8,3)."""
    d = np.empty((8, 3))
    for a in range(8):
        xa, ya, za = _NAT[a]
        d[a, 0] = 0.125 * xa * (1 + ya * eta) * (1 + za * zeta)
        d[a, 1] = 0.125 * (1 + xa * xi) * ya * (1 + za * zeta)
        d[a, 2] = 0.125 * (1 + xa * xi) * (1 + ya * eta) * za
    return d


def quality(nodes, elems):
    GP = _gp_coords()
    detmin = np.empty(len(elems))
    detmax = np.empty(len(elems))
    vol = np.empty(len(elems))
    aspect = np.empty(len(elems))
    for e, el in enumerate(elems):
        X = nodes[el]
        dets = []
        for g in range(8):
            J = X.T @ _dN(*GP[g])      # ∂x/∂ξ
            dets.append(np.linalg.det(J))
        dets = np.array(dets)
        detmin[e] = dets.min()
        detmax[e] = dets.max()
        vol[e] = dets.mean() * 8.0     # ~ integral of detJ over [-1,1]^3
        # edge-length aspect ratio (12 hex edges)
        E = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]
        L = np.array([np.linalg.norm(X[i] - X[j]) for i, j in E])
        aspect[e] = L.max() / (L.min() + 1e-30)
    return detmin, detmax, vol, aspect


def main():
    for nt, nmu, nth in [(2, 12, 16), (3, 14, 18), (3, 16, 20)]:
        m = geom.build_mesh(nt, nmu, nth, apex_offset=0.0)
        dmin, dmax, vol, asp = quality(m.nodes, m.elems)
        scaled = dmin / (dmax + 1e-30)         # min/max Jacobian (1=uniform)
        cx = np.array([m.nodes[el][:, 0].mean() for el in m.elems])
        apex = cx < -0.06
        print(f"--- mesh nt={nt} nmu={nmu} ntheta={nth}  ({m.n_elem} hexes) ---")
        print(f"  min Jac det (per-elem, over GPs): overall min={dmin.min():.2e}  "
              f"#inverted(<=0)={int(np.sum(dmin <= 0))}")
        print(f"  scaled Jacobian (min/max): min={scaled.min():.3f} "
              f"median={np.median(scaled):.3f}  #poor(<0.2)={int(np.sum(scaled < 0.2))}")
        print(f"     of poor: {int(np.sum((scaled < 0.2) & apex))} at apex, "
              f"{int(np.sum((scaled < 0.2) & ~apex))} in bulk")
        print(f"  aspect ratio: median={np.median(asp):.2f} max={asp.max():.2f}")
        print(f"  bulk(x>-0.06) scaled-Jac: min={scaled[~apex].min():.3f} "
              f"median={np.median(scaled[~apex]):.3f}")


if __name__ == "__main__":
    main()
