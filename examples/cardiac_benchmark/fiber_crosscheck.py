"""Cross-check our fiber field against the authors' reference generator.

Adapts `cardiac_benchmark_toolkit/ellipsoid_fiber_generation.py`
(`FiberExpression.eval`) to NumPy and compares, at matched points, to our
`geometry.fiber_frame` for both helix orientations. A dot product ≈ +1 means
the two local formulas agree; this script does not execute the upstream code.

This file is one of the distribution's two incorporated third-party source
adaptations and remains MIT-licensed. The pinned source, copyright notice, and
the exact provenance and MIT license record are in `THIRD_PARTY_NOTICES.md` and
`LICENSES/Reidmen-cardiac_benchmark_toolkit-MIT.txt`.

SPDX-License-Identifier: MIT
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
import geometry as geom

# authors' DEFAULTS (data.py) — identical to ours
RSE, RSP, RLE, RLP = 2.5e-2, 3.5e-2, 9.0e-2, 9.7e-2
AE, AP = -60.0, 60.0


def authors_fiber(x, td):
    """NumPy adaptation of the pinned FiberExpression.eval formulas."""
    return authors_frame(x, td)[0]


def authors_frame(x, td):
    """Independent NumPy oracle for the pinned toolkit's complete frame."""
    r_s = RSE + (RSP - RSE) * td
    r_l = RLE + (RLP - RLE) * td
    a = np.sqrt(x[1] * x[1] + x[2] * x[2]) / r_s
    b = x[0] / r_l
    u = np.arctan2(a, b)
    v = 0.0 if u < 1e-7 else np.pi - np.arctan2(x[2], -x[1])
    cos_u, sin_v = np.cos(u), np.sin(v)
    sin_u, cos_v = np.sin(u), np.cos(v)
    e_2 = np.array([-r_l * sin_u, r_s * cos_u * cos_v, r_s * cos_u * sin_v])
    e_3 = np.array([0.0, -r_s * sin_u * sin_v, r_s * sin_u * cos_v])
    e_2 = e_2 / np.linalg.norm(e_2)
    e_3 = e_3 / np.linalg.norm(e_3)
    alpha = (AE + (AP - AE) * td) * np.pi / 180.0
    f = np.sin(alpha) * e_2 + np.cos(alpha) * e_3
    f /= np.linalg.norm(f)
    n = np.cross(e_2, e_3)
    n /= np.linalg.norm(n)
    s = np.cross(f, n)
    s /= np.linalg.norm(s)
    return f, s, n


def compare_fiber_formulas():
    """Return dot-product arrays for the two local helix conventions."""
    ts = [0.0, 0.25, 0.5, 0.75, 1.0]
    mus = np.linspace(geom.MU_BASE_ENDO, -np.pi + 0.15, 8)
    ths = np.linspace(-np.pi + 0.1, np.pi - 0.1, 12)
    dot_unflip, dot_flip = [], []
    for t in ts:
        for mu in mus:
            for th in ths:
                x = geom.point(t, mu, th)
                auth = authors_fiber(x, t)
                mine = geom.fiber_frame(t, mu, th, asign=+1.0)[0]
                mine_f = geom.fiber_frame(t, mu, th, asign=-1.0)[0]
                dot_unflip.append(float(mine @ auth))
                dot_flip.append(float(mine_f @ auth))
    return np.array(dot_unflip), np.array(dot_flip)


def main():
    du, df = compare_fiber_formulas()
    print(f"samples: {len(du)}")
    print(f"  unflipped (asign=+1): mean dot={du.mean():+.4f}  "
          f"min|dot|={np.abs(du).min():.4f}  frac aligned(+)={np.mean(du > 0.99):.2f}")
    print(f"  flipped   (asign=-1): mean dot={df.mean():+.4f}  "
          f"min|dot|={np.abs(df).min():.4f}  frac aligned(+)={np.mean(df > 0.99):.2f}")
    aligned_fraction = float(np.mean(df > 0.999))
    if aligned_fraction <= 0.99:
        raise RuntimeError(
            "flipped helix does not match the adapted upstream fiber formula: "
            f"aligned fraction={aligned_fraction:.6f}"
        )
    print("\n=> FLIPPED helix matches the adapted upstream fiber formula.")


if __name__ == "__main__":
    main()
