"""Benchmark-toolkit physical-coordinate structural directions.

This module is a NumPy adaptation of ``FiberExpression.eval``,
``SheetNormalExpression.eval``, and ``SheetExpression.eval`` from the pinned
``Reidmen/cardiac_benchmark_toolkit`` commit recorded in
``examples/REFERENCES.md``.  In particular, it reconstructs the ellipsoidal
coordinates from the physical point and the supplied transmural coordinate;
it does not assume that a point lies on the ellipsoid indexed by a mesh
generator's stored parametric coordinates.

SPDX-License-Identifier: MIT
"""
from __future__ import annotations

import numpy as np


R_SHORT_ENDO = 2.5e-2
R_SHORT_EPI = 3.5e-2
R_LONG_ENDO = 9.0e-2
R_LONG_EPI = 9.7e-2
ALPHA_ENDO = np.deg2rad(-60.0)
ALPHA_EPI = np.deg2rad(+60.0)


def ellipsoid_parameters(tbar, physical_point):
    """Return the toolkit's reconstructed ``(u, v)`` at a physical point.

    The ``v=0`` branch is the apex rule used verbatim by the pinned toolkit.
    ``tbar`` controls the radii of the ellipsoid used for the inverse map.
    """
    tbar = float(tbar)
    point = np.asarray(physical_point, dtype=float)
    if point.shape != (3,) or not np.all(np.isfinite(point)):
        raise ValueError("physical_point must contain three finite coordinates")
    if not np.isfinite(tbar):
        raise ValueError("tbar must be finite")

    r_short = R_SHORT_ENDO + (R_SHORT_EPI - R_SHORT_ENDO) * tbar
    r_long = R_LONG_ENDO + (R_LONG_EPI - R_LONG_ENDO) * tbar
    a = np.hypot(point[1], point[2]) / r_short
    b = point[0] / r_long
    u = float(np.arctan2(a, b))
    v = 0.0 if u < 1.0e-7 else float(
        np.pi - np.arctan2(point[2], -point[1])
    )
    return u, v


def ellipsoid_frame(tbar, physical_point, *, alpha_scale=1.0):
    """Return the pinned toolkit's ``(fiber, sheet, sheet_normal)`` frame.

    ``alpha_scale=1`` is the upstream helix prescription.  The application
    passes ``-asign`` so its historical ``asign=-1``/``flip_helix=True`` CLI
    convention remains source matched.
    """
    tbar = float(tbar)
    point = np.asarray(physical_point, dtype=float)
    alpha_scale = float(alpha_scale)
    if not np.isfinite(alpha_scale):
        raise ValueError("alpha_scale must be finite")
    u, v = ellipsoid_parameters(tbar, point)

    r_short = R_SHORT_ENDO + (R_SHORT_EPI - R_SHORT_ENDO) * tbar
    r_long = R_LONG_ENDO + (R_LONG_EPI - R_LONG_ENDO) * tbar
    cos_u, sin_u = np.cos(u), np.sin(u)
    cos_v, sin_v = np.cos(v), np.sin(v)
    e_u = np.array(
        [
            -r_long * sin_u,
            r_short * cos_u * cos_v,
            r_short * cos_u * sin_v,
        ]
    )
    e_v = np.array(
        [
            0.0,
            -r_short * sin_u * sin_v,
            r_short * sin_u * cos_v,
        ]
    )
    norm_u = float(np.linalg.norm(e_u))
    norm_v = float(np.linalg.norm(e_v))
    if (
        not np.isfinite(norm_u)
        or not np.isfinite(norm_v)
        or norm_u == 0.0
        or norm_v == 0.0
    ):
        raise ValueError("physical point gives a degenerate ellipsoidal basis")
    e_u /= norm_u
    e_v /= norm_v

    helix_angle = alpha_scale * (
        ALPHA_ENDO + (ALPHA_EPI - ALPHA_ENDO) * tbar
    )
    fiber = np.sin(helix_angle) * e_u + np.cos(helix_angle) * e_v
    fiber /= np.linalg.norm(fiber)
    sheet_normal = np.cross(e_u, e_v)
    sheet_normal /= np.linalg.norm(sheet_normal)
    sheet = np.cross(fiber, sheet_normal)
    sheet /= np.linalg.norm(sheet)
    return fiber, sheet, sheet_normal
