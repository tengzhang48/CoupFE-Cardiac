"""Monoventricular truncated-ellipsoid geometry + analytic fibers.

Adapts the geometry and fiber equations used in Benchmark 1 of Aróstica et al.
(CMAME 2025), "A software benchmark for cardiac elastodynamics". All lengths
are in metres (SI).

Two application-owned Hex8 topologies are available.  The historical polar-ring
mesh uses parametric coordinates (t̃, μ, θ) and either collapses its last ring
or truncates it before the apex.  The benchmark-reproduction topology maps a
structured square bijectively to the endocardial/epicardial surface disks and
extrudes it through the wall.  The square center is one ordinary mesh vertex at
the apex, so this topology is closed without collapsed Hex8 cells or an added
tip boundary.

The mesh coordinates remain:

    t̃ in [0,1]   transmural   (0 = endocardium, 1 = epicardium)
    μ            longitudinal (μ_base ... -π apex)
    θ in [-π,π]  circumferential

The historical polar-ring mesh uses the ellipsoid map below directly.  The
closed benchmark mesh instead follows the toolkit's planar meridional Gmsh
surface: corresponding points on the endocardial and epicardial ellipses are
joined by straight Cartesian segments.  Its stored parametric coordinates are
retained for topology provenance and the historical polar rule.  Closed-mesh
structural directions reconstruct the toolkit's local coordinates from each
physical point and can replace the analytic transmural value with the
prescribed Laplace field.

Surface parametrisation (paper Eq. 12):
    x(μ,θ,t̃) = ( r_long(t̃) cos μ,
                 r_short(t̃) sin μ cos θ,
                 r_short(t̃) sin μ sin θ )

Endo/epi radii (Eqs. 10-11) and base truncation angles.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

if __package__:
    from . import structural_directions
else:
    import structural_directions

# --- geometry constants (metres / radians) --------------------------------
RS_ENDO = structural_directions.R_SHORT_ENDO
RL_ENDO = structural_directions.R_LONG_ENDO
RS_EPI = structural_directions.R_SHORT_EPI
RL_EPI = structural_directions.R_LONG_EPI
MU_BASE_ENDO = -np.arccos(5.0 / 17.0)      # endo base truncation angle
MU_BASE_EPI = -np.arccos(5.0 / 20.0)       # epi  base truncation angle
ALPHA_ENDO = structural_directions.ALPHA_ENDO
ALPHA_EPI = structural_directions.ALPHA_EPI

# reference output points (paper, metres)
P0 = np.array([0.025, 0.030, 0.0])
P1 = np.array([0.000, 0.030, 0.0])

PHYSICAL_DIRECTION_RECONSTRUCTION = "toolkit-physical-coordinate-u-v-v1"
HISTORICAL_DIRECTION_RECONSTRUCTION = "historical-parametric-mu-theta-v1"


def r_long(t):
    return RL_ENDO + (RL_EPI - RL_ENDO) * t


def r_short(t):
    return RS_ENDO + (RS_EPI - RS_ENDO) * t


def mu_base(t):
    return MU_BASE_ENDO + (MU_BASE_EPI - MU_BASE_ENDO) * t


def alpha(t):
    return ALPHA_ENDO + (ALPHA_EPI - ALPHA_ENDO) * t


def point(t, mu, theta):
    """Reference position for parametric (t̃, μ, θ)."""
    rl, rs = r_long(t), r_short(t)
    return np.array([rl * np.cos(mu),
                     rs * np.sin(mu) * np.cos(theta),
                     rs * np.sin(mu) * np.sin(theta)])


def closed_wall_point(t, rho, theta):
    """Map closed-mesh coordinates using the toolkit's Cartesian ruling.

    ``rho`` runs from the apex (0) to the base (1).  At fixed ``(rho,
    theta)``, the endocardial and epicardial locations use their own radii and
    base angles.  The interior wall point is the affine Cartesian
    interpolation of those two locations.  In particular, both the apex and
    base transmural curves are straight, matching the ``Line`` boundaries of
    the Benchmark 1 Gmsh meridional surface.
    """
    mu_endo = -np.pi + rho * (np.pi + MU_BASE_ENDO)
    mu_epi = -np.pi + rho * (np.pi + MU_BASE_EPI)
    endo = point(0.0, mu_endo, theta)
    epi = point(1.0, mu_epi, theta)
    return (1.0 - t) * endo + t * epi


def fiber_frame(t, mu, theta, asign=1.0):
    """Return (f, s, n) orthonormal fiber/sheet/normal frame (paper Eq. 14-15).

    f = sin(α) e_μ + cos(α) e_θ ;  n = (e_μ × e_θ)/|·| ;  s = (f × n)/|·|
    ``asign=-1`` reverses the transmural helix gradient (chirality flip): it
    flips the longitudinal component sin(α) while leaving cos(α) unchanged,
    reversing the prescribed helix handedness. This change does not by itself
    predict the sign of a ventricular displacement or twist response.
    Raises ValueError near the apex (μ→-π) where e_θ degenerates.
    """
    rl, rs = r_long(t), r_short(t)
    smu, cmu = np.sin(mu), np.cos(mu)
    sth, cth = np.sin(theta), np.cos(theta)
    e_mu = np.array([-rl * smu, rs * cmu * cth, rs * cmu * sth])
    e_th = np.array([0.0, -rs * smu * sth, rs * smu * cth])
    nmu, nth = np.linalg.norm(e_mu), np.linalg.norm(e_th)
    if nmu < 1e-12 or nth < 1e-12:
        raise ValueError("degenerate basis (apex)")
    e_mu /= nmu
    e_th /= nth
    a = asign * alpha(t)
    f = np.sin(a) * e_mu + np.cos(a) * e_th
    f /= np.linalg.norm(f)
    n = np.cross(e_mu, e_th)
    n /= np.linalg.norm(n)
    s = np.cross(f, n)
    s /= np.linalg.norm(s)
    return f, s, n


def toolkit_fiber_frame(tbar, physical_point, *, asign=-1.0):
    """Return the pinned toolkit frame reconstructed from physical coordinates.

    The application's historical ``asign=-1`` convention is the source-matched
    orientation, so the sign is reversed when it is passed to the toolkit's
    helix-angle prescription.
    """
    return structural_directions.ellipsoid_frame(
        tbar,
        physical_point,
        alpha_scale=-float(asign),
    )


def structural_direction_reconstruction(topology):
    """Return the provenance identifier for a mesh topology's frame rule."""
    if topology == "closed_multiblock_disk":
        return PHYSICAL_DIRECTION_RECONSTRUCTION
    if topology == "polar_ring":
        return HISTORICAL_DIRECTION_RECONSTRUCTION
    raise ValueError(f"unsupported structural-frame topology {topology!r}")


def _circular_shape_average(theta, weights):
    """Interpolate a periodic angle with shape-function weights."""
    theta = np.asarray(theta, dtype=float)
    weights = np.asarray(weights, dtype=float)
    finite = np.isfinite(theta)
    if theta.shape != weights.shape or not np.any(finite):
        raise ValueError("circular interpolation needs finite nodal angles")
    scale = float(np.sum(weights[finite]))
    if not np.isfinite(scale) or abs(scale) <= 1.0e-14:
        raise ValueError("circular interpolation has zero finite-node weight")
    phasor = np.sum(weights[finite] * np.exp(1j * theta[finite])) / scale
    if not np.isfinite(phasor) or abs(phasor) <= 1.0e-12:
        raise ValueError("circular interpolation is directionally degenerate")
    return float(np.angle(phasor))


def sample_structural_frame(mesh, element_nodes, weights, tbar, *, asign=-1.0):
    """Evaluate one structural frame using the topology's declared rule.

    Closed benchmark meshes use the pinned toolkit inverse map at the physical
    Q1 point ``X = N @ nodes``.  Historical polar meshes retain their original
    parametric interpolation and near-apex fallback.
    """
    element_nodes = np.asarray(element_nodes, dtype=int)
    weights = np.asarray(weights, dtype=float)
    tbar = np.asarray(tbar, dtype=float)
    if element_nodes.shape != weights.shape:
        raise ValueError("element nodes and shape weights must have equal length")
    if tbar.shape != (mesh.n_node,):
        raise ValueError("fiber transmural coordinates must be per-node data")
    element_tbar = tbar[element_nodes]
    if not np.all(np.isfinite(element_tbar)):
        raise ValueError("element fiber transmural coordinates must be finite")
    t_gp = float(weights @ element_tbar)

    if mesh.topology == "closed_multiblock_disk":
        physical_point = weights @ mesh.nodes[element_nodes]
        return toolkit_fiber_frame(t_gp, physical_point, asign=asign)
    if mesh.topology != "polar_ring":
        raise ValueError(f"unsupported structural-frame topology {mesh.topology!r}")

    parametric = mesh.param[element_nodes]
    mu_gp = float(weights @ parametric[:, 1])
    try:
        theta_gp = _circular_shape_average(parametric[:, 2], weights)
        return fiber_frame(t_gp, mu_gp, theta_gp, asign=asign)
    except ValueError:
        return fiber_frame(t_gp, -np.pi + 0.05, 0.0, asign=asign)


def _hex_volume(p):
    """Signed volume of an 8-node hex (p: (8,3)), CoupFE node order.

    Order: bottom face [0,1,2,3] CCW, top face [4,5,6,7] above them.
    Decompose into 5 tets (verified to tile the unit cube to volume 1.0).
    """
    # tet decomposition of a hexahedron (consistent with VTK_HEXAHEDRON)
    tets = [(0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6),
            (1, 4, 5, 6), (3, 4, 6, 7)]
    vol = 0.0
    for a, b, c, d in tets:
        vol += np.dot(np.cross(p[b] - p[a], p[c] - p[a]), p[d] - p[a]) / 6.0
    return vol


@dataclass
class CardiacMesh:
    nodes: np.ndarray          # (Nnode, 3) reference coords
    elems: np.ndarray          # (Nelem, 8) Hex8 connectivity
    param: np.ndarray          # (Nnode, 3) parametric (t̃, μ, θ); apex θ=nan
    fiber: np.ndarray          # (Nelem, 3) per-element centroid fiber f0
    sheet: np.ndarray          # (Nelem, 3) per-element centroid sheet s0
    normal: np.ndarray         # (Nelem, 3) per-element centroid normal n0
    fiber_node: np.ndarray     # (Nnode, 3) per-node fiber f0 (for per-GP interp)
    sheet_node: np.ndarray     # (Nnode, 3) per-node sheet s0
    facets_endo: np.ndarray    # (M,4) node ids on endocardium (t̃=0)
    facets_endo_elem: np.ndarray  # (M,) parent element index per endo facet
    facets_epi: np.ndarray     # (M,4) node ids on epicardium  (t̃=1)
    facets_base: np.ndarray    # (M,4) node ids on base        (μ=μ_base)
    n_t: int = 0
    n_mu: int = 0
    n_theta: int = 0
    topology: str = "polar_ring"
    n_side: int = 0
    n_core: int = 0
    n_radial: int = 0
    core_half_width: float = 0.0
    tip_refine: float = 1.0

    @property
    def n_node(self):
        return len(self.nodes)

    @property
    def n_elem(self):
        return len(self.elems)


def update_structural_frames(mesh, tbar, *, asign=-1.0):
    """Rebuild stored nodal and element frames for a transmural field.

    This is the common path used both for the mesh's initial analytic field and
    for an injected Laplace field.  Only the closed benchmark topology adopts
    physical-coordinate reconstruction; the historical polar topology remains
    bit-for-bit on its parametric convention, including its apex fallbacks.
    """
    tbar = np.asarray(tbar, dtype=float)
    if tbar.shape != (mesh.n_node,) or not np.all(np.isfinite(tbar)):
        raise ValueError("fiber transmural coordinates must be finite per-node data")
    asign = float(asign)

    if mesh.topology == "closed_multiblock_disk":
        for node in range(mesh.n_node):
            fiber, sheet, _normal = toolkit_fiber_frame(
                tbar[node], mesh.nodes[node], asign=asign
            )
            mesh.fiber_node[node] = fiber
            mesh.sheet_node[node] = sheet
        center_weights = np.full(8, 1.0 / 8.0)
        for element, connectivity in enumerate(mesh.elems):
            fiber, sheet, normal = sample_structural_frame(
                mesh,
                connectivity,
                center_weights,
                tbar,
                asign=asign,
            )
            mesh.fiber[element] = fiber
            mesh.sheet[element] = sheet
            mesh.normal[element] = normal
        return

    if mesh.topology != "polar_ring":
        raise ValueError(f"unsupported structural-frame topology {mesh.topology!r}")

    for node in range(mesh.n_node):
        mu, theta = mesh.param[node, 1:]
        try:
            if not np.isfinite(theta):
                raise ValueError
            fiber, sheet, _normal = fiber_frame(
                tbar[node], mu, theta, asign=asign
            )
        except ValueError:
            fiber, sheet, _normal = fiber_frame(
                tbar[node], -np.pi + 0.05, 0.0, asign=asign
            )
        mesh.fiber_node[node] = fiber
        mesh.sheet_node[node] = sheet

    for element, connectivity in enumerate(mesh.elems):
        parametric = mesh.param[connectivity]
        t_center = float(np.mean(tbar[connectivity]))
        mu_center = float(np.nanmean(parametric[:, 1]))
        theta = parametric[:, 2]
        theta = theta[np.isfinite(theta)]
        theta_center = (
            float(np.angle(np.mean(np.exp(1j * theta)))) if len(theta) else 0.0
        )
        try:
            fiber, sheet, normal = fiber_frame(
                t_center, mu_center, theta_center, asign=asign
            )
        except ValueError:
            fiber = np.array([1.0, 0.0, 0.0])
            sheet = np.array([0.0, 1.0, 0.0])
            normal = np.array([0.0, 0.0, 1.0])
        mesh.fiber[element] = fiber
        mesh.sheet[element] = sheet
        mesh.normal[element] = normal


def build_mesh(n_t: int = 3, n_mu: int = 16, n_theta: int = 24,
               flip_helix: bool = True, apex_offset: float = 0.2) -> CardiacMesh:
    """Structured Hex8 LV-wall mesh.

    n_t    : transmural element layers
    n_mu   : longitudinal element layers (base -> final ring/apex)
    n_theta: circumferential element layers (periodic)
    flip_helix: use the chirality that matches the pinned NumPy adaptation of the
        upstream reference formula (DEFAULT; checked by fiber_crosscheck.py). Our (μ,θ)
        parametrisation differs from their inverse-map (u,v) by a circumferential
        offset that flips the helix; negating the helix angle realigns it.
    apex_offset: angular distance of the final ring from the apex. The
        noncollapsed demonstration default is 0.2; zero collapses the ring and creates
        degenerate apex cells/facets.
    """
    if n_t < 1 or n_mu < 1 or n_theta < 3:
        raise ValueError("n_t and n_mu must be positive; n_theta must be at least 3")
    apex_limit = np.pi + min(MU_BASE_ENDO, MU_BASE_EPI)
    if (
        not np.isfinite(apex_offset)
        or apex_offset < 0.0
        or apex_offset >= apex_limit
    ):
        raise ValueError(
            f"apex_offset must satisfy 0 <= value < {apex_limit:.6g} radians"
        )
    asign = -1.0 if flip_helix else 1.0
    nodes: list[list[float]] = []
    param: list[list[float]] = []
    nid: dict = {}

    def add(key, t, mu, theta):
        if key in nid:
            return nid[key]
        nid[key] = len(nodes)
        nodes.append(point(t, mu, theta).tolist())
        param.append([t, mu, theta if key[1] != "apex" else np.nan])
        return nid[key]

    n_t_node = n_t + 1
    # Apex handling.  apex_offset == 0: the i_mu = n_mu ring collapses to a
    # single apex node -> degenerate hexes + degenerate endo/epi facets at the apex.
    # With apex_offset > 0 the last ring sits at mu = -pi + apex_offset (a
    # nonzero-radius full theta ring), removing repeated apex nodes in the
    # checked configurations. This leaves an open apex boundary and therefore
    # defines a different geometry; whole-case sensitivity to that choice
    # remains a research question.
    collapse = apex_offset <= 0.0
    mu_apex = -np.pi + apex_offset
    n_mu_ring = n_mu if collapse else n_mu + 1   # number of full theta-rings generated
    for it in range(n_t_node):
        t = it / n_t
        for imu in range(n_mu_ring):
            s = imu / n_mu
            mu = mu_base(t) * (1.0 - s) + mu_apex * s
            for ith in range(n_theta):
                theta = -np.pi + 2.0 * np.pi * ith / n_theta
                add((it, imu, ith), t, mu, theta)
        if collapse:
            add((it, "apex"), t, -np.pi, 0.0)    # single apex node per t-layer

    def node(it, imu, ith):
        if imu == n_mu and collapse:
            return nid[(it, "apex")]
        return nid[(it, imu, ith % n_theta)]

    elems: list[list[int]] = []
    for it in range(n_t):                    # transmural element layer
        for imu in range(n_mu):              # longitudinal element layer
            for ith in range(n_theta):       # circumferential
                j = ith
                jp = (ith + 1) % n_theta
                # bottom face = inner transmural layer it, top = it+1
                b = [node(it, imu, j), node(it, imu, jp),
                     node(it, imu + 1, jp), node(it, imu + 1, j)]
                t_ = [node(it + 1, imu, j), node(it + 1, imu, jp),
                      node(it + 1, imu + 1, jp), node(it + 1, imu + 1, j)]
                elems.append(b + t_)

    nodes = np.asarray(nodes, dtype=float)
    param = np.asarray(param, dtype=float)
    elems = np.asarray(elems, dtype=int)

    # fix orientation per element (positive Jacobian)
    for e in range(len(elems)):
        p = nodes[elems[e]]
        if _hex_volume(p) < 0:
            elems[e] = elems[e][[4, 5, 6, 7, 0, 1, 2, 3]]

    # per-element centroid fiber frame
    fiber = np.zeros((len(elems), 3))
    sheet = np.zeros((len(elems), 3))
    normal = np.zeros((len(elems), 3))
    for e in range(len(elems)):
        pr = param[elems[e]]
        tc = np.nanmean(pr[:, 0])
        muc = np.nanmean(pr[:, 1])
        # circumferential mean via complex average (handle wrap + apex nan)
        th = pr[:, 2]
        th = th[~np.isnan(th)]
        thc = np.angle(np.mean(np.exp(1j * th))) if len(th) else 0.0
        try:
            f, s, n = fiber_frame(tc, muc, thc, asign=asign)
        except ValueError:
            f, s, n = np.array([1, 0, 0.]), np.array([0, 1, 0.]), np.array([0, 0, 1.])
        fiber[e], sheet[e], normal[e] = f, s, n

    # per-node fiber/sheet (for per-Gauss-point interpolation, matching the
    # reference's CG1 fiber field).  Apex nodes (degenerate e_θ) fall back to a
    # near-apex evaluation.
    fiber_node = np.zeros((len(nodes), 3))
    sheet_node = np.zeros((len(nodes), 3))
    for nd in range(len(nodes)):
        t, mu, th = param[nd]
        try:
            if np.isnan(th):
                raise ValueError
            f, s, _ = fiber_frame(t, mu, th, asign=asign)
        except ValueError:
            f, s, _ = fiber_frame(t, -np.pi + 0.05, 0.0, asign=asign)
        fiber_node[nd], sheet_node[nd] = f, s

    # boundary facets (node-id quads)
    f_endo, f_epi, f_base, f_endo_elem = [], [], [], []
    for imu in range(n_mu):
        for ith in range(n_theta):
            j, jp = ith, (ith + 1) % n_theta
            # endo (it=0) and epi (it=n_t) outward quads on transmural surfaces
            f_endo.append([node(0, imu, j), node(0, imu, jp),
                           node(0, imu + 1, jp), node(0, imu + 1, j)])
            f_endo_elem.append(imu * n_theta + ith)   # parent elem (it=0, imu, ith)
            f_epi.append([node(n_t, imu, j), node(n_t, imu, jp),
                          node(n_t, imu + 1, jp), node(n_t, imu + 1, j)])
    for it in range(n_t):
        for ith in range(n_theta):
            j, jp = ith, (ith + 1) % n_theta
            f_base.append([node(it, 0, j), node(it, 0, jp),
                           node(it + 1, 0, jp), node(it + 1, 0, j)])

    return CardiacMesh(
        nodes=nodes, elems=elems, param=param,
        fiber=fiber, sheet=sheet, normal=normal,
        fiber_node=fiber_node, sheet_node=sheet_node,
        facets_endo=np.asarray(f_endo, int),
        facets_endo_elem=np.asarray(f_endo_elem, int),
        facets_epi=np.asarray(f_epi, int),
        facets_base=np.asarray(f_base, int),
        n_t=n_t, n_mu=n_mu, n_theta=n_theta,
        topology="polar_ring", n_side=0,
        n_core=0, n_radial=0, core_half_width=0.0,
    )


def _cluster_meridian_to_apex(rho, strength):
    """Remap the meridian coordinate ``rho`` to cluster elements at the apex.

    ``g(rho) = rho * (a + (1 - a) * rho)`` with ``a = 1 / strength`` keeps
    the apex (0) and the base rim (1) fixed, stays monotone with slope in
    ``[1/strength, 2 - 1/strength]``, and scales the meridian element size
    adjacent to the apex by ``1 / strength``.  ``strength = 1`` is the
    identity.  Only the physical placement of the unchanged disk lattice
    moves; connectivity, boundary labels, and the base rim are untouched.
    """
    a = 1.0 / strength
    return rho * (a + (1.0 - a) * rho)


def build_closed_mesh(
    n_t: int = 2,
    n_core: int = 20,
    n_radial: int = 17,
    *,
    core_half_width: float = 0.36,
    flip_helix: bool = True,
    tip_refine: float = 1.0,
) -> CardiacMesh:
    """Build a closed, noncollapsed five-block Hex8 LV-wall mesh.

    The continuous domain and its three boundary classes match Benchmark 1:
    endocardium, epicardium, and base.  There is no fourth/free apex surface.
    A regular central square contains the apex. Four surrounding structured
    blocks map its sides to quarter-circle base arcs. This avoids both a
    collapsed polar ring and the near-singular corner cells of a single square
    mapped directly onto a smooth disk.

    ``tip_refine`` changes only the local element size near the apex: the
    meridian coordinate is remapped so the physical spacing adjacent to the
    apex is scaled by ``1 / tip_refine``, growing smoothly to at most
    ``2 - 1 / tip_refine`` at the base rim.  The disk lattice, element
    counts, connectivity, boundary labels, wall ruling, and the base rim
    are unchanged; ``tip_refine = 1`` reproduces the uniform benchmark mesh
    bit-for-bit.
    """
    if n_t < 1 or n_core < 4 or n_core % 2 or n_radial < 1:
        raise ValueError(
            "n_t >= 1, even n_core >= 4, and n_radial >= 1 are required"
        )
    if not 0.1 < core_half_width < 1.0 / np.sqrt(2.0):
        raise ValueError("core_half_width must lie in (0.1, 1/sqrt(2))")
    if not np.isfinite(tip_refine) or not 1.0 <= tip_refine <= 8.0:
        raise ValueError("tip_refine must lie in [1, 8]")
    asign = -1.0 if flip_helix else 1.0
    disk = []
    disk_lookup = {}

    def add_disk(u, v):
        # Round-key dedup with a tolerance fallback: adjacent blocks compute
        # shared seam columns through different floating-point expressions
        # whose results can differ in the last bit, so the bare rounded key
        # can split a shared point in two.  The 3x3 neighborhood of the
        # rounded key is searched and true coordinates are compared within
        # an absolute tolerance far below the smallest legitimate spacing.
        # Meshes whose keys never split keep the first writer's bytes, so
        # existing resolutions are reproduced bit-for-bit.
        ku = round(float(u), 14)
        kv = round(float(v), 14)
        for du in (-1.0e-14, 0.0, 1.0e-14):
            for dv in (-1.0e-14, 0.0, 1.0e-14):
                existing = disk_lookup.get((ku + du, kv + dv))
                if existing is not None:
                    eu, ev = disk[existing]
                    if abs(eu - u) < 1.0e-12 and abs(ev - v) < 1.0e-12:
                        return existing
        key = (ku, kv)
        disk_lookup[key] = len(disk)
        disk.append([float(u), float(v)])
        return disk_lookup[key]

    surface_quads = []
    core = np.linspace(-core_half_width, core_half_width, n_core + 1)
    core_ids = np.empty((n_core + 1, n_core + 1), dtype=int)
    for j, v in enumerate(core):
        for i, u in enumerate(core):
            core_ids[j, i] = add_disk(u, v)
    for j in range(n_core):
        for i in range(n_core):
            surface_quads.append(
                [core_ids[j, i], core_ids[j, i + 1],
                 core_ids[j + 1, i + 1], core_ids[j + 1, i]]
            )

    tangential = np.linspace(-1.0, 1.0, n_core + 1)
    radial = np.linspace(0.0, 1.0, n_radial + 1)

    def inner_point(block, s):
        h = core_half_width
        return (
            (h, h * s) if block == 0
            else (-h * s, h) if block == 1
            else (-h, -h * s) if block == 2
            else (h * s, -h)
        )

    for block in range(4):
        start_angle = -0.25 * np.pi + block * 0.5 * np.pi
        block_ids = np.empty((n_radial + 1, n_core + 1), dtype=int)
        for j, q in enumerate(radial):
            for i, s in enumerate(tangential):
                inner = np.asarray(inner_point(block, s))
                theta = start_angle + 0.25 * np.pi * (s + 1.0)
                outer = np.asarray([np.cos(theta), np.sin(theta)])
                u, v = (1.0 - q) * inner + q * outer
                block_ids[j, i] = add_disk(u, v)
        for j in range(n_radial):
            for i in range(n_core):
                surface_quads.append(
                    [block_ids[j, i], block_ids[j + 1, i],
                     block_ids[j + 1, i + 1], block_ids[j, i + 1]]
                )

    disk = np.asarray(disk, dtype=float)
    surface_quads = np.asarray(surface_quads, dtype=int)
    edge_owners = {}
    for quad in surface_quads:
        for a, b in zip(quad, np.roll(quad, -1)):
            edge = tuple(sorted((int(a), int(b))))
            edge_owners[edge] = edge_owners.get(edge, 0) + 1
    boundary_edges = np.asarray(
        [edge for edge, count in edge_owners.items() if count == 1], dtype=int
    )
    if len(boundary_edges) != 4 * n_core:
        raise RuntimeError(
            "closed-mesh disk rim has an invalid boundary edge count "
            f"{len(boundary_edges)} (expected {4 * n_core}): the five-block "
            "disk seams are not conforming"
        )
    surface_node_count = len(disk)

    def node(it, surface_node):
        return it * surface_node_count + surface_node

    nodes = []
    param = []
    for it in range(n_t + 1):
        t = it / n_t
        for u, v in disk:
            rho = min(1.0, float(np.hypot(u, v)))
            if tip_refine != 1.0:
                rho = min(
                    1.0, float(_cluster_meridian_to_apex(rho, tip_refine))
                )
            theta = 0.0 if rho < 1.0e-14 else float(np.arctan2(v, u))
            mu = -np.pi + rho * (np.pi + mu_base(t))
            nodes.append(closed_wall_point(t, rho, theta))
            param.append([t, mu, np.nan if rho < 1.0e-14 else theta])

    nodes = np.asarray(nodes, dtype=float)
    param = np.asarray(param, dtype=float)
    elems = []
    for it in range(n_t):
        for quad in surface_quads:
            inner = [node(it, surface_node) for surface_node in quad]
            outer = [node(it + 1, surface_node) for surface_node in quad]
            elems.append(outer + inner)
    elems = np.asarray(elems, dtype=int)

    fiber = np.zeros((len(elems), 3))
    sheet = np.zeros((len(elems), 3))
    normal = np.zeros((len(elems), 3))
    fiber_node = np.zeros((len(nodes), 3))
    sheet_node = np.zeros((len(nodes), 3))

    f_endo = [
        [node(0, surface_node) for surface_node in quad]
        for quad in surface_quads
    ]
    f_epi = [
        [node(n_t, surface_node) for surface_node in quad]
        for quad in surface_quads
    ]
    f_endo_elem = np.arange(len(surface_quads), dtype=int)
    f_base = []
    for it in range(n_t):
        for a, b in boundary_edges:
            f_base.append(
                [node(it, a), node(it, b),
                 node(it + 1, b), node(it + 1, a)]
            )

    mesh = CardiacMesh(
        nodes=nodes,
        elems=elems,
        param=param,
        fiber=fiber,
        sheet=sheet,
        normal=normal,
        fiber_node=fiber_node,
        sheet_node=sheet_node,
        facets_endo=np.asarray(f_endo, dtype=int),
        facets_endo_elem=np.asarray(f_endo_elem, dtype=int),
        facets_epi=np.asarray(f_epi, dtype=int),
        facets_base=np.asarray(f_base, dtype=int),
        n_t=n_t,
        n_mu=0,
        n_theta=0,
        topology="closed_multiblock_disk",
        n_side=0,
        n_core=n_core,
        n_radial=n_radial,
        core_half_width=float(core_half_width),
        tip_refine=float(tip_refine),
    )
    update_structural_frames(mesh, param[:, 0], asign=asign)
    return mesh


if __name__ == "__main__":
    m = build_mesh(n_t=3, n_mu=16, n_theta=24, apex_offset=0.2)
    print(f"nodes={m.n_node}  elems={m.n_elem}")
    print(f"facets endo/epi/base = "
          f"{len(m.facets_endo)}/{len(m.facets_epi)}/{len(m.facets_base)}")

    # --- verification ----------------------------------------------------
    vols = np.array([_hex_volume(m.nodes[e]) for e in m.elems])
    print(f"element volumes: min={vols.min():.3e} max={vols.max():.3e}  "
          f"(all positive: {bool((vols > 0).all())})")
    total_vol = vols.sum()
    print(f"total wall volume = {total_vol*1e6:.2f} cm^3")

    # fiber frame orthonormality
    on_err = 0.0
    for e in range(m.n_elem):
        f, s, n = m.fiber[e], m.sheet[e], m.normal[e]
        on_err = max(on_err, abs(np.dot(f, s)), abs(np.dot(f, n)),
                     abs(np.dot(s, n)),
                     abs(np.linalg.norm(f) - 1), abs(np.linalg.norm(s) - 1),
                     abs(np.linalg.norm(n) - 1))
    print(f"max fiber-frame orthonormality error = {on_err:.2e}")

    # transmural helix-angle sweep at a fixed (mu,theta): expect -60..+60 deg
    print("helix angle vs t̃ (deg):")
    for t in (0.0, 0.5, 1.0):
        f, s, n = fiber_frame(t, -1.0, 0.3)
        # angle of fiber w.r.t. circumferential direction e_theta in tangent plane
        print(f"  t̃={t:.1f}: alpha={np.rad2deg(alpha(t)):+.1f}")

    # endo/epi radii sanity at equator-ish point
    pe = point(0.0, -np.pi / 2, 0.0)
    pp = point(1.0, -np.pi / 2, 0.0)
    print(f"endo |y| at mu=-pi/2: {abs(pe[1]):.4f} (expect {RS_ENDO})")
    print(f"epi  |y| at mu=-pi/2: {abs(pp[1]):.4f} (expect {RS_EPI})")

    if not (vols > 0).all():
        raise RuntimeError("non-positive element volume detected")
    if on_err >= 1e-10:
        raise RuntimeError(f"fiber frame is not orthonormal: error={on_err:.3e}")
    print("\nGEOMETRY SELF-CHECK PASSED")
