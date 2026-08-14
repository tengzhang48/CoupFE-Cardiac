# Case B mesh error decomposition and resolution plan

**Date:** 2026-08-06; **updated:** 2026-08-07
**Branch:** `case-b-tip-refine-hex8-20260805`
**Status:** analysis, plan, and completed results through the full-cycle
`2x20x17`/`4x20x17 tip 6.0` pair and its clean isolated replay;
`6x64x72` deferred

## Answer first

The controls identify **at least two contributors** to the mesh sensitivity:

1. **Boundary-operator sensitivity.** The discrete epicardial normal changes
   the long-axis rotational restraint and materially changes the trajectory.
2. **Bulk spatial resolution.** Q1 Hex8 interpolation through a curved,
   snapping wall remains sensitive to surface and through-wall resolution.

These controls do not prove that the two contributors are exhaustive, do not
uniquely allocate the remaining error, and do not establish distinct solution
branches.

**Working interpretation:** the observed trajectory split is numerical
sensitivity near snap-through. There is no evidence here for physically
distinct branches.

Through-wall resolution: **two Q1 layers is a recognized risk for a
bending/snap-dominated wall** and four layers belong in a future convergence
ladder. The completed pair is nearly identical before snap, separates in the
snap window, and nearly reunites late in the cycle. That phase-dependent
effect is numerical-sensitivity evidence, not a convergence estimate.

Nothing in the current evidence indicts the Hex8 local-pressure formulation
itself; no element switch is required by the evidence so far.

## Evidence distinguishing the measured contributors

### Robin-normal 2x2 gate (produced on the `p2p1-hex20-hex8` line at `59324a6`)

Current-frame 0.32 s trajectories, `{facet, ellipsoid-gradient} x
{2x20x17, 2x36x32}`. The raw gate metrics are retained off-repository; the
values below are the record:

| Comparison | `p0` endpoint | `p0` max | `p0` RMS | `p1` endpoint | `p1` max | `p1` RMS |
|---|---:|---:|---:|---:|---:|---:|
| facet, coarse-minus-surface (mm) | 9.962 | 14.428 | 9.651 | 8.742 | 13.365 | 8.987 |
| ellipsoid-gradient, coarse-minus-surface (mm) | 3.695 | 3.695 | 1.325 | 3.383 | 3.383 | 1.219 |
| ratio EG/facet | 0.371 | 0.256 | 0.137 | 0.387 | 0.253 | 0.136 |

Snap timing (first downward `u_z = -5 mm` crossing): the facet pair is
separated by `0.0129 s` (`p0`) and `0.0130 s` (`p1`); the
ellipsoid-gradient pair by `0.0008 s` at both points. Against the
decision rules pre-registered before the runs (all six
ratios `<= 0.25` and `D_snap <= 0.002 s` for "strongly causal"; endpoint and
RMSE ratios `>= 0.75` for "insufficient"), the outcome is **mixed**: the
normal representation is associated with most of the measured split
(including the snap-timing shift) but not all of it. The remainder cannot be
uniquely assigned by this gate.

This gate also **confirms the current-frame facet split**: 9.962/8.742 mm,
matching the pre-frame 9.994/8.772 mm measurement; the physical-frame
correction does not explain the split.

### Six-rigid-mode reduction of the assembled Robin matrices (independent review)

Reference-configuration reduction of the exact assembled `K` and `C` onto
the six rigid modes (unit-norm columns, symmetrized, rotations about the
ellipsoid center), reproduced from the released sources:

| Surface | `rx` (long-axis) total | base part | epi faceting part | `ry` / `rz` total |
|---|---:|---:|---:|---:|
| CoupFE `2x20x17` | 1.4446 | 0.1595 | 1.2851 | 1382.34 / 1382.34 |
| CoupFE `2x36x32` | 0.5558 | 0.1597 | 0.3962 | 1383.42 / 1383.42 |
| local FEniCS P1 mesh | 1.8570 | 0.1590 | 1.6980 | 1382.22 / 1382.23 |

Units N m/rad. The reduced 6x6 matrices are rank 6/6 and PSD for both `K`
and `C` on all three tested meshes: the compliant support is complete. Within
these reference-configuration reductions, the measured cross-mesh variation
is concentrated in the long-axis rotation mode (`ry`/`rz` agree to <=0.09%).
On every tested mesh the discrete epicardial contribution exceeds the physical
base restraint for that mode
(2.5x on the surface-only mesh, 8.1x on the coarse mesh, 10.7x on the
FEniCS mesh). The retained 0.32 s states carry about -32 degrees of
epicardial long-axis twist, so the trajectories exercise a mode with a large,
surface-dependent discrete epicardial contribution.

### Four-way mesh split (pre-frame, corrected geometry)

| Change | `p0` endpoint (mm) | `p1` endpoint (mm) |
|---|---:|---:|
| through-wall `2x20x17 -> 4x20x17` | 0.669 | 0.635 |
| through-wall at fine surface `2x36x32 -> 4x36x32` | 0.881 | 1.217 |
| surface `2x20x17 -> 2x36x32` | 9.994 | 8.772 |

## Why refinement looked non-monotone

Naive convergence reasoning says the finer mesh should move toward the true
solution. Instead the surface-refined facet trajectory moved *away* from the
local FEniCS result. Surface refinement simultaneously changes Q1 resolution
and the discrete boundary operator: the measured epicardial long-axis
restraint changes from 1.285 to 0.396 N m/rad on the compared `2x20x17` and
`2x36x32` surfaces, while the local FEniCS mesh measures 1.698 N m/rad. These
are different discrete problems, so their errors need not order monotonically.
FEniCS itself integrates with `FacetNormal` on planar P1 triangles, so it
is a finite-mesh comparison with its own surface-discretization sensitivity,
not a continuum
oracle; its producing revision was not retained.

## Through-wall assessment: two layers of a linear element

A Q1 Hex8 has linear displacement through the wall thickness, so the
through-thickness strain gradient of a bending wall is piecewise-constant
per element layer. The condensed Q1/P0 local pressure addresses volumetric
locking; it does not add strain-gradient resolution. For bending-dominated
response the standard minimum is about four linear layers (or quadratic
elements such as the Hex20 line under investigation on the `p2p1` branch).

Measured at the 0.32 s prefix: `<= 0.67 mm` endpoint effect at the coarse
surface and `<= 1.22 mm` at the fine surface - smaller than the measured
surface/normal effects at this horizon. Snap-through is bending-dominated
and the largest curvatures develop later in the cycle, so this prefix bound
must not be treated as a full-cycle clearance.

**Suggestion:** keep 2 layers for the current 0.32 s diagnostic gates
(comparability with the retained clean archives), and put 4 through-wall
layers on the declared convergence ladder below. Do not run a generic
through-wall sweep; the two through-wall controls above already exist.

## The tip-refine control (this branch) — results

`--tip-refine F` (default `1.0`, validated to `[1, 8]`) remaps the meridian
coordinate `rho -> rho(a + (1-a) rho)` with `a = 1/F`, clustering physical
elements toward the apex - the highest-curvature region of the wall. This is
**global meridional node grading (r-adaptation), not local h-refinement**: it
keeps element counts, connectivity, boundary labels, wall ruling, and the base
rim fixed, but moves almost every non-base meridional node and coarsens cells
toward the base. The default reproduces the uniform benchmark mesh
bit-for-bit. On `2x20x17` at `F=2.5`,
apex-adjacent epicardial meridian edges shrink `2.31 -> 0.97 mm` while
base-adjacent edges grow by at most `1.59x`; the graded mesh passes the
production pre-solve geometry audit (scaled Jacobian 0.239 vs 0.258
uniform, condition 8.2 vs 7.9, reference measures within 0.06% vs the 1%
gate). The value is recorded in result archives and Laplace-field sidecars;
stale fields are rejected fail-closed.

### Tip-refinement trend (2x20x17, facet operator, current frame, 0.32 s)

Endpoint landmark gap against the local FEniCS result, all components:

| Run | `|gap|` p0 (mm) | `|gap|` p1 (mm) | `Delta x` p0 (mm) | `Delta y/z` p0 (mm) | RMSE 0-0.32 p0 (mm) | snap p0 (s) |
|---|---:|---:|---:|---:|---:|---:|
| uniform `1.0` | 6.154 | 6.279 | -6.09 | +0.17 / +0.88 | 2.593 | 0.2487 |
| `tip 2.5` | 3.608 | 3.923 | -3.53 | -0.06 / +0.77 | 1.535 | 0.2450 |
| `tip 4.0` | 2.587 | 2.945 | -2.48 | -0.14 / +0.73 | 1.421 | 0.2456 |
| `tip 6.0` | 1.998 | 2.366 | -1.86 | -0.18 / +0.71 | 1.424 | 0.2461 |

Reading:

- **Monotone endpoint sensitivity, not a convergence sequence.** The endpoint
  gap shrinks 6.15 -> 3.61 -> 2.59 -> 2.00 mm (p0) as the tested grading
  strengthens. The full-history RMSE instead levels at 1.421/1.424 mm for
  `F=4/6`, and the maximum history error grows 4.19 -> 4.93 -> 5.55 mm for
  `F=2.5/4/6`. Because this is fixed-count global grading, the data do not
  establish an error order, apex ownership, or convergence. At `F=6.0`, x is
  still ~0.8 mm outside the ten-team band and the z residual remains ~0.7 mm.
- **The 0.32 s endpoint error is ~90% an x (long-axis) overshoot** in every
  run; the y/z endpoint residuals remain sub-millimetre in this particular
  grading series. This does not establish full-history or full-cycle
  insensitivity of those components.
- **Snap timing clusters after `F=2.5`**: 0.2487 -> 0.2450 -> 0.2456 ->
  0.2461 s. The separate normal-representation gate motivates another
  contributor, but the available controls do not causally separate it from
  bulk discretization.
- **Direction check**: uniform surface refinement moves the endpoint
  *away* from FEniCS (-22.2 -> -32.1 x); tip refinement moves it *toward*
  (-22.2 -> -19.6 -> -18.5). The opposite trends are consistent with multiple
  coupled numerical sensitivities. These controls do not provide a causal
  decomposition of the historical non-monotone behavior.

### Ten-team benchmark envelope at 0.32 s

Published participant curves (Zenodo dataset, canonical 0.01 s grid, ten
teams; duplicate SimVascular alias excluded per the gate policy):

| Landmark | x band (mm) | y band (mm) | z band (mm) |
|---|---|---|---|
| p0 | [-17.15, -15.37] (median -16.33) | [+1.11, +1.29] | [-16.40, -15.09] |
| p1 | [-18.03, -15.80] (median -17.23) | [+0.32, +0.59] | [-17.06, -15.62] |

The local FEniCS endpoint sits essentially inside (p1 fully; p0 within
0.02 mm). CoupFE endpoints at 0.32 s: uniform coarse is 5.0 mm outside the
x band, `tip 2.5` is 2.4 mm outside, `tip 4.0` is ~1.4 mm outside, and the
surface-refined `2x36x32` facet run is 15 mm outside. The team envelope is
a comparison frame with a wide spread of participant methods, not a
converged answer; it corroborates the direction of the tip effect.

### Provenance of the new trajectories

All facet-operator, current-frame, eight MPI ranks, 320/320 accepted
increments, clean source trees, produced on this branch:

- `tip 2.5` `2x20x17`: NPZ SHA-256
  `3b47ad2cc6d5a6c5ff8ee0dff50e027ed243ea83df7e9fcf7a473dc2b661ca8b`
  (app `2f5c633`);
- `tip 4.0` `2x20x17`: NPZ SHA-256
  `2144bc6e67794181cafb1c42d687d18a35da4767825fc258d2f199f83b56bcf9`
  (app `c965103`);
- `tip 6.0` `2x20x17`: NPZ SHA-256
  `b3937a41f53481944fa4dd70a6722479944770ada7a581b90855fd487c2b67d8`
  (app `fbccf35`);
- `4x32x48` uniform: NPZ SHA-256
  `9e99eeb9b1d1e1441f3b5401077db9479e4aebf1b0399cc95e0bdf1fdac1c6c4`;
- `4x32x48` `tip 2.5`: NPZ SHA-256
  `0db2423b8c2f428f0807ac0dfee081031f619fd1fc9a91ba433cf1c795270de8`
  (both app `fbccf35`).
- `4x20x17` `tip 6.0`: NPZ SHA-256
  `774a7dc5fc970bb744ff0188f0f428cff54c532b401a365643f4b626584d7acf`
  (app `a5824d3`, Core `454f73c`; 320/320 increments, eight ranks).

### Tip effect at four wall layers and 1.9 mm surface (4x32x48 pair)

A one-variable pair at `4x32x48` (28,672 elements, 1.91 mm mean
epicardial edge, ~2.3 mm wall layers), facet operator, current frame:

| Run | p0 x at 0.32 s (mm) | `|gap|` vs FEniCS p0 (mm) | snap p0 (s) |
|---|---:|---:|---:|
| `4x32x48` uniform | -32.09 | 16.04 | 0.2368 |
| `4x32x48` `tip 2.5` | -31.87 | 15.82 | 0.2352 |

The endpoint difference between this pair is 0.23 mm at p0 (0.095 mm at
p1), but the maximum transient difference is 2.39/2.10 mm. Against FEniCS,
the p0 full-prefix RMSE changes 5.638 -> 5.866 mm and the p1 RMSE changes
5.108 -> 5.348 mm. Thus the two runs have similar endpoints but not
identical histories. One pair does not establish a vanishing grading effect,
an apex-converged mesh, or a distinct solution branch.

### Four-wall-layer graded control (`4x20x17 tip 6.0`, full cycle)

This predeclared control changes only the wall-layer count relative to the
retained `2x20x17 tip 6.0` run. It completed 1,000/1,000 increments on eight
ranks. The geometry, pressure, Robin, Laplace-field, solver, consistent-mass,
generalized-alpha, Q1/P0 mean-`log(J)` local-pressure, `eta=100 Pa s`, and
physical-landmark audits passed. Its independently completed 0.32 s prefix
agrees with the full-run prefix to roundoff.

| Pair metric, four minus two layers | p0 | p1 |
|---|---:|---:|
| pre-snap history RMSE (mm) | 0.009831 | 0.002622 |
| maximum transient difference vs two-layer graded run (mm) | 1.633 at 0.250 s | 1.406 at 0.247 s |
| full-cycle history RMSE (mm) | 0.407604 | 0.361832 |
| late 0.75--0.999 s history RMSE (mm) | 0.051556 | 0.062702 |
| cycle-end separation (mm) | 0.034598 | 0.051776 |
| normalized late / endpoint ratio | 0.031570 / 0.021186 | 0.044600 / 0.036829 |
| FEniCS full-shared-history RMSE, two -> four layers (mm) | 1.092675 -> 1.016944 | 1.168502 -> 1.072951 |
| maximum FEniCS gap, two -> four layers (mm) | 5.552 -> 4.188 | 5.024 -> 3.753 |

The pair is essentially coincident before snap, separates transiently during
snap, and nearly reunites late in the cycle. Downward/upward `u_z=-5 mm`
crossings are 0.242327/0.561007 s and 0.239063/0.558177 s for FEniCS,
0.246148/0.567678 s and 0.242905/0.564664 s for two layers, and
0.244986/0.568785 s and 0.241707/0.565612 s for four layers (`p0`, then
`p1`). Four layers improve the downward timing and full-shared-history FEniCS
RMSE,
but relaxation RMSE worsens by 16.2%/14.9%; ten-team containment fractions
are also mixed. This supports numerical transient/timing sensitivity, not
physically distinct branches or mesh convergence.

At 0.32 s the new p0 displacement is
`(-18.093, +0.904, -15.904) mm` and p1 is
`(-19.052, +0.411, -16.629) mm`. Relative to the ten-team envelope, p0 z and
p1 y/z are inside while both x components and p0 y remain outside. The axial
endpoint therefore does not by itself improve monotonically with every mesh
change even though the full-prefix transient comparison improves.

Prefix provenance: stdout SHA-256
`0b30a93d01f0e5419f8b5acace8c17931355cfab6bd486e114aac2c2045d5e5a`;
Laplace NPY SHA-256
`1578362593495b6fe48d6a2fd2e1332150121be4d6b361915d04f3980d78da8f`;
Laplace metadata SHA-256
`ecdadd335e41922eab459f4d0d6a17cf7fd4a3add496ddc0b179ca9f25daceeb`.
The first full-cycle candidate is an independently audited but non-retained
dirty-tree provenance diagnostic. A clean isolated replay completed under Core
`454f73c` and runtime-source SHA-256
`f8d9469a101709d11460a0803b6c031001192ccce61c472d032355b01070da05`;
it completed with states matching to roundoff. Retained provenance is
**application `2458e7c`, NPZ
`1e333b29b05f01dedce9272b32b82ce6ccfda56036c1c8f57eb395b9b4494800`,
stdout `0810a9d0c944345464f48917c707a7d3df7187b543c74efd4b0384542b95aed7`,
elapsed `1778.4 s (29.6 min)`**.

### Current trajectory groups at 0.32 s (p0 x, mm)

| Configuration | x (mm) | Observed group |
|---|---:|---|
| ten-team published band | [-17.15, -15.37] | published range |
| local FEniCS | -16.06 | within published range |
| `2x20x17` facet + tip `1.0 / 2.5 / 4.0 / 6.0` | -22.15 / -19.59 / -18.54 / -17.92 | approaches published range |
| `4x20x17` facet + tip `6.0` | -18.09 | close to two-layer graded result |
| fine facet (`2x36x32`, `4x32x48`, uniform or tip) | -32.09 | more-negative response |
| smooth-normal (ellipsoid-gradient) coarse / surface (Robin-normal gate) | -33.71 / -37.40 | more-negative response |

The boundary-normal operator is the single largest lever measured (coarse
facet -22.15 vs coarse ellipsoid-gradient -33.71: 11.6 mm of endpoint x
from the operator alone), and the smooth-normal mode moves CoupFE *away*
from the ten-team consensus, not toward it. Faceting explains much of the
coarse-vs-fine *split* (the gate's pre-registered ratios), but neither
normal choice explains the CoupFE-vs-consensus *gap*. The remaining open
question is why fine/smooth CoupFE meshes develop a more-negative post-snap x
response while the coarse graded run, local FEniCS, and the ten published
teams cluster closer together. These are descriptive trajectory groups, not
evidence for distinct mathematical solution branches. The `6x64x72` uniform
run is **deferred**:
another large uniform run would not isolate the through-wall or grading
contribution.

## Snap-window amplification (observation; cause unresolved)

The measured landmark histories remain close through the pre-snap loading and
their separation grows rapidly during the instability: `4x32x48` matches the
clean coarse trajectory to 0.16 mm at t=0.20 s and
first separates by 0.5/1/2 mm at t=0.226/0.229/0.232 s - exactly the snap
window (snap at 0.236-0.249 s across runs). The `tip 6.0` p0 history first
separates from the same reference by 0.5 mm at 0.237 s. This localizes the
visible amplification in time, but close landmark values do not prove
identical fields or physics. Snap-through can amplify small discretization
and operator differences; no continuation or stability evidence establishes
distinct solution branches here.

One candidate contributor is the epicardial faceting representation's
discrete restraint of the long-axis rotation (the smooth ellipsoid value is
zero; the
physical restraint is the base spring 0.159 N m/rad plus elastic
coupling). The measured association is suggestive but confounded:

| Configuration | epi faceting twist restraint (N m/rad) | snap (s) | p0 x at 0.32 s (mm) |
|---|---:|---:|---:|
| local FEniCS | 1.698 | -- | -16.06 (near) |
| coarse `2x20x17` | 1.285 | 0.2487 | -22.15 (near) |
| fine `4x32x48` | 0.492 | 0.2368 | -32.09 (far) |
| ellipsoid-gradient (GPT gate) | ~0 | 0.223 | -33.7 / -37.4 (farthest) |

Lower measured twist restraint co-occurs with earlier snap and a more-negative
p0-x endpoint in these configurations. This is not a causal proof: confounders
include facet
pattern (quad vs triangle), within-facet normal variation (CoupFE's
bilinear quads are non-planar and carry four Gauss-point normals per
facet; FEniCS P1 triangles are planar with one), Q1-Hex8 vs P2-tet
kinematics, and any pre-snap stiffness differences the restraint itself
introduces. The available controls establish association, not causality.

**Load and constraint areas (factual).** The follower pressure loads the
endocardium per deformed area (`-p J F^-T N`, Nanson; resultant and
moment gated against the analytic closed-cavity projection at reference).
The Robin springs/dashpots constrain per unit *reference* area
(`alpha` Pa/m, `beta` Pa s/m at facet Gauss points) on the epicardium
(normal-only) and base (full-vector). Both are distributed per-area
operators using the same facet geometry; neither is a point constraint.

**Discriminating tests (no law changes).** (a) Full-cycle comparison: does the
axial agreement of `2x20x17 tip 6.0` at the 0.32 s endpoint continue through
the peak and unloading phases against the ten-team full-cycle envelope -
**done**, see the next section (mostly, but not continuously, for x).
(b) Manufactured twist loading on retained states: apply a prescribed
long-axis twist field to each archived state and compare the total
(spring + elastic + dashpot) generalized moment - separates the Robin
restraint from the elastic response without new trajectories.
(c) Pre-snap stiffness tracking: monitor the twist-mode generalized
stiffness along each accepted trajectory through 0.20-0.25 s and check
whether the more-negative trajectories soften earlier.

### Full-cycle comparison (1 s, `2x20x17`/`4x20x17 tip 6.0`, facet operator)

The two-layer full-cycle run on this branch (1000/1000 increments, 8 ranks,
clean tree; NPZ SHA-256
`5bb152c47b693af1dc2c0d650dde8b07ba28ef51d44e45db8c493cfe9a339375`).
Declared question: does the favorable 0.32 s axial comparison of the
coarse+graded configuration continue through the peak and relaxation phases?

**The x component is inside the team envelope for most, but not all, of the
peak and relaxation; the three components are rarely inside simultaneously.**
Same-time comparisons below are computed on the corrected
time alignment: CoupFE retains states at t = 0.000...1.000 while the
FEniCS record spans t = 0.001...0.999, so FEniCS comparisons use matched
timestamps and the cycle-end gap is the last shared time (t = 0.999 s),
not a literal t = 1.000 s comparison. An earlier version of this section
reported scalars from a one-step-misaligned comparison; the corrected
numbers are slightly better.

Axial (x) component position relative to the ten-team full-cycle
envelope (per-component interpolation of the published curves):

| Phase | ours p0 x (mm) | team band x (mm) | x inside |
|---|---:|---:|---|
| snap window t=0.25 | +10.05 | [+0.36, +7.67] | no (pre-snap rise overshoots) |
| t=0.32 | -17.92 | [-17.15, -15.37] | no (~0.8 mm out) |
| peak t=0.44-0.48 | -22.26 / -22.89 | [-29.13, -19.28] | yes |
| t=0.60 | +2.42 | [+1.67, +2.78] | yes |
| t=0.80 | -0.28 | [-0.50, -0.22] | yes |
| end t=1.00 | -0.17 | [-0.37, -0.13] | yes |

Component-specific inside-envelope fractions over the full cycle:

| Component | p0 | p1 | Where it leaves the envelope |
|---|---:|---:|---|
| x (axial) | 86.4% | 84.0% | snap/early systole (0.23-0.34 s) and ~0.54-0.57 s |
| y | 34.5% | 39.5% | outside throughout the peak; repeatedly outside in unloading |
| z | 59.7% | 68.2% | snap window and parts of unloading |

All three components are inside simultaneously for **21.0%** of the full
cycle at p0 and **22.2%** at p1. In the declared peak window
`[0.350, 0.484] s`, that fraction is 0% at p0 (its y component is outside)
and 99.3% at p1. During relaxation `[0.484, 1.000] s`, it is 40.4% at p0
and 2.9% at p1. These phase-specific results are why a single statement that
the trajectory is "inside" is misleading.

Against the retained FEniCS arrays at shared times: last
shared-time difference (t = 0.999 s) **0.031 mm** (p0) and **0.082 mm**
(p1); shared-time vector RMSE 1.093/1.169 mm; maximum separation 5.552 mm
(p0, t = 0.248 s) and 5.024 mm (p1, t = 0.246 s), both in the snap
window; snap timing 0.2461 s, unchanged from the prefix.

Interpretation, with boundaries: after the snap window the two-layer axial component
approaches the FEniCS curve and is usually inside the published team envelope,
but it leaves the envelope again around 0.535-0.574 s. The transverse
components repeatedly leave it. The largest same-time FEniCS differences are
in the snap window; smaller but nonzero differences remain during unloading
(post-0.4 s vector gaps reach about 2.07 mm at p0 and 2.81 mm at p1). This is
one coarse, globally graded, facet-normal configuration and one landmark pair;
it does not establish mesh convergence or explain the more-negative response
of the fine and smooth-normal trajectories.

The source-bound metrics are in
`examples/cardiac_benchmark/results/step0b_tip6p0_full_cycle_comparison.report.json`;
the corresponding six-component figure is
`docs/figures/step0b_tip_refine_full_cycle.svg`. Both are regenerated by the
fail-closed comparison script and carry exact CoupFE, FEniCS, and ten-team
input identities. The report and figure now also carry the four-layer
full-cycle trajectory. Its full-history pairwise RMSE against two layers is
0.407604/0.361832 mm at `p0`/`p1`; its cycle-end separation is only
0.034598/0.051776 mm. Against FEniCS, full-shared-history RMSE improves by
6.9%/8.2% and the maximum snap gap by 24.6%/25.3%, while relaxation RMSE
worsens by 16.2%/14.9%. The improvement is therefore phase-specific rather
than monotone.



## Disk seam dedup defect found by the 6x64x72 gate (fixed)

The first `6x64x72` launch failed the pre-solve geometry audit
(`base_area_cm2` +3.62% vs reference). Root cause: adjacent disk blocks
compute shared seam-corner columns through different floating-point
expressions; at `64x72` their last-bit differences split the bare
`round(u, 14)` dedup key, duplicating six seam points and creating twelve
overlapping full-size base quads. The audit refused the mesh before any
solve - the fail-closed gate working as designed. Fix: tolerance-based
neighborhood fallback in the dedup (first writer's bytes preserved, so
`2x20x17` and `2x36x32` remain bit-identical to the retained archives,
verified array-equal), plus a build-time rim-edge-count invariant and a
64x72 regression test. Post-fix `6x64x72` passes the full audit (base
+0.033%, scaled Jacobian 0.250).

## Global mesh sizes (epicardial surface)

| Mesh | epi edge mean (mm) | wall layers (mm) | elements |
|---|---:|---:|---:|
| `2x20x17` (coarse gate) | 4.00 | ~4.5 | 3,520 |
| `4x32x48` | 1.91 | ~2.3 | 28,672 |
| `2x36x32` (surface) | 2.17 | ~4.5 | 11,808 |
| `6x64x72` (global fine) | 1.09 | ~1.5 | 135,168 |
| local FEniCS (P1 tets, P2 field) | 3.88 (effective ~1.9) | -- | -- |

## Recommended ordered program

1. **Graded-node 0.32 s sensitivity control** (this branch, facet operator,
   current frame): **done** - endpoint gaps decrease for the tested
   `F = 2.5, 4.0, 6.0` values, while history metrics plateau or worsen. Treat
   this as sensitivity, not an apex convergence sequence.
2. **Tip effect at four wall layers and mid surface** (`4x32x48` uniform
   vs `tip 2.5`, one-variable pair): **done** - only 0.23/0.095 mm endpoint
   change at p0/p1, but 2.39/2.10 mm maximum transient change and slightly
   worse FEniCS RMSE. One pair does not close apex or grading sensitivity.
3. **Global fine ladder point** (`6x64x72`, facet operator, current
   frame): **deferred**; another large uniform run does not isolate the
   surface-grading and through-wall contributions.
4. **Four-wall-layer graded control (`4x20x17 tip 6.0`, 0.32 s): done.**
   The exact prefix is retained as the qualification and continuity record.
5. **Full-cycle extension of the four-layer graded control: done.** The pair
   is nearly identical before snap and late in the cycle but separates by up
   to 1.41--1.63 mm during snap. FEniCS full-shared-history and maximum-gap
   metrics improve, while relaxation worsens. Retain the phase-mixed result as
   numerical-sensitivity evidence, not convergence.
6. **Grading x ellipsoid-gradient diagnostic pair**, if needed, to measure
   interaction between node placement and the alternate normal operator.
   It is a diagnostic, not the benchmark replacement: the reference uses
   finite facet normals.
7. **Declared convergence ladder under one preselected operator:** surface
   `{2x20x17, 2x36x32}` x wall `{2, 4}`, 0.32 s, with pre-declared success
   bands on landmark splits and snap timing. Keep facet normals as the
   benchmark-comparable primary operator; use ellipsoid-gradient only as a
   separately labelled diagnostic. This mesh-independence evidence does not
   yet exist for CoupFE.
8. **Element-order decision (Hex20/P2P1) only after the ladder**, on
   evidence; higher order and additional spatial-resolution controls probe
   related bulk-discretization sensitivity and should not run as competing
   explanations in parallel.

## Do not

- Do not rerun the completed `2x20x17 tip 6.0` full cycle merely to reconfirm
  it; its exact retained archive is source- and input-bound.
- Do not tune Robin coefficients or replace the normal-only physical law.
- Do not call ellipsoid-gradient/analytic normals "the benchmark
  correction"; the reference code itself uses finite facet normals.
- Do not extrapolate the 0.32 s through-wall bound to the full cycle.
- Do not treat any retained trajectory (CoupFE coarse/fine or local
  FEniCS) as a converged reference; none exists yet.

## Claim boundaries and provenance

- The 2x2 gate numbers come from the locally retained
  `robin_normal_gate_metrics.json` cited above; its producing application
  revision is `59324a6` on the `p2p1-hex20-hex8` line. The facet-mode
  coarse archive is the reviewed clean gate (`efc7e42a...`).
- The six-rigid-mode numbers are reproducible from the released sources
  alone: build the closed mesh, assemble the Robin spring and dashpot, reduce
  onto the six unit-normalized rigid modes about the ellipsoid centre with
  rotations scaled by `RL_EPI`, and symmetrize before the eigenvalue query.
  The reduction cross-checks the published long-axis audit totals to 1e-12
  before reporting.
- The four-way split numbers come from
  `corrected_mesh_split_prefix_metrics.json` (pre-frame campaign).
- The local FEniCS arrays' producing source revision was not retained; all
  FEniCS comparisons carry that evidence limitation.
- The runs on this branch (tip series, `4x32x48` pair, full cycle) were
  produced with the environment's CoupFE Core at `454f73c`, while the
  clean gate used the approved Core `e2f42ed`; the two revisions are on
  different lineages and the numeric core differs. A dedicated anchor
  rerun of the clean-gate configuration (application `056c02d`, Core
  `454f73c`, identical command) reproduces the retained clean archive to
  1.3e-16 m (roundoff). Thus the Core difference is numerically inert for
  this audited `2x20x17`, 0.32 s anchor trajectory; it is not a blanket
  equivalence claim for every cardiac workload. The comparisons in this
  document use that anchor as their Core-lineage check;
  the revision difference is recorded here for provenance.
- FEniCS comparisons use matched timestamps: the CoupFE archives retain
  t = 0.000...1.000 s and the FEniCS record spans t = 0.001...0.999 s;
  an earlier one-step-misaligned comparison overstated the full-cycle
  RMSE (1.204/1.288 mm) relative to the corrected same-time values
  (1.093/1.169 mm).

## Lessons learned

- **Call the mesh operation what it is.** `tip_refine` moves nodes throughout
  the meridian at fixed element count. It is global r-adaptation, not local
  h-refinement. A favorable endpoint trend cannot be attributed only to the
  apex without a genuinely local control.
- **An endpoint trend is not convergence.** The p0 endpoint gap improves with
  stronger grading, but the prefix RMSE plateaus and its maximum error grows.
  Report endpoint, history, phase, and transient extrema together.
- **Similar curves around a snap do not prove distinct branches.** Snap-through
  amplifies small discretization differences. Branch language requires
  continuation or stability evidence that these runs do not provide.
- **Envelope membership must state its quantifier.** Per-component membership,
  all-components-simultaneous membership, and phase windows answer different
  questions. The full-cycle record now reports all three.
- **Regeneration must fail closed.** The comparison tool pins the exact retained
  CoupFE archive, all four FEniCS inputs, and the selected ten upstream team
  files; it validates timestamp alignment and writes deterministic SVG. This
  prevents a plausible-looking figure from silently mixing campaigns.
- **Preserve physical-law boundaries.** The pressure and Robin operators are
  distributed area laws. Mesh-dependent discrete normals can affect their
  response without implying that their physical coefficients should be tuned.
