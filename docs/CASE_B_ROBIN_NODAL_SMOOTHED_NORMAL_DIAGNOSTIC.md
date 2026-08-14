# Case B — nodal-smoothed epicardial Robin normal diagnostic

This document records an opt-in **mechanism diagnostic** for the Step 0B snap-window
response family. It is not a benchmark-faithful operator and does not change any
retained run.

## Scientific question

Fine and smooth-normal CoupFE meshes select a *far* post-snap response family at
0.32 s (`p0 x ≈ −32 to −37 mm`), while the coarse tip-graded facet mesh, the local
FEniCS model, and all ten published teams sit on the *near* family (`−15 to −19 mm`).
The candidate selector is the epicardial faceting artifact's spurious restraint of
long-axis rotation, which shrinks with surface refinement. The question is whether
the selector is normal **direction** error (fixable by smoothing the normal field)
versus facet **size/pattern/quadrature** (not fixed by smoothing).

## The opt-in mode

`examples/cardiac_benchmark/robin.py` adds a third Robin projection mode,
`normal-smoothed`, alongside the default `facet` (`normal`) and `full` modes:

- epicardial nodal normals are computed by **area-weighted averaging** of the
  adjacent facet normals at each surface node, then normalized;
- at each 2×2 Gauss point the active normal is `N_q = normalize(Σ_a N_a(ξ,η) n_a)`
  with the Q1 shape functions, and the projection is `P = N_q ⊗ N_q`;
- the facet area weights (`detA`), quadrature, coefficients, and assembly structure
  are identical to the default mode, so **the default `facet` mode is bit-identical**
  to the historical operator.

Both drivers expose `--epicardial-normal-mode {facet,nodal-smoothed}` (default
`facet`); the mode is recorded in the result archive and solver configuration.
The operator's built-in rigid-spring diagnostic (`_rigid_spring_reduction`) also
supports the smoothed mode.

## Operator-level evidence (verified, current source)

On the coarse `2×20×17` mesh, the operator's own rigid-spring reduction reports the
epicardial long-axis rotation restraint (N·m/rad, alpha-weighted quadrature):

| mode | epi long-axis rotation restraint |
|---|---|
| `facet` (default) | 23.9381 |
| `nodal-smoothed` | **−0.0000** |

The smoothed normal field collapses the long-axis twist restraint to numerical zero,
while the base translation block and the other rotation components are preserved.
This is the predicted operator-level signature: normal-direction error is the source
of the spurious `rx` restraint. Gates recorded with this change:

- default `facet` mode bit-identical on every non-epicardial DOF;
- `K` symmetric (`max|K−Kᵀ| < 1e-9`) and PSD;
- `rx_smoothed < 0.05·rx_facet`; `ry`/`rz`/translations preserved;
- smoothed mode bit-identical to facet mode on any flat facet.

These are properties of the assembled operator only.

## Not established (trajectory-level verdict pending)

The pre-registered decision rules concern the **trajectory** response — whether the
`nodal-smoothed` runs move toward the *far* family relative to the `facet` pair
(endpoint/RMS ratios against the facet pair, 0.25/0.75 bands). That verdict requires
the controlled `{facet, nodal-smoothed} × {2×20×17, 2×36×32}` experiment through
`t = 0.32 s` (8 MPI ranks, `dt = 0.001 s`, generalized-alpha, consistent mass,
Q1/P0 local pressure, follower pressure, identical everything else) under the pinned
runtime of `docs/CONTROLLED_BENCHMARK_RUNS.md`. **That experiment is not run in this
change.** Until it is, this document states only the operator-level mechanism result;
it does not claim that smoothing selects the response family in the full solve.

This mode is a diagnostic. FEniCS itself uses finite facet normals; smoothed normals
are **not** "the benchmark correction." No Robin coefficients or the normal-only
physical law were changed, and no retained run was relabeled or altered.

## Source

Implementation: `examples/cardiac_benchmark/robin.py` (`_nodal_smoothed_normals`,
`_assemble` and `_rigid_spring_reduction` mode `normal-smoothed`); plumbed through
`run.py` and `run_mpi.py`. Tests: `tests/test_cardiac_fast.py` (bit-identity,
symmetry/PSD, `rx` collapse, CLI validation, unknown-mode rejection).

## Trajectory-level verdict (2026-08-13, experiment now run)

The controlled `{facet, nodal-smoothed} x {2x20x17, 2x36x32}` experiment was
run on 2026-08-13 through `t=0.32 s` (8 MPI ranks, `dt=0.001 s`,
generalized-alpha, consistent mass, Q1/P0 local pressure, follower pressure,
pinned petsc4py 3.18.4 runtime, clean tree at `b28ebc9`). The coarse facet
run reproduces the retained clean gate archive to 5.2e-17 m, so the default
mode is bit-compatible with the clean gate.

Endpoints (`p0` x, mm, at 0.32 s):

| Configuration | facet | nodal-smoothed |
|---|---:|---:|
| coarse `2x20x17` | -22.15 (near family) | -33.78 (far) |
| surface `2x36x32` | -32.10 (far) | -37.46 (far) |

Pre-registered decision-rule ratios (smoothed pair split / facet pair split):
`p0` endpoint 0.371, maximum 0.256, RMS 0.140; `p1` endpoint 0.379, maximum
0.248, RMS 0.135; snap timing within the smoothed pair is 0.9-1.0 ms
(0.2229/0.2239 s) versus 12.9/13.0 ms within the facet pair. These ratios
are numerically identical (to 0.001-0.01) to the ellipsoid-gradient 2x2
gate's ratios, and the smoothed trajectories land 0.06-0.07 mm from the
ellipsoid-gradient trajectories on both meshes. Two independent
smooth-normal constructions (analytic ellipsoid gradient, mesh-based nodal
smoothing) produce the same far-family trajectory.

**Verdict: normal-direction error is the response-family selector.**
Smoothing the normal field (either way) removes the faceting artifact's
spurious twist restraint and moves the trajectory to the far family; the
move is mesh-robust (coarse and surface agree to 3.3-3.7 mm under smoothing
and snap timing converges to ~1 ms). Per the pre-registered bands the
verdict is *mixed* on whether faceting explains the whole coarse-vs-surface
split (ratios 0.14-0.38, between 0.25 and 0.75), but unambiguous on
direction: the normal-direction error drives the family selection.

Corollary recorded without overclaim: within the benchmark's discrete
definition (the reference code itself uses finite facet normals), more
accurate surface normals move the trajectory *away* from the ten-team
consensus, not toward it. The near-family agreement of the coarse facet
mesh is partly a cancellation between its faceting artifact and other
error sources; FEniCS's own artifact (1.698 N m/rad, the largest measured)
holds its trajectory the same way. The smooth-normal family (-33.7 to
-37.5 mm) is the mesh-consistent CoupFE answer for this operator; whether
it or the consensus is the better continuum estimate is not established by
this experiment and remains the open question.

Archive SHA-256 (all 320/320 steps, 8 ranks, clean tree):
- coarse facet: `fc4e49963e42568a812a68bcb75a90bc357497286701fe15d56ede8eadf2a5da`
- coarse nodal-smoothed: `65b2274bde453d9b974ceee3ba91eb49f0603dd3b2c014322699a95dd742d114`
- surface facet: `446944855f79dfc2f71eb4ea56dccf2001b9d4f59d1ac0589a83235a9c4a3a82`
- surface nodal-smoothed: `e1c4461ebbf2fce681b9a7a3240bc074672b81f71ea1e2849ef1d2ab5da8f496`

Note: these runs' environment resolves CoupFE Core `454f73c`, not the
approved `e2f42ed`; a dedicated anchor rerun (application `056c02d`, Core
`454f73c`, identical command) reproduces the retained clean gate to
1.3e-16 m, so the Core difference is numerically inert for this workload.

## Context note: the participant field is entirely tetrahedral

For interpretation of the family comparison: the benchmark distributed a
Gmsh-created tetrahedral domain to all participants (paper, Section 3), and
the solver table (Table 10) lists P2 or stabilized P1-P1 tetrahedral
discretizations with GL/KL quadrature for all ten teams. No participant used
hexahedral elements. The CoupFE Hex8 entry is therefore the only
hexahedral-discretization comparison in this study, and its facet artifact
pattern differs systematically from every participant's: at comparable
surface edge length the long-axis-rotation faceting stiffness is 1.285
N m/rad for the coarse CoupFE Q1 quadrilateral surface versus 1.698 for the
FEniCS P1 triangular surface. Quad-versus-triangle surface pattern is one of
the named confounders in the correlation table above.
