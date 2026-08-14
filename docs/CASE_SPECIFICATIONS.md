# Case specifications

What each retained case actually solves: geometry, boundary conditions,
material model, element type, and time integration. Rows identify the bound
report for each case and explicitly mark the later diagnostic whose compact
report does not bind its complete configuration and execution provenance.

Two records that share a case name can differ in formulation. The selected
Step 0 Case A and Step 0 Case B records use the condensed Q1/P0 mean-`log(J)`
volumetric response. The source-bound Step 2 Case B development record uses the
paper's pointwise `kappa` penalty, while a later provenance-incomplete Step 2
diagnostic uses Q1/P0 mean-`log(J)`. Read the formulation and provenance rows
before comparing any two numbers.

Numerical constants are not duplicated here. The single source of truth is
[`examples/cardiac_benchmark/benchmark_parameters.py`](../examples/cardiac_benchmark/benchmark_parameters.py),
which every driver reads and every report records.

---

## Intended current-driver physical contract

The current drivers implement these shared Step 0 Case A, Step 0 Case B, and
Step 2 Case B properties. Historical retained records below are explicit
exceptions when their source predates a geometry or frame correction; they are
kept under their recorded identity rather than relabeled as current output.

**Geometry.** A truncated-ellipsoid left ventricle, built as a **closed
five-block straight-wall Hex8 domain** (`closed_multiblock_disk`). Endocardial
and epicardial surfaces are confocal ellipsoids with
`R_SHORT_ENDO = 2.5e-2`, `R_SHORT_EPI = 3.5e-2`, `R_LONG_ENDO = 9.0e-2`,
`R_LONG_EPI = 9.7e-2` m. The wall is ruled in straight Cartesian segments
between them, matching the reference toolkit. There are exactly **three
boundary classes** — endocardium, epicardium, base — and **no free apex
surface**; every exterior face carries exactly one label. Landmarks are the
physical points `p0 = [0.025, 0.03, 0]` m and `p1 = [0, 0.03, 0]` m, sampled by
Hex8 inverse-isoparametric interpolation (reconstruction error below
`1.5e-16` m).

**Fibers.** Rule-based, reconstructed at Gauss points from physical coordinates
and a Q1 Laplace transmural field, using the MIT-licensed port in
[`structural_directions.py`](../examples/cardiac_benchmark/structural_directions.py).
The helix angle runs `-60°` at the endocardium to `+60°` at the epicardium.
Fibers are carried as 3x3 structural tensors in Gauss-point state, never as
vectors.

**Boundary conditions.** There are **no Dirichlet conditions and no clamped
nodes**. All support is compliant, in the reference configuration:

| Surface | Condition | Coefficients |
|---|---|---|
| Endocardium | deformation-dependent follower pressure `p J F^-T N` (consistent Nanson load) | load only; not a support |
| Epicardium | **normal-only** spring + dashpot, projector `N (x) N` | `alpha = 1e8` Pa/m, `beta = 5e3` Pa s/m |
| Base annular cut | **full-vector** spring + dashpot, identity projector | `alpha = 1e5` Pa/m, `beta = 5e3` Pa s/m |

The base full-vector term is the intended rigid-body support. All six rigid
modes are restrained (rank 6/6 for both spring and dashpot); see
[`MESH_REFINEMENT_GUIDE.md`](MESH_REFINEMENT_GUIDE.md) for the caveat that most
of the long-axis rotational restraint is a facet-normal discretization artifact
rather than prescribed physics.

**Material.** Holzapfel--Ogden transversely-isotropic passive myocardium with
an isochoric split: isotropic `a, b`; fiber `a_f, b_f` and sheet `a_s, b_s`
terms with a Heaviside tension-only switch; and a fiber--sheet `a_fs, b_fs`
coupling. Active stress is **additive**, `S += Ta (f0 (x) f0)`, driven by a
Bestel--Clement--Sorine activation ODE integrated in Python and applied as a
per-step property. Viscous damping is a pseudo-potential
`eta = 100` Pa s evaluated as a velocity-consistent Green--Lagrange rate at the
`alpha_f` stage.

**Element.** Trilinear **Q1 Hex8** displacement, 3 DOF per node,
`2x2x2` Gauss quadrature, with **consistent** (not lumped) mass
(`consistent_q1_hex8`).

**Time integration.** Source-matched **generalized-alpha**, contract
`simula-source-matched-v1`, with `alpha_m = 0.2`, `alpha_f = 0.4`,
`gamma = 0.7`, `beta = 0.36` — the Chung--Hulbert optimal-damping member at
spectral radius `rho_inf = 2/3`. Convention:
`x_{n+1-alpha_f} = alpha_f x_n + (1-alpha_f) x_{n+1}`. The Robin spring and
dashpot, the follower pressure (**including its load time**), and the condensed
volumetric block are all evaluated at the `alpha_f` stage; inertia uses the
`alpha_m` acceleration stage. Time step `dt = 1e-3` s.

**Solver.** PETSc SNES; the distributed runs use an FGMRES/GAMG profile with a
rigid-body near-nullspace rebuild, on eight MPI ranks.

---

## Step 0 Case A — active contraction

| Property | Value |
|---|---|
| Load | active contraction only; **cavity pressure identically zero** |
| Mesh | closed `4x36x32`, uniform (no tip refinement), 23,616 elements |
| Volumetric response | condensed Q1/P0 **mean-`log(J)`** local pressure, `K = 1.0` MPa (`hex8_local_pressure_p0_condensed_logj`) |
| Element / mass | Q1 Hex8 / consistent |
| Integration | generalized-alpha (0.2, 0.4, 0.7, 0.36), `dt = 1e-3` s, `t_end = 1.0` s, 1000/1000 steps |
| Bound report | [`case_a_local_pressure_4x36x32_dt0p001.report.json`](../examples/cardiac_benchmark/results/case_a_local_pressure_4x36x32_dt0p001.report.json) |

**Status: historical approximate-agreement evidence, not a current-setup
reproduction.** Vector-history relative L2 against official Simula is 8.57% at
`p0` and 12.26% at `p1`. The retained run **predates the straight-wall geometry
and physical-coordinate-frame corrections**, so it must not be quoted as a
current result. Its Step 0A identity is labelled `legacy-inferred` because the
archive predates explicit benchmark-identity fields.

## Step 0 Case B — passive pressure loading

| Property | Value |
|---|---|
| Load | endocardial follower pressure only; **active tension disabled** |
| Mesh | closed `2x20x17` and `4x20x17`, **tip refinement strength 6.0** for the full-cycle pair; the 0.32 s clean gate is uniform |
| Volumetric response | condensed Q1/P0 **mean-`log(J)`** local pressure, `K = 1.0` MPa |
| Element / mass | Q1 Hex8 / consistent |
| Integration | generalized-alpha (0.2, 0.4, 0.7, 0.36), `dt = 1e-3` s; full cycle `t_end = 1.0` s (1000/1000), clean gate `t_end = 0.32` s (320/320) |
| Bound reports | [full cycle](../examples/cardiac_benchmark/results/step0b_tip6p0_full_cycle_comparison.report.json), [0.32 s clean gate](../examples/cardiac_benchmark/results/step0b_case_b_clean_frame_0p32.report.json) |

**Status: partial quantitative reproduction** — the headline case. Relative L2
against the exact ten-team mean is 8.98% / 9.20% at `p0` / `p1`. The response
exhibits pressure-driven **snap-through and recovery**, which is why this case
is sensitive to mesh and boundary discretization. The volumetric response is a
documented near-incompressibility variant, **not the paper's pointwise
volumetric law**; that bounds the claim to method comparison rather than exact
discrete-equation identity. See
[`BENCHMARK_REPRODUCTION_STATUS.md`](BENCHMARK_REPRODUCTION_STATUS.md).

## Step 2 Case B — combined activation and pressure

The source-bound development record applies active contraction and endocardial
follower pressure on a closed `2x20x17` Q1 Hex8 mesh with consistent mass and
source-matched generalized-alpha at `dt = 1e-3` s. A later artifact is described
by its renderer and reproduction log as a corrected-setup Q1/P0 run, but its
compact report does not bind that complete configuration. The records are not
interchangeable:

| Record | Volumetric response | Provenance and status |
|---|---|---|
| [Source-bound development comparison](../examples/cardiac_benchmark/results/step2_case_b_std_kappa_2x20x17_dt0p001.report.json) | paper pointwise `kappa` penalty (`hex8_standard_pointwise_kappa`) | Exact dirty-tree runtime-source manifest; predates the straight-wall geometry and physical-frame corrections. Global full-history relative L2 is 9.80%, but `p1-z` has the opposite sign from all ten official curves. **Development evidence only.** |
| [Later corrected-setup diagnostic](../examples/cardiac_benchmark/results/step2b_current_rerun_comparison.report.json) | condensed Q1/P0 mean-`log(J)` local pressure | The compact report binds the result hash and comparison metrics, but not the complete source, command/environment, solver/deformation, and ten-team role/hash provenance required for release-grade evidence. **Provenance-incomplete diagnostic only.** |

**Status: no Step 2 reproduction claim.** The later trajectory is promising,
but it must not replace the source-bound development record or be promoted to
headline evidence until its full provenance is reconstructed and guarded.

---

## CLI geometry default

Both serial and MPI drivers now default to `closed-multiblock`. Selecting
`polar-ring` with a positive apex offset is an explicit historical choice and
emits a not-the-benchmark-geometry warning.

A bare `--case B` command selects the closed geometry, but it is not by itself a
controlled reproduction command: formulation, mass, time integration, fiber
sampling, and the Laplace field must still be selected explicitly. Use the
commands in [`CONTROLLED_BENCHMARK_RUNS.md`](CONTROLLED_BENCHMARK_RUNS.md).

## Archived open-tip campaign — not part of the comparison

The `truncated_polar` records under
[`results/archive/truncated_polar/`](../examples/cardiac_benchmark/results/archive/truncated_polar/)
and
[`docs/figures/archive/truncated_polar/`](figures/archive/truncated_polar/)
are retained **for comparison of methods and for lessons only**. They are not
benchmark evidence and are excluded from the current comparison documents.

They differ from every current case in all five respects:

| | Current cases | Archived open-tip campaign |
|---|---|---|
| Geometry | closed five-block, straight wall, three boundary classes | `polar_ring` with `apex_offset = 0.2`: **the apex tip is cut off**, leaving a fourth traction-free annular surface |
| Boundary conditions | endo pressure + epi normal Robin + base full Robin | same three, **plus an unintended free tip**; the cut removed pressure and Robin facets and lowered the undeformed pressure resultant by 4.59% |
| Material | Holzapfel--Ogden, as above | same law, but F-bar or Q1/P0 depending on the record |
| Element | Q1 Hex8, consistent mass | Q1 Hex8, **lumped** mass in several records |
| Time integration | generalized-alpha, `dt = 1e-3` s | **backward Euler**, `dt = 2e-3` or `4e-3` s |
| Mesh labels | `n_t x n_core x n_radial` | `n_t x n_mu x n_theta` — a different convention that looks identical |

That last row is a live hazard: `2x36x48` (archived, open apex) and `2x36x32`
(current, closed) differ by one character and denote different domains.

---

## Reading a mesh label

`A x B x C` is `n_t x n_core x n_radial`: through-wall element layers, elements
per side of the central apex square, and radial layers in each outer block.

A trailing **`tip F`** means tip refinement at **strength `F`** — a clustering
factor, not an element count. `F` scales apex-adjacent meridian spacing by
`1/F` (so `tip 6.0` gives one sixth of uniform spacing at the apex) while base
spacing grows by up to `2 - 1/F`. It relocates nodes at a fixed element count,
so **it adds no resolution**. Absence of a `tip` token means strength 1.0,
which is the uniform mesh exactly.

Always write the strength with its meaning — "tip refinement strength 6.0
(apex spacing 1/6 of uniform)" — rather than a bare `tip_refine=6.0`, which
does not say what 6.0 measures. Full detail, including why tip refinement and
uniform refinement move the solution in opposite directions, is in
[`MESH_REFINEMENT_GUIDE.md`](MESH_REFINEMENT_GUIDE.md).
