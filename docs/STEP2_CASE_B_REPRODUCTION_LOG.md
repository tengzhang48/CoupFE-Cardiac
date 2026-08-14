# Benchmark 1 Step 2 Case B reproduction log

This is the dated evidence and decision log for reproducing Benchmark 1,
Step 2, Case B from Aróstica et al., *A software benchmark for cardiac
elastodynamics*. It distinguishes verified inputs, numerical observations,
failed hypotheses, and pivots. A short run that does not traverse the snap is
only a software-plumbing check and cannot qualify the reproduction.

## Target identity

- Paper DOI: `10.1016/j.cma.2024.117485`.
- Audited open-access final PDF SHA-256:
  `d72cf68d41b20ab00694d99077a18ea9b9e18bd5b3c37b538fe8ec801d09a4ae`.
- Official data DOI: `10.5281/zenodo.14260459` (CC BY 4.0).
- Official archive identity: 23,180,741,494 bytes; MD5
  `75602be4777c4ca2262c2bcfd2134b15`.
- Application starting point: commit
  `e9b7d9084b24f7098170a221061eb159d0b090c1` on the local branch
  `step2-case-b-reproduction`.
- Target loads: the Step 1 combination of time-dependent active stress and
  endocardial follower pressure, with the Step 2 Case B material changes.
- Target time step: 0.001 s. Published displacement files contain 101 samples
  on the 0.01 s output grid over the 1 s cycle (Ambit has a slightly different
  first sample).

The Step 2 Case B changes are `a=295 Pa`, `a_f=92360 Pa`,
`a_s=12405 Pa`, `a_fs=1080 Pa`, and `sigma_0=100000 Pa`. The unchanged
dimensionless coefficients are `b=8.023`, `b_f=16.026`, `b_s=11.12`, and
`b_fs=11.436`. Other paper inputs remain fixed, including `rho=1000 kg/m^3`,
`eta=100 Pa s`, and `kappa=1e6 Pa`.

## Official spatial and output oracle

The selectively extracted official XDMF/HDF5 mesh has 4,577 nodes and 17,625
linear tetrahedra. Its coordinate extents are
`[-0.097, 0.02647058823529411] x [-0.0349998942, 0.0349998942] x
[-0.0349998942, 0.0349998942] m`. The HDF5 boundary mesh function contains
2,500 tag-1 triangles, 3,474 tag-2 triangles, and 256 tag-3 triangles; its
other 32,135 entries are unmarked interior facets. The extracted files and
their SHA-256 identities are:

- `ellipsoid_0.005.h5`:
  `b3a894a1616bae81d1aae66cb4e369d4cab1c3a0cf614390bcaa1d761adb4847`
- `ellipsoid_0.005.xdmf`:
  `007590d7ea1414dae1678c1e069fed1021bdc1f52a10af15d059044a98d43891`

The publisher plotting script selects exactly ten Step 2 Case B teams:
CARPentry, Ambit, 4C, Simula, CHimeRA, CHeart, lifeX, SimVascular P1,
SimVascular P2, and COMSOL. The unsuffixed SimVascular file is byte-identical
to the P2 file and is excluded rather than counted as an eleventh result.
The official pickle files contain usable displacement curves only. All stress
entries are placeholder mappings whose values are `None`. Nine volume entries
contain only `None`; CARPentry instead contains an all-NaN array. These fields
must be rejected as unusable, so no official stress or volume oracle is
available from this dataset.

The benchmark observation points are `p0=(0.025, 0.030, 0) m` and
`p1=(0, 0.030, 0) m`. This case contracts mainly in the x direction. The old
Step 0 Case B `p1 u_z=-5 mm` snap definition is therefore invalid for Step 2
Case B; the official `p1 u_z` does not reach -5 mm.

## Snap-aware acceptance measurements

All discrepancy experiments must include the contraction transition and a
settled post-snap interval. The initial event window is 0.16--0.32 s. For both
`p0` and `p1`, report normalized 10%, 50%, and 90% x-displacement drop times,
the full vector trajectory error, and the post-snap branch. Snap onset is a
diagnostic, not the principal success criterion: Hex8 is known to snap with
some onset offset, while its displacement difference becomes much larger on
the published post-snap branch. The
0.24--0.48 s component/vector histories, branch levels, slopes, and extrema
are therefore first-class comparison outputs. The official ten-team timing
ranges are:

| Point | 10% time (s) | 50% time (s) | 90% time (s) |
| --- | ---: | ---: | ---: |
| p0 | 0.1845--0.1943 | 0.2256--0.2374 | 0.2953--0.3155 |
| p1 | 0.1848--0.1942 | 0.2230--0.2357 | 0.2869--0.3092 |

A full-cycle comparison is also required because the rebound has a second
rapid transition around 0.48--0.58 s, with peak rebound speed near
0.505--0.515 s. Comparisons report each team, the all-team mean and spread,
and CoupFE-to-ensemble errors; the ensemble is not treated as an exact single
trajectory.

For phase reporting, use 0--0.16 s as pre-transition, 0.16--0.32 s as the
forward transition, 0.32--0.48 s as the loaded post-snap plateau, 0.48--0.56 s
as unloading/rebound, and 0.56--1.00 s as late relaxation. Across every
official sample in the 0.32--0.48 s plateau, the component ranges in mm are:

| Point | ux | uy | uz | Unanimous branch sign |
| --- | ---: | ---: | ---: | --- |
| p0 | -31.877 to -27.564 | +2.468 to +2.721 | +4.277 to +5.240 | (-,+,+) |
| p1 | -25.230 to -21.895 | -1.718 to -1.487 | +0.326 to +0.787 | (-,-,+) |

An opposite transverse sign or a different plateau branch is outside all ten
published results, not merely away from their mean.

## Decisions and pivots

### 2026-08-05: correct the benchmark identity before running

The starting driver labeled its public Case B path as Step 0: it set active
stress to zero and used the Step 0 material and activation magnitudes. Running
that path would not test the requested target. The first code change is an
explicit, fail-closed benchmark-step selector whose Step 2 Case B identity
requires both active stress and pressure and records the exact parameters and
loads in the result archive. Existing Step 0 behavior remains separately
named and is not relabeled.

### 2026-08-05: snap traversal is the diagnostic gate

The user noted that the difference appears only when snap occurs. Tiny serial
and MPI runs are retained solely as setup and load-plumbing checks. No element,
solver, or rank-equivalence conclusion will be drawn unless the simulated
interval crosses snap and includes post-snap behavior.

### 2026-08-05: understand Q1--P0 Hex8 before using Tet4/VMS

Q1--P0 Hex8 remains the primary formulation and is expected to work. The
investigation will first isolate snap-sensitive differences in active-stress
state/tangent handling, follower-pressure residual and tangent, local-pressure
condensation, time integration and continuation, geometry/fibers, quadrature,
and nonlinear globalization. The official Tet4 mesh with a VMS formulation is
a later controlled cross-formulation check if needed; it is not a substitute
for explaining a Hex8 branch discrepancy.

### 2026-08-05: gross geometry measures agree closely

A read-only comparison of the official linear-tetrahedron mesh against the
closed multiblock Hex8 geometry found close wall volume and boundary-area
agreement. The official Tet4 values are 177.6945 cm3 wall volume and
155.2639/228.5060/18.5773 cm2 for tags 1/2/3. The `(2,20,17)` Hex8 mesh gives
177.5475 cm3 and 155.3752/228.5308/18.5742 cm2 for endocardium,
epicardium, and base; the `(4,36,32)` Hex8 mesh gives 177.8217 cm3 and
155.4800/228.6825/18.5888 cm2. Thus a missing apex or a gross size/area error
is not the leading explanation for a post-snap displacement mismatch on the
closed mesh. Element topology, local geometry interpolation, fiber evaluation,
and discretized boundary/load operators remain open and must not be conflated
with these integral checks.

### 2026-08-05: post-snap branch is the primary failure signal

The user clarified that Q1--P0 Hex8 does undergo snap, with some onset
difference, but its displacement difference becomes substantially larger
afterward. Experiments will therefore separate (1) snap timing from (2)
post-snap error amplification, branch selection, and evolution. Active-stress/follower-load
tangents, local-pressure condensation, time integration, and nonlinear
globalization will be checked using accepted-state histories and force/energy
diagnostics where available. The current archive has full output-point, load,
and nonlinear histories, but only peak-time whole-field `det(F)` and condensed
pressure summaries; per-step energy, cavity volume, whole-field deformation,
and stability-eigenvalue diagnostics are a known instrumentation gap.

### 2026-08-05: no published participant used Q1--P0 Hex8

Table 10 of the paper reports tetrahedral discretizations for all published
monoventricle curves: either quadratic displacement with volumetric
penalization (`P2`) or stabilized linear displacement--pressure (`P1--P1`).
No official curve is an element-for-element Q1--P0 Hex8 oracle. This does not
invalidate the expectation that Hex8 should converge to the same response,
but it makes spatial/wall refinement and post-snap convergence evidence
mandatory. The current fixed generalized-alpha parameters
`alpha_m=0.2, alpha_f=0.4, gamma=0.7, beta=0.36` match Simula's reported
time discretization at 1 ms, making the Simula curve the cleanest named
temporal comparator; the full ten-team envelope remains the overall published
comparison. CARPentry and SimVascular P1 are the closest published mixed-linear
spatial comparators, but neither is a Q1--P0 analog: both use continuous
linear displacement and pressure with stabilization and KGen-alpha. CARPentry
reports Keast--Lyness degree-6 quadrature; SimVascular P1 reports
Gauss--Legendre degree 4. The paper does not publish per-team mesh or DOF
counts, and it does not give CARPentry's stabilization coefficient.

The official SimVascular P1/P2 pair gives a useful controlled scale for one
published discretization change. At 0.48 s their vector gaps are about
1.055 mm at p0 and 0.937 mm at p1, and both retain the same post-snap component
signs. CARPentry's stabilized P1--P1 curve selects that same branch. Therefore
a gross post-snap Hex8 branch or sign mismatch cannot be explained merely by
pointing to the published linear-versus-quadratic tetrahedral difference;
Hex8 convergence and its discrete operators still need direct investigation.

### 2026-08-05: interpret the small onset offset separately from branch error

A source-level audit found that the current source-matched generalized-alpha
path uses the same staging convention as Simula: inertia at the
`alpha_m=0.2` stage, constitutive/Robin/follower-pressure terms at the
`alpha_f=0.4` force stage, and both prescribed loads at
`t_(n+1)-0.4*dt`. The follower-pressure residual sign and its `0.6`
displacement-chain factor also agree with the upstream weak form. Therefore no
load-time shift or pressure-sign change will be introduced merely to align the
snap visually.

The named Simula curve itself differs from the ten-team mean fixed-threshold
crossings by as much as about 1.9 ms, while the published ten-team normalized
onset ranges span roughly 10--18 ms depending on point and fraction. A small
CoupFE onset shift can therefore be within inter-method spread. Each run will
report the signed timing bias in milliseconds on both the native 1 ms grid and
the publisher's 10 ms grid. Tuning is not considered successful unless the
0.24--0.48 s shape and the 0.32--0.48 s component levels/signs improve as
well.

The audit also identified a controlled spatial/fiber difference from Simula.
Simula uses a P2 tetrahedral displacement field and constructs its fibers at
P2 degrees of freedom from a CG1 Laplace field before component-wise
interpolation to quadrature. The current closed Hex8 path interpolates the Q1
Laplace/transfinite coordinates to each 2x2x2 Gauss point, reevaluates the
analytic frame there, and applies Gram--Schmidt normalization. This is a
plausible critical-load perturbation near a bifurcation, not proof that Hex8
is incapable of reproducing the branch. The diagnostic order is therefore:
pointwise `std-kappa` first, wall refinement, in-plane refinement, then
same-mesh switches to the element-mean `local-pressure-paper` and logarithmic
Q1--P0 laws. Solver-tolerance/globalization controls follow those spatial and
formulation tests.

### 2026-08-05: harden comparison identity and phase separation

An adversarial integration review found that a merely self-consistent custom
hash manifest could previously be labeled as publisher data. The Step 2B
comparator now anchors the byte identity of the reviewed manifest as well as
its official DOI and `figures.py::TEAMS_DATASETS_B` selection source. It also
recomputes the exact staged activation and pressure ODE histories rather than
accepting broad peak ranges.

The same review found that the first branch-shape metric began at 0.24 s,
overlapping the defined transition through 0.32 s. That could turn a small
phase shift into an apparent branch-shape error. The transition error now uses
0.16--0.32 s and the onset-independent settled plateau uses 0.32--0.48 s,
aligned at 0.32 s. The comparator reports the named Simula trajectory and
signed CoupFE-minus-Simula crossing-time errors alongside ensemble
mean/range statistics. Absolute plateau component signs and the vector octant
are reported separately from the direction of drift across the plateau.

### 2026-08-05: content-identify synchronized dirty-tree runs

The remote VM receives reviewed files before they are committed, so Git HEAD
plus `tree_state=dirty` was not enough to distinguish evolving experiments.
Every explicit benchmark archive now contains a canonical path-to-SHA-256
manifest of all application-owned result-producing sources and an aggregate
manifest SHA-256. Both the general reporter and the Step 2B comparator
recompute and verify that identity. The comparator additionally requires the
reviewed CoupFE Core revision, MPI source-matched full-cycle profile, exact
density/Robin/material/mass/fiber/Laplace/geometry setup metadata, successful
pre-solve audits, and one accepted nonlinear diagnostic per time step. A dirty
application worktree is accepted only when this exact content identity is
present; it is never asserted to be clean.

### 2026-08-05: the first cross-snap prefix shows a rate error, not a constant phase shift

The first scientifically meaningful run used four MPI ranks, the audited
`2x20x17` closed Hex8 mesh, pointwise `std-kappa`, consistent mass, the
injected Laplace field with GP-direct fibers, the source-matched
generalized-alpha method, and `dt=1 ms`. It completed 360/360 steps through
0.36 s with no invalid-domain rejection and at most three Newton iterations.
The retained archive SHA-256 is
`b04fe1ffbb418fae53a30a3f0f9bc07f7879b18ceecce050696b086920024279`;
its runtime-source aggregate SHA-256 is
`6b96395761dd3203f0e9ffab90a77d6389dca13cdad43490a1deac95073480f1`.

Linear-interpolated native-grid x-displacement crossings give:

| Point | Landmark | CoupFE (ms) | Simula (ms) | 10-team mean (ms) | CoupFE - Simula (ms) | CoupFE - mean (ms) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| p0 | -5 mm | 199.556 | 201.324 | 199.447 | -1.769 | +0.109 |
| p0 | -15 mm | 228.879 | 232.665 | 231.449 | -3.786 | -2.570 |
| p1 | -5 mm | 204.193 | 204.464 | 203.162 | -0.272 | +1.031 |
| p1 | -15 mm | 236.778 | 240.844 | 240.759 | -4.066 | -3.981 |

Thus the first -5 mm landmark is close to both published comparators, but the
-15 mm landmark becomes about 3.8--4.1 ms early relative to Simula. The
-5-to--15 mm travel time is 29.323 ms at p0 versus 31.341 ms for Simula and
32.585 ms at p1 versus 36.380 ms for Simula. This is evidence of a faster
computed transition, not a uniform load-clock offset. The prescribed load is
still evaluated at the verified `t_(n+1)-0.4*dt` stage; no timing correction
will be introduced to hide this response-rate difference.

The post-transition vector mismatch is already much larger than the onset
shift can explain. At 0.32 s the CoupFE-minus-Simula vector errors are
2.683 mm at p0 and 2.439 mm at p1. Delaying the CoupFE histories by their
respective -15 mm crossing offsets reduces these only to 2.534 mm and
2.299 mm. At 0.36 s the raw errors remain 2.555 mm and 2.341 mm, and the same
phase correction reduces them only to 2.510 mm and 2.310 mm. Component-wise,
at 0.36 s CoupFE gives p0 `(-32.553,+2.749,+3.597) mm` versus the contemporaneous
ten-team envelope `[-30.862,-28.977] x [+2.522,+2.712] x [+4.373,+5.008] mm`.
For p1 it gives `(-26.319,-1.472,-0.580) mm` versus
`[-24.619,-22.948] x [-1.589,-1.572] x [+0.338,+0.656] mm`; notably, p1 z has
the opposite sign from every published solution. This supports the user's
observation: onset is only modestly shifted, while the computed post-snap
branch is materially different.

The nonlinear record does not implicate a failed solve or branch jump caused
by rejected Newton states. There were 243 two-iteration steps and 117
three-iteration steps, the latter forming one continuous 0.166--0.282 s
interval; there were zero function-domain rejections. The largest tracked
point increments are smooth clusters near 0.233--0.234 s rather than an
isolated accepted-step discontinuity. This justified completing the identical
configuration through 1.0 s before changing any physical or numerical factor.

### 2026-08-05: the full-cycle run completes with a bounded global error and a branch-sign discrepancy

The canonical `mpi4_ga_full` run completed all 1,000/1,000 requested
increments through `t=1.0 s` on four MPI ranks. It retained the prefix's
2×20×17 closed Hex8 mesh, pointwise `std-kappa`, consistent mass, injected
Laplace field with GP-direct fibers, source-matched generalized-alpha method,
and `dt=0.001 s` Step 2 active-stress-plus-pressure contract. The external NPZ
has SHA-256
`23312a5e0147544eb9a4e6de004a166ada2722b70d3d39742f93aacd8a0fa0e6`.

The strict publisher comparison passed its input and physical-identity
plumbing gates. Against the official all-team mean on the 101-point full-cycle
grid, global relative L2 is 9.8038%, `p0`/`p1` vector relative L2 is
9.0555%/10.9271%, aggregate RMSE is 0.829875 mm, and maximum component error
is 2.846225 mm. The benchmark paper Eq. 21 RED is 28.3004% at `p0` and
35.2774% at `p1`. These percentages use different aggregation: relative L2
divides norms accumulated over the full history, while RED averages the
pointwise relative vector error over time. Neither metric has a paper-defined
acceptance threshold.

The full cycle preserves the prefix's principal branch concern. Throughout the
settled 0.32--0.48 s plateau, CoupFE `p1-z` is negative while all ten official
curves are positive. The bounded global norms therefore do not erase the
component-level wrong-sign result. This run is development evidence, not a
validation, convergence, or rank-independence claim.

The archive records application revision
`e9b7d9084b24f7098170a221061eb159d0b090c1` with a dirty tree and Core
`e2f42ed5772850a0a23a2ce434f430c287eae5c8`. It is content-identified, not
relabeled as clean: the exact result-producing runtime-source manifest has
aggregate SHA-256
`6b96395761dd3203f0e9ffab90a77d6389dca13cdad43490a1deac95073480f1`.
The public derived records are:

- [comparison JSON](../examples/cardiac_benchmark/results/step2_case_b_std_kappa_2x20x17_dt0p001.report.json),
  SHA-256 `098e316daaea369a2a595cf43829d28597e53d2ff5a38cf32388e01c8dfa74aa`;
- [normalized raw stdout](../examples/cardiac_benchmark/results/step2_case_b_std_kappa_2x20x17_dt0p001.raw.stdout.txt),
  SHA-256 `2bc7ddf633e1b905a1e8b42551ccf22eda4702f122cccd4370dd4c561a2a381c`;
- [reviewed runtime-source manifest](../examples/cardiac_benchmark/step2b_case_b_runtime_source_hashes.json),
  SHA-256 `d39ccbd6e67d7517b31a536f0c34472afb770bfd76916a3572e20d52acf39a41`;
  it binds the dirty result-producing file map to repository snapshot
  `d06c3e9` without relabeling later compatibility changes;
- [six-panel official comparison figure](figures/step2_case_b_comparison.svg),
  SHA-256
  `6a4c1b16fb2af6098d923380401fa1768ef85bd6062d670eee0114fdee694935`.

The report contains the derived 101-point CoupFE, ten-team range/mean, and
named Simula curves, so the figure can be regenerated without the external
NPZ or publisher pickle files:

```bash
python examples/cardiac_benchmark/plot_step2b_case_b.py --canonical
```

## Controlled-experiment ledger

Every retained run will record: source commit and dirty-tree diff identity;
Core commit; environment; mesh topology and hash; formulation; quadrature;
fiber construction and sampling; mass and time integrator; time step and load
horizon; rank/thread mapping; solver tolerances; accepted-step residual and
line-search history; invalid-deformation rejections; deformation-Jacobian and
energy summaries; complete `p0`/`p1`, active-stress, and pressure histories;
output hash; and comparison-report hash.

| Date | Run ID | One changed factor | Snap crossed? | Result | Interpretation / next action |
| --- | --- | --- | --- | --- | --- |
| 2026-08-05 | environment gates | none (setup only) | no | historical pre-run state: pending before execution | This row preserves the initial gate state; the later VM fast/native/MPI rows below record the completed passes. |
| 2026-08-05 | official-vs-closed geometry audit | mesh representation | not applicable | integral volume/areas agree within about 0.1% | Deprioritize gross geometry error; retain topology/interpolation/fiber/load-discretization checks. |
| 2026-08-05 | corrected Step 2B focused/full Python gates | implementation only | no | initial 41 focused tests pass; full suite 371 pass with only two intentionally stale current-source hash guards | Adversarial review identified and prompted official-manifest, exact-load, source-content, canonical-run, and phase-separation hardening before VM execution. Historical result hashes remain immutable. |
| 2026-08-05 | hardened local integration gate | implementation only | no | 379 fast tests pass, 29 slow/MPI tests deselected | Official manifest and exact-load/source/canonical-profile guards, named Simula metrics, and onset/plateau separation are integrated. |
| 2026-08-05 | VM fast/native/MPI gates | environment and parallel implementation only | no | 378 fast pass, 1 reference-data test skipped; 5 native compiled pass; 24 MPI pass across 1/2/4 ranks | The Step 2B source set is byte-identical on the VM, the generalized-alpha pressure smoke regression is repaired, and compiled/MPI plumbing is ready for minimal runtime checks. |
| 2026-08-05 | tiny `1x4x1` closed-mesh attempt | mesh coarsening for plumbing | no | rejected before compilation: wall volume -20.476%, endo/epi areas about -10.5%, base area -2.418% | Preserve the 1% geometry audit; do not use a geometrically invalid tiny mesh. Pivot plumbing to the audited `2x20x17` mesh. |
| 2026-08-05 | `2x20x17` Laplace field | geometry-faithful coarse mesh | not applicable | 5,403 nodes, 3,520 Hex8; field SHA-256 `135dad96e20572140b36fd1e4857d6e6551e33d2f88c956c2d9915f6ed282e6c`; metadata SHA-256 `9012fd453734e8a70348002fe0410adc40684d43fb297c806a4195f12e3994f2` | Geometry audit passes within 0.083%; Laplace field differs from analytic layer coordinate by at most 0.1616, corresponding to a reported 19.39-degree maximum fiber-angle change. |
| 2026-08-05 | serial BE plus MPI generalized-alpha 1/2/4 rank, two-step plumbing | execution/rank only | no | all archives complete with exact Step 2B/source identities; MPI max output-point difference versus rank 1 is `1.10e-19 m`; staged loads bitwise equal | Plumbing only. No snap or reproduction claim; proceed with 4-rank canonical prefix through the plateau. |
| 2026-08-05 | `mpi4_ga_t0p36` | first canonical 4-rank cross-snap prefix | yes, through early plateau | -5 mm crossings are within 0.27--1.77 ms of Simula, but -15 mm crossings are 3.79--4.07 ms early; 0.32 s vector errors are 2.683/2.439 mm and p1 z selects the opposite published sign; zero rejected states | Treat onset-rate and plateau-branch errors separately. Complete the same configuration through 1.0 s before changing one factor. |
| 2026-08-05 | `mpi4_ga_full` | extend the same canonical configuration from 0.36 to 1.0 s | yes, full cycle | completed 1,000/1,000; global relative L2 9.8038%, aggregate RMSE 0.829875 mm, paper RED 28.3004%/35.2774%; `p1-z` plateau has the opposite official sign | Retain as content-identified dirty-tree development evidence; use the component-level branch discrepancy when choosing the first one-factor refinement. |
| 2026-08-11 | corrected-setup full-cycle diagnostic (`mpi8`, local pressure, `tip_refine=6.0`) | rerun Step 2 Case B with the corrected straight-wall geometry, physical-coordinate frame reconstruction, current generalized-alpha Q1/P0 path, consistent mass, and tip grading | reported full cycle | The retained compact comparison reports 1,000/1,000 steps, relative L2 of 3.76% (`p0`) / 5.94% (`p1`), and a positive `p1-z` plateau (+0.24 to +0.50 mm vs the official band +0.33 to +0.79 mm). It binds NPZ SHA-256 `63a8de59b7b8b9ab309896ff69989d6ff89f6dfe2532151605486ad67967dd41` but not the complete release-grade source, command/environment, solver/deformation, and ten-team role/hash provenance. | **Promising provenance-incomplete diagnostic; not publisher or reproduction evidence.** The formulation also differs from the old run (Q1/P0 local pressure vs pointwise `kappa`), so the change bundles geometry/frame and formulation. The separately recorded stdout SHA-256 is `b05197b3471bdec7da9d23e6d4c6486387d8a88755c21f142e01d047ff0ed7d2`; reported elapsed time is 1114.8 s. |

## Background execution and recovery boundary

The Compute Engine VM is a 16-vCPU standard instance (8 physical cores with
SMT) and production commands are limited to 4 or 8 MPI ranks with one thread
per rank. A detached `tmux` session named `coupfe-step2b` survives a browser or
SSH disconnect. It does not survive a guest reboot. Restartable checkpoints
must therefore be implemented and tested before the fine run; VM automatic
restart alone is not process restart.
