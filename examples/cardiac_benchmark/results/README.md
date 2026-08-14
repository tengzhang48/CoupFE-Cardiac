# Cardiac benchmark result index

This directory keeps reviewable output with the runnable cardiac examples. Its
root contains a selected closed Case A report and a closed-multiblock Step 2
Case B development report. Both producing sources predate the toolkit-matched
straight-wall geometry mapping and physical-coordinate structural-frame
reconstruction, so they remain qualified under their exact recorded source
identities rather than being relabeled as current-source output. The root also
contains a current-source Step 0B **dual-run full-cycle comparison report** for
the closed `2x20x17` and `4x20x17`, `tip_refine=6.0` runs and the earlier **prefix diagnostic**
that binds two clean 0.32 s gates. The prefix report's pause decision is a
historical decision for that gate and is superseded for the separately
configured tip-refined runs. The earlier external `4x20x17`,
`tip_refine=6.0` controlled prefix remains a prefix-continuity record; its
shared states agree with the full-cycle result to roundoff. The directory also
includes a later corrected-setup Step 2 Q1/P0 diagnostic whose compact
provenance is explicitly insufficient for a reproduction claim. All checked-in
Case A and Case B records made with
the historical `polar_ring` mesh and `apex_offset=0.2` are
under the
[`archive/truncated_polar/`](archive/truncated_polar/) hierarchy.

The archived open-tip construction is not Benchmark 1 geometry. Its reports,
logs, and reduced record remain reproducible project history and regression
evidence for their exact configurations, not current Benchmark 1 validation
evidence. Nothing in the move relabels or changes their numerical content.

The generated CoupFE NPZ archives are not committed. The external benchmark
pickle files are not committed either; they remain part of the separately
distributed CC BY 4.0 Zenodo dataset. The historical Step 0 v2 reports retain
complete result histories and accepted-step diagnostics. The selected fine
Case A record is deliberately compact: it retains the reviewed 101-point
curves, exact configuration and sampler metadata, external result/manifest/log
identities, and full-precision RED without copying multi-megabyte accepted-step
diagnostics. The Step 2 development report uses a distinct compact publisher
schema. Retained paper-comparison reports include source identities,
external-file hashes, derived benchmark curves, and full-precision metrics
appropriate to their declared schema. The distinct Step 0B prefix-diagnostic
schema instead binds the clean short runs, source identities, and diagnostic
artifacts without embedding the ten-team curves.

For a completed, explicit recorded Step 0 run, first create the full forensic
generic-v2 report with `post.py`, then distill it with
`publish_step0_comparison.py`. The publisher verifies the external campaign
manifest and raw stdout identities, completion, closed geometry and boundary
audits, selected local-pressure/generalized-alpha method, official reference
selection, and RED recomputation. It emits the case-neutral
`coupfe-cardiac-step0-retained-comparison-v1` shape without copying accepted-
step/KSP histories; exact public runs remain hash-pinned by the release guard.

The comparison reports use the exact 10 files selected by upstream
`results_time_curves/figures.py`. The additional base-name SimVascular file in
the archive is byte-identical to selected SimVascular P2 and is recorded as an
excluded duplicate alias. Reports that correct the earlier 11-file wildcard
record its exact predecessor SHA-256. The separately named corrected-law Case
A report is a new physical-law record, not a relabeling or replacement of the
`62ad760` result.

## Files

| Record | JSON | Console output |
|---|---|---|
| Step 0 Case B current dual-run full cycle, closed 2×20×17 and 4×20×17, `tip_refine=6.0`, generalized-alpha, `dt=0.001 s`, eight ranks | [`step0b_tip6p0_full_cycle_comparison.report.json`](step0b_tip6p0_full_cycle_comparison.report.json) | Not distributed; both external NPZ identities are retained in the comparison report; clean four-layer replay metadata is **`1e333b29b05f01dedce9272b32b82ce6ccfda56036c1c8f57eb395b9b4494800` / `0810a9d0c944345464f48917c707a7d3df7187b543c74efd4b0384542b95aed7` / `2458e7c` / `1778.4 s (29.6 min)`** |
| Step 0 Case B completed through-wall prefix, closed 4×20×17, `tip_refine=6.0`, generalized-alpha, `dt=0.001 s`, eight ranks, through 0.32 s | Not packaged; external NPZ SHA-256 `774a7dc5fc970bb744ff0188f0f428cff54c532b401a365643f4b626584d7acf` | Not distributed; retained as a prefix-continuity record |
| Step 0 Case B current-frame diagnostic, closed 2×20×17 and wall-only 4×20×17, generalized-alpha, `dt=0.001 s`, through 0.32 s | [`step0b_case_b_clean_frame_0p32.report.json`](step0b_case_b_clean_frame_0p32.report.json) | Not distributed; both external stdout identities are retained in the compact report |
| Selected fine Step 0 Case A, Q1/P0, pre-straight-wall closed 4×36×32, generalized-alpha, `dt=0.001 s` | [`case_a_local_pressure_4x36x32_dt0p001.report.json`](case_a_local_pressure_4x36x32_dt0p001.report.json) | Not distributed; external stdout identity is retained in the compact report |
| Step 2 Case B development, pointwise `kappa`, pre-straight-wall closed 2×20×17, generalized-alpha, `dt=0.001 s` | [`step2_case_b_std_kappa_2x20x17_dt0p001.report.json`](step2_case_b_std_kappa_2x20x17_dt0p001.report.json) | [`step2_case_b_std_kappa_2x20x17_dt0p001.raw.stdout.txt`](step2_case_b_std_kappa_2x20x17_dt0p001.raw.stdout.txt) |
| Step 2 Case B corrected-setup Q1/P0 comparison — provenance-incomplete diagnostic only | [`step2b_current_rerun_comparison.report.json`](step2b_current_rerun_comparison.report.json) | Not distributed; the compact report binds the NPZ hash but not complete release-grade execution provenance |
| Archived truncated-polar Case A and Case B records | [archive inventory](archive/truncated_polar/README.md) | [archive policy](archive/truncated_polar/README.md) |

## Step 0 Case B current tip-refined full-cycle pair

As of 2026-08-07, the closed `2x20x17` and `4x20x17`, `tip_refine=6.0`
Q1/P0 runs have completed 1,000/1,000 increments on eight MPI ranks with
consistent mass, source-matched generalized-alpha, `eta=100 Pa s`, and
`dt=0.001 s`. The report above binds both external archives plus the exact
local FEniCS and ten-team inputs used for the same-time comparison. The
corresponding figure is
[`docs/figures/step0b_tip_refine_full_cycle.svg`](../../../docs/figures/step0b_tip_refine_full_cycle.svg).

Relative to the two-layer result, the four-layer trajectory has
0.009831/0.002622 mm pre-snap pairwise RMSE, 1.633083/1.405866 mm maximum
snap-window separation at 0.250/0.247 s, 0.407604/0.361832 mm full-cycle RMSE,
and 0.034598/0.051776 mm cycle-end separation (`p0`/`p1`). Over
0.75--0.999 s its pairwise RMSE is only 0.051556/0.062702 mm. Against FEniCS,
full-shared-history RMSE changes from 1.092675/1.168502 mm for two layers to
1.016944/1.072951 mm for four, and the maximum gap changes from
5.551997/5.023942 to 4.188388/3.753416 mm. Relaxation RMSE nevertheless
worsens by 16.2%/14.9%, and the ten-team containment fractions are mixed.
This phase-resolved result supports numerical transient/timing sensitivity;
it does not prove convergence or physically distinct solution branches.

The first full-cycle four-layer candidate was produced from a dirty
application tree and is explicitly non-retained provenance
diagnostic evidence. A clean isolated replay completed under Core `454f73c`
and runtime-source SHA-256
`f8d9469a101709d11460a0803b6c031001192ccce61c472d032355b01070da05`;
its states match the candidate to roundoff. Clean retained identity is
**application `2458e7c`, NPZ SHA-256
`1e333b29b05f01dedce9272b32b82ce6ccfda56036c1c8f57eb395b9b4494800`,
stdout SHA-256
`0810a9d0c944345464f48917c707a7d3df7187b543c74efd4b0384542b95aed7`,
elapsed `1778.4 s (29.6 min)`**.

## Historical Step 0 Case B current-frame prefix diagnostic

> Superseded status: this artifact correctly preserves the 2026-08-05 pause
> decision for its two clean prefix gates. The separate tip-refined full-cycle
> run above later completed; the prefix report itself remains unchanged.

The compact Step 0B report binds two clean eight-rank runs at application/Core
revisions `056c02d`/`e2f42ed`: closed 2×20×17 and through-wall-refined
4×20×17. Both use the toolkit-matched straight-wall geometry,
physical-coordinate structural-frame reconstruction, condensed Q1/P0
`log(J)` local pressure, consistent mass, source-matched generalized-alpha,
the paper viscosity, `dt=0.001 s`, and the canonical 1 s load schedule. They
complete only the 0--0.32 s diagnostic prefix.

The report retains exact identities for the external archives, logs, local
FEniCS reference inputs, the four-way mesh-split analysis, and the Robin
rigid-mode audit. It records the measured physical-frame change (at most
0.0357 mm) and the strong surface-refinement/Robin-faceting association. The
full 1 s run is explicitly paused because the snap window already exposes the unresolved
local axial difference. This stop is a diagnostic decision, not a solver
failure; the prefix is not convergence or full-cycle validation evidence.

## Selected fine Case A comparison record (pre-straight-wall mapping)

The compact report is the public figure input for the selected fine Case A
trajectory. Its source predates the toolkit-matched straight-wall mapping and
physical-coordinate structural-frame reconstruction; the report remains bound
to that recorded geometry/frame implementation and is not relabeled as output
from current source. The run completed 1,000/1,000 increments on eight MPI
ranks with a closed 4×36×32 Hex8 mesh (23,616 elements, 29,885 nodes, 89,655
displacement DOFs), Q1/P0 condensed local pressure, consistent mass,
source-matched generalized-alpha, and `dt=0.001 s`. It records clean
application/Core revisions `016a4f9`/`e2f42ed`.

The external NPZ is 16,232,720 bytes with SHA-256
`ba9b31ec533398be1f39fc9a898e72f77d9587c90f9b7d9e00ce91e4d2ae6a6c`.
The external campaign manifest is 4,080 bytes with SHA-256
`db769328c9ba13079311cecaa33f1bbea0c1b1d9b33c5681118cf116b829b938`;
the external stdout is 19,467 bytes with SHA-256
`beaef006e462bdfad4ce2e827620488aed633471d5907fd4debb3fab23187331`.
The raw transcript is not copied here because it contains machine-local
host/path text.

Against the official ten-team mean, benchmark-paper Eq. 21 RED is
0.3337402/0.5024615 at `p0`/`p1`. The compact report retains the exact ten-team
file manifest and the 101-point CoupFE and mean curves needed to reproduce the
figure. The archive predates explicit `BENCHMARK_ARCHIVE_FIELDS`, so Step 0A is
clearly labeled `legacy-inferred`, based on `case=A`, the source-identified
closed Case A implementation, canonical active tension, and identically zero
cavity pressure. It is not represented as an identity recorded directly in
the archive.

## Step 2 Case B full-cycle development record (pre-straight-wall mapping)

The Step 2 run completed all 1,000/1,000 requested increments on four MPI
ranks. Its content-identified source predates the toolkit-matched straight-wall
mapping and physical-coordinate structural-frame reconstruction, so it remains
development evidence for that recorded implementation. It uses the
topologically closed 2×20×17 Hex8 mesh, pointwise `kappa`, consistent
mass, GP-direct fibers from an injected Laplace field, source-matched
generalized-alpha at `dt=0.001 s`, and the Step 2 active-stress-plus-pressure
physical contract. The external `mpi4_ga_full.npz` is 3,835,766 bytes with
SHA-256
`23312a5e0147544eb9a4e6de004a166ada2722b70d3d39742f93aacd8a0fa0e6`.

Against the official all-team mean on the full 101-point cycle, aggregate
relative L2 is 9.8038%, aggregate RMSE is 0.829875 mm, and maximum component
error is 2.846225 mm. Vector relative L2 is 9.0555% at `p0` and 10.9271% at
`p1`. The benchmark paper Eq. 21 RED values are 28.3004%/35.2774% at
`p0`/`p1`. The relative-L2 metrics divide norms accumulated over the full
history; RED instead averages the pointwise relative vector error. The values
are therefore not interchangeable, and the paper provides no acceptance
threshold. The `p1-z` plateau has the opposite sign from all ten official
curves.

Application revision `e9b7d90` reported `tree_state=dirty`; it is not relabeled
as clean. The exact result-producing source content is identified by runtime
manifest SHA-256
`6b96395761dd3203f0e9ffab90a77d6389dca13cdad43490a1deac95073480f1`,
whose exact file map is frozen in the
[`reviewed runtime-source manifest`](../step2b_case_b_runtime_source_hashes.json)
and repository snapshot `d06c3e9`. Current Python 3.9 compatibility changes
are not relabeled as result-producing content. Core is `e2f42ed`. This is
development evidence, not validation,
convergence, or rank-independence evidence. The corresponding
[`six-panel figure`](../../../docs/figures/step2_case_b_comparison.svg) can be
regenerated without the external NPZ or publisher pickles:

```bash
python examples/cardiac_benchmark/plot_step2b_case_b.py --canonical
```

## Step 2 Case B corrected-setup diagnostic

The separate
[`step2b_current_rerun_comparison.report.json`](step2b_current_rerun_comparison.report.json)
describes a later closed, corrected-setup Q1/P0 full-cycle trajectory. It binds
the corrected NPZ hash, the legacy report hash, the dataset DOI, and derived
comparison metrics. It does not bind the complete application/Core identities
and tree states, command/environment and stdout, solver completion and
configuration, deformation audit, or exact ten-team role/hash manifest.

The corresponding
[`diagnostic figure`](../../../docs/figures/step2b_current_rerun_comparison.svg)
is retained so the observation remains reviewable, but neither artifact is
release-grade Step 2 reproduction evidence. Do not substitute it for the
source-bound development record above.

## Archived truncated-polar source-identified records

The earlier reports are retained in three application-checkpoint groups rather
than being rewritten as if one source produced all of them:

| Report group | Application revision | CoupFE Core revision |
|---|---|---|
| Historical-law Case A and the five checkpoint-`62ad760` Case B reports listed above | `62ad760d2a1731bb9668897863ac026d3768194e` | `454f73ce2de284262b214a2b37bd676c6aca3c0a` |
| Q1/P0 2×24×32, `dt=0.004 s`, and 2×36×48, `dt=0.002 s` | `e07993bcf1166bd20eb87370c0b458552753e7ee` | `454f73ce2de284262b214a2b37bd676c6aca3c0a` |
| Corrected-law Case A | `6839c13b5bc80ec06c897684c51f503e80bd4b19` | `e2f42ed5772850a0a23a2ce434f430c287eae5c8` |

All three groups used Python 3.10.8, NumPy 1.26.4, SciPy 1.15.2, CoupFE 0.0.1,
the nondegenerate structured mesh with `apex_offset=0.2`,
`cg1_gram_schmidt` fibers, and `hex8_reference_isoparametric` output sampling.
Each report identifies exact files from the verified external dataset.

These six checkpoint-`62ad760` reports predate the later
`function_domain_rejections` and `last_function_domain_error` diagnostic
fields. Their absence must not be rewritten as evidence that the newer
invalid-trial recovery policy was active. Runs made with the newer policy
record those fields per step; they still require a valid accepted state and
the unchanged final residual rule.

The historical Case A run exercises the 1 s active transient with F-bar, eight
Hex8 elements, backward Euler, and guarded Core Newton. It completed 500/500
requested steps. At peak activation (`t=0.48 s`), `p0` was
(-20.901, -2.637, -0.142) mm and `p1` was
(-14.571, -1.886, -0.042) mm. The 64 retained peak-load Gauss-point
deformation Jacobians ranged from 0.827202 to 1.314427. RED against the
`step_0A` all-team mean is 0.5218981 for `p0` and 0.6816066 for `p1`.

Commit `6839c13` corrected the fiber/sheet stress to differentiate the complete
smooth-switch energy. Its clean rerun keeps the same mesh, loads, time step,
integrator, and solver and completed 500/500 steps. Peak displacements remain
(-20.901, -2.637, -0.142) mm at `p0` and
(-14.571, -1.886, -0.042) mm at `p1`; its peak `det(F)` range is
0.827206--1.314409 and its RED is 0.5061198/0.6790790. The corrected report
records material model identity
`holzapfel-ogden-smooth-switch-complete-energy-derivative-v1`. The old report
remains unchanged and is not used as the corrected-law numerical oracle.

The archived truncated-polar Case B records are:

| App | Formulation | Mesh and time | Accepted steps | Peak time | Peak `p0` (mm) | Peak `p1` (mm) | Peak Gauss-point `det(F)` | RED (`p0`, `p1`) |
|---|---|---|---:|---:|---|---|---|---|
| `62ad760` | Q1/P0 local pressure | 2×12×16, `dt=0.004 s` | 250/250 | 0.480 s | (+17.917, -0.405, -0.601) | (+13.380, +1.099, -1.026) | 0.894270–1.058401 | 0.8831859, 0.8657258 |
| `62ad760` | Q1/P0 local pressure | 2×12×16, `dt=0.002 s` | 500/500 | 0.482 s | (+17.923, -0.406, -0.602) | (+13.384, +1.100, -1.028) | 0.894372–1.058326 | 0.8853305, 0.8699802 |
| `62ad760` | Q1/P0 local pressure | 2×24×32, `dt=0.002 s` | 500/500 | 0.482 s | (-7.886, +2.269, +12.725) | (-10.525, +1.039, +13.704) | 0.739792–1.197484 | 1.0742661, 1.1626890 |
| `e07993b` | Q1/P0 local pressure | 2×24×32, `dt=0.004 s` | 250/250 | 0.480 s | (-7.700, +2.261, +12.715) | (-10.358, +1.045, +13.688) | 0.741052–1.196233 | 1.0741823, 1.1607469 |
| `e07993b` | Q1/P0 local pressure | 2×36×48, `dt=0.002 s` | 500/500 | 0.482 s | (-36.229, +2.037, -14.009) | (-37.018, +0.985, -15.030) | 0.523357–1.475752 | 0.5511448, 0.6538962 |
| `62ad760` | F-bar | 2×24×32, `dt=0.002 s` | 500/500 | 0.482 s | (-12.230, +2.689, +13.088) | (-14.457, +0.767, +14.332) | 0.699566–1.185359 | 1.0709970, 1.2200494 |
| `62ad760` | F-bar | 2×36×48, `dt=0.002 s` | 500/500 | 0.482 s | (-36.956, +1.415, -14.022) | (-37.205, +1.089, -14.517) | 0.589525–1.316852 | 0.5953528, 0.6971874 |

For Q1/P0 on 2×12×16, the retained condensed element pressure ranged from
-15,705.99 Pa to -2,431.77 Pa at `dt=0.004 s` and from -15,707.11 Pa to
-2,428.72 Pa at `dt=0.002 s`. On 2×24×32, the 1,536 values ranged from
-44,954.29 Pa to -3,965.68 Pa at `dt=0.002 s` and from -44,441.69 Pa to
-4,043.97 Pa at `dt=0.004 s`. On the 2×36×48 Q1/P0 run, 3,456 values ranged
from -160,254.02 Pa to +11,309.58 Pa, with mean -21,286.60 Pa. All seven Case
B runs used persistent PETSc SNES, and their reports contain one independently
checked acceptance record for every requested increment. “Accepted steps”
describes those solver records; it does not mean the curves are converged with
respect to mesh or time step.

The `e07993b` 2×24×32, `dt=0.004 s` run exercised the dedicated invalid-trial
recovery. PETSc rejected 46 Q1/P0 residual evaluations: 44 during step 132
(`t=0.528 s`) and two during step 133 (`t=0.532 s`). Backtracking then found
valid trials, and all 250 accepted/final states met the unchanged residual and
domain checks. The rejected evaluations were not committed states. This record
shows that the recovery path operated in this execution; it is not a tolerance
change, an accuracy result, or evidence of mesh/time convergence.

The `e07993b` 2×36×48, `dt=0.002 s` run reached peak index 241 at
`t=0.482 s`; its peak Gauss-point `det(F)` mean was 0.985462950484. It rejected
168 invalid residual trials: 83 at step 277 (`t=0.554 s`) and 85 at step 279
(`t=0.558 s`). Backtracking again found valid trials. All 500 accepted/final
states met the unchanged checks; the largest
final-residual-to-acceptance-threshold ratio was 0.973921, below one. This is a
retained solver diagnostic, not a relaxed tolerance or an accuracy,
convergence, validation, or bifurcation result.

For the 2×12×16 Q1/P0 time-step pair, every `dt=0.004 s` loading sample is
identical to the corresponding sample from the `dt=0.002 s` run. On those 251
common times, the Euclidean vector-history difference has maximum 0.144513 mm
and RMS 0.024943 mm at `p0`; at `p1` the maximum is 0.115627 mm and RMS is
0.020255 mm. These are bounded two-step sensitivity measurements, not an
asymptotic time-convergence result.

Across all seven retained Case B configurations, CoupFE RED ranges from
0.5511448 to 1.0742661 at `p0` and from 0.6538962 to 1.2200494 at `p1`. Each is
above the maximum participant-team RED retained in the fine report,
0.2875376036 at `p0` and 0.3701682113 at `p1`. Representative Q1/P0 values are
about 0.88/0.87 on 2×12×16, about 1.07/1.16 on 2×24×32, and 0.551/0.654 on
2×36×48. The close coarse time-step histories therefore coexist with
substantial spatial sensitivity in signed point components. RED against the
all-team mean combines amplitude and mesh-dependent signed point-component
differences; it is not an amplitude-only error measure.

For a signed-component-aware check, the canonical arrays already retained in
[`case_b_local_pressure_2x36x48_dt0p002.report.json`](archive/truncated_polar/case_b/case_b_local_pressure_2x36x48_dt0p002.report.json)
give component-history correlations with the all-team mean of
(0.939462, 0.837752, 0.969152) at `p0` and
(0.944736, 0.825150, 0.964841) at `p1`. Component RMSE is
(5.432790, 0.358965, 2.161431) mm at `p0` and
(5.613262, 0.244329, 2.245171) mm at `p1`. At the common `t=0.48 s` sample,
CoupFE `p0` is (-36.078865, +2.027750, -14.014764) mm and the all-team mean is
(-24.009661, +1.554442, -16.480099) mm; the exact Zenodo participant curves identified
by the report all have signs (-, +, -) there. Thus the 2×36×48 configuration
has similar curve shapes and signs but a material amplitude difference. The
2×24×32 `p0.u_z` history instead has correlation -0.814 and the opposite sign
at peak loading, while the 2×12×16 `p0.u_x` correlation is negative. These are
descriptive post hoc diagnostics, not a validation criterion or bifurcation
test.

The checkpoint-`62ad760` 2×24×32, `dt=0.002 s` pair holds mesh, time step,
geometry, fiber convention, output sampler, source, and nonlinear solver fixed.
Relative to F-bar, Q1/P0 changes the peak
`p0` vector by (+4.345, -0.420, -0.362) mm and `p1` by
(+3.932, +0.272, -0.628) mm, with vector magnitudes 4.380 and 3.991 mm. The
maximum vector-history differences are 4.382 and 3.991 mm. Q1/P0 RED is 0.00327
higher for `p0` and 0.05736 lower for `p1`. These mixed component and RED
differences are reported without declaring either formulation more accurate.

The two Q1/P0 2×24×32 records have nested `dt=0.004 s` and `dt=0.002 s`
time grids. On their 251 common times, the maximum/RMS vector-history
differences are 0.519200/0.132138 mm at `p0` and 0.549237/0.138402 mm at `p1`.
They were produced at different application checkpoints, so this is a bounded
cross-checkpoint sensitivity record, not a controlled one-factor time-step or
temporal-convergence study.

The F-bar 2×24×32 peak has the same component signs as the legacy-reported
2×24×32 observation. The source, open-apex treatment, output sampler, and
retained evidence differ, so this is qualitative consistency rather than a
claim of reproduction or equivalence. Its positive `p0.u_z` is opposite the
negative value in the recorded benchmark comparison row.

The archived F-bar 2×36×48 peak `p0` is within 1.11593 mm, or 2.89%, of the
rounded legacy-reported (-35.86, +1.49, -14.22) mm vector. This is not an exact
or byte-for-byte reproduction because the legacy raw archive, log, and
complete configuration evidence are absent, and the archived run uses an open
apex and the Hex8 output sampler. This 2×36×48 vector also has the same
component signs as the recorded benchmark comparison row. Thus the retained
2×24×32 and 2×36×48 rows have different mesh-dependent signed point-component
responses; the evidence does not establish a bifurcation, global twist
direction, or solution branch.

The archived Q1/P0 and F-bar 2×36×48 rows share the nominal mesh, time
step, geometry, fiber convention, output sampler, and solver family, but their
application checkpoints and invalid-trial policies differ. Their peak-vector
differences are 0.956714 mm at `p0` and 0.555621 mm at `p1`; maximum
vector-history differences are 9.336698 and 8.242373 mm. Q1/P0 RED is
0.0442080 lower at `p0` and 0.0432912 lower at `p1`. These are side-by-side,
source-identified output differences, not a controlled one-factor formulation
comparison or an accuracy ranking. The Q1/P0 row is not a reproduction of the
legacy F-bar observation.

Between the archived F-bar 2×24×32 and 2×36×48 records at the same time step,
the peak-vector differences are 36.713 mm at `p0` and 36.740 mm at `p1`; the
maximum vector-history differences are 50.677 and 48.298 mm. RED changes from
1.0709970 to 0.5953528 at `p0` and from 1.2200494 to 0.6971874 at `p1`. This is
clear spatial sensitivity between two meshes, not evidence of a converged mesh
sequence or an accuracy ranking.

RED is reported as the benchmark defines it. No repository-defined pass/fail
threshold is assigned. The checkpoint-`62ad760`, `dt=0.002 s` matched pair
isolates the volumetric formulation only at one discretization. The two
same-checkpoint 2×12×16 time steps, the cross-checkpoint 2×24×32 time steps,
and the two F-bar spatial resolutions are bounded sensitivity records, not
convergence studies. The current tip-refined full-cycle comparison now closes
the earlier execution gap for a controlled `2x20x17`/`4x20x17` pair. Their
near-coincident pre-snap and late-cycle histories, snap-window separation, and
mixed phase-specific FEniCS metrics support numerical timing sensitivity but
do not establish a spatial limit. Controlled mesh-ladder, temporal,
formulation, and rank
comparisons on the current closed mesh, plus real-device or clinical
validation, remain open. Signed changes between
configurations do not establish a bifurcation, bistability, or unique twist
direction.

The serial driver, archived NPZ identities, and public logs remain unchanged as
historical example and regression evidence; the v2 JSON reports make the
reporting-only reference correction explicit. On 2026-08-01, a detached
clean-checkpoint-`e07993b`
Q1/P0 3×36×48, `dt=0.001 s` rerun was stopped by user direction after printing
accepted steps 1 and 2 so development could prioritize an MPI companion. It
produced no completed NPZ or report. The external partial log has SHA-256
`d47ee0ce04a312db3808a2a3b373bb02a0c7f71ff5192f9dd88999e53b50554a`.
That interruption is neither a numerical failure nor a retained result. MPI
rank equivalence is a separate claim: a distributed companion must establish
same-configuration agreement before it can support rank-independence evidence,
while a named-rank run may still provide configuration-specific diagnostics.

Configuration-equivalent driver commands are:

```bash
python examples/cardiac_benchmark/run.py \
  --case A --formulation fbar --integrator be \
  --nonlinear-solver core-newton \
  --nt 1 --nmu 2 --ntheta 4 --apex-offset 0.2 \
  --mass lumped --fiber-sampling cg1 \
  --dt 0.002 --tend 1.0 \
  --out case_a_fbar_1x2x4_dt0p002.npz

python examples/cardiac_benchmark/run.py \
  --case B --formulation local-pressure --integrator be \
  --nonlinear-solver petsc-snes \
  --nt 2 --nmu 12 --ntheta 16 --apex-offset 0.2 \
  --mass lumped --fiber-sampling cg1 \
  --dt 0.004 --tend 1.0 \
  --out case_b_local_pressure_2x12x16_dt0p004.npz

python examples/cardiac_benchmark/run.py \
  --case B --formulation local-pressure --integrator be \
  --nonlinear-solver petsc-snes \
  --nt 2 --nmu 12 --ntheta 16 --apex-offset 0.2 \
  --mass lumped --fiber-sampling cg1 \
  --dt 0.002 --tend 1.0 \
  --out case_b_local_pressure_2x12x16_dt0p002.npz

python examples/cardiac_benchmark/run.py \
  --case B --formulation local-pressure --integrator be \
  --nonlinear-solver petsc-snes \
  --nt 2 --nmu 24 --ntheta 32 --apex-offset 0.2 \
  --mass lumped --fiber-sampling cg1 \
  --dt 0.002 --tend 1.0 \
  --out case_b_local_pressure_2x24x32_dt0p002.npz

python examples/cardiac_benchmark/run.py \
  --case B --formulation local-pressure --integrator be \
  --nonlinear-solver petsc-snes \
  --nt 2 --nmu 24 --ntheta 32 --apex-offset 0.2 \
  --mass lumped --fiber-sampling cg1 \
  --dt 0.004 --tend 1.0 \
  --out case_b_local_pressure_2x24x32_dt0p004.npz

python examples/cardiac_benchmark/run.py \
  --case B --formulation fbar --integrator be \
  --nonlinear-solver petsc-snes \
  --nt 2 --nmu 24 --ntheta 32 --apex-offset 0.2 \
  --mass lumped --fiber-sampling cg1 \
  --dt 0.002 --tend 1.0 \
  --out case_b_fbar_2x24x32_dt0p002.npz

python examples/cardiac_benchmark/run.py \
  --case B --formulation fbar --integrator be \
  --nonlinear-solver petsc-snes \
  --nt 2 --nmu 36 --ntheta 48 --apex-offset 0.2 \
  --mass lumped --fiber-sampling cg1 \
  --dt 0.002 --tend 1.0 \
  --out case_b_fbar_2x36x48_dt0p002.npz

python examples/cardiac_benchmark/run.py \
  --case B --formulation local-pressure --integrator be \
  --nonlinear-solver petsc-snes \
  --nt 2 --nmu 36 --ntheta 48 --apex-offset 0.2 \
  --mass lumped --fiber-sampling cg1 \
  --dt 0.002 --tend 1.0 \
  --out case_b_local_pressure_2x36x48_dt0p002.npz
```

Use `post.py` with `--json` and `--run-log` plus a verified external reference
directory to recreate a retained report. Exact build-directory names and
process placement are operational details; they are not numerical parameters
or performance claims.

## Historical reduced Case A record

The original retained Case A run remains unchanged as
[`case_a_reduced.json`](archive/truncated_polar/case_a/case_a_reduced.json) and
[`case_a_reduced_stdout.txt`](archive/truncated_polar/case_a/case_a_reduced_stdout.txt).
It used application
revision `44cbfed9e09d4150203faae3087f2e4617d1fc47`, the same Core revision
`454f73ce2de284262b214a2b37bd676c6aca3c0a`, 24 nodes, eight Hex8 elements,
backward Euler with `dt=0.002 s`, and 500/500 requested steps. Its output points
were sampled with the earlier global Delaunay-tetra policy.

At its recorded peak, `p0` was (-20.9006, -2.6377, -0.1411) mm and `p1` was
(-14.3984, -1.9592, +0.1221) mm. Those values remain useful implementation
history. They are not relabeled as current Hex8-sampler output and are not a
paper-curve comparison.

The model citation, external dataset identity, and reproducibility boundary
are documented in the [benchmark comparison
guide](../../../docs/BENCHMARK_COMPARISON.md), [Case B status
record](../../../docs/CASE_B_STATUS.md), and [example reference
map](../../REFERENCES.md).

Repository-authored evidence records are licensed [CC BY
4.0](../../../docs/LICENSE.md). Archived report JSON also retains derived
all-team mean curves and RED values from the cited CC BY 4.0 benchmark dataset;
the raw archive and team pickle files are not redistributed.
