# Closed Case A stopping record

Status snapshot: 2026-08-04; retained-publication update: 2026-08-05.
Campaign state: stopped by user direction.

This record closes the current Benchmark 1, Step 0A investigation. No further
Case A mesh refinement or production simulation is scheduled. The completed
fine generalized-alpha local-pressure trajectory remains the selected CoupFE
result. Stopping here is a practical decision; it is not a claim of exact
agreement, formal spatial or temporal convergence, or a uniquely identified
cause for the remaining displacement difference.

Here, "selected" means the scientific result of the closed Case A campaign.
Its compact 101-point comparison report and six-panel figure are now checked
in, while the full NPZ remains external and hash-bound. The checked-in
eight-element F-bar reports remain release/regression records for their own
historical configurations; they were not overwritten or relabeled as the
selected fine closed-geometry trajectory.

## Retained result

The selected trajectory uses the corrected closed five-block geometry and
benchmark boundary roles: zero cavity pressure, normal-only epicardial Robin
support, and full-vector base Robin support. Its key configuration is:

| Item | Retained value |
|---|---|
| Benchmark | Benchmark 1, Step 0A: active contraction with zero cavity pressure |
| Mesh | `4x36x32`; 23,616 Hex8 cells, 29,885 nodes, 89,655 displacement DOFs |
| Spatial method | Q1 Hex8 displacement with one condensed P0 pressure per element |
| Volume law | `m_e=<log(J)>_e`, `p_e=K*m_e`, `K=1.0e6 Pa`; pointwise material `kappa=0` |
| Material fields | Laplace transmural field, GP-direct fibers, `eta=100 Pa s` |
| Mass and time | Consistent Q1 mass; generalized alpha; `dt=0.001 s`; 1,000 steps |
| Time parameters | `alpha_m=0.2`, `alpha_f=0.4`, `gamma=0.7`, `beta=0.36` |
| Parallel solver | 8 MPI ranks, one thread per rank, FGMRES/GAMG profile |
| Completion | 1,000/1,000 steps; 4,149.6 s solve elapsed; zero deformation-domain rejections |
| Peak-state `det(F)` | `0.696001` to `1.129621` at the retained diagnostic state |

The clean source identities are application
`016a4f9eec6f2a4c74d10c734ddff3e24cf343de` and CoupFE Core
`e2f42ed5772850a0a23a2ce434f430c287eae5c8`. The external result archive is
`caseA_ga_local_pressure_rank8_t100.npz`, SHA-256
`ba9b31ec533398be1f39fc9a898e72f77d9587c90f9b7d9e00ce91e4d2ae6a6c`.

A loaded coarse-mesh generalized-alpha gate also passed on 1, 2, and 4 MPI
ranks. The largest `p0`/`p1` history differences from the 4-rank archive were
`5.20e-18 m` and `3.47e-18 m`. This establishes rank invariance to roundoff
for that named direct-solver prefix; it does not establish spatial convergence
or turn the fine 8-rank iterative trajectory into the same solver-profile
experiment.

The repository retains the compact self-contained derived report
`case_a_local_pressure_4x36x32_dt0p001.report.json` (62,674 bytes, SHA-256
`bbd26f3b30819ff2b67ffb48c9ad52cc9825c7fa0486e3984673c9e349bf82b1`).
It preserves the external NPZ, 4,080-byte campaign manifest, and 19,467-byte
stdout identities; the raw stdout is not distributed because it contains
machine-local host/path text. It also retains the clean source identities,
closed-mesh/method configuration, exact physical-point sampler metadata,
official ten-team file manifest, and the CoupFE and all-team-mean curves on the
101-point canonical grid.

The NPZ was produced before explicit `BENCHMARK_ARCHIVE_FIELDS` were added. Its
Step 0A identity is therefore labeled `legacy-inferred`, not recorded: the
reviewed inference uses the archive's `case=A`, the source-identified closed
Case A generalized-alpha implementation, the canonical active-tension cycle,
and identically zero cavity pressure.

## Comparison with official Simula

The comparison uses all three displacement components on the canonical,
resampled `0:0.01:1 s` comparison grid. Relative L2 is one norm of the complete
history, not the maximum percentage error at an individual time.

| Physical point | Vector RMSE | Full-history relative L2 | Maximum vector gap |
|---|---:|---:|---:|
| `p0` | 1.51128 mm | 8.56568% | 3.53428 mm at 0.66 s |
| `p1` | 1.55098 mm | 12.26311% | 3.77550 mm at 0.67 s |

The result follows the overall contraction and recovery response. Agreement
is close to the project's approximate 10% working level at `p0` and modestly
outside it at `p1`. The largest vector differences are temporary and occur
during unloading. The benchmark teams do not give identical curves, so these
numbers are a named comparison with Simula rather than a pass/fail statement
about one universal discrete answer.

The official Simula displacement file is
`monoventricular_nonblinded_step_0A_group_simula.pickle`, SHA-256
`e5804d9f7bb9f99690512a55fbab30c45cffcf3fa812e014113fb4779f386575`.

The public figure instead uses the official ten-team mean. Against that mean,
the benchmark-paper Eq. 21 RED is 0.3337402 at `p0` and 0.5024615 at `p1`.
Those values should not be confused with the Simula full-history relative L2
percentages above: RED averages pointwise relative vector errors, while
relative L2 divides norms accumulated over the whole history.

## What interpolation can and cannot explain

The benchmark landmarks are physical coordinates, not nearby mesh nodes:

```text
p0 = (0.025, 0.030, 0.000) m
p1 = (0.000, 0.030, 0.000) m
```

The CoupFE sampler inverts the reference trilinear Hex8 map and archives the
selected element, natural coordinates, and eight shape weights. On the
retained fine mesh, the coordinate-reconstruction errors are
`4.13e-20 m` at `p0` and `1.04e-17 m` at `p1`. Both points lie on shared
element faces; candidate elements give the same continuous Q1 value to
roundoff and the selected element is deterministic. There is therefore no
evidence of a wrong landmark, nearest-node substitution, or element-choice
discontinuity.

The displacement returned at each landmark is the trilinear interpolation of
the computed Q1 nodal displacement. Evaluating those shape functions gives the
discrete Hex8 field value up to roundoff; no additional point-sampler error has
been identified. The Q1 field is still only an approximation to the unknown
continuum displacement, while Simula uses quadratic P2 tetrahedral
displacement. Consequently, the measured CoupFE-versus-Simula difference can
include:

- global spatial-discretization error in the two solved fields;
- Q1-versus-P2 approximation-space and element differences at the landmarks;
  and
- formulation and discrete structural-field differences.

It is therefore reasonable that what is informally called interpolation error
-- the Q1 spatial-approximation error relative to the continuum field --
contributes to the remaining gap. The current data do not isolate its
magnitude, so the gap must not be described as a point-sampling bug or as being
caused solely by interpolation.

## Mesh evidence and claim boundary

The retained local-pressure generalized-alpha mesh prefixes reach only
0.20 s on the coarser meshes. The Simula metrics in this table use the common
0.01 s comparison grid:

| Mesh | Ranks / linear solver | `p0` relative L2 vs Simula | `p1` relative L2 vs Simula |
|---|---|---:|---:|
| `1x12x8` | 4 / direct | 10.879% | 5.051% |
| `2x20x17` | 8 / iterative | 4.386% | 6.956% |
| `4x36x32` | 8 / iterative | 7.103% | 10.033% |

Only `2x20x17` to `4x36x32` is a strict controlled pair. Over 0--0.20 s,
their trajectories differ by 5.331% at `p0` and 2.821% at `p1` relative to
the finer trajectory on their native 0.001 s grid, with vector RMSE of
0.1512 mm and 0.0564 mm. This controlled refinement moved both Simula prefix
metrics in the unfavorable direction. The retained levels do not show
monotonic improvement toward Simula and cannot support a convergence rate.

No coarser local-pressure Case A archive reaches the decisive 0.65--0.67 s
unloading interval. The present evidence therefore does not determine how the
largest gap changes with mesh size. The external mesh-evidence record is
`mesh_evidence.json`, SHA-256
`2d70558023f4a51d42a6afcdff02a33c536f113708ef5f8bb7f625daf3178e8f`;
its technical report is `mesh_refinement_decision.html`, SHA-256
`2e11fce8bb35f459ab2fe96a075728fa5a39581052ea194ecd84f9938dd77f64`.

## Controlled volume-law result

Changing only the condensed scalar volume response from `p_e=K*m_e` to
`p_e=K*(exp(2*m_e)-1)/2`, with its exact tangent, did not improve the fine
trajectory through 0.70 s. The Simula relative L2 values changed from 8.273%
to 8.348% at `p0` and from 11.480% to 11.578% at `p1`; the maximum-gap region
was unchanged to less than 0.001 mm. This is a retained null result, not a
reason to tune the volume law further.

The candidate archive is `caseA_ga_local_pressure_paper_rank8_t070.npz`,
SHA-256
`a3c19a6bfa040a6780466eae20a86ee653982015dc7d1fd049baf9ca1a2d21a9`.

## Closeout decision

Case A is closed at this point with the following boundaries:

- retain the full fine generalized-alpha log-law result as the selected
  CoupFE trajectory;
- retain the paper-law change as a controlled null experiment;
- do not claim spatial or temporal convergence or assign the remaining gap to
  one cause;
- do not launch the proposed `3x28x24` or `5x44x40` refinement now; and
- reopen Case A only for a defined new need, such as a same-approximation-space
  cross-code test, a changed element/interpolation method, or a formal mesh-
  convergence requirement.

This preserves useful benchmark-comparison evidence for the new
CoupFE-Cardiac package without spending more time on hypotheses that the
existing evidence does not support.

The cross-case answer-first
[`benchmark reproduction status`](BENCHMARK_REPRODUCTION_STATUS.md) classifies
this retained trajectory as historical approximate trajectory-agreement
evidence and records why it is not a current-setup reproduction or convergence
claim.
