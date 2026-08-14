# Benchmark reproduction status

Status date: 2026-08-07.

## Answer

**For the current Step 0B `p0`/`p1` displacement histories: partially. For the
whole paper benchmark: no.** The clean current Step 0B trajectory is a credible
partial reproduction: it captures pressure-driven snap-through and recovery
with close event timing and approximately 9% relative L2 against the ten-team
mean, but its snap amplitudes, unloading history, and small displacement
components do not agree closely enough to claim full quantitative
reproduction. The retained Step 0 Case A trajectory shows approximate
historical agreement at the project's roughly 10% working level, but it
predates the later straight-wall geometry and physical-frame corrections and
is not a current-setup reproduction.

The paper reports contributions from ten teams, not one unique discrete
curve, and defines no numerical pass threshold. Consequently, "reproduced"
must not mean byte-for-byte agreement with one FEniCS/Simula trajectory. It also
must not be inferred from one time, one component, or one landmark.

## Evidence by case

| Case | Strongest retained comparison | Status |
|---|---|---|
| Step 0 Case A, active contraction | Against official Simula on the canonical full-cycle grid, vector-history relative L2 is 8.57% at `p0` and 12.26% at `p1`; vector RMSE is 1.51/1.55 mm. The largest gaps, 3.53/3.78 mm, are temporary and occur during unloading. | **Historical approximate trajectory-agreement evidence.** The overall contraction and recovery agree at `p0` near and `p1` modestly outside the project's approximate 10% working level. The retained run predates the later straight-wall and physical-frame corrections, so it is not a current-setup reproduction or convergence result. |
| Step 0 Case B, passive pressure loading | Against the exact ten-team mean on the canonical 101-point grid, relative L2 is 8.98%/9.20% and vector RMSE is 1.19/1.25 mm at `p0`/`p1`. Against the local hash-pinned FEniCS reference, shared-history RMSE is 1.017/1.073 mm; snap and recovery crossings are 2.64--2.66 ms and 7.43--7.78 ms late, with 4.19/3.75 mm maximum gaps near snap. | **Partial quantitative reproduction.** The overall full-history comparison reaches the project's approximate 10% working level and the correct qualitative snap/recovery is present, but transient amplitude, relaxation, and small-component agreement remain mixed. The current Q1/P0 mean-`log(J)` volume response is not the paper's pointwise scalar volumetric law. |
| Step 2 Case B, combined activation and pressure | The source-bound development comparison has 9.80% global full-history relative L2 against the all-team mean, but its `p1-z` plateau has the opposite sign from all ten official curves; its producing content predates later geometry/frame corrections. A later corrected-setup Q1/P0 run reported a promising trajectory, but its compact comparison record does not bind the complete source, execution, solver/deformation, and ten-team role/hash provenance required for release-grade evidence. | **Development and provenance-incomplete diagnostic evidence only; no reproduction claim.** |

The complete Case A evidence and its limits are in
[`CASE_A_STATUS.md`](CASE_A_STATUS.md). The current Step 0B numbers come from
the source-bound
[`step0b_tip6p0_full_cycle_comparison.report.json`](../examples/cardiac_benchmark/results/step0b_tip6p0_full_cycle_comparison.report.json)
and its
[`six-component comparison figure`](figures/step0b_tip_refine_full_cycle.svg).

The Step 0B ten-team-mean values above were independently recomputed from the
report's exact hash-gated ten-file manifest on `t=0:0.01:1 s`. Relative L2 is
`||U-M||/||M||` over the complete 101-by-3 history, where `M` is the arithmetic
team mean. The corresponding benchmark Eq. 21 RED values are 0.175895 at `p0`
and 0.311732 at `p1`; RED weights each time separately and is not
interchangeable with full-history relative L2. The loader, team-selection
policy, and file hashes are in
[`compare_tip_refine_full_cycle.py`](../examples/cardiac_benchmark/compare_tip_refine_full_cycle.py).

## Why the Step 0B verdict is partial

The current clean Step 0B run establishes all of the following:

- the paper's pressure-only physical case is selected on the current closed
  straight-wall geometry;
- the benchmark boundary roles, physical points, parameters, fiber rule,
  viscosity, consistent mass, and source-matched generalized-alpha staging are
  recorded;
- the Q1/P0 local-pressure run completed 1,000/1,000 increments on eight MPI
  ranks with no accepted invalid deformation; and
- the full trajectory shows the same principal snap-through and recovery as
  the FEniCS reference.

The run intentionally uses the application-owned condensed Q1/P0
mean-`log(J)` pressure response. That is a documented near-incompressibility
variant, not the paper's pointwise scalar volumetric response. It is therefore
method-comparison evidence rather than exact discrete-equation identity with
the paper or FEniCS implementation. This distinction bounds the claim; the
current Step 0B evidence does not isolate the volume law as the cause of the
remaining gap. Separately, the paper-law change was a null only through 0.70 s
on one named fine, pre-correction Case A configuration and must not be
generalized to Step 0B.

The approximately 9% relative-L2 result establishes useful overall trajectory
agreement under the project's working criterion, but it does not establish
exact component-wise agreement. The four-layer mesh
reduces full-shared-history FEniCS RMSE by 6.9%/8.2% and the maximum gap by
24.6%/25.3% relative to the two-layer mesh, but relaxation RMSE becomes
16.2%/14.9% worse. The trajectory lies inside the official ten-team envelope
in all three components simultaneously at only 15.3% of timestamps for `p0`
and 26.6% for `p1`. This all-component statistic is intentionally strict and
is strongly affected by the small components; it is a diagnostic, not a
paper-defined acceptance rule. Together, the phase-dependent changes show
numerical transient sensitivity rather than monotonic convergence.

## Claim boundary

The defensible package statement is:

> For current Step 0B `p0`/`p1` displacement histories, CoupFE-Cardiac
> reproduces the displacement scale and principal snap/recovery response. It
> matches the physical geometry/load/boundary/time contract but uses the
> documented Q1/P0 mean-`log(J)` volume variant; this is a partial quantitative
> reproduction whose largest differences occur near snap-through and unloading.
> The retained Case A trajectory is historical approximate agreement evidence,
> not a current-setup reproduction.

Do not shorten this to "the paper benchmark is reproduced" or "the benchmark
is validated." No retained study yet establishes spatial and temporal
convergence of the current Step 0B method, equivalence to another team's
discretization, or experimental/clinical validation.
