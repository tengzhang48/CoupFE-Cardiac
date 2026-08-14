# Benchmark comparison figures

These figures answer one narrow question: how do selected CoupFE-Cardiac
displacement histories compare with the official benchmark curves? They expose
component-level agreement and discrepancies. The Case A figure renders the
selected fine closed-multiblock result. The generic Step 0 Case B figure is an
archived open-tip `polar_ring`, `apex_offset=0.2` record, not Benchmark 1
geometry. As of 2026-08-07, the separately labeled current Step 0B figure shows
the completed eight-rank `2x20x17` and `4x20x17`, `tip_refine=6.0` full
cycles. The figure is now a dual-run comparison with FEniCS and the published
ten-team envelope. The clean isolated four-layer replay completed and matched
the independently audited non-retained candidate to roundoff. The source-bound
Step 2 Case B figure is development evidence; a second corrected-setup Step 2
figure is retained only as a provenance-incomplete diagnostic. None is a pass
or validation claim.

The displaced open-tip Case A and Step 0 Case B figures are preserved in the
[`archive/truncated_polar/`](archive/truncated_polar/) directory with their
claim boundary and source-report links.

## Case A

- Figure: [`case_a_comparison.svg`](case_a_comparison.svg)
- Selected report:
  `examples/cardiac_benchmark/results/case_a_local_pressure_4x36x32_dt0p001.report.json`
- Retained source identity: application `016a4f9`, Core `e2f42ed`; the plotting
  script pins the full revisions, external result identity, and report digest.
- Source boundary: this run predates the toolkit-matched straight-wall mapping
  and physical-coordinate structural-frame reconstruction.
- Execution: 1,000/1,000 increments on eight MPI ranks; closed
  t/core/radial 4×36×32 Hex8; Q1/P0 condensed local pressure with the log-J
  volume law; consistent mass; source-matched generalized-alpha; `dt=0.001 s`.
- Identity boundary: the NPZ predates explicit benchmark-identity metadata.
  The report and figure therefore label Step 0A as `legacy-inferred`, based on
  the recorded Case A label, source-identified closed Case A implementation,
  canonical active load, and identically zero pressure history.
- Historical boundary: the archived
  [`6839c13`](../../examples/cardiac_benchmark/results/archive/truncated_polar/case_a/case_a_fbar_1x2x4_dt0p002_corrected_switch.report.json)
  corrected-law and
  [`62ad760`](../../examples/cardiac_benchmark/results/archive/truncated_polar/case_a/case_a_fbar_1x2x4_dt0p002.report.json)
  historical-law open-tip records remain unchanged and do not drive this
  figure.
- Data shown: `comparison.ours_on_canonical_grid_m` and
  `reference.mean_curves_m`, at p0 and p1 for x, y, and z, on the report's 101
  canonical points from 0 to 1 s.
- Encoding: solid blue is CoupFE-Cardiac; dashed charcoal is the benchmark
  all-team mean. Line pattern and labels distinguish the series without color.
- Supported takeaway: this selected fine closed Hex8 Q1/P0 run can be compared
  component by component with the Case A mean, including the visible
  discrepancies and RED 0.3337402/0.5024615 at `p0`/`p1`. The figure does not
  establish validation or convergence.

## Current Step 0 Case B tip-refined full-cycle comparison

- Figure: [`step0b_tip_refine_full_cycle.svg`](step0b_tip_refine_full_cycle.svg)
- Selected report:
  [`step0b_tip6p0_full_cycle_comparison.report.json`](../../examples/cardiac_benchmark/results/step0b_tip6p0_full_cycle_comparison.report.json)
- Execution: 1,000/1,000 increments on eight MPI ranks for both current closed
  `2x20x17` and `4x20x17` Hex8 meshes with `tip_refine=6.0`; Q1/P0 condensed
  mean-`log(J)` local pressure; consistent mass; source-matched generalized
  alpha; `eta=100 Pa s`; `dt=0.001 s`.
- Data shown: both CoupFE trajectories and local FEniCS p0/p1 displacement
  components at matched timestamps plus the exact published ten-team min-max
  envelope.
- Paired result: the four-layer/two-layer full-cycle RMSE is
  0.4076/0.3618 mm at `p0`/`p1`. The pair is almost identical before snap
  (0.0098/0.0026 mm RMSE), differs by at most 1.6331/1.4059 mm in the snap
  window, and returns to 0.0346/0.0518 mm separation at cycle end.
  Four layers improve FEniCS full-shared-history RMSE and maximum snap gap, but
  worsen relaxation RMSE by 16.2%/14.9%.
- Provenance gate: the current v2 figure uses the clean four-layer replay. The
  first dirty-tree candidate is non-retained diagnostic evidence; the replay
  matches it to roundoff and records
  **application `2458e7c`, NPZ
  `1e333b29b05f01dedce9272b32b82ce6ccfda56036c1c8f57eb395b9b4494800`,
  stdout `0810a9d0c944345464f48917c707a7d3df7187b543c74efd4b0384542b95aed7`,
  elapsed `1778.4 s (29.6 min)`**. Runtime source is
  `f8d9469a101709d11460a0803b6c031001192ccce61c472d032355b01070da05`;
  Core is `454f73c`.
- Claim boundary: this controlled pair supports a numerical transient/timing
  interpretation. It does not establish mesh convergence, rank equivalence,
  physically distinct branches, or validation.

## Historical Step 0 Case B comparison

This subsection describes the older open-tip full-cycle figure below, not the
current tip-refined closed-geometry result above.

- Figure: [`case_b_comparison.svg`](archive/truncated_polar/case_b_comparison.svg)
- Selected report:
  `examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_local_pressure_2x36x48_dt0p002.report.json`
- Retained source identity: application `e07993b`, Core `454f73c`; the plotting
  script pins the full revisions, result identity, and report digest.
- Source boundary: this historical open-tip run also predates the
  physical-coordinate structural-frame reconstruction.
- Data shown: `comparison.ours_on_canonical_grid_m` and
  `reference.mean_curves_m`, at p0 and p1 for x, y, and z, on the report's 101
  canonical points from 0 to 1 s.
- Encoding: solid blue is CoupFE-Cardiac; dashed charcoal is the benchmark
  all-team mean. Line pattern and labels distinguish the series without color.
- Supported takeaway: this retained fine Hex8 Q1/P0 local-pressure run can be
  compared component by component with the Case B mean, including the visible
  discrepancies and the report's RED values. The figure does not establish
  validation or convergence.

## Step 2 Case B development comparison

- Figure: [`step2_case_b_comparison.svg`](step2_case_b_comparison.svg)
- Selected report:
  [`step2_case_b_std_kappa_2x20x17_dt0p001.report.json`](../../examples/cardiac_benchmark/results/step2_case_b_std_kappa_2x20x17_dt0p001.report.json)
- Console record:
  [`step2_case_b_std_kappa_2x20x17_dt0p001.raw.stdout.txt`](../../examples/cardiac_benchmark/results/step2_case_b_std_kappa_2x20x17_dt0p001.raw.stdout.txt)
- Execution: 1,000/1,000 full-cycle increments on four MPI ranks, closed
  2×20×17 Hex8, pointwise `std-kappa`, consistent mass, source-matched
  generalized-alpha, `dt=0.001 s`, and Step 2 active stress plus ventricular
  pressure.
- Source identity: application `e9b7d90` reported a dirty tree, Core is
  `e2f42ed`, and the exact result-producing application content is identified
  by runtime-source SHA-256
  `6b96395761dd3203f0e9ffab90a77d6389dca13cdad43490a1deac95073480f1`.
- Source boundary: the producing content predates the toolkit-matched
  straight-wall mapping and physical-coordinate structural-frame
  reconstruction.
- Data shown: the checked-in report's 101 samples from 0 to 1 s for
  CoupFE-Cardiac, the official ten-team range and mean, and the named Simula
  curve, at `p0` and `p1` for x, y, and z.
- Full-cycle metrics against the all-team mean: global relative L2 9.8038%;
  `p0`/`p1` vector relative L2 9.0555%/10.9271%; aggregate RMSE 0.829875 mm;
  maximum component error 2.846225 mm; and paper Eq. 21 RED
  28.3004%/35.2774% at `p0`/`p1`.
- Interpretation: global and pointwise-vector relative L2 divide accumulated
  history norms, whereas paper RED averages the pointwise relative vector
  error over time. The metrics therefore have different values, and the paper
  defines no acceptance threshold. The `p1-z` plateau has the wrong sign
  relative to all ten official curves. This is development evidence, not a
  validation, convergence, or rank-independence claim.

## Step 2 Case B corrected-setup diagnostic

- Figure:
  [`step2b_current_rerun_comparison.svg`](step2b_current_rerun_comparison.svg)
- Compact report:
  [`step2b_current_rerun_comparison.report.json`](../../examples/cardiac_benchmark/results/step2b_current_rerun_comparison.report.json)
- Configuration described by the artifact: closed 2×20×17 Hex8 with
  `tip_refine=6.0`, Q1/P0 mean-`log(J)` local pressure, consistent mass,
  generalized-alpha, and Step 2 active stress plus pressure.
- Provenance boundary: the compact report binds the corrected NPZ and legacy
  report hashes plus the dataset DOI, but not the complete application/Core
  identities and tree states, command/environment, solver and deformation
  audits, or exact ten-team role/hash manifest.
- Status: promising provenance-incomplete diagnostic only. It is not the
  release's Step 2 reproduction evidence and must not replace the source-bound
  development record above.

The benchmark means in the source-identified reports are derived from the
separately distributed benchmark dataset, DOI
[`10.5281/zenodo.14260459`](https://doi.org/10.5281/zenodo.14260459), licensed
under CC BY 4.0. The CoupFE-Cardiac reports retain the source identities and
team-file hashes used to form those means.

Regenerate the Case A and Step 0 Case B figures from any working directory; no
external pickle files are needed. Their canonical checked-in bytes use Python
3.10.8, Matplotlib 3.10.9, FreeType 2.6.1, and the bundled DejaVu Sans file
with SHA-256
`3fdf69cabf06049ea70a00b5919340e2ce1e6d02b0cc3c4b44fb6801bd1e0d22`:

```bash
python examples/cardiac_benchmark/plot_retained_comparisons.py --canonical
```

Regenerate the Step 2 development figure from its checked-in report with:

```bash
python examples/cardiac_benchmark/plot_step2b_case_b.py --canonical
```

Its canonical bytes use Python 3.12.3, Matplotlib 3.11.1, FreeType 2.14.3,
and the same pinned DejaVu Sans file.

Without `--canonical`, another supported Matplotlib stack renders the same
validated report data, but library-dependent text coordinates or internal SVG
IDs can change the file hash.
