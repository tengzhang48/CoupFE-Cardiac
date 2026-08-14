# Historical truncated-polar comparison figures

These figures preserve the presentations that were generated from the old
open-tip `polar_ring`, `apex_offset=0.2` mesh. That domain is not the closed
Benchmark 1 ventricular geometry, so these files are historical artifacts and
must not be cited as current benchmark-validation evidence.

- [`case_a_comparison.svg`](case_a_comparison.svg) is bound to the archived
  corrected-law eight-element F-bar report.
- [`case_b_comparison.svg`](case_b_comparison.svg) is bound to the archived
  2×36×48 Q1/P0 local-pressure report.

The corresponding reports and raw logs live in
[`examples/cardiac_benchmark/results/archive/truncated_polar/`](../../../../examples/cardiac_benchmark/results/archive/truncated_polar/).

**Resolved 2026-08-10.** A byte-identical copy of `case_b_comparison.svg`
(SHA-256 `1814d84f...`) previously also sat at the active path
`docs/figures/case_b_comparison.svg`, inlined on the README front page between two
closed-geometry figures. Case A had been migrated to the closed geometry and Case B
had not, so the two active files looked like a matched pair and were not one. That
duplicate has been removed; every reference now points here, and
`docs/figures/` contains only closed-geometry results.

These archived figures are **comparison and lessons material only**. Note that
their bytes carry no open-tip marker — the Case B `<title>` reads
"Benchmark 1, Case B — displacement comparison" — so extracted or screen-read on
their own they are indistinguishable from current benchmark results. Always
reproduce the surrounding caption when reusing them.
