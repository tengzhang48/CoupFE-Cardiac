# Truncated-polar result archive

This directory preserves reproducible records made with the historical
`polar_ring` mesh and `apex_offset=0.2`. That construction leaves a truncated,
open tip and is a **non-benchmark geometry**: these files are project history
and regression evidence for their exact configurations. They are not current Benchmark 1 validation evidence.

The 20 files below retain their original bytes and basenames. Each comparison
report binds its adjacent normalized console log by filename, size, and
SHA-256. The reports also bind the uncommitted generated NPZ result and the
external benchmark data used for comparison. Moving the records here does not
change any recorded numerical value or provenance identity.

## Case A

| Record | JSON | Console output |
|---|---|---|
| Historical material law, F-bar, 1×2×4, `dt=0.002 s` | [`case_a_fbar_1x2x4_dt0p002.report.json`](case_a/case_a_fbar_1x2x4_dt0p002.report.json) | [`case_a_fbar_1x2x4_dt0p002.raw.stdout.txt`](case_a/case_a_fbar_1x2x4_dt0p002.raw.stdout.txt) |
| Corrected complete-switch-energy derivative, same run contract | [`case_a_fbar_1x2x4_dt0p002_corrected_switch.report.json`](case_a/case_a_fbar_1x2x4_dt0p002_corrected_switch.report.json) | [`case_a_fbar_1x2x4_dt0p002_corrected_switch.raw.stdout.txt`](case_a/case_a_fbar_1x2x4_dt0p002_corrected_switch.raw.stdout.txt) |
| Historical reduced Case A execution | [`case_a_reduced.json`](case_a/case_a_reduced.json) | [`case_a_reduced_stdout.txt`](case_a/case_a_reduced_stdout.txt) |

The Historical reduced Case A record used the older global Delaunay-tetra policy
for output sampling. It is an executable-pipeline record, not a
paper-curve comparison.

## Case B

| Formulation and mesh | `dt` | JSON | Console output |
|---|---:|---|---|
| Q1/P0, 2×12×16 | 0.004 s | [`report`](case_b/case_b_local_pressure_2x12x16_dt0p004.report.json) | [`log`](case_b/case_b_local_pressure_2x12x16_dt0p004.raw.stdout.txt) |
| Q1/P0, 2×12×16 | 0.002 s | [`report`](case_b/case_b_local_pressure_2x12x16_dt0p002.report.json) | [`log`](case_b/case_b_local_pressure_2x12x16_dt0p002.raw.stdout.txt) |
| Q1/P0, 2×24×32 | 0.002 s | [`report`](case_b/case_b_local_pressure_2x24x32_dt0p002.report.json) | [`log`](case_b/case_b_local_pressure_2x24x32_dt0p002.raw.stdout.txt) |
| Q1/P0, 2×24×32 | 0.004 s | [`report`](case_b/case_b_local_pressure_2x24x32_dt0p004.report.json) | [`log`](case_b/case_b_local_pressure_2x24x32_dt0p004.raw.stdout.txt) |
| Q1/P0, 2×36×48 | 0.002 s | [`report`](case_b/case_b_local_pressure_2x36x48_dt0p002.report.json) | [`log`](case_b/case_b_local_pressure_2x36x48_dt0p002.raw.stdout.txt) |
| F-bar, 2×24×32 | 0.002 s | [`report`](case_b/case_b_fbar_2x24x32_dt0p002.report.json) | [`log`](case_b/case_b_fbar_2x24x32_dt0p002.raw.stdout.txt) |
| F-bar, 2×36×48 | 0.002 s | [`report`](case_b/case_b_fbar_2x36x48_dt0p002.report.json) | [`log`](case_b/case_b_fbar_2x36x48_dt0p002.raw.stdout.txt) |

The generated CoupFE NPZ archives are not committed. External benchmark
pickles are also not committed and remain part of the separately distributed
CC BY 4.0 dataset. The JSON and text evidence records are licensed under CC BY
4.0 as described by the repository license files.

The closed-multiblock Step 2 Case B record is intentionally not in this
archive; it remains in the parent [`results/`](../..) directory.
