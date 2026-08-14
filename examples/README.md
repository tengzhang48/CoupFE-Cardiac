# Code, tests, examples, and retained outputs

The repository contains runnable source applications, focused checks, and small
text/JSON result records so readers can inspect both the implementation and
what a documented run produced.

| Area | Runnable code | Tests | Retained output |
|---|---|---|---|
| Cardiac Benchmark 1, Case A | serial F-bar and closed distributed Q1/P0/generalized-alpha paths in [`cardiac_benchmark/`](cardiac_benchmark/) | fast component checks, compiled `slow` run, and opt-in MPI gates | [selected pre-straight-wall closed report and archived truncated-polar history](cardiac_benchmark/results/) |
| Cardiac Benchmark 1, Case B | serial and distributed F-bar, pointwise-`kappa`, and application-owned Q1/P0 local-pressure paths on explicitly named historical and closed Hex8 topologies in [`cardiac_benchmark/`](cardiac_benchmark/) | geometry/boundary, load, tangent, mass, local-pressure, determinant-domain step/recovery, Hex8-sampling, result, and solver checks under `tests/` | [archived truncated-polar reports and qualified closed-domain development status](../docs/CASE_B_STATUS.md) |
| Distributed implementation | five scripts under [`mpi_smoke/`](mpi_smoke/) | opt-in PETSc/MPI tests at 1, 2, and 4 ranks | [five-script rank-comparison record](mpi_smoke/README.md#retained-release-configuration-output) |
| Scientific and software provenance | [`REFERENCES.md`](REFERENCES.md) | source, license, and comparison-contract checks | citations, external-data identities, and stated numerical boundaries |

The truncated-polar archive preserves the reduced Case A history and both
source-identified 500-step reports. These are exact-configuration regression
records for the historical open-tip path, not current Benchmark 1 geometry or
validation evidence. The reduced record uses a coarse open-apex mesh and the
Delaunay-tetra output sampler. The corrected-law regression record was rerun at
checkpoints `6839c13`/`e2f42ed` with the complete smooth-switch-energy
derivative and reference-Hex8 isoparametric sampling. The
`62ad760`/`454f73` historical-law report remains separately labeled; their
provenance is not merged.

The selected Case A and Step 2 Case B records at the result-index root are
topologically closed, but their producing sources predate the toolkit-matched
straight-wall mapping and physical-coordinate structural-frame reconstruction.
They remain qualified development/history evidence for their recorded source.
The result index also contains a separately identified, corrected-geometry
Step 0B dual-run full-cycle comparison and its earlier 0.32 s prefix
diagnostic; neither is relabeled as evidence from a different source state.

Archived truncated-polar Case B outputs include paired 2×12×16 Q1/P0 runs at
`dt=0.004 s` and `0.002 s`, plus matched Q1/P0 and F-bar runs on the 2×24×32
mesh at `dt=0.002 s`. Each has a full-history report, all-step solver
diagnostics, positive peak-load Gauss-point deformation Jacobians, result/log
hashes, and a separately verified comparison with the external curves. The
checkpoint-`62ad760` 2×24×32 pair isolates the volumetric formulation at one
named configuration. The 2×12×16 pair gives a bounded two-step
time-sensitivity check on nested loading grids. A further 2×36×48 F-bar run
retains a second spatial resolution;
its peak `p0` is within 1.11593 mm, or 2.89%, of the rounded legacy-reported
vector under an archived, different apex/sampling record. None of these
comparisons establishes accuracy or mesh/time convergence.

A newer checkpoint-`e07993b` Q1/P0 2×24×32 run at `dt=0.004 s` completed
250/250 increments. PETSc rejected 46 invalid trial residuals at steps 132–133
and backtracked to valid trials; the rejected trials were not committed, and
all accepted/final states met the unchanged checks. This is execution evidence
for the recovery path, not a tolerance, convergence, accuracy, validation, or
bifurcation claim. Its earlier `dt=0.002 s` comparison record has a different
application checkpoint, so the pair is described as cross-checkpoint
sensitivity rather than a controlled time-step study.

The same checkpoint also retains a completed 2×36×48, `dt=0.002 s` Q1/P0 run.
PETSc rejected 168 invalid trial residuals at steps 277 and 279 before finding
valid trials; all 500 accepted/final states met the unchanged checks. The
existing F-bar row at the same nominal mesh and time step was produced at
`62ad760`, so the two outputs are a source-identified side-by-side comparison,
not a controlled formulation, convergence, accuracy, reproduction,
validation, or bifurcation result.

The 2026-06-27 Case B table remains development history from the F-bar path.
Its raw archives and logs are absent, so it is not substituted for the archived
source-identified records.

No external Zenodo pickles or generated NPZ archives are committed. The JSON
reports, normalized console logs, and small historical comparison SVGs are
retained so readers can inspect what the documented runs produced. The
truncated-polar files are exact-configuration regression history, not current
Benchmark 1 validation evidence.
