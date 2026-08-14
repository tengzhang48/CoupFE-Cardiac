# Example references and provenance

This map links each public example to its scientific source, local checks,
external comparison data, and reproducibility limits. A citation defines a
model or method; it does not by itself establish that a particular numerical
result is correct.

## Primary benchmark

Reidmen Aróstica et al., “A software benchmark for cardiac elastodynamics,”
*Computer Methods in Applied Mechanics and Engineering* **435** (2025), 117485,
<https://doi.org/10.1016/j.cma.2024.117485>.

The repository implements parts of Benchmark 1, Step 0:

| Files | Benchmark source | Local implementation and check |
|---|---|---|
| `cardiac_benchmark/geometry.py` and `boundary_audit.py` | geometry, boundaries, output points, and fiber formulas in Eqs. 10–16 | preserves the historical truncated polar mesh and provides a noncollapsed five-block closed mesh with the toolkit-matched straight Cartesian wall ruling; checks exact mesh identity, exterior-face ownership, the declared historical free tip, extended Jacobians, measures, signed pressure resultant/moment, Robin matrices, and rigid-mode Robin reductions before integration |
| `cardiac_benchmark/tbar_laplace.py` | Laplace transmural coordinate in Eq. 13 | solves the Q1 Hex8 scalar problem with endocardial/epicardial Dirichlet values and natural base flux; checks bounds, boundary values, Jacobians, and linear residual; the dynamics driver verifies the generated NPY and native sidecar before use and records both hashes |
| `cardiac_benchmark/material.py` | passive, active, and viscous laws in Eqs. 2–4 and Tables 1–3 | compiled F-bar or standard Q1 material kernel, application-owned Gauss-point Gram–Schmidt, constitutive-tangent checks, and a zero-reference-stress check |
| `cardiac_benchmark/structural_directions.py` | complete fiber, sheet, and sheet-normal construction in Eqs. 14–16 and the pinned toolkit implementation | reconstructs the local coordinates from each physical point and Laplace value, including the toolkit apex branch; closed serial and MPI paths share this implementation and record its reconstruction identity |
| `cardiac_benchmark/activation.py` | activation and pressure histories in Eqs. 5–8 and Tables 3–4 | CC-BY-4.0 adaptation of the exact Finsberg/Sundnes/van den Brink Zenodo v1.0.0 deposit; the file header indicates CoupFE-Cardiac's changes, and peak-scale and finite-history checks are local tests |
| `cardiac_benchmark/robin.py` and `pressure.py` | boundary terms in Eq. 1 | analytic spring reaction, rigid-translation and long-axis-rotation reductions, orientation reversal, scalar/batched equality, and finite-difference follower-load tangent checks; rigid-mode values diagnose discrete facet geometry without changing the boundary condition or defining an accuracy threshold |
| `cardiac_benchmark/newmark.py` | Newmark's method | optional constant-average-acceleration inertia and Robin velocity map; material viscosity remains a backward strain difference |
| `cardiac_benchmark/consistent_mass.py` | consistent-mass convention used by the local Simula/FEniCS reference | standard 2×2×2 Q1 Hex8 consistent mass with symmetry, total-mass, row-sum, residual/tangent, and state-commit checks; the historical row-summed driver path remains selectable |
| `cardiac_benchmark/distributed_mass.py` and `distributed_material.py` | distributed operators used by separately identified closed formulations | assemble complete consistent-mass rows from every Hex8 touching the rank-owned rows and evaluate the selected generated material on rank-local elements; archive implementation identities distinguish pointwise-`kappa` and condensed local-pressure paths, and serial/partition component checks precede any claimed rank gate |
| `cardiac_benchmark/local_pressure.py` | application discretization experiment, not a formula copied from the benchmark implementation | one P0 pressure per Q1 Hex8 is algebraically eliminated; the `local-pressure` CLI variant uses the volume-average of `log(det(F))`, while the separately identified paper-law variant applies its scalar response to the element mean log-volume; affine, isochoric, symmetry, finite-difference tangent, inverted-deformation, degenerate-reference, and Core `max_step` controls check the operators |
| `cardiac_benchmark/sampling.py` | material-point output at benchmark points `p0` and `p1` | checked inverse of the reference Hex8 trilinear map followed by the eight Hex8 shape weights; outside, degenerate, distorted, and shared-boundary controls exercise point location and interpolation |
| `cardiac_benchmark/solver.py` | application nonlinear acceptance policy | guarded Core Newton and optional persistent PETSc SNES; only `InvalidDeformationError` in a PETSc residual trial becomes IEEE positive infinity for `bt`, rejection diagnostics are serialized, and the final residual is reassembled before any physical-state commit |
| `cardiac_benchmark/viscous_evidence.py` | retained historical forensic implementation, not a benchmark path | documents how an accepted-state eta split was guarded during debugging; the public driver no longer invokes it because the reproduction keeps the paper physical parameters fixed and runs no parameter sweeps |
| `cardiac_benchmark/run.py` | composition of Benchmark 1, Step 0, Cases A and B | selects F-bar, application-owned Q1/P0, or paper pointwise-kappa Q1; consistent or historical lumped mass; fiber/tbar policy; backward Euler or experimental Newmark; guarded Core Newton or PETSc SNES; and writes only completed archives with discretization and source metadata |
| `cardiac_benchmark/run_mpi.py` | distributed composition of supported Benchmark 1 cases | keeps historical backward-Euler and source-matched generalized-alpha identities separate; supports closed Step 0 Cases A/B and Step 2 Case B under the consistent-mass, GP-direct/Laplace contract and records formulation, source, mesh, frame, solver, and rank provenance |
| `cardiac_benchmark/result_io.py` | local result-archive contract | atomic replacement after all requested steps complete; this protects artifact integrity but is not an independent physics oracle |
| `cardiac_benchmark/post.py` | RED metric in Eq. 21 | compares a qualified completed archive with user-supplied participating-code curves and reports measurements without assigning a pass threshold |
| `cardiac_benchmark/compare_fenics_case_b.py` | direct landmark comparison with the locally completed Simula/FEniCS Case B output | reads only the explicitly supplied, hash-checked parameters/time/p0/p1 files and one qualified CoupFE archive; reports fixed-grid RMSE, relative L2, component error, and fixed snap-onset metrics without treating the P2-tet and Q1-Hex8 discretizations as equivalent |
| `cardiac_benchmark/results/archive/truncated_polar/` | repository-authored historical execution and comparison records | preserves the reduced Case A record and archived Case A F-bar, Case B F-bar, and Case B Q1/P0 reports with full histories, accepted-step diagnostics, exact source/Core and external-input identities, normalized-log hashes, and full-precision RED values; these open-tip records are not current Benchmark 1 validation evidence |

Supporting method references:

- G. A. Holzapfel and R. W. Ogden, “Constitutive modelling of passive
  myocardium: a structurally based framework for material characterization,”
  *Philosophical Transactions of the Royal Society A* **367** (2009),
  <https://doi.org/10.1098/rsta.2009.0091>.
- J. Bestel, F. Clément, and M. Sorine, “A Biomechanical Model of Muscle
  Contraction,” MICCAI 2001,
  <https://doi.org/10.1007/3-540-45468-3_143>.
- N. M. Newmark, “A Method of Computation for Structural Dynamics,” *Journal
  of the Engineering Mechanics Division* **85** (1959),
  <https://doi.org/10.1061/JMCEA3.0000098>.

## External comparison data

R. A. Arostica Barrera and Cristobal Bertoglio, “A software benchmark for
cardiac elastodynamics,” Zenodo,
<https://doi.org/10.5281/zenodo.14260459>.

- License: CC BY 4.0.
- File: `benchmark_article_data.zip`.
- Published size: 23,180,741,494 bytes.
- Zenodo checksum: `md5:75602be4777c4ca2262c2bcfd2134b15`.
- Independently computed SHA-256:
  `134951af5e38d147b0223f0a83666eb3fe1b75acb5bfa9f1b9aa30f255f8f1f5`.
- Repository policy: the archive and raw pickles are not vendored; load pickle
  files only from a verified copy of the trusted Zenodo archive.

The benchmark article and external dataset receive credit here, in `NOTICE`,
and in `THIRD_PARTY_NOTICES.md`. Checked-in reports and SVGs contain modified
CC-BY-4.0 material: CoupFE-Cardiac selects hash-pinned team files, maps curves
to documented common time grids, and computes team means, standard deviations,
minimum/maximum envelopes, selected Simula curves, metrics, and graphics. The
23.2 GB archive and raw pickle files are not redistributed, and transformed
publisher material is not presented as project-authored source data.

A retained comparison should record the archive identity, exact curve-file
identities used, loaded team labels, result hash, application and Core
revisions, environment, command, solver configuration, and complete RED output.
The archived source-identified and legacy Case B records are summarized in
[`docs/CASE_B_STATUS.md`](../docs/CASE_B_STATUS.md). The archived JSON reports
name and hash the exact `step_0B` inputs consumed; the external pickle files
and generated CoupFE NPZ archives are not redistributed.

## Source-derived formula checks

`cardiac_benchmark/fiber_crosscheck.py` and
`cardiac_benchmark/structural_directions.py` are NumPy adaptations of the
fiber/sheet expressions from `Reidmen/cardiac_benchmark_toolkit`, pinned to
commit `e8d47553cfc83eb274eba3e177de33148e7f441c`. The pinned project is
MIT-licensed; see `THIRD_PARTY_NOTICES.md`. The cross-check remains an
independent oracle; the runtime module reconstructs `(u,v)` from the physical
point and Laplace transmural coordinate without executing upstream FEniCS.

The Finsberg/Sundnes activation sources are pinned and credited in
`THIRD_PARTY_NOTICES.md`. The port conservatively uses the exact v1.0.0 Zenodo
deposit's CC BY 4.0 record and also credits Jonas van den Brink, the third
record creator. The file header states how CoupFE-Cardiac combined, renamed,
and hardened the two source modules. The archived README says MIT but its
referenced `LICENSE` file is absent, so this project does not manufacture or
rely on a downstream MIT copyright notice for that source.

## Diagnostics

These scripts and fields are implementation diagnostics, not independent
literature benchmarks:

| Script or field | Purpose | Interpretation |
|---|---|---|
| `cardiac_benchmark/diagnose.py` | reports deformation Jacobian, fiber stretch, and endocardial radial motion from a completed archive | finite values and interpreted signs; no universal pass threshold |
| `det_f_gauss_peak` and `element_pressure_peak_pa` result fields | retain Gauss-point `det(F)` for both formulations and eliminated P0 pressure for the Q1/P0 path at peak load | positive finite Jacobians are a necessary validity condition; pressure and mesh sensitivity remain problem-level questions |
| `cardiac_benchmark/mesh_quality.py` | reports reference Hex8 Jacobians, volume, and edge aspect ratio | detects inverted, degenerate, or highly distorted cells; no universal accuracy threshold |
| `cardiac_benchmark/fiber_crosscheck.py` and `structural_directions.py` | provide an independent oracle and the shared runtime adaptation of the pinned complete-frame formula | pytest checks fiber/sheet/normal at physical points and the apex branch; these checks establish formula identity, not ventricular twist or a whole trajectory |
| `solver_configuration_json` and `nonlinear_step_diagnostics_json` | retain nonlinear parameter values and one accepted-step record per increment | allow convergence to be audited separately from the displacement curve |
| `pre_solve_audit_json` | retains geometry, boundary, signed pressure-resultant/moment, Robin checks, and rigid translation/rotation reductions | establishes the recorded setup passed its declared pre-solve gates and exposes mesh-dependent discrete Robin stiffness; it does not establish curve agreement or impose a mesh-dependent validation target |
| `cardiac_benchmark/viscous_evidence.py` | retained historical accepted-state eta-split implementation | preserves the debugging method and its unit tests; it is not callable from the public benchmark driver and is not evidence for the final fixed-parameter comparison |
| module `__main__` checks in `activation.py`, `geometry.py`, `material.py`, `pressure.py`, and `robin.py` | focused scale, geometry, tangent, orientation, and analytic-reaction checks | explicit local failures or reported diagnostics; pytest carries the automated gates |

Run a diagnostic from the repository root, for example:

```bash
python examples/cardiac_benchmark/fiber_crosscheck.py
python examples/cardiac_benchmark/mesh_quality.py
python examples/cardiac_benchmark/diagnose.py caseB_full.npz
```

For an open-apex mesh, a cavity volume is not uniquely defined without an
explicit cap construction. Do not report a closed-cavity volume unless the cap
policy, orientation, and numerical check are stated.

## Distributed implementation examples

The five `mpi_smoke/distributed_cardiac_*.py` scripts cover passive assembly,
transient Robin dynamics, viscous state, follower pressure, and a
timing-oriented larger-mesh path. The first four compare with the serial
solution of the same local problem; the scaling script compares each
multi-rank result with its distributed rank-1 result. These checks can detect
partition, scatter, boundary-row, or state-commit errors; they cannot establish
physiological validity, benchmark-curve agreement, or a general scaling law.
All five use the open-apex smoke geometry and a softened volumetric penalty
(`kappa=1e3`), so they are not benchmark-parameter comparisons.

The opt-in `mpi` pytest selection also includes an application PETSc SNES
callback/reuse check. With a matched PETSc/MPI installation:

```bash
python -m pip install -e ".[dev,mpi]"
python -m pytest -q -m mpi
```

## Reproducibility boundary

The repository makes the driver, formulations, solver policy, diagnostics,
source adaptations, and distributed scripts available. Component checks show
what their named controls evaluate. The existing reduced Case A output remains
a completed historical execution, while archived source-identified Case A and
Case B reports retain their exact source checkpoint, configuration, all accepted-step
records, result and log hashes, external-data identities, complete histories,
and RED output together. The legacy Case B table records observations whose raw
archives are absent and remains labeled accordingly.

The selected closed Case A report and the checked-in Step 2 Case B development
report were produced before the toolkit-matched straight-wall mapping and
physical-coordinate structural-frame reconstruction were introduced. They are
retained under their exact source identities and are not relabeled as output
from later source. A separately identified corrected-geometry Step 0B record
binds the closed tip-refined dual-run full cycle; the earlier 0.32 s report
remains a prefix diagnostic for its own recorded configuration.

The matched 2×24×32 Case B reports hold mesh, time step, geometry, sampler,
source, and solver fixed at checkpoint `62ad760` while changing the volumetric
formulation, so their output difference is a controlled formulation comparison
at that configuration. It does not establish which path is more accurate. The
2×12×16 Q1/P0 reports hold the spatial problem fixed and halve the time step;
their common-grid history difference is a two-step sensitivity measurement,
not a time-convergence claim. The archived 2×36×48 F-bar peak is close to its
legacy-reported vector, but absent legacy artifacts and different
apex/sampling policies prevent an exact-reproduction claim.

The `e07993b` Q1/P0 2×24×32, `dt=0.004 s` record preserves 46 rejected trial
residual evaluations at steps 132–133 followed by valid accepted/final states
under the unchanged solver rule. This establishes that the documented
backtracking recovery operated in that run only. Its `dt=0.002 s` counterpart
was produced at `62ad760`, so the two histories provide cross-checkpoint
sensitivity rather than a controlled time-step comparison. None of these
records establishes bifurcation, qualifies the open-apex geometry, or provides
clinical or real-device validation.

The `e07993b` Q1/P0 2×36×48, `dt=0.002 s` record similarly preserves 168
rejected trial residual evaluations at steps 277 and 279 followed by 500 valid
accepted/final states under the unchanged rule. Its nominally matched F-bar
record comes from `62ad760`; the source and invalid-trial policy therefore
differ. The side-by-side outputs do not isolate formulation, rank accuracy,
establish convergence, reproduce the legacy F-bar observation, or support a
bifurcation claim.
