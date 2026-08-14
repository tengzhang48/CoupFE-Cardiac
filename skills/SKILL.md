# CoupFE cardiac benchmark skill

Use this repository as a cardiac application over the exact CoupFE revision
pinned in `pyproject.toml`. Start with the fast tests, then select the smallest
relevant optional check: `slow` for the compiled reduced serial driver, or
`mpi` for the PETSc callback check and documented 1/2/4-rank scripts.

Before changing or interpreting Case B, read `docs/CASE_B_STATUS.md`,
`docs/CASE_B_DEBUGGING_POSTMORTEM.md`, and `docs/lessons_learned.md`. Preserve
the published `A_EPI=1e8 Pa/m`; do not tune
it merely to fit a curve. Preserve the benchmark `kappa=1e6 Pa` and
`eta=100 Pa s` for every remaining run. Earlier changed-parameter controls are
historical forensic records; the public driver rejects new nonpaper `eta` and
isotropic requests, and its former accepted-state eta-split option, before
setup. Keep ventricular mesh adapters, facet semantics, output sampling, and
volumetric experiments in this application, not Core.

Treat geometry and boundary-condition identity as a pre-solve gate. The
historical `polar_ring` mesh with `apex_offset=0.2` is a truncated-domain
record with an additional free tip surface; do not call it the paper's closed
geometry. A closed Case B candidate must use the noncollapsed
`closed_multiblock_disk` topology and pass the retained exterior-label,
extended-Jacobian, reference-measure, pressure-resultant/moment, and Robin
audits before dynamics. Passing those audits qualifies the setup, not the
trajectory or its agreement.

For every run intended as evidence, report the application and Core revisions
and tree states, formulation (`fbar`, `local-pressure`, or `std-kappa`), mesh
topology and topology-specific counts, element/DOF counts, apex treatment,
mass, viscosity, transmural-coordinate definition, fiber sampling, point
sampler, time step and integrator, nonlinear solver configuration, pre-solve
audit, and all-step acceptance status. Current output points must use the
checked reference Hex8 isoparametric sampler and retain their elements, natural
coordinates, weights, and reconstruction errors. A dirty-tree result remains
development evidence; do not relabel it after a later commit.

For PETSc SNES, require finite iterates, acceptable KSP/SNES reasons, and the
independently reassembled final residual before state commit. Retain every
accepted-step diagnostic rather than inferring convergence from an iteration
count or the existence of an NPZ file.

Treat a nonpositive or non-finite trial `det(F)` through the narrow recovery
contract. Core Newton uses the Q1/P0 operator's `max_step` determinant-domain
bound. The PETSc residual callback catches only `InvalidDeformationError`,
returns IEEE positive infinity so `bt` shortens the trial, and records the
rejection count and last detail. Never catch unrelated residual exceptions or
Jacobian exceptions as domain rejections. Invalid initial, accepted, or final
states remain fail-closed, and no state or completed archive may survive them.
Do not change nonlinear tolerances to obtain this recovery.

Treat the 2026-06-27 F-bar table as legacy-reported history: its raw archives
and logs are absent. Do not relabel it as Q1/P0 output or current-source
validation. Retain completed source-identified results with explicit
qualification boundaries, result and transcript hashes,
an explicit constitutive `material_model_id`,
complete `p0`/`p1` histories, positive deformation Jacobians, all-step solver
records, the verified CC BY 4.0 Zenodo archive and exact curve-file identities,
and full-precision RED output. State that it is a result at one named
configuration; require controlled spatial/time sensitivity before claiming
mesh/time convergence.

Keep the public surface legible as five separate things: application code,
automated tests, runnable example commands, benchmark-comparison reports, and
retained console output. Raw Zenodo pickle files and generated NPZ files remain
external. A retained report may include the cited CC BY 4.0 derived all-team
mean curves and RED values. State when no repository-defined RED pass/fail
threshold exists.

Generate a closed-mesh Laplace transmural field with the application-owned
`tbar_laplace.py` utility, and retain its mesh parameters, numerical checks,
and SHA-256. Keep the native `.meta.json` beside the NPY; require the driver to
validate and record both portable identities before compilation. Record
whether structural directions use nodal CG1 interpolation plus Gram--Schmidt
or direct analytic evaluation at Hex8 Gauss points. These are numerical
representation choices, not physical-parameter tuning.

Before using MPI for a closed Case B result, require the same five-block mesh,
pointwise-`kappa` material, assembled consistent mass, Laplace field, and
GP-direct fiber policy as the serial candidate. Retain a serial-versus-1/2/4-
rank comparison across the snap window before any wall-only, intermediate, or
fine production run. Do not substitute the historical Q1/P0/lumped/CG1 MPI
companion as if rank count were the only changed variable. Use
`--load-horizon 1.0` for the shortened gate, retain it in the archive/report,
and verify its load arrays are the exact prefix of the full production
schedule. Use `compare_mpi_rank_gate.py` for the hard-pinned 2×20×17 serial,
one-, two-, and four-rank comparison; do not replace it with visual log checks.

Use `compare_fenics_case_b.py --retained` for a public direct comparison with
the local FEniCS displacement output. Supply and verify all five role hashes,
including the clean CoupFE result. Require the corrected constitutive model
ID, clean public app/Core identities, closed setup audits, and the fixed paper
parameters. The snap window and `u_z=-5 mm` onset rule are declared metrics,
not tuning controls. Reject both POSIX and Windows machine-local paths and any
nonpublic Core URL.

Do not use the locally retained FEniCS point-stress arrays as a quantitative
oracle. Their supplied postprocessor omits the accepted velocity/acceleration
state, projects to DG1, and contains a dimensionally inconsistent von Mises
expression. A future stress comparison must reconstruct a corrected tensor
from accepted `u/v/a` and retain matched element/quadrature locations and
physical separation.

When F-bar and Q1/P0 runs share mesh, time step, geometry, sampling, and solver,
report their vector and RED differences as a controlled formulation comparison
at that configuration. Do not turn a lower RED at one point or component into
an accuracy ranking, and do not call one matched pair a mesh/time study.
For a time-step pair, compare vector histories on the exact common time grid
and report both maximum and RMS differences; two steps are sensitivity evidence,
not convergence. A current result close to a legacy-reported vector is a
recovery only within the surviving record, not an exact reproduction when the
legacy artifacts and complete configuration are absent.
