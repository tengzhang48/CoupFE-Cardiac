# Cardiac implementation guidance

Read `docs/CASE_B_DEBUGGING_POSTMORTEM.md` before starting a Case B diagnostic
campaign. Establish the benchmark evidence inventory, geometry/boundary gate,
and discrete-method manifest before allowing a full trajectory run. Treat a
null experiment as configuration-specific, especially when it was performed
on the historical open domain.

Keep the historical polar mesh and the benchmark-target closed multiblock mesh
as explicitly named topologies. A positive historical `apex_offset` truncates
the domain; zero collapses cells. Neither is interchangeable with the
noncollapsed closed topology. Before time integration, require exactly one
boundary label per exterior face, positive Gauss and extended Jacobians,
declared reference measures/extents, the expected pressure resultant and
moment, and Robin matrix checks.

The geometry provides analytic material-coordinate data. At Gauss points,
record whether the run interpolates nodal fiber/sheet directions and applies
Gram--Schmidt or evaluates the analytic frame directly. When using the paper's
Laplace transmural coordinate, generate it with `tbar_laplace.py`, verify its
boundary values, bounds, residual, mesh identity, field hash, and native
sidecar hash, then pass the exact field and adjacent sidecar to the driver.

For Case A, preserve the helix convention checked against the pinned upstream
formula adaptation and the structural-tensor state path. For Case B, orient
pressure normals away from each parent-element interior. Check vectorized and
scalar pressure residuals, finite-difference tangents, and an
opposite-orientation broken control.

Keep all three volumetric paths explicit:

- F-bar is the historical single-displacement-field path.
- Standard-kappa uses the standard Q1 kernel with the paper's pointwise
  `kappa=1e6 Pa` volumetric penalty.
- The application-owned Q1/P0 path uses the standard cardiac kernel with its
  material bulk penalty set to zero and composes one algebraically eliminated
  pressure per Hex8. Check affine and isochoric behavior, the full condensed
  tangent, symmetry, and invalid reference/deformation controls.

Sample `p0` and `p1` through the checked reference Hex8 isoparametric map. Do
not substitute global node-tetrahedralization weights or relabel older
Delaunay-sampled histories as current Hex8 output.

Use backward Euler for the coherent first-order path. Treat Newmark as an
experimental option because material viscosity remains a backward strain
difference. Before committing state, require the application residual rule;
for PETSc also retain SNES/KSP reasons, iterations, norms, thresholds, and
residual histories for every step.

For Q1/P0, reject a nonpositive or non-finite trial `det(F)` with the dedicated
`InvalidDeformationError`. Core Newton consumes the operator's `max_step`
determinant-domain bound. In the PETSc residual callback only, translate that
exact exception into an IEEE positive-infinity residual so `bt` can shorten
the trial; retain the rejection count and last detail. Let Jacobian and
unrelated residual exceptions abort. Invalid initial, accepted, and final
states remain fail-closed, and this recovery must not change the solver
tolerances or final residual rule.

For a claimed comparison, record the mesh and apex treatment, formulation,
fiber and output sampling, time step, integrator, solver, application and Core
commits, result hash, all-step nonlinear acceptance status, Gauss-point
deformation Jacobians, complete `p0`/`p1` curves, and verified external-data
identity. Require the closed setup audit for a closed-domain claim. An
open-apex surface has no unique cavity volume without an explicit, checked cap,
and a capped volume does not restore missing myocardial material or boundary
support.
Do not infer bifurcation or a unique twist direction from a sign change between
two discrete runs.

Keep `eta=100 Pa s` for every remaining benchmark run. The public driver must
stop before setup for nonpaper `eta`, isotropic, or accepted-state eta-split
requests. A stopped or failed historical changed-eta trajectory is not a
result and does not authorize tuning. Declare
the trilinear Q1 Hex8/backward-Euler versus quadratic P2
tetrahedral-displacement/generalized-alpha numerical
gap with any comparison.
Require every new archive/report to identify the constitutive semantics with
`material_model_id`. Preserve pre-`6839c13` smooth-switch-stress records as
historical; never use them as numerical oracles for the corrected complete
smooth-switch-energy derivative, and never relax a tolerance to bridge that
law change.
For a Case B spatial study, refine through the wall as well as in-surface,
retain at least three qualified Q1 levels, and include a wall-only control.
Comparable displacement-DOF counts do not make the approximation orders
equivalent.
Port and validate MPI before the wall-only or fine production runs. The MPI
gate must compare the same closed pointwise-`kappa`, consistent-mass,
Laplace/GP-direct configuration against serial at one, two, and four ranks
across the snap window. The historical Q1/P0/lumped/CG1 companion is not that
gate. Generate shortened-gate loads with `--load-horizon 1.0`, retain that
field in every archive/report, and prove the shortened schedule is the exact
prefix of the production schedule before comparing ranks.
Run `compare_mpi_rank_gate.py` on the four completed 2×20×17 archives. Its
mesh, time interval, and tolerances are fixed; do not introduce a smaller
acceptance mesh or a caller-selected tolerance.

For the direct local-FEniCS landmark comparison, use
`compare_fenics_case_b.py` in retained mode. Require expected SHA-256 values
for all five inputs: the clean CoupFE archive, FEniCS parameters, time grid,
and both displacement histories. Require the corrected
`material_model_id`, the public Core URL and clean source identities, and the
closed fixed-parameter setup. Keep the fixed 0.20--0.32 s snap window and
`u_z=-5 mm` onset rule; do not expose them as fitting parameters. Do not copy
machine-local or Windows paths from input metadata into the report.

Do not compare against the existing FEniCS point-stress arrays. Reconstruct a
corrected stress from retained accepted `u/v/a` and compare declared nearby
element-interior or quadrature samples only after geometry and displacement
histories are controlled.

Distinguish nonlinear increment acceptance from mesh/time convergence. A run
whose every requested step met the recorded residual rule can be retained and
reported even when resolution sensitivity is still open. Preserve its JSON
report and console log, leave generated NPZ and raw external pickles unbundled,
and state the unresolved qualification boundary. State when no
repository-defined RED pass/fail threshold exists. Compare formulations as such
only when mesh, time step, geometry, and other model choices are controlled.
For a matched pair, report all `p0`/`p1`
component and RED changes, including mixed directions, and limit the conclusion
to that discretization. For nested time grids, compare full vector histories on
the shared samples and report maximum and RMS differences. Two time steps or
two meshes can expose sensitivity but do not establish convergence. When a
current peak is close to a legacy-reported vector, retain the quantitative
difference and the missing-artifact/configuration caveat rather than claiming
byte-for-byte reproduction.
