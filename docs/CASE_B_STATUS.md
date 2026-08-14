# Case B implementation and result record

This document keeps the implemented Case B path, legacy development history,
and current evidence record separate. That distinction matters because the
2026-06-27 observations survive as a table, while their raw simulation archives
and logs do not.

For the detailed account of why geometry and mass were checked too late, what
the failed and null experiments establish, and which FEniCS/CoupFE differences
remain unresolved, see the separate
[`CASE_B_DEBUGGING_POSTMORTEM.md`](CASE_B_DEBUGGING_POSTMORTEM.md).
For the measured numerical-sensitivity contributors (boundary-normal
representation and bulk spatial resolution), the through-wall layer assessment,
and the ordered resolution plan, see
[`CASE_B_MESH_ERROR_LAYERS.md`](CASE_B_MESH_ERROR_LAYERS.md).
For what each case actually solves — geometry, boundary conditions, material
model, element type and time integration — see
[`CASE_SPECIFICATIONS.md`](CASE_SPECIFICATIONS.md); for the mesh axes and their
naming, [`MESH_REFINEMENT_GUIDE.md`](MESH_REFINEMENT_GUIDE.md).

Private coordination notes and raw simulation archives are outside this
repository. Their release-relevant conclusions and the identities of retained
evidence are carried here, in the documents listed above.

## 2026-08-07 current status: two- and four-layer full cycles complete

The current closed `2x20x17` and `4x20x17`, `tip_refine=6.0` Step 0B runs
completed all 1,000/1,000 increments on eight MPI ranks with `dt=0.001 s`.
Both use the facet Robin operator, consistent mass, source-matched
generalized-alpha, condensed Q1/P0 mean-`log(J)` pressure, and `eta=100 Pa s`.
The retained two-layer archive records clean application/Core revisions
`ae2c2eb`/`454f73c`. The checked-in dual-run
[`comparison report`](../examples/cardiac_benchmark/results/step0b_tip6p0_full_cycle_comparison.report.json)
and [`six-component figure`](figures/step0b_tip_refine_full_cycle.svg) use
matched CoupFE/FEniCS timestamps and the exact ten-team manifest. This is one
controlled comparison pair, not convergence, rank equivalence, or validation.

Relative to the two-layer trajectory, the four-layer run is nearly identical
before snap: vector-history RMSE is 0.009831/0.002622 mm at `p0`/`p1`.
Maximum separation is 1.633083 mm at 0.250 s and 1.405866 mm at 0.247 s;
full-cycle RMSE is 0.407604/0.361832 mm; late-cycle (0.75--0.999 s) RMSE is
0.051556/0.062702 mm; and cycle-end separation is 0.034598/0.051776 mm.
The normalized late/cycle-end ratios are 0.031570/0.044600 and
0.021186/0.036829. The shared 0--0.32 s states agree with the independently
completed prefix to roundoff.

Against FEniCS, full-shared-history RMSE improves from 1.092675/1.168502 mm
(two layers) to 1.016944/1.072951 mm (four layers), and maximum gap improves
from 5.551997/5.023942 to 4.188388/3.753416 mm. The downward/upward
`u_z=-5 mm` crossings are 0.242327/0.561007 s and 0.239063/0.558177 s for
FEniCS, 0.246148/0.567678 s and 0.242905/0.564664 s for two layers, and
0.244986/0.568785 s and 0.241707/0.565612 s for four layers (`p0`, then
`p1`). Four layers move the downward crossing toward FEniCS, but the upward
crossing remains late and relaxation RMSE worsens by 16.2%/14.9%. The
ten-team containment fractions are likewise mixed. This phase-resolved
evidence supports numerical transient/timing sensitivity, not a physically
distinct branch or a convergence claim.

The first full-cycle four-layer candidate was produced from a dirty
application tree and is explicitly non-retained. A clean isolated
replay completed under Core `454f73c` and runtime-source SHA-256
`f8d9469a101709d11460a0803b6c031001192ccce61c472d032355b01070da05`;
its states match the candidate to roundoff. Retained provenance is
**application `2458e7c`, NPZ SHA-256
`1e333b29b05f01dedce9272b32b82ce6ccfda56036c1c8f57eb395b9b4494800`,
stdout SHA-256
`0810a9d0c944345464f48917c707a7d3df7187b543c74efd4b0384542b95aed7`,
elapsed `1778.4 s (29.6 min)`**.

### Reproduction verdict

The answer is **partial reproduction, not exact reproduction**. The current
run reproduces the pressure-only physical geometry/load/boundary/time contract
and the principal snap-through/recovery response. It uses the documented
Q1/P0 mean-`log(J)` volume variant rather than the paper's pointwise scalar
volumetric response, so the discrete equations are not identical. This is a
claim boundary, not evidence that the volume law causes the remaining gap.
Separately, the paper-law change was a null only through 0.70 s on one named
fine, pre-correction Case A configuration; it is not a Step 0B null result.
Against the local hash-pinned FEniCS trajectory, the four-layer shared-history
vector RMSE is 1.016944/1.072951 mm at `p0`/`p1`.
Its downward `u_z=-5 mm` crossings are 2.659/2.644 ms late and its upward
crossings are 7.779/7.435 ms late. The largest vector gaps, 4.188/3.753 mm,
occur near snap.

An independent verification against the exact hash-gated ten-team mean on the
canonical `0:0.01:1 s` grid gives full-history relative L2 of 8.9797%/9.1970%,
vector RMSE of 1.1924/1.2484 mm, and benchmark Eq. 21 RED of
0.175895/0.311732 at `p0`/`p1`. Thus the aggregate histories meet the project's
approximate 10% working level at both landmarks; this is not a paper-defined
pass threshold and does not erase phase/component-specific differences.

That agreement is quantitatively useful but incomplete. The four-layer result
improves full-shared-history FEniCS RMSE and the maximum snap gap relative to
two layers while worsening relaxation RMSE. Agreement for small components and
the ten-team envelope also remains phase dependent. Therefore this result
supports a credible qualitative and partial quantitative reproduction, not
exact curve agreement, spatial convergence, exact method equivalence, or
validation. See the answer-first
[`benchmark reproduction status`](BENCHMARK_REPRODUCTION_STATUS.md) for the
cross-case claim boundary.

## Superseded 2026-08-05 status: 0.32 s gate complete; 1 s paused

> Historical status: the pause decision below preceded the completed
> `2x20x17`, `tip_refine=6.0` full-cycle run reported above. It remains here to
> preserve the decision chronology and the claim boundary of the separate
> clean prefix artifact.

The current clean Step 0B gate is complete through the snap window on eight MPI
ranks. Both the coarse `2x20x17` mesh and the through-wall-only `4x20x17` mesh
accepted all 320 requested increments through `t=0.32 s`. They used the current
straight-wall closed geometry, physical-coordinate toolkit frame, condensed
Q1/P0 mean-`log(J)` pressure with `K=1.0 MPa`, consistent Q1 Hex8 mass,
`eta=100 Pa s`, and source-matched generalized alpha
`(alpha_m,alpha_f,gamma,beta)=(0.2,0.4,0.7,0.36)`. The time step was `0.001 s`
and the load was evaluated from the unchanged 1 s schedule. PETSc used the
reviewed FGMRES/GAMG rigid-body-rebuild profile.

The endpoint comparison is:

| Trajectory | `p0` displacement at 0.32 s (mm) | `p1` displacement at 0.32 s (mm) | 0--0.32 s vector RMSE versus FEniCS, `p0`/`p1` (mm) |
|---|---|---|---:|
| CoupFE `2x20x17` | (-22.153406, +1.258442, -15.162644) | (-23.301088, +0.557504, -16.143093) | 2.58885 / 2.53042 |
| CoupFE wall-only `4x20x17` | (-22.276570, +1.234342, -15.821064) | (-22.877477, +0.419922, -16.575841) | 2.34986 / 2.21569 |
| local FEniCS P2 reference | (-16.064352, +1.088629, -16.039858) | (-17.042406, +0.394483, -16.626946) | -- |

The physical-frame correction is a controlled null result for these two
trajectories: relative to their otherwise matched stored-parametric-frame
prefixes, the maximum pointwise vector-history change over both landmarks and
both meshes is `0.03566 mm`, and the largest vector RMSE is `0.01077 mm`.
Therefore the frame bug needed correction, but it does not explain the present
post-snap displacement gap.

The full 1 s extension is intentionally paused. The 0.32 s prefix already
contains the snap-through interval and distinguishes the dominant surface-mesh
effect described below; continuing through unloading and recovery would not
resolve which finite boundary-surface discretization should be compared. The
stop rule is: first declare the benchmark surface target, rerun that exact
choice through this short gate, and inspect the complete displacement and
boundary-response histories. Extend that same configuration to 1 s only if the
short gate is informative enough to justify a full-cycle record, then verify
that its retained 0.32 s prefix agrees with the gated prefix. Until then these
are diagnostic prefixes, not full retained results, RED values, a recovery
comparison, rank-equivalence evidence, or a spatial-convergence proof.

The CoupFE archives identify clean application revision
`056c02df7c2a56bbc36e41973b7c8a8d8c917e2a`, clean Core revision
`e2f42ed5772850a0a23a2ce434f430c287eae5c8`, driver
`examples/cardiac_benchmark/run_mpi.py`, and embedded runtime-source manifest
`bec13f9ab1bbd9e50116e05b7342501e6ff8992c5710f62b83bcd02a089b6cf1`.
They were produced from an isolated detached checkout at the recorded clean
revision and retained below the external campaign root as:

- `clean/coarse/step0b_local_pressure_ga_nt2_core20_rad17_t0p32_frame_clean.npz`,
  SHA-256 `efc7e42a60218ab275df2250ae4383f081554e25c46beee1c25cff3248c47785`;
- `clean/wall/step0b_local_pressure_ga_nt4_core20_rad17_t0p32_frame_clean.npz`,
  SHA-256 `dc79faa158f04007c2430592024c69030a0894bde02967caac6dc338194b706f`;
- `diagnostics/clean_physical_frame_prefix_metrics.json`, SHA-256
  `8350e74322daa2b663f3d69511b83c9d9e770d137401f3eae730fa1893cacfaa`;
- `diagnostics/clean_physical_frame_prefix_comparison.svg`, SHA-256
  `fa13b940f0d3ca2b08c21a193f0391000fed0bb1d6d91f6a05e545b1bda2a2b5`.

The reviewed compact package record is
[`step0b_case_b_clean_frame_0p32.report.json`](../examples/cardiac_benchmark/results/step0b_case_b_clean_frame_0p32.report.json),
SHA-256 `faef788ede2b42a175ed422292f64d916a0a30088faabff1e9252408388c0a3f`.
It binds the external archives, logs, comparison metrics, and stop decision
without copying multi-megabyte solver histories into the repository.

The FEniCS inputs are hash-gated as `componentwise_displacement_up0.npy`
`4344a4f599a6eabb16159682339a735bff572eaa18eedd1fe2a97ebd3ee7f4a0`,
`componentwise_displacement_up1.npy`
`88a679de2189bc137de5d64186c698f1702e9df20333849377cdfd01aac8bf1e`,
`time_stamps.npy`
`ddba330b1c8f8c1bb61282e187047f3aa99d0df37b2c4ed2139ea1b0e0ff0f0c`,
and `parameters.json`
`c1cd4c8d2521fd6c28774975843740a8af12568edd1240f5daa133d469e6fb76`.
The producing FEniCS revision was not retained contemporaneously; the nearby
clean clone is contextual and is not claimed as the producer of those arrays.

## Implemented model pieces

The Case B driver includes the published passive material parameters,
`A_EPI=1e8 Pa/m` epicardial Robin stiffness, pressure history, and a
deformation-dependent endocardial follower load. Automated component checks
cover pressure-history scale, facet orientation reversal, scalar/batched
residual equality, and finite-difference follower-load tangent agreement.

The driver offers four explicit volumetric discretizations:

- `fbar`: the historical single-displacement-field compiled Hex8 path.
- `std-kappa`: the standard Q1 Hex8 material kernel with the pointwise
  `kappa` penalty used in the paper's volumetric energy.
- `local-pressure`: a standard Q1 material kernel with its bulk penalty
  disabled, composed with an application-owned operator that computes and
  eliminates one P0 pressure per element from the volume-average of
  `log(det(F))`.
- `local-pressure-paper`: the same condensed Q1/P0 pattern, but with the
  paper-law response evaluated at the geometric element-mean dilatation,
  `Jbar=exp(<log(det(F))>)` and
  `p=kappa*(Jbar**2-1)/2`, including the matching condensed tangent.

The Q1/P0 operator has affine, isochoric, tangent, symmetry, and invalid-cell
controls. Those tests establish the named element-level properties. They do not
by themselves show that a ventricular Case B curve is spatially or temporally
converged.

## Output sampling and nonlinear acceptance record

Current runs locate `p0` and `p1` by inverting the reference trilinear Hex8 map
and applying the eight Hex8 shape weights. The archive records the selected
element, natural coordinates, weights, and reconstruction error under
`point_sampling=hex8_reference_isoparametric`. This is the displacement field's
own interpolation, not the global Delaunay-tetra convenience used by the
existing historical Case A record.

The application can use guarded Core Newton or persistent PETSc SNES. The PETSc
parameter values match the recovered 2026-06-27 adapter: `newtonls`, `bt`,
`preonly`, `lu`, `rtol=1e-9`, `atol=1e-10`, `stol=1e-12`, and `max_it=60`.
The new implementation additionally requires a valid KSP status, finite
iterate, positive SNES reason, and an independently reassembled final residual
no larger than `max(atol, rtol*|R_initial|)` before one state commit.

Successful PETSc archives retain the solver configuration and one diagnostic
record per accepted increment, including SNES/KSP reasons, nonlinear and linear
iterations, residual norms and history, acceptance threshold, and measured
assembly/solve time. Timing fields diagnose a run; they are not a performance
claim without a controlled hardware and load policy. Core-Newton archives
retain the smaller diagnostic set available from that path.

For Q1/P0 invalid trial deformations, Core Newton uses the operator's
`max_step` determinant-domain bound. The PETSc residual callback catches only
`InvalidDeformationError`, supplies an IEEE positive-infinity trial residual so
`bt` can shorten the step, and records the rejection count and last detail.
Jacobian and unrelated residual exceptions still abort. An invalid initial,
accepted, or final state cannot be committed, and the recovery changes no
solver tolerances or final residual rule. Result records remain source-specific:
reports produced before these diagnostic fields were introduced legitimately
omit them.

The distributed pointwise-`std-kappa` batch does not yet have the equivalent
explicit `det(F)>0` trial-state guard. Its kinematic probe is used for retained
output, while the generated material evaluates the pointwise paper penalty.
Consequently, `function_domain_rejections=0` on a `std-kappa` archive does not
prove that every rejected or line-searched trial remained orientation
preserving. The retained fine Case B pressure-peak state is positive at every
stored Gauss point and independently reconstructs exactly, so this is a
formulation-coverage gap in the admissibility guard, not evidence that it
caused the accepted trajectory.

## Fiber convention and signed response

In this repository's local coordinate basis, the default `flip_helix=True`
negates the prescribed helix-angle field so that its analytic `f0` agrees with
the pinned NumPy adaptation of the upstream fiber formula. Toggling the option
changes prescribed helix handedness; it is not a harmless global
`f0 -> -f0` sign change.

Closed-mesh runs now reconstruct the complete fiber/sheet/normal frame from
the physical point and Laplace transmural coordinate using the pinned toolkit
formula. Serial and MPI Gauss-point paths share that implementation. The
historical polar topology retains its parametric convention. See
`PHYSICAL_FRAME_RECONSTRUCTION.md` for the source, tests, and coarse/fine field
audit. This local-frame check still does not validate ventricular twist or the
sign of every output displacement; mesh resolution, apex treatment,
formulation, and nonlinear path can also affect the signed response.

## Legacy-reported 2026-06-27 F-bar observations

A development campaign dated 2026-06-27 used the F-bar formulation, a
temporary PETSc SNES adapter, a collapsed structured apex, and in-plane mesh
refinement. A surviving development record says that the listed runs completed
and that their archives were read to create the following table. The associated
implementation and narrative were committed in the legacy repository as
`26fe8e21b84fa2ad8d94b2d56a4b75f49e37056c`, followed by merge commit
`2711547552173c9371c65c5aae290456e616ad03`.

Legacy-reported displacement of `p0` at peak pressure, in millimetres:

| mesh (`nt×nmu×ntheta`), `dt` | `u_x` | `u_y` | `u_z` |
|---|---:|---:|---:|
| 2×12×16, 1e-3 | +17.59 | -0.18 | -0.51 |
| 2×24×32, 2e-3 | -6.91 | +2.10 | +12.47 |
| 2×28×36, 2e-3 | -14.22 | +2.02 | +11.53 |
| 2×36×48, 2e-3 | -35.86 | +1.49 | -14.22 |
| 2×44×60, 2e-3 | -38.74 | +1.34 | -14.35 |
| 3×36×48, 1e-3 | -31.26 | +1.62 | +14.54 |
| comparison value recorded in the legacy note | -24.68 | +1.57 | -16.55 |

The raw CoupFE NPZ archives, complete commands, environments, normalized logs,
and exact external-curve inputs behind this table were not pushed and have not
been recovered. The comparison row therefore remains a value recorded in the
legacy note, not a newly recomputed all-team mean. The table documents project
history and motivates rerunning Case B; it is not current-source quantitative
evidence.

The observations also do not show a monotone signed limit: `u_z` changes sign,
and the magnitude of `u_x` varies materially with discretization. Earlier
development notes described this as a chiral near-bifurcation. Sign changes
between discrete runs do not by themselves establish a bifurcation,
bistability, a converged signed response, or a unique global twist direction.

## 2026-08-02 topology and 2026-08-05 straight-wall corrections

The paper reference mesh is closed at the ventricular apex. The seven
checked-in records in the archived source-identified results table instead
used the historical polar-ring mesh with `apex_offset=0.2`. That choice stops
1.9335 mm before the apex and creates 96 unclassified exterior quads with
total area 2.672330 cm². They act as a traction-free annular cut; they are not
part of the paper's boundary. The older clean fine external record documented
below used the first topologically closed replacement, but predates the current
straight-wall mapping and physical-coordinate frame. The current clean prefixes
in the answer-first section use both corrections.

The reference/truncated-polar-mesh measures are, respectively: wall volume
177.694656/176.483537 cm³, endocardial area 155.263738/154.489983 cm²,
epicardial area 228.505027/226.864499 cm², and base area
18.577290/18.541517 cm². Although the corresponding differences are all below
0.72%, the fine truncated polar mesh's undeformed pressure resultant is
4.5937% below the analytic closed-cavity projection, about 1.324 N at its
retained peak pressure.
The open-apex runs remain valid records of exactly what that truncated model
computed, but they are not executions of the paper's closed Case B domain.

The current replacement mesh is an application-owned, noncollapsed five-block
square-to-disk Hex8 topology with the corrected straight-wall mapping. At
`n_t=2`, `n_core=20`, `n_radial=17`, and `core_half_width=0.36`, it has 5,403
nodes, 3,520 elements, and 16,209 displacement DOFs. All 3,680 exterior faces
are assigned exactly once to endocardium, epicardium, or base. Extended
sampling has no nonpositive Jacobians, minimum `det(J)=1.236455937e-9 m³`,
maximum condition `7.876486664`, and minimum scaled Jacobian `0.258134961`.
Its wall volume and endocardial, epicardial, and base areas differ from the
retained reference by -0.122926%, +0.069468%, +0.009025%, and -0.056501%; its
unit-pressure resultant error is 0.102777%, with negligible transverse
resultant and moment. The current clean coarse and wall-only prefixes passed
these stored pre-solve audits. The historical full closed runs below retain
their earlier curved mapping and frame identities rather than being relabeled
after the correction.

Geometry and boundary-condition identity is now a pre-solve gate. A candidate
must pass units/extents, exterior-face ownership and exactly-one labeling,
extended reference-Jacobian quality, measures by boundary marker, pressure
orientation/resultant/moment, and Robin checks before time integration. This
ordering prevents solver or material experiments from compensating for a
different domain.

The intended reproduction keeps the paper physical parameters, including
`kappa=1.0e6 Pa`, `eta=100 Pa s`, the published pressure history, density, and
Robin coefficients. A controlled `eta=0` diagnostic failed closed at 0.026 s;
an `eta=50 Pa s` run was then stopped by user direction after its last printed
accepted state at 0.256 s and wrote no completed result. No viscosity tuning
was accepted, and neither diagnostic is included as benchmark output. Those
are historical forensic records from the development investigation. The
current public driver rejects nonpaper `eta`, isotropic-control, and
accepted-state eta-split requests before mesh or solver setup; all remaining
studies vary discrete choices only.

Numerical-discretization differences remain explicit. The current clean short
gate uses Q1 Hex8/Q1-P0 with source-matched generalized alpha
(`alpha_m=0.2`, `alpha_f=0.4`), while the local FEniCS displacement reference
uses P2 tetrahedra with those generalized-alpha parameters. The historical
full closed runs below used backward Euler and different volumetric paths.
Fiber-field representation also remains source-recorded. Matching the domain,
time integrator, and physical parameters does not erase the spatial and
boundary-discretization gaps.

The locally retained FEniCS point-stress arrays are quarantined as quantitative
oracles. Their postprocessor loads displacement only into a newly constructed
problem, leaves old `u/v/a` fields at zero, reconstructs a velocity of about
`1944.44*u` per second instead of loading accepted velocity, projects to DG1,
and uses a von Mises expression with unsquared shear terms. A future stress
comparison must reload accepted `u/v/a`, use a corrected Cauchy-stress
definition, and compare declared nearby element-interior or Gauss points. This
quarantine concerns the derived stress archive, not the reference displacement
solution.

### Historical pre-straight-wall/pre-frame coarse backward-Euler run

On 2026-08-02, the then-current closed 2×20×17 mesh completed all 1,000 Case B
increments with `dt=0.001 s`, the paper values `kappa=1.0e6 Pa` and
`eta=100 Pa s`, consistent mass, the Laplace transmural field, GP-direct fiber
evaluation, backward Euler, and PETSc SNES. All 1,000 diagnostics report SNES
reason +2 and KSP reason +4, with zero function-domain rejections. The stored
peak-state Gauss-point `det(F)` range is 0.660073--1.152459. No tolerance,
viscosity, or bulk modulus was changed.

Against the local FEniCS quadratic-P2 tetrahedral displacement trajectory on
its native 1 ms grid,
the full-history vector RMSE is 1.5803 mm at `p0` and 1.7903 mm at `p1`;
relative L2 error is 0.1176 and 0.1299. The first downward `u_z=-5 mm`
crossing is 5.3/5.6 ms earlier than FEniCS. Pressure peaks are
16,070.954/16,073.142 Pa at the same 0.482 s sample. Against the exact
ten-team mean on the canonical 10 ms grid, vector RMSE is 1.4835/1.7621 mm.
The dominant x/z histories capture the qualitative snap shape and recovery in
this development comparison; the transverse y histories and parts of unloading
remain visibly different. Only 31.35%/40.59% of all p0/p1 component samples
lie inside the pointwise ten-team min/max envelope.

This execution is deliberately classified as historical development evidence.
It predates the straight-wall and physical-frame corrections and was run
from application revision `6839c13` with a dirty tree while the closed-mesh and
reporting changes were being reviewed. The result archive SHA-256 is
`96e28ea247503a94ee95a2513be32628934f3b201f88ebc51907f99a4eaa31cf`;
the guarded comparison JSON and displacement figure SHA-256 values are
`441df228c88883fec2e599980b0fa2cd1c51aaf813a9f207a952955e52d01486`
and `16a3e1fd44a744ae899350de8d3057534a2c08f13ace5fffa8a524a9ff80caf3`.
They are not placed in the public retained-result table. A public record
requires reviewed committed source and an exact clean-tree rerun; the metadata
will not be relabeled after the fact.

### Historical pre-straight-wall/pre-frame clean fine MPI-4 run

On 2026-08-03, the then-current closed
`(n_t,n_core,n_radial)=(4,36,32)` mesh completed all 1,000 Case B increments on
four MPI ranks with `dt=0.001 s`. It has 29,885
nodes, 23,616 Q1 Hex8 elements, 89,655 displacement DOFs, and four elements
through the wall. The run used the paper pointwise `std-kappa` energy with
`kappa=1.0e6 Pa`, `eta=100 Pa s`, consistent mass, the cleanly regenerated
Laplace transmural field, the then-current GP-direct structural directions,
backward Euler,
and PETSc SNES with `SuperLU_DIST`. It did **not** use condensed local pressure
or generalized alpha.

This run predates the straight-wall and physical-coordinate-frame corrections;
it is retained under its exact historical identity. The application revision
was clean
`1ec27165f9e471404db8c5efb46da04da75c0250`; the public CoupFE Core revision
was clean `e2f42ed5772850a0a23a2ce434f430c287eae5c8`. All geometry, pressure,
and Robin pre-solve audits passed. Every increment has positive SNES and KSP
reasons, all independently checked final residuals satisfy the retained rule,
and no function-domain rejection was recorded. The completed archive SHA-256
is `63d41a1f69dceaa8c1fe7f3c7d46a6de4e40c270a35977214de741e06fc580a3`.
The full archive agrees exactly with the independently completed 0.32 s prefix
in their stored time, pressure, and landmark-displacement arrays, which rules
out restart or prefix-file drift for this execution.

A hash-gated full-history displacement comparison uses the 999 shared native
samples from 0.001 through 0.999 s without interpolation. Against the retained
local FEniCS P2 trajectory, its vector RMSE is 7.5955/7.2363 mm and its
relative L2 error is 0.5654/0.5250 at `p0`/`p1`. Before 0.20 s, vector RMSE is
only 0.1689/0.0551 mm. The difference begins at snap: the first downward
`u_z=-5 mm` crossing is 14.293/14.506 ms earlier than FEniCS, and the later
error is dominated by excessive negative `x` displacement rather than the
visually conspicuous `y` component. At 0.320 s, for example, CoupFE/FEniCS
`u_x` is -27.728/-16.064 mm at `p0` and -27.790/-17.042 mm at `p1`.

The sanitized external full-comparison JSON and its CoupFE/FEniCS figure have
SHA-256 values
`602a1e904973d9700b6cd99f1b922fa81651815c7722c157c9a571bd9d7640f5` and
`f365ec930110741bfba553a0d95fabfaafdc65b80d6984918b207697d0f40343`.
They identify the exact CoupFE archive and all four FEniCS array/parameter
inputs by hash. They are not bundled in this repository, so the row is a clean
external result record rather than a bundled regression oracle. The completed
MPI-4 execution also does not replace the still-missing closed
serial/1/2/4-rank equivalence gate.

### Peak-volume and deformation-Jacobian audit

The fine archive retains the eight material-point deformation Jacobians at
peak pressure, `t=0.482 s`. Its minimum is `0.3272872317`, but that scalar is
not representative of the wall. Only 48 of 188,928 Gauss samples have
`det(F)<0.4`; they occupy about 0.0100% of the reference-volume quadrature
weight and lie in 24 elements of the outermost transmural layer near the
epicardial apex. The reference-volume-weighted mean is `0.9913340` and the
median is `0.9915595`.

The worst Gauss point lies in a well-shaped reference element: its reference
Jacobian condition is 1.468, versus a mesh-wide Gauss maximum of 7.404. Its
principal stretches are 1.0661, 0.7676, and 0.4000. The same element has
inner-side `det(F)` values 0.882--0.980 and epicardial-side values
0.327--0.450, so the minimum is a steep localized through-wall compression,
not an inverted element or a poor undeformed Jacobian. Four symmetry-related
values lie within `1.4e-12` of the one exact minimum, and an independent Q1
kinematic reconstruction matches the archive exactly. This rules out a single
corrupt MPI element or an output-calculation error.

Most importantly, a 2026-08-04 local post-hoc reconstruction from the retained
FEniCS P2 state gives minimum `det(F)=0.332702` at its actual degree-6
tetrahedral quadrature at the same 0.482 s state. At 0.320 s the fine
CoupFE/FEniCS minima are 0.570442/0.527761, and at peak pressure their
reconstructed wall-volume ratios are 0.991334/0.992210. The FEniCS input
`result.h5` used for this check has SHA-256
`68fdbb222bf948aebcc6cdc287450601efaed9eec0aadaf4abbdbc9eaccb74d3`;
the P2 local-DOF map was checked against the saved visualization field and the
24 quadrature weights against the retained compiled FFC form. The
reconstruction script and derived JSON are not yet packaged, so this is a
hash-identified local diagnostic rather than a solver-retained field or
checked-in regression gate.

The CoupFE fine run and the audited nearby FEniCS source implement the paper
energy `kappa*(J**2-1-2*log(J))/4`, a finite penalty rather than an exact
incompressibility constraint. The exact producing FEniCS source revision was
not retained with the output. At the CoupFE minimum the energy produces about
-1.364 MPa hydrostatic Cauchy stress while the strongly compressed fiber and
sheet invariants leave their tension-only reinforcement nearly inactive.
Localized post-snap compression can therefore coexist with nearly conserved
global wall volume. The coarse closed run's `det(F)` minimum of 0.660 could be
consistent with smearing this apex layer; the fine value is not evidence that
refinement made the physical solution worse or that Q1 alone caused the large
landmark-trajectory gap. CoupFE does resolve a broader `det(F)<0.5` tail than
the reconstructed FEniCS field, so distribution and trajectory differences remain
worth measuring. The localized pocket is observed at the post-snap pressure-
peak state; its benchmark-like minimum is not evidence that it caused the
earlier onset discrepancy.

## Controlled closed-mesh prefixes and archived source-identified runs

### Prior clean four-way mesh split and Robin mechanism (2026-08-05)

Before the physical-frame correction, four clean eight-rank prefixes from
application revision `a2006b78104109c625ea3c502753b5cff15452d4` and Core
revision `e2f42ed5772850a0a23a2ce434f430c287eae5c8` completed all 320 increments.
They used the corrected straight-wall geometry and the same local-pressure,
consistent-mass, generalized-alpha, `dt=0.001 s`, `eta=100 Pa s`, and 1 s load-
horizon contract as the current gate, but retained the older stored-parametric
Gauss-point frame. This two-factor split isolates through-wall refinement from
surface refinement:

| Mesh | `p0` displacement at 0.32 s (mm) | `p1` displacement at 0.32 s (mm) |
|---|---|---|
| coarse `2x20x17` | (-22.187842, +1.259134, -15.162859) | (-23.331800, +0.557480, -16.142457) |
| wall-only `4x20x17` | (-22.292914, +1.232925, -15.823229) | (-22.888464, +0.419490, -16.575826) |
| surface-only `2x36x32` | (-32.162179, +1.868219, -15.317890) | (-32.102131, +0.560974, -16.322544) |
| fine `4x36x32` | (-31.775744, +1.735410, -16.098150) | (-31.077364, +0.364589, -16.949891) |

At the coarse surface, adding two wall layers changes the endpoint vector by
only `0.66919/0.63514 mm` at `p0/p1`. Holding two wall layers and refining the
surface changes it by `9.99412/8.77218 mm`; the endpoint effect is
`14.93/13.81` times larger. Thus coarse and wall-only form one response pair,
while surface-only and fine form another. This is a controlled surface-mesh
split, not evidence that the prefix is a full retained result or that the four
meshes establish a converged sequence.

The prescribed normal-only epicardial Robin boundary provides a strong
discrete mechanism for that split. A rigid rotation about the long `x` axis is
tangent to the smooth axisymmetric epicardium, so its continuum epicardial
spring energy is zero. Piecewise-planar facet normals are not exactly
orthogonal to that rotation field. The reference-configuration audit gives
total discrete long-axis stiffness `1.444555 N m/rad` on both coarse-surface
meshes and `0.555828 N m/rad` on both fine-surface meshes. Their epicardial
parts are `1.285103` and `0.396153 N m/rad`; the full-vector base part stays
near `0.1595 N m/rad`. For context, the retained FEniCS mesh gives
`1.857048 N m/rad` total (`1.698036` epicardial plus `0.159012` base). The mean
epicardial edge lengths are `3.996 mm`, `2.173 mm`, and `3.880 mm` for the
coarse CoupFE, fine CoupFE, and FEniCS surfaces, respectively; the paper's
nominal target is approximately `5 mm`. The participant trajectories are
therefore finite-mesh comparisons, not a continuum oracle. The benchmark Robin
law and coefficient must not be changed to compensate for its facet-normal
discretization.

The actual nonlinear states retain the same pairing:

| Mesh | Epicardial/total spring energy at 0.32 s (J) | Epicardial/total long-axis internal moment at 0.32 s (N m) |
|---|---:|---:|
| coarse `2x20x17` | 0.308407 / 0.376541 | -0.771763 / -0.853824 |
| wall-only `4x20x17` | 0.305358 / 0.375250 | -0.774477 / -0.860245 |
| surface-only `2x36x32` | 0.116757 / 0.235046 | -0.248906 / -0.332284 |
| fine `4x36x32` | 0.117304 / 0.235185 | -0.250147 / -0.338022 |

These elastic spring energies and internal moments describe materially
non-rigid states. The archived state does not retain enough velocity data to
reconstruct the dashpot contribution exactly. The pairing is strong mechanism
evidence, but neither this response audit nor the rigid-rotation quadratic form
proves that Robin faceting is the sole cause of the CoupFE/FEniCS trajectory
gap.

The four archives and derived records remain in the external corrected-
geometry campaign directory. Their NPZ SHA-256 values are
`16f2f67032ce35f09a22d2d3ab6accf693861c545f17ad958fbae13c99c933aa`
(coarse),
`1ac0f0d84c12ec7d90b347c8d8d4c532c42728b3abf3bccd92082be202d2bf4d`
(wall-only),
`0a6a034200298928884547aa9b89dadd36b49753de515f22f04f83ec6841f16b`
(surface-only), and
`94279e58d99f497638e48c46bd00510600475faa0260a79b241c97c44e8897b5`
(fine); the same identities and source manifests are in
`diagnostics/corrected_mesh_split_prefix_metrics.json`, SHA-256
`53f77f98fde598c31c151c33f1a16cf9aae18601767a8f5dcda5bb5c5d2241fc`.
The Robin audit is `diagnostics/robin_rigid_rotation_mesh_audit.json`, SHA-256
`406c4226cd1f6cf6879acad2f46291c7d78786706f37f63707f0109716cfdd79`.
Its nearby FEniCS source clone confirms the documented Robin form but is not
claimed as the producer of the retained FEniCS displacement arrays.

All retained `apex_offset=0.2` polar-ring reports describe a fourth free tip
boundary and are therefore historical non-benchmark-geometry evidence. They
are hash-preserved in the dedicated
[`truncated-polar archive`](../examples/cardiac_benchmark/results/archive/truncated_polar/)
rather than in the active closed-geometry comparison set.

Source-identified Case B records are retained under the archive's
[`case_b/`](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/)
directory.
Five earlier records identify clean application checkpoint
`62ad760d2a1731bb9668897863ac026d3768194e`; the newer Q1/P0 2×24×32,
`dt=0.004 s` and 2×36×48, `dt=0.002 s` records identify clean checkpoint
`e07993bcf1166bd20eb87370c0b458552753e7ee`. Both groups identify public Core
checkpoint `454f73ce2de284262b214a2b37bd676c6aca3c0a`, use the nondegenerate
open-apex geometry (`apex_offset=0.2`), the documented fiber convention,
reference-Hex8 output sampling, backward Euler, and persistent PETSc SNES. The
earlier records remain under their original source identity. In light of the
domain correction above, all seven rows are historical truncated-domain
comparisons rather than closed-domain benchmark reproductions; none is deleted
or relabeled as if it came from the new mesh.

| App | Formulation | Mesh and time | Accepted steps | Peak `p0` (mm) | Peak `p1` (mm) | Peak Gauss-point `det(F)` | RED (`p0`, `p1`) |
|---|---|---|---:|---|---|---|---|
| `62ad760` | Q1/P0 local pressure | 2×12×16, 384 elements, 1,872 DOFs, `dt=0.004 s` | 250/250 | (+17.917, -0.405, -0.601) at 0.480 s | (+13.380, +1.099, -1.026) | 0.894270–1.058401 | 0.8831859, 0.8657258 |
| `62ad760` | Q1/P0 local pressure | 2×12×16, 384 elements, 1,872 DOFs, `dt=0.002 s` | 500/500 | (+17.923, -0.406, -0.602) at 0.482 s | (+13.384, +1.100, -1.028) | 0.894372–1.058326 | 0.8853305, 0.8699802 |
| `62ad760` | Q1/P0 local pressure | 2×24×32, 1,536 elements, 7,200 DOFs, `dt=0.002 s` | 500/500 | (-7.886, +2.269, +12.725) at 0.482 s | (-10.525, +1.039, +13.704) | 0.739792–1.197484 | 1.0742661, 1.1626890 |
| `e07993b` | Q1/P0 local pressure | 2×24×32, 1,536 elements, 7,200 DOFs, `dt=0.004 s` | 250/250 | (-7.700, +2.261, +12.715) at 0.480 s | (-10.358, +1.045, +13.688) | 0.741052–1.196233 | 1.0741823, 1.1607469 |
| `e07993b` | Q1/P0 local pressure | 2×36×48, 3,456 elements, 15,984 DOFs, `dt=0.002 s` | 500/500 | (-36.229, +2.037, -14.009) at 0.482 s | (-37.018, +0.985, -15.030) | 0.523357–1.475752 | 0.5511448, 0.6538962 |
| `62ad760` | F-bar | 2×24×32, 1,536 elements, 7,200 DOFs, `dt=0.002 s` | 500/500 | (-12.230, +2.689, +13.088) at 0.482 s | (-14.457, +0.767, +14.332) | 0.699566–1.185359 | 1.0709970, 1.2200494 |
| `62ad760` | F-bar | 2×36×48, 3,456 elements, 15,984 DOFs, `dt=0.002 s` | 500/500 | (-36.956, +1.415, -14.022) at 0.482 s | (-37.205, +1.089, -14.517) | 0.589525–1.316852 | 0.5953528, 0.6971874 |

The Q1/P0 material kernel used `kappa=0` and the application operator used
bulk modulus `K=1e6 Pa`. Peak-load condensed element pressure on 2×12×16 ranged
from -15,705.99 to -2,431.77 Pa at `dt=0.004 s` and from -15,707.11 to
-2,428.72 Pa at `dt=0.002 s`; on 2×24×32 it ranged from -44,954.29 to
-3,965.68 Pa at `dt=0.002 s` and from -44,441.69 to -4,043.97 Pa at
`dt=0.004 s`; on Q1/P0 2×36×48 it ranged from -160,254.02 to +11,309.58 Pa,
with mean -21,286.60 Pa. The F-bar material kernel used `kappa=1e6 Pa`. All
retained peak-load Gauss-point deformation Jacobians are finite and positive.

The machine-readable records are:

- [`case_b_local_pressure_2x12x16_dt0p004.report.json`](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_local_pressure_2x12x16_dt0p004.report.json)
  and its [console log](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_local_pressure_2x12x16_dt0p004.raw.stdout.txt);
- [`case_b_local_pressure_2x12x16_dt0p002.report.json`](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_local_pressure_2x12x16_dt0p002.report.json)
  and its [console log](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_local_pressure_2x12x16_dt0p002.raw.stdout.txt);
- [`case_b_local_pressure_2x24x32_dt0p002.report.json`](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_local_pressure_2x24x32_dt0p002.report.json)
  and its [console log](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_local_pressure_2x24x32_dt0p002.raw.stdout.txt);
- [`case_b_local_pressure_2x24x32_dt0p004.report.json`](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_local_pressure_2x24x32_dt0p004.report.json)
  and its [console log](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_local_pressure_2x24x32_dt0p004.raw.stdout.txt);
- [`case_b_local_pressure_2x36x48_dt0p002.report.json`](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_local_pressure_2x36x48_dt0p002.report.json)
  and its [console log](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_local_pressure_2x36x48_dt0p002.raw.stdout.txt);
- [`case_b_fbar_2x24x32_dt0p002.report.json`](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_fbar_2x24x32_dt0p002.report.json)
  and its [console log](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_fbar_2x24x32_dt0p002.raw.stdout.txt);
- [`case_b_fbar_2x36x48_dt0p002.report.json`](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_fbar_2x36x48_dt0p002.report.json)
  and its [console log](../examples/cardiac_benchmark/results/archive/truncated_polar/case_b/case_b_fbar_2x36x48_dt0p002.raw.stdout.txt).

Each JSON report retains the complete `p0`/`p1` histories, formulation and
mesh parameters, output-point locations, all accepted-step SNES/KSP and
residual diagnostics, NPZ and normalized-log hashes, exact verified `step_0B`
file identities, all-team mean curves, and full-precision RED values. The NPZ
archives and external CC BY 4.0 pickle files are not bundled.

“250/250” and “500/500” mean that every requested increment met the recorded
nonlinear acceptance rule before state commit. They do not assert spatial or
temporal convergence. RED is reported without a repository-defined cutoff.
Across the seven retained Case B configurations, CoupFE RED ranges from
0.5511448 to 1.0742661 at `p0` and from 0.6538962 to 1.2200494 at `p1`. Those
values are above the maximum participant-team RED retained in the fine report,
0.2875376036 at `p0` and 0.3701682113 at `p1`. This is a descriptive comparison,
not a pass/fail classification: RED against the all-team mean combines
amplitude differences with mesh-dependent signed point-component differences.
The checkpoint-`62ad760` 2×24×32, `dt=0.002 s` rows hold mesh, time step,
geometry, fiber and output sampling, source, and nonlinear solver fixed.
Relative to F-bar, Q1/P0 changes peak `p0` by
(+4.345, -0.420, -0.362) mm and peak `p1` by
(+3.932, +0.272, -0.628) mm. The peak-vector magnitudes are 4.380 and
3.991 mm, and the maximum vector-history differences are 4.382 and 3.991 mm.
Its RED is 0.00327 higher for `p0` and 0.05736 lower for `p1`. This controlled
formulation comparison reports what both paths produce at one discretization;
the mixed RED changes do not rank their accuracy.

The two 2×12×16 Q1/P0 records differ only in time step. The `dt=0.004 s`
loading grid is an exact subset of the `dt=0.002 s` grid. On those shared
times, the maximum/RMS vector-history differences are 0.144513/0.024943 mm at
`p0` and 0.115627/0.020255 mm at `p1`. This is a bounded two-step sensitivity
measurement, not a time-convergence claim. Thus the coarse pair is internally
close under time-step halving—its RED is about 0.88/0.87—while the retained
Q1/P0 2×24×32 rows are about 1.07/1.16 and the 2×36×48 row is 0.551/0.654.
Together with the signed-component changes, this shows coarse-grid temporal
repeatability but substantial spatial sensitivity in the signed point
components. It does not establish quantitative reproduction, spatial
convergence, a bifurcation, a global twist direction, or a solution branch.

A signed-component-aware recomputation from the canonical arrays retained in the
Q1/P0 2×36×48 report gives full-history component correlations with the
all-team mean of (0.939462, 0.837752, 0.969152) at `p0` and
(0.944736, 0.825150, 0.964841) at `p1`; component RMSE is
(5.432790, 0.358965, 2.161431) mm and
(5.613262, 0.244329, 2.245171) mm, respectively. At the common `t=0.48 s`
sample, its `p0` is (-36.078865, +2.027750, -14.014764) mm versus the mean
(-24.009661, +1.554442, -16.480099) mm, and the exact Zenodo participant
curves identified by that report all have `p0` signs (-, +, -). This supports
similar curve shape and signs in the retained 2×36×48 configuration while also
showing an amplitude gap. By contrast, the 2×24×32 `p0.u_z` history has
correlation -0.814 with the mean and the 2×12×16 `p0.u_x` correlation is
negative. These descriptive, post hoc checks explain why neither blanket
agreement nor blanket inconsistency is justified.

## Numerical execution record

Historical fine Q1/P0 runs can propose trials with nonpositive `det(F)`. The current
application catches only `InvalidDeformationError` from a residual evaluation
and lets the existing PETSc `bt` line search reduce the proposed nonlinear
trial correction. These events are correction overshoots outside the
admissible deformation domain, not mesh perturbations: all retained runs here
use `perturb=0`. Jacobian and unrelated errors still abort, and neither the
nonlinear tolerances nor the independently checked final residual rule changes.

The `e07993b` 2×24×32, `dt=0.004 s` run completed all 250 requested
increments under the same PETSc tolerance settings and final residual rule.
Its residual callback rejected 46 invalid Q1/P0 trial evaluations: 44 at step
132 (`t=0.528 s`) and two at step 133 (`t=0.532 s`). Backtracking found valid
trials, and every accepted/final state remained valid; none of the rejected
trials was committed. This is retained evidence that the narrow recovery path
operated, not evidence of accuracy, convergence, validation, or a bifurcation.

The Q1/P0 2×24×32 `dt=0.004 s` and `dt=0.002 s` grids are nested. On their
251 common times, maximum/RMS vector-history differences are
0.519200/0.132138 mm at `p0` and 0.549237/0.138402 mm at `p1`. Because the runs
identify different application checkpoints, this is a cross-checkpoint
sensitivity record, not a controlled one-factor time-step or temporal
convergence study.

The `e07993b` Q1/P0 2×36×48, `dt=0.002 s` run reached peak index 241 at
`t=0.482 s`, where mean Gauss-point `det(F)` was 0.985462950484. It rejected
168 invalid trial residual evaluations: 83 at step 277 (`t=0.554 s`) and 85
at step 279 (`t=0.558 s`). Backtracking found valid trials, and all 500
accepted/final states passed the unchanged checks. The largest
final-residual-to-threshold
ratio was 0.973921, below one. This retains a successful named execution and
its solver behavior without treating a rejected trial as accepted or changing
a tolerance.

A separate Q1/P0 3×36×48, `dt=0.002 s` attempt ended fail-closed with SNES
reason `-5` before state commit because backtracking did not yield a converged
increment. It produced no completed NPZ or comparison report and is therefore
absent from the retained-result table; the failed attempt is an execution
record, not numerical result evidence.

On 2026-08-01, the detached clean-checkpoint-`e07993b` Q1/P0 3×36×48 rerun
with `dt=0.001 s` was stopped by user direction after it printed accepted steps
1 and 2, so work could prioritize an MPI companion. It produced no completed
NPZ or comparison report. Its external partial console log has SHA-256
`d47ee0ce04a312db3808a2a3b373bb02a0c7f71ff5192f9dd88999e53b50554a`.
This is a user-directed interruption record, not a numerical failure or result.
The current serial driver and serial example were left unchanged. The archived
truncated-polar serial results remain regression evidence for their exact
historical configurations, not current Benchmark 1 validation evidence.

MPI support is additive, not a replacement for that serial reference. The
distributed companion preserves the historical truncated-polar Q1/P0/lumped/
CG1 contract and has separate closed contracts for the historical pointwise-
`kappa`/consistent-mass/backward-Euler path and the current condensed-Q1/P0/
consistent-mass/generalized-alpha path. The latter includes the physical-
coordinate Laplace/GP-direct frame and stages its velocity-dependent material
and Robin terms consistently. The closed path assembles complete owned mass
rows, retains rank and mass-partition provenance, and shares the benchmark load
history.

The historical pre-straight-wall/pre-frame fine 4×36×32 MPI-4 execution
completed 1,000 backward-Euler increments. The current clean `2x20x17` and
wall-only `4x20x17` executions each completed 320 generalized-alpha increments
on eight ranks. These establish successful executions of their named contracts;
neither is rank-equivalence or scaling evidence. Reduced historical paths keep
serial-versus-1/2/4-rank automated checks, but no completed current closed
serial/1/2/4/8-rank equivalence record is bundled.

The archived F-bar 2×24×32 peak has the same component signs as the
legacy-reported 2×24×32 row. Its recorded source, open-apex treatment, Hex8
sampler, and evidence record differ, so this qualitative consistency is not a
claim of reproduction or equivalence. Its positive `p0.u_z` is opposite the
negative `p0.u_z` in the recorded benchmark comparison row, so it occupies a
different mesh-dependent signed point-component response.

For 2×36×48, the archived F-bar peak `p0` is within 1.11593 mm, or 2.89%,
of the rounded legacy-reported (-35.86, +1.49, -14.22) mm vector. This is not an
exact or byte-for-byte reproduction. The old archive, log, and complete
configuration evidence are absent, and the archived run uses an open apex and
the Hex8 output sampler. The archived 2×36×48 peak has the same component signs
as the recorded benchmark comparison row, so it has a narrow qualitative sign
consistency that the 2×24×32 row does not. The archived 2×24×32 and 2×36×48
F-bar records are two spatial resolutions. Their peak-vector differences are
36.713 mm at `p0` and
36.740 mm at `p1`, while the maximum history differences are 50.677 and
48.298 mm. RED
changes from 1.0709970 to 0.5953528 at `p0` and from 1.2200494 to 0.6971874 at
`p1`. This is clear spatial sensitivity, not a spatial-convergence or accuracy
claim.

The `e07993b` Q1/P0 and `62ad760` F-bar 2×36×48 records can be inspected side
by side, but the application checkpoint and invalid-trial policy differ. Their
peak-vector differences are 0.956714 mm at `p0` and 0.555621 mm at `p1`, and
their maximum vector-history differences are 9.336698 and 8.242373 mm. Q1/P0
RED is 0.0442080 lower at `p0` and 0.0432912 lower at `p1`. These measurements
do not isolate formulation, rank accuracy, establish convergence, reproduce
the legacy F-bar row, or support a bifurcation claim.

## Remaining numerical questions

- Decide which finite epicardial surface is the intended benchmark target
  before another production run. The clean four-way split has already executed
  the coarse, wall-only, surface-only, and fine controls; the response shift is
  strongly associated with surface resolution rather than the count through
  the wall, and the Robin audit supplies a quantified, consistent mechanism.
  It does not select a continuum limit or prove convergence. The completed
  `2x20x17`/`4x20x17`, `tip_refine=6.0` full-cycle pair does not resolve that
  fine-surface-mesh question. It supplies a controlled through-wall
  comparison: the pair separates primarily in the snap window, then recovers
  to 0.0346/0.0518 mm at cycle end. Its mixed FEniCS phase metrics establish
  numerical timing sensitivity rather than convergence.
- The compact current-source prefix report is now packaged and hash-binds the
  external archives, metrics, and figure. The multi-megabyte archives and plot
  remain external. Its recorded pause decision is historical and superseded
  for the separately configured tip-refined full-cycle run. Historical
  polar-ring records remain explicitly
  truncated-domain results; the collapsed polar alternative remains
  degenerate.
- Any future claim that the current trajectory is rank independent requires an
  exact same-configuration rank gate through the snap window. The present
  eight-rank runs prove completion on eight ranks, not serial/1/2/4/8
  equivalence or parallel scaling.
- Do not use the fine native-quadrature `det(F)` minimum as evidence that Q1
  caused the trajectory gap. FEniCS has the same peak-pressure minimum to
  within 0.0055 and nearly the same wall-volume ratio. Future spatial studies
  should compare low-`J` reference-volume fractions, fields on a declared
  common sampling set, and landmark histories rather than ranking meshes by a
  single extremum.
- F-bar and Q1/P0 are different volumetric discretizations and should be
  compared as such; one should not be silently substituted for the other in a
  historical claim.
- The retained Q1/P0 parameterization uses `K=1.0e6 Pa` with material
  `kappa=0`; F-bar uses material `kappa=1.0e6 Pa`. A separate external
  2×12×16, `dt=0.002 s`, `t_end=0.6 s` diagnostic compared Q1/P0 at
  `K=0.1`, `0.3`, and `1.0 MPa`. All three completed 300/300 increments with
  zero invalid trials; total nonlinear iterations were 613, 614, and 645, but
  peak reference-volume-weighted `V/V0` was 0.889889, 0.958346, and 0.987143.
  Lowering `K` therefore changed volume far more than nonlinear effort. It did
  not show that `K` caused the fine-run failure, so the benchmark comparison
  retains `1.0 MPa`. The changed-`K` runs are diagnostics, not retained
  benchmark results; their external summary has SHA-256
  `13cec75a0134a8599197ed1f4c371fb32509c36817d4bb3d3b911663ad1108d3`.
- The optional Newmark path is not the source-matched generalized-alpha path
  for all rate-dependent terms and must remain separately named. The current
  closed, consistent-mass, physical-frame MPI generalized-alpha Step 0B path
  has completed the two clean 0.32 s gates and the separate eight-rank
  `2x20x17`, `tip_refine=6.0` full-cycle trajectory reported above.
- An open endocardial surface does not define a unique closed-cavity volume.
  Any reported volume needs an explicit, oriented, tested cap policy.
- Complete displacement curves, deformation-Jacobian ranges, fiber convention,
  twist measures, and external reference identity should be considered
  together rather than reducing Case B to one signed peak component.
