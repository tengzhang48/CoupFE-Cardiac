# Case B debugging postmortem: benchmark identity before solver tuning

This document records the failed hypotheses, null experiments, corrections,
and process mistakes from the CoupFE-Cardiac Benchmark 1, Step 0, Case B
campaign. The purpose is not to turn a completed run into a stronger result.
It is to make the investigation reproducible and to prevent the same ordering
mistake in future benchmark work.

> **2026-08-07 update:** later graded-mesh controls support numerical
> sensitivity near snap-through, not distinct mathematical or physical
> branches. The current-result section at the end records the completed
> two-/four-wall-layer full-cycle evidence; earlier pause language is retained
> only as investigation chronology. The clean four-layer archive replay is
> complete and its exact provenance is retained below.

The most important conclusion is procedural:

> A multi-code benchmark reproduction must first establish the identity of the
> physical problem: domain, boundary roles, loads, material parameters, and
> compared outputs. Discrete-method differences must then be declared and
> qualified before a solver trajectory is interpreted.

We did not follow that order. We investigated the volumetric formulation,
transmural coordinate, time step, time integrator, mass matrix, wall
resolution, fiber sampling, and viscosity before completing a face-by-face
audit of the geometry and boundary partition. The mass representation and the
geometry should both have been explicit at the beginning. Consistent mass was
eventually implemented and retained as the reference-aligned default, but it was not the
cause of the observed snap delay in its controlled test. The geometry audit,
performed still later, exposed a different computational domain and boundary
value problem.

That distinction is central to this postmortem: a design can be wrong without
causing the discrepancy currently under investigation, and a run can converge
without solving the benchmark problem.

## Executive findings

1. **The historical mesh was numerically usable but was not the paper's
   domain.** Its positive `apex_offset` removed 1.9335 mm of the ventricular
   tip, left a 2.672330 cm² annular free surface, removed terminal
   endocardial pressure and epicardial Robin facets, and reduced the
   undeformed pressure resultant by 4.5937%. Calling this only an “apex
   regularization” hid a topology and boundary-condition change.

2. **The mass matrix was an undeclared reference-method mismatch, not the
   measured root cause.** The inherited driver used row-summed mass although
   the audited nearby Simula/FEniCS implementation uses the consistent
   variational mass form. Checked Q1 Hex8 consistent mass therefore became the
   reference-aligned default; row-summed mass is not asserted to be
   mathematically invalid in general. On the named historical truncated-polar,
   2 ms comparison,
   consistent and lumped trajectories differed by no more than about 0.03 mm
   through the snap window and crossed the declared threshold at the same
   time.

3. **Two interim conclusions were false because log lines were compared at
   different times.** Newmark first appeared to advance the snap, and
   consistent mass first appeared to remove half the delay. Equal-time array
   comparisons reversed both conclusions. Console output is an execution
   trace, not a comparison grid.

4. **The retained FEniCS run should have been audited before changing the
   driver.** A nearby source clone, the completed displacement and pressure
   histories, parameters, and geometry were already on the machine, although
   the result directory was misleadingly named as Case A. Inspecting those
   assets earlier would have exposed the closed apex and the nearby source's
   consistent-mass, generalized-alpha, and P2 implementation before expensive
   campaigns. The exact producing source revision remains a provenance gap.

5. **The first historical closed-apex run improved the comparison, but that
   result was not an isolated geometry experiment.** Relative to the earlier
   open-domain, standard-`kappa`, 1 ms result, the dirty-tree coarse closed run
   reduced full-shared-history vector RMSE against local FEniCS from
   3.502/3.282 mm
   to 1.580/1.790 mm at `p0`/`p1`. Its onset was about 5--6 ms early rather
   than the open run's 13--14 ms late. Topology, mesh, mass, transmural field,
   and fiber evaluation changed together, so geometry remains a necessary
   correction but not the isolated cause of that improvement. This run also
   predates the straight Cartesian wall ruling and physical-coordinate frame
   reconstruction documented in the 2026-08-05 update below.

6. **The historical clean fine run did not show monotonic improvement.** The
   completed clean `(4,36,32)` MPI-4 run has 89,655 Q1 displacement DOFs and
   four Hex8 layers through the wall, nearly the same algebraic displacement
   size as the 89,700-DOF FEniCS P2 field. Its full-shared-history FEniCS RMSE
   is
   instead 7.596/7.236 mm, with relative L2 error 56.54%/52.50% and snap onset
   about 14.3/14.5 ms early. The historical dirty coarse and clean fine runs are not a
   controlled refinement pair because their source and execution provenance
   differ. Both full runs predate the straight-wall and physical-frame
   corrections. Their nonmonotonic errors therefore do not identify mesh
   refinement as the cause and must not be combined with the newer prefixes
   as a convergence sequence.

7. **The small peak-pressure deformation Jacobian is benchmark-like, not a
   Q1-only failure.** At 0.482 s the historical clean fine CoupFE minimum is
   0.327287,
   while a local reconstruction of the accepted FEniCS P2 state at its native
   degree-6 quadrature gives 0.332702. Their wall-volume ratios are 0.991334
   and 0.992210. The CoupFE minimum occupies a sharply localized epicardial-
   apex layer, not a globally collapsed or initially bad cell. CoupFE has a
   broader low-`J` tail, but the minimum alone cannot explain the landmark
   trajectory gap or establish a Q1 pathology.

   The retained fine Q1 Hex8 Case A result has only 8.566%/12.263%
   full-history relative L2 error at its two landmarks, compared with
   56.54%/52.50% here. Different benchmark cases cannot isolate an element-
   family effect, but this contrast is further evidence against “Q1 Hex8
   alone” as an adequate explanation for the Case B gap.

8. **The latest snap-window controls identify a strong mesh-dependent Robin
   mechanism, not an exclusive cause.** Corrected straight-wall,
   source-matched generalized-alpha runs through 0.32 s show that refining the
   surface changes the trajectory far more than refining only through the
   wall. The discrete normal-only epicardial Robin operator also loses most of
   its artificial long-axis rotational stiffness with surface refinement.
   This controlled association is diagnostic evidence; it is not a full-cycle
   validation, a convergence result, or proof that the Robin discretization
   explains the entire CoupFE/FEniCS difference.

The final curve is therefore not the most valuable product of the campaign.
The more durable product is a benchmark-identity gate, an experiment ledger,
and a clearer distinction among correctness, causality, convergence, and
validation.

## 2026-08-05 answer-first update: stop at the snap-window gate

**Do not extend the current Case B calculation to 1 s yet.** Four completed
eight-rank, 1 ms prefixes using the corrected straight Cartesian wall ruling,
consistent mass, condensed Q1/P0 mean-`log(J)` volume response, and the
source-matched generalized-alpha path isolate the large mesh effect to surface
resolution rather than through-wall resolution. The evidence points to a
strong discrete boundary-operator mechanism: planar facet normals give the
normal-only epicardial Robin support an artificial resistance to rigid
rotation about the ellipsoid's long axis. The mechanism is not established as
the only source of the cross-code trajectory gap, so changing the Robin law or
tuning its coefficients would be unjustified.

The controlled `a2006b7` mesh split completed through 0.32 s on eight MPI
ranks for `(n_t,n_core,n_radial)=(2,20,17)`, `(4,20,17)`, `(2,36,32)`, and
`(4,36,32)`. A controlled old-versus-corrected coarse comparison, validated
against the exact historical source checkout, produced `p0`/`p1` vector RMSE
of 0.039/0.037 mm and maximum differences of 0.145/0.131 mm over the prefix.
At 0.32 s, adding wall layers on the coarse surface changed the
`p0`/`p1` displacement vectors by only 0.669/0.635 mm. Refining only the
surface changed them by 9.994/8.772 mm. The same grouping appears when wall
layers are held at four: the surface change is 9.500/8.198 mm, while the
through-wall change on the fine surface is 0.881/1.217 mm. These are endpoint
sensitivities over a truncated prefix, not spatial error estimates.

For the rigid long-axis rotation field, the reference-configuration Robin
quadratic form is 1.44455 N m/rad on both the coarse and wall-only CoupFE
surfaces and 0.555828 N m/rad on both the surface-only and fine surfaces. The
local FEniCS P1 tetrahedral boundary mesh gives 1.85705 N m/rad under the same
read-only operator audit. In the smooth axisymmetric continuum, the
epicardial normal spring contributes zero stiffness to that rigid tangential
rotation; the nonzero values arise from faceted normals. The FEniCS number is
not an infinite-resolution target. The paper's nominal mesh target is about
`h=5 mm`, the local FEniCS mean epicardial edge is 3.880 mm, the coarse CoupFE
mean is 3.996 mm, and the fine CoupFE mean is 2.173 mm. The ten-team curves are
likewise finite-mesh results, not a continuum oracle.

The actual 0.32 s states reproduce the operator grouping. The elastic
epicardial Robin energy and internal long-axis moment are:

| Mesh | Elastic epicardial energy (J) | Internal long-axis moment (N m) |
|---|---:|---:|
| Coarse `(2,20,17)` | 0.308407 | -0.771763 |
| Wall-only `(4,20,17)` | 0.305358 | -0.774477 |
| Surface-only `(2,36,32)` | 0.116757 | -0.248906 |
| Fine `(4,36,32)` | 0.117304 | -0.250147 |

These state descriptors include the elastic spring contribution only, not the
dashpot contribution, and the motion is materially non-rigid. They support the
mechanism but do not prove sole causality.

The subsequent physical-frame correction reconstructs the toolkit's
structural directions from each physical point instead of treating mesh-
construction parameters as an inverse coordinate map. Its Gauss-point fiber
change is small over nearly all of the domain: the coarse/fine median changes
are 0.0197°/0.0205°, the 99th percentiles are 0.3860°/0.3217°, and the maxima
are localized near the apex at 2.8229°/2.6916°. Clean eight-rank coarse and
wall-only 0.32 s gates from application `056c02d` and Core `e2f42ed` completed
with this correction. Against their `a2006b7` counterparts, the maximum
landmark trajectory-vector change is 0.03566 mm. This is a required source-
fidelity correction, but it does not explain the large surface-refinement
response.

Unless a later paragraph explicitly says “corrected straight-wall,” the older
full “closed” results in this postmortem refer to the pre-`a2006b7` curved
through-wall mapping and pre-`9f5e0f1` frame reconstruction. They remain valid
records of those named executions, but they are not full-cycle evidence for
the current geometry and frame. At this point in the investigation the 1 s
extension was paused pending a targeted short diagnostic; the later graded
full-cycle result is recorded in the 2026-08-07 update above and the current
section below.

## What “reproduce Case B” has to mean

Curve similarity is the last check, not the definition of reproduction. A
candidate calculation has several independently reviewable layers:

| Layer | Question that must be answered before a full run |
|---|---|
| Source evidence | Do we have the paper, parameter record, reference code or input deck, mesh, output definitions, and comparison data? |
| Geometry | Is the physical domain closed in the same places, with the same units, axes, extents, and material-point coordinates? |
| Boundary partition | Is every exterior face assigned exactly once to the same physical role? Are there new free surfaces or missing loaded/supported surfaces? |
| Loads and supports | Do pressure orientation, resultant, moment, time evaluation, and Robin operators match? |
| Physical model | Are material parameters, volumetric energy, active/passive setting, density, and viscosity unchanged? |
| Spatial discretization | What element family, order, quadrature, mesh resolution, direction field, and mass representation are used? |
| Temporal discretization | What time grid, load-evaluation time, acceleration/velocity update, and rate definition are used? |
| Solver acceptance | Did each accepted increment meet a declared residual and admissibility rule without committing rejected state? |
| Output comparison | Are the same material points and components compared at the same physical times with a declared metric? |
| Evidence status | Was the run made from clean reviewed source, and are its commands, environment, diagnostics, and hashes retained? |

These layers support different statements:

- **setup identity** means the intended domain, boundaries, loads, and physical
  parameters have been matched;
- **a completed numerical execution** means the selected discrete equations
  met the recorded solver and deformation checks;
- **trajectory agreement** is a measured comparison for the selected outputs;
- **spatial or temporal convergence** requires a controlled refinement study;
- **benchmark equivalence** requires the remaining discrete-method differences
  to be either matched or bounded; and
- **physical validation** requires experimental or real-device evidence.

Passing one layer must not be used as evidence for another.

## The investigation in the order it actually happened

The following table records the main rounds. “Null” means that the named
change did not explain the observed gap in the tested configuration. It does
not mean the choice is universally irrelevant.

| Round | Hypothesis or action | What happened | Lesson |
|---|---|---|---|
| Volumetric formulation | The application-owned element-mean `log(J)` Q1/P0 path might represent the paper model. | The paper instead uses a pointwise `kappa` energy. In the named open-domain comparison, restoring that form removed the earlier axial overshoot; other refined and apex-offset runs still encountered failures. | Match the equation before tuning the algorithm. A stable alternative formulation is still a different model. |
| Transmural coordinate | Analytic layer coordinates might rotate fibers relative to the paper's Laplace field. | The difference was real and locally several degrees. The Laplace field changed parts of the response, but did not remove the snap-timing gap on the historical truncated polar mesh. | A real specification difference need not be the dominant cause. Adopt the correct definition anyway. |
| Time-step size | Backward Euler at 2 ms might delay the instability. | Halving to 1 ms left the historical truncated-polar onset and trajectory nearly unchanged. | A small step is not evidence that the temporal formulation matches; it is one controlled sensitivity result. |
| Integrator family | Backward-Euler damping might cause the delay. | An experimental constant-average-acceleration Newmark run was nearly coincident and marginally worse on the declared grid. It did not implement Simula's full generalized-alpha force staging or velocity-consistent rate. The first contrary reading came from mismatched log times. | Compare arrays at common physical times, and distinguish a kinematic probe from a matched reference integrator. |
| Reference-code forensics | Determine what the locally reproduced Simula/FEniCS contribution actually solves. | The retained run was Case B despite its directory name. It uses a closed tetrahedral domain, P2 displacement, consistent mass, generalized alpha, and a velocity-consistent viscous rate. | Inventory and inspect reference assets before launching candidate sweeps. Filenames are not provenance. |
| Mass representation | Row-summed mass might delay the inertially triggered snap. | Consistent Q1 Hex8 mass was implemented, unit-checked, and made the default. Its controlled historical truncated-polar trajectory differed from lumped mass by at most about 0.03 mm in the snap window; onset was unchanged. | Correct the design, but do not claim it caused a result change when the controlled experiment is null. |
| Through-wall resolution | Two Q1 layers might underresolve bending or twist. | Three layers did not move onset. They deepened the post-snap response, introduced axial overshoot, and first encountered post-snap nonlinear failures; a 1 ms run later completed. | Refinement can change different response features in different directions. Two meshes do not establish convergence. |
| Fiber evaluation | Historical nodal CG1 direction interpolation plus Gram--Schmidt might differ materially from analytic-rule evaluation at Hex8 Gauss points. | The direction change was measurable, but a controlled historical truncated-polar run did not move onset and slightly worsened the selected discrepancy metrics. The GP-direct option removes that historical interpolation; it does not reproduce Simula's rule-based P2 structural field. | Verify that a code change actually changes the intended field, then accept a null trajectory result without overstating method identity. |
| Viscosity | The backward-difference viscous term might dominate the snap. | Accepted-state decomposition measured a roughly 3--7% assembled material-force contribution in the snap window. `eta=0` accepted steps 1--12, then failed before committing step 13 at 0.026 s; `eta=50` was stopped before completion. | A force fraction does not predict the nonlinear trajectory change. Failed and interrupted controls are not trajectories and do not authorize parameter tuning. |
| Small apex-offset sensitivity | Reducing the truncation might test the apex hypothesis. | The smaller-offset mesh still had an open boundary and later failed. It did not test the paper's closed topology. | A smaller version of the wrong topology is not a controlled test of the right topology. |
| Detailed FEniCS comparison | Compare pressure, displacement phases, geometry, and boundary measures directly. | The pressure histories were close, while the displacement error concentrated in the snap window. The face audit then exposed the open tip and missing pressure/Robin surfaces. | When a curve gap persists, return to the problem statement before adding another numerical knob. |
| First closed Hex8 prototype | Close the apex without collapsed cells. | Eight Gauss points looked acceptable, but vertex/edge/face/center checks found near-singular square-corner mappings. The mesh was rejected before dynamics. | A mesh check must search the mapping's vulnerable locations, not only its current integration points. |
| Historical five-block closed mesh | Build a noncollapsed all-Hex8 domain with only endocardium, epicardium, and base boundaries. | The geometry, boundary, pressure, and Robin audits passed. A 1,000-step development run completed and substantially improved trajectory agreement. It predates the straight-wall and physical-frame corrections. | Pre-solve qualification makes the subsequent run interpretable; it does not by itself validate the curve. |
| Historical clean fine MPI-4 run | Test the pre-straight-wall closed geometry at four wall layers and approximately the FEniCS algebraic displacement size. | All 1,000 increments completed from clean source. The pre-snap curves agree closely, but CoupFE enters snap about 14.3--14.5 ms early and develops excessive post-snap negative `x` displacement. It also predates the physical-frame correction. | The numerical trajectory differs, but the cause is not isolated. Only a same-source sequence can attribute the change to resolution. |
| Straight-wall geometry correction | Match the toolkit's meridional construction instead of curving interior wall lines along interpolated ellipsoids. | `a2006b7` joins corresponding endocardial and epicardial points by straight Cartesian segments. A read-only old-versus-new coarse prefix comparison found the geometry correction itself had a small trajectory effect, but the correction is required for setup identity. | A source-fidelity correction can be causal-null for the current discrepancy and must still be retained. |
| Source-matched generalized-alpha mesh split | Run corrected-geometry `(2,20,17)`, `(4,20,17)`, `(2,36,32)`, and `(4,36,32)` controls to 0.32 s on eight ranks. | Through-wall-only endpoint changes were 0.669/0.635 mm at `p0`/`p1`; surface-only changes were 9.994/8.772 mm. The surface-only and fine curves formed one pair, while coarse and wall-only formed another. | The previously missing Step 0B generalized-alpha prefix now exists. Surface resolution, not wall-layer count, controls the large observed split over this window. |
| Faceted-Robin rotation audit | Apply an exact rigid long-axis rotation field to each reference-configuration Robin matrix. | Total stiffness paired exactly by surface mesh: 1.44455 N m/rad for coarse/wall-only and 0.555828 N m/rad for surface-only/fine; the local FEniCS boundary mesh gave 1.85705 N m/rad. Actual 0.32 s spring energies and moments showed the same pairing. | This is a strong discrete mechanism, not proof of exclusive causality. Preserve the benchmark Robin law and diagnose the discretization; do not tune the physical coefficients. |
| Physical-coordinate frame correction and clean gate | Replace stored-parametric frame reconstruction with the pinned toolkit's physical-point inverse and rerun benchmark-scale surfaces. | Field rotations were small and localized. Clean `056c02d`/`e2f42ed` coarse and wall-only eight-rank prefixes completed, and the maximum landmark trajectory change was 0.03566 mm. | Correct the structural field for source fidelity, record its identity, and accept the null trajectory result. A 1 s continuation is not warranted by this gate. |

This chronology should remain visible. Deleting the failed and null experiments
would make the final configuration look obvious when it was not, and would
invite later developers to repeat the same sweeps.

### Evidence-status ledger

The campaign produced artifacts with different evidentiary meanings. They must
not be collapsed into a folder count or described collectively as “runs.”

| Status | Records | Permitted use |
|---|---|---|
| Earlier completed development result | Dirty-tree coarse closed five-block Case B, 1,000/1,000 increments; historical truncated-polar trajectories | Compare only the named numerical configurations. The coarse closed result predates the straight-wall and physical-frame corrections, is not public-retained evidence, and the polar-ring results describe a truncated domain. |
| Historical clean full external result record | Pre-straight-wall, pre-physical-frame fine `(4,36,32)` Case B on MPI-4, 1,000/1,000 increments, archive SHA-256 `63d41a1f69dceaa8c1fe7f3c7d46a6de4e40c270a35977214de741e06fc580a3` | Supports only the named historical execution and its hash-gated full-trajectory comparison. It is not a current-geometry full result and does not establish rank equivalence, mesh convergence, or a bundled regression oracle. |
| Completed corrected-geometry mechanism controls | `a2006b7`/`e2f42ed` four-way coarse, wall-only, surface-only, and fine 0.32 s prefixes on eight ranks | Supports the surface-versus-wall sensitivity and its association with discrete Robin rotational stiffness over the retained prefix. It is not a full-cycle comparison, convergence proof, or exclusive causal attribution. |
| Completed current-source snap-window gate | Clean `056c02d`/`e2f42ed` coarse and wall-only 0.32 s generalized-alpha prefixes on eight ranks | Supports the physical-frame null result on the benchmark-scale surface: maximum landmark trajectory change no greater than 0.03566 mm. It does not validate the trajectory beyond 0.32 s. |
| Paused full-horizon extension | Pre-correction fine generalized-alpha execution stopped by direction at 0.600 s and retained as diagnostic; corrected current-frame 1 s continuation not launched | Records an explicit evidence-based stop decision. No current-geometry, current-frame 1 s trajectory exists. |
| Completed mechanism control | Coarse isotropic control | Shows that the selected coarse response and twist depend on anisotropy. It is not a Case B candidate trajectory. |
| Configuration-specific null sensitivity | Named open-domain mass, `dt`, experimental Newmark, and CG1-versus-GP-direct comparisons | Constrains those one-variable hypotheses only on the tested open configuration. |
| Failed before a completed result | `eta=0` before step 13 commit; first three-layer and reduced-offset-apex attempts | Records the exact solvability boundary. No completed trajectory or RED is available. |
| User-interrupted | `eta=50` after its last printed accepted state at 0.256 s; the named three-layer 1 ms diagnostic stopped by direction | Records that execution stopped without a completed result; it is neither convergence failure nor comparison evidence. |
| Rejected comparison or analysis | First viscous-evidence run with `t_end=0.32 s`; first viscous-work identity | Preserves a reproducibility or definition failure after the numerical execution itself passed its checks. It must not enter the accepted comparison. |
| Quarantined derived data | Existing local FEniCS point-stress arrays | May identify files for later repair, but cannot support a stress-agreement claim. |
| Incomplete reference rerun | Separate FEniCS duplicate ending at 0.124 s after the process died | Must not replace the completed 0.001--0.999 s FEniCS Case B record. |

The rejected `t_end=0.32 s` viscous-evidence comparison exposed a subtle input
reproducibility problem: the adaptive Radau pressure-history generator changed
shared-time pressure values by as much as 0.4921 Pa when only its requested
final horizon changed from 0.40 to 0.32 s. The material-force diagnostic was
not at fault, but the attempted one-variable comparison no longer held the load
byte-identical. The serial and MPI drivers now expose `--load-horizon`: a
shortened run can integrate the canonical 1 s schedule once and retain its
exact prefix through `t_end`. The archive and reports record the horizon, and
shared-time load samples must be byte-identical before a sensitivity or rank
comparison is accepted.

## Why geometry and mass were checked late

The late ordering was not caused by one missing test. It came from several
reinforcing assumptions and workflow biases.

### 1. A completed inherited run created false confidence in the setup

The historical polar-ring mesh had positive Jacobians, produced completed
trajectories, and resembled the ventricular shape. Those facts established
that it was executable. They did not establish that it was the benchmark
domain. We implicitly promoted “runs successfully” to “setup is substantially
correct.”

The positive `apex_offset` was treated as a harmless anti-degeneracy device.
Geometrically, however, it changed a point closure into an annular cut. The
change was topological before it was quantitative.

### 2. Small scalar measure errors hid a large semantic change

The historical truncated polar mesh differed from the retained closed reference by less than 0.72%
in wall volume and the three main boundary areas. Those small percentages made
the geometry look close. But the missing endocardial surface reduced the
   pressure resultant by 4.5937%, and the new annulus changed the boundary
   partition in a region that can affect the instability.

Volume and area totals cannot detect whether the same faces carry pressure,
Robin support, or zero traction. Boundary identity needs topology and marker
checks, not only scalar geometry summaries.

### 3. We audited numerical validity before benchmark identity

The early mesh questions were “Are the cells nondegenerate?” and “Can the
solver advance?” The correct prior questions were “Is the apex closed?” and
“Does every exterior face have the benchmark role?” A valid finite-element
mesh can discretize the wrong boundary-value problem perfectly well.

The same inversion affected mass. The inherited diagonal operator was
convenient and stable, so mass representation was treated as a later accuracy
choice. In an implicit dynamic benchmark it is part of the discrete governing
equations and belongs in the initial method manifest.

### 4. Code-proximity bias favored easy experiments

Changing `dt`, selecting an integrator, switching a fiber option, or changing
a material parameter is easy to express as a flag and easy to run. Rebuilding
the apex, enumerating every exterior facet, tracing the reference mesh, and
checking resultant and moment require more source forensics and implementation
work. The experiment sequence therefore followed what was easiest to vary,
not what had the highest diagnostic priority.

This bias is stronger in agent-assisted work because automation can generate
and launch plausible sweeps quickly. More throughput does not improve the
ordering of the questions. Without an explicit gate, it can make a weak
diagnostic strategy look productive.

### 5. Solver symptoms focused attention on solver mechanisms

The visible discrepancy occurred during a dynamic snap, and some fine runs
encountered line-search or linear-solver failures. That naturally directed
attention toward damping, time integration, viscosity, bulk modulus, and mass.
Those were plausible hypotheses, but the symptom did not prove the defect was
inside the solver. A different pressure surface and support surface can create
the same kind of timing and numerical-trajectory difference.

### 6. Reference assets were discovered rather than inventoried

A nearby Simula/FEniCS source clone, a complete local FEniCS execution, the
benchmark toolkit, and the paper data were present locally. They were not
placed in one evidence inventory at the start. The completed Case B result was
also stored under a directory name containing “CaseA”; only its parameter
record established the actual case. The clone's revision was not retained as
the contemporaneous producer identity for that output.

This made reference-code inspection appear like a later research task when it
should have been step zero. A filesystem inventory and content-based identity
check would have been cheaper than any full trajectory.

### 7. Early success criteria emphasized curves instead of invariants

The campaign initially focused on reducing RED and moving the snap time into
the team band. Those are useful outcome measures, but they do not tell us
whether a closer curve came from a faithful setup or compensating errors. The
geometry, boundary, mass, load-time, and field-representation invariants should
have been acceptance criteria before any curve metric was computed.

### 8. We allowed each plausible null result to nominate another run

The single-variable discipline was useful once an experiment began, but the
hypothesis queue itself was not gated. A null `dt` result led to integrator; a
null integrator result led to mass; a null mass result led to resolution and
fibers. We accumulated careful answers to lower-priority questions while the
highest-priority setup question remained open.

The correction is not to abandon controlled experiments. It is to control the
order in which hypotheses are allowed to consume full-run time.

## What the FEniCS run establishes

The local reference is valuable because it is a completed execution of one of
the benchmark contributions, with full parameter and displacement state
retained. Its `parameters.json` declares Case B even though the containing
directory is misnamed. The trusted comparison inputs are its parameters,
pressure history, time stamps, and `p0`/`p1` displacement histories. A nearby
audited FEniCS source clone establishes the method implementation described
below, but the exact producing source revision was not retained
contemporaneously and must not be inferred from that clone.

The retained output plus the audited nearby source establish the following
method comparison, subject to that producer-revision provenance gap:

| Choice | Local Simula/FEniCS run | Historical clean fine CoupFE run |
|---|---|---|
| Domain | Closed ventricular wall with exactly endocardial, epicardial, and base boundary roles | Same intended continuum domain and roles, represented by a five-block Hex8 mesh |
| Mesh and spatial field | Quadratic Lagrange P2 vector displacement on 17,648 tetrahedral cells; 4,575 mesh vertices and 89,700 vector DOFs | Trilinear Q1 vector displacement on 23,616 Hex8 elements; 29,885 nodes, 89,655 vector DOFs, and four elements through the wall |
| Quadrature | UFL form-compiler quadrature degree 6 | Eight-point `2×2×2` Hex8 material integration |
| Volumetric energy | Paper pointwise `kappa=1.0e6 Pa` form | Paper pointwise `kappa=1.0e6 Pa` form |
| Mass | Consistent variational mass | Consistent Q1 Hex8 mass |
| Time method | Generalized alpha, `alpha_m=0.2`, `alpha_f=0.4` | Backward Euler |
| Load evaluation | Prescribed load evaluated at `t-alpha_f*dt` | Prescribed load stored/evaluated on the backward-Euler time level |
| Viscous rate | Rate from generalized-alpha velocity | Backward strain difference |
| Transmural/fiber field | Laplace coordinate solved in CG1 and interpolated to P2; rule-based directions normalized at P2 DOFs and used as P2 structural fields | Q1 Laplace coordinate plus analytic-rule frame evaluated and re-orthonormalized at Hex8 Gauss points |
| Pressure and Robin forms | Follower pressure; full-vector base and normal epicardial Robin terms | Same continuum forms, independently audited |
| Compared outputs | P2 displacement at `p0` and `p1` | Q1 Hex8 displacement at the same material points |

The CoupFE column describes the pre-straight-wall, pre-physical-frame full
execution retained for historical comparison; it is not the current setup.
Both completed calculations are displacement-only penalty formulations. The
Simula element is quadratic P2 displacement on tetrahedral cells, not a mixed
Q2/P1 displacement-pressure element. The closed CoupFE calculation uses the
standard pointwise-`kappa` Q1 Hex8 path, not the application-owned Q1/P0
`local_pressure` formulation used by some historical truncated-polar runs.

The corrected CoupFE geometry is close in integral measures, but it is not the
same node-for-node mesh. The previously reported -0.0828%, +0.0695%, +0.0090%,
and -0.0164% wall-volume and boundary-area differences, and the 0.1028%
analytic pressure-resultant error, belong to the earlier pre-straight-wall
coarse closed mesh.
The fine mesh independently passed the same geometry, face-ownership,
pressure, and Robin gates, with all retained-reference integral differences
within 0.139%. Those are setup checks and declared discrete-geometry/load
differences, not proof that the two spaces have identical instability spectra.

The two fiber representations are also related but not identical. A read-only
cross-mesh diagnostic on the earlier coarse mesh evaluated the retained Simula
P2 fields at all 28,160 CoupFE material Gauss-point coordinates. The fiber
axial-angle difference had
0.834° RMS, 0.727° volume-weighted RMS, and 2.09° at the 99th percentile. A
localized near-apex point reached 43.89°, where the interpolated Simula P2
fiber norm was also unusually small; the dynamic effect of that outlier is not
isolated. The CoupFE GP-direct option therefore removes the historical nodal
CG1 direction interpolation, but it must not be described as matching the
Simula P2 structural-field definition.

The FEniCS pressure peak is 16,073.142 Pa at the 0.482 s label; CoupFE records
16,070.954 Pa at the same label. Their same-label pressure RMSE is 27.813 Pa.
Shifting the sampled FEniCS pressure array by its known
`alpha_f*dt=0.4 ms` evaluation offset reduces that pressure RMSE to 0.591 Pa.
This explains the array-label discrepancy. Its effect on the nonlinear
trajectory has not been isolated, so this observation cannot exclude it as a
contributor to the displacement difference.

The existing FEniCS point-stress arrays are not trusted comparison oracles.
Their supplied postprocessor does not reload the accepted velocity and
acceleration state into the reconstructed problem, projects to DG1, and uses a
dimensionally inconsistent von Mises expression. This does not invalidate the
FEniCS displacement solve. It means that the current comparison must stop at
pressure and displacement until the stress path is corrected and matched at
declared element-interior or quadrature locations.

## Why differences remain after correcting the apex

We can identify known differences, but we cannot yet assign the remaining
trajectory gap to one proven cause.

### What is measured

The full-history numbers in this subsection describe the historical
pre-straight-wall, pre-physical-frame fine archive. They remain useful for
auditing that execution, but they must not be reported as current-setup
accuracy. That archive and the local FEniCS trajectory have 999 exactly shared
native samples from 0.001 through 0.999 s. No time interpolation is used. The
fine CoupFE vector RMSE is 7.5955/7.2363 mm at `p0`/`p1`, and full-history
relative L2 error is 56.54%/52.50%. Before 0.20 s the RMSE is only
0.1689/0.0551 mm; through the 0.20--0.32 s snap window it grows to
10.0531/9.3429 mm. The first downward `u_z=-5 mm` crossing is
14.293/14.506 ms earlier in CoupFE. At 0.320 s the `p0` and `p1` `x`
components are -27.728/-27.790 mm in CoupFE versus -16.064/-17.042 mm in
FEniCS. The dominant discrepancy is therefore the snap transition and
post-snap negative `x` response, not the sign of the smaller `y` component.

Against the exact ten-team mean on the canonical 10 ms grid, the historical
clean fine vector RMSE is 7.5678/7.2656 mm and 21.33%/33.00% of component
samples fall
inside the pointwise participant envelope. These are coverage fractions, not
percentage displacement errors, and no official pass/fail threshold is
assigned. The hash-gated full comparison report has SHA-256
`602a1e904973d9700b6cd99f1b922fa81651815c7722c157c9a571bd9d7640f5`;
the earlier 0.32 s prefix is a different artifact and is not a full-cycle
comparison candidate.

For historical context, the dirty coarse closed run had FEniCS RMSE
1.5803/1.7903 mm, relative L2 error 11.76%/12.99%, and onset 5.3/5.6 ms early.
The fine clean result is much worse after snap. This is informative but not a
controlled convergence observation: the application revision, clean-tree
status, serial/MPI execution, and every mesh direction differ. A clean coarse
rerun at the fine source revision is required before assigning the change to
resolution.

The peak-pressure deformation audit also has to be interpreted as a field,
not a single minimum. At 0.482 s CoupFE has `det(F)=0.327287` at its worst
native Hex8 Gauss point. Only 48 of 188,928 samples are below 0.4, spanning 24
outer-layer elements near the epicardial apex and about 0.01003% of reference
volume. The worst reference Jacobian condition is 1.468, its principal
stretches are 1.0661, 0.7676, and 0.4000, and the global wall-volume ratio is
0.991334. This is localized compression, not inversion or global collapse.

A 2026-08-04 read-only reconstruction from the hash-identified accepted
FEniCS P2 state gives `det(F)_min=0.332702` and wall-volume ratio 0.992210 at
its native degree-6, 24-point tetrahedral quadrature. The HDF5 input SHA-256 is
`68fdbb222bf948aebcc6cdc287450601efaed9eec0aadaf4abbdbc9eaccb74d3`.
The native extrema are in the same apex region but are not the same physical
point. CoupFE has a broader low-`J` tail: the reference-volume fractions below
0.4 and 0.5 are 0.01003% and 0.03431%, versus reconstructed FEniCS values
0.000247% and 0.001727%. That distribution may accompany the different
post-snap numerical response; no causality has been isolated. The reconstruction script
and derived JSON are not yet packaged, so this is a reproducible-input local
diagnostic rather than a checked-in regression gate.

CoupFE and the audited nearby FEniCS source implement the same paper penalty
energy `kappa*(J**2-1-2*log(J))/4` with `kappa=1.0e6 Pa`; the bulk modulus is
not doubled in the CoupFE fine run. The penalty strongly resists volume change
but does not enforce `J=1` pointwise. The near-equal cross-code minima and
global volumes rule out the minimum itself as an explanation for the much
larger landmark-trajectory difference.

### What has a supported causal interpretation

- The old open domain and boundary partition were wrong for the benchmark.
  The earlier pre-straight-wall coarse closed configuration had substantially lower error than
  the open configuration. Because topology/mesh, mass, transmural field, and
  fiber evaluation changed together, neither the geometry contribution nor
  any interaction is isolated.
- The FEniCS generalized-alpha evaluation offset explains nearly all of the
  pressure-array label difference after a 0.4 ms shift. Its effect on the
  nonlinear displacement trajectory has not been tested independently.
- The historical clean fine CoupFE and reconstructed FEniCS states have nearly the same
  peak-pressure minimum `det(F)` and global wall-volume ratio. A low native-
  quadrature minimum is therefore not a supported explanation for the curve
  gap, although their low-`J` distributions differ.
- On the old open domain, mass choice, time-step halving, Newmark versus
  backward Euler, and CG1 versus direct Gauss-point fibers did not materially
  move the declared onset. These are configuration-specific null results.
  They do not prove those choices are irrelevant on the corrected closed
  domain or in a coupled interaction.
- On the corrected closed domain, the old curved-through-wall and straight-
  wall coarse generalized-alpha prefixes differ little. This is a controlled
  null for the trajectory effect of that geometry correction over 0.32 s, not
  permission to retain the source-inconsistent geometry.
- In the four-way straight-wall split, the large trajectory change follows
  surface refinement, while the discrete Robin rigid-rotation stiffness and
  actual elastic epicardial energy/moment pair in the same way. This supports
  a mesh-dependent faceted-boundary mechanism. It does not prove that the
  mechanism is the only cause of cross-code disagreement.
- Reconstructing structural frames from physical coordinates changes the
  clean coarse/wall-only landmark trajectories by no more than 0.03566 mm
  through 0.32 s. This bounds that correction on the tested surfaces; it does
  not establish a full-cycle null.

### What remains unresolved

The source-matched generalized-alpha and velocity-consistent rate path is no
longer an unrun control: it completed on all four `a2006b7` split meshes
through 0.32 s and on the clean current coarse and wall-only meshes. This
removes the stale temporal-method experiment from the immediate plan. It does
not make the CoupFE and FEniCS discretizations identical. Condensed Q1/P0
Hex8 versus P2 tetrahedral displacement, surface topology, quadrature, and
structural-field representation remain declared differences, and the local
FEniCS producing revision remains a provenance gap.

The controlled surface split and Robin audit provide a quantified mechanism
consistent with the finer surface moving farther from the retained FEniCS curve
without implying that the finite-element solution became intrinsically less
accurate. Both codes discretize
the benchmark's normal-only support with facet normals, and their different
surface meshes generate different artificial rotational stiffness. The local
FEniCS and ten-team curves are finite-mesh benchmark evidence at approximately
the paper's `h=5 mm` scale, not a continuum target against which every finer
CoupFE surface must move monotonically. This mechanism is strong but not
exclusive; it does not bound the remaining Q1/P0-versus-P2, quadrature, or
nonlinear-path sensitivities.

The four-way split is a controlled prefix study, not a convergence sequence.
The clean rerun discussed in this historical section intentionally retained
the benchmark-scale surface and tested only the physical-frame correction. At
that time no corrected-straight-wall, current-frame full-cycle result existed,
and the old full archives could not fill that role. The next run at that stage
was required to be another predeclared short
diagnostic that isolates one remaining spatial or boundary representation
without changing the physical Robin law or its coefficients. Resume a 1 s
execution only if such a gate changes the interpretation enough to make the
full horizon decision-useful.

## Deeper lessons from the failures

### Correctness, causality, and accuracy are different questions

Consistent mass is the clearest example. Row-summed mass was not the reference
choice and had no documented justification as the default for this implicit
dynamic application. Making consistent mass the reference-aligned default
removed an undeclared method mismatch; it was not a claim that lumped mass is
universally wrong. The controlled comparison then showed that mass choice was
not the cause of the measured open-domain snap lag in that configuration.
Neither statement weakens the other.

The same reasoning applies to the Laplace transmural field and GP-direct fiber
option. They remove avoidable historical approximations, but GP-direct Hex8
fibers still do not match Simula's P2 structural-field representation. A
closer implementation can remain a declared method difference. We should not
keep an unjustified setup merely because a compensating error happens to fit a
curve.

### A converged solve can provide strong evidence for the wrong question

The open-domain runs met their nonlinear residual and deformation rules. They
are legitimate records of that truncated model. Their convergence does not
make them Case B executions. Solver verification answers “Did we solve these
discrete equations?” Benchmark verification first asks “Are these the intended
equations on the intended domain?”

### Null results need narrow language

“Mass was eliminated” is too broad. The evidence says that changing from the
historical lumped mass to the checked consistent mass did not materially alter
one named open-domain trajectory through the comparison window. The same
narrow wording is required for `dt`, Newmark, fibers, and transmural
coordinates. Nonlinear interactions and the corrected domain prevent a
universal conclusion.

### Failed runs are useful only when their failure boundary is exact

The `eta=0` run accepted steps 1--12, failed before committing step 13 at
0.026 s, and produced no completed result.
The `eta=50` run was interrupted and produced no completed result. The first
three-layer attempt and the reduced-offset apex attempt reached parts of the
snap but later failed. These records constrain solver behavior and guide the
next experiment. They are not partial evidence of a completed trajectory and
must not be plotted as if they were.

### A diagnostic can be mathematically wrong even when its code runs

The first viscous-work audit omitted a finite-strain geometric term and failed
its own identity. The check was corrected rather than its tolerance being
weakened. The FEniCS stress postprocessor provides a second example: files and
plots exist, but the reconstructed state and invariant are not suitable for a
quantitative comparison. Evidence code needs the same review as solver code.

### Visual and log agreement are hypotheses, not measurements

Twice, nearby console lines created a plausible but false time-shift story.
The correction was simple: load the completed arrays, define a common grid,
interpolate only under a declared policy, and compute onset and vector errors
programmatically. A future campaign should generate this comparison
automatically as soon as a result closes.

### Mesh quality must be adversarial

The rejected single-block closed-apex mesh passed checks at the eight material
Gauss points but was nearly singular at mapped square corners. Testing only the
locations already used by the current operator can miss fragile geometry that
may be exposed by changed quadrature, interpolation, refinement, or boundary
evaluation. Extended vertex, edge, face, center, and Gauss sampling is now a
pre-solve requirement, while still being a finite sample rather than a proof
over every natural coordinate.

### Interpret deformation Jacobians as fields, not verdicts

A minimum `det(F)` is meaningful only with its physical time, sampling rule,
location, affected measure, and global volume. In the historical clean fine
Case B run,
`det(F)_min=0.327287` initially looked severe. The small values are confined to
24 outer-layer elements near the epicardial apex and only 0.01003% of
reference volume is below 0.4. The corresponding reconstructed FEniCS P2 state
has `det(F)_min=0.332702` and nearly the same global wall-volume ratio. Native
Hex8 and tetrahedral quadrature extrema are not values at a shared physical
point, so they support a scale-and-localization comparison rather than
pointwise equality.

The coarser CoupFE minimum of about 0.66 is consistent with smearing this thin
layer, but the present evidence does not prove that interpretation. Conversely,
the broader fine CoupFE low-`J` volume fraction may accompany the different
post-snap trajectory without causing it. Report a distribution, reference-volume
fractions, a global volume measure, and a common-time field comparison before
using one extremum to diagnose inversion, locking, an element pathology, or a
trajectory error.

### More experiments do not compensate for a missing specification

The campaign was disciplined about one variable per run, yet still used the
wrong global order. Controlled experimentation is powerful only after the
baseline has been proven to represent the target. Otherwise it produces
precise knowledge about a surrogate problem.

### Agent assistance needs explicit epistemic roles

AI agents made source inspection, implementation, testing, and long-run
monitoring faster. They also made it easy to continue from one plausible
numerical hypothesis to the next. The corrective design is to separate roles:

- a **specification auditor** maps every paper statement and reference-code
  choice to the candidate setup;
- a **geometry and boundary auditor** proves topology, labels, measures,
  resultants, moments, and support operators before execution;
- an **experiment owner** declares one hypothesis, baseline, changed variable,
  falsifying observation, and stop rule;
- a **result auditor** compares completed arrays at common physical times and
  distinguishes fact from inference; and
- a **release auditor** verifies clean source, retained provenance, licenses,
  and claim boundaries.

The same agent may perform several roles, but the artifacts and review points
must remain separate. A count of completed runs is not a progress metric.
Reduction in uncertainty and closure of specification gaps are.

### Owner memory and disagreement are evidence triggers

Several decisive corrections began when the owner challenged an apparently
reasonable conclusion: whether Case B had already been run, whether the tip
was actually the paper geometry, and whether a design choice was being called
a difficult numerical feature. Such disagreement should trigger a source and
artifact audit, not a rhetorical defense of the current model. The reference
files, not confidence on either side, decide the issue.

## The setup-first protocol adopted from this case

Future benchmark work should use the following ordered gates.

### Gate 0 — evidence inventory

Before editing the solver, locate and hash the paper, supplemental data,
reference source or input decks, meshes, parameter records, completed local
runs, output definitions, licenses, and prior notes. Determine identity from
file contents rather than directory names. Record what is missing.

### Gate 1 — executable specification

Create a machine-readable manifest of units, coordinate frames, geometry
parameters, boundary roles, loads, material parameters, mass, spatial and
temporal methods, quadrature, field representations, solver acceptance, and
output points. Mark every entry as matched, intentionally different, or
unknown.

### Gate 2 — geometry and boundary qualification

Before kernel compilation or time integration:

- check units, axes, extents, topology, and reference point positions;
- enumerate all exterior faces and require exactly one physical role per face;
- reject unexpected free, duplicate, nonmanifold, or interior-labeled faces;
- check Jacobians at Gauss and extended natural-coordinate samples while
  stating the finite-sampling limitation;
- compare pointwise distance to the reference implicit surface, not only
  extents and integrated measures;
- compare wall volume and each marked surface measure;
- define the same virtual base/cavity closure when a cavity or wall volume is
  reported;
- check signed pressure resultant, transverse force, and moment; and
- check Robin support/damping symmetry and active degrees of freedom.

No full benchmark run should start while this gate is unresolved.

### Gate 3 — discrete-method manifest and gap qualification

Record element family/order, mesh resolution, quadrature, volumetric form,
mass representation, time integrator, load-evaluation time, viscous-rate
definition, transmural field, fiber interpolation, and output interpolation.
For a multi-code comparison, classify differences as matched, bounded by a
study, or unresolved; exact method identity is required only when it is the
declared objective of a direct code-to-code test. Correct unjustified
differences before using them as sensitivity variables.

### Gate 4 — component and broken-control checks

Check material stress/tangent, volumetric response, mass symmetry, positive
definiteness, row sums, integrated total mass, rigid-acceleration response, and
partition consistency. Check pressure tangent and reversed normal, Robin
matrices, state commit, output-point reconstruction, and malformed geometry.
Every important guard needs a broken control that demonstrates it can fail.

### Gate 5 — short execution and automatic comparison

Run a two-step smoke to prove composition and metadata, then a short window
that crosses the feature under study. Compare arrays at exact common times;
do not draw conclusions from log tails. Confirm that the changed option is
recorded and that invariant inputs stayed identical. Set an explicit canonical
full horizon, retain it with the result, and require byte-identical shared-time
samples when shorter trajectories are sliced for comparison.

### Gate 6 — controlled full run

Declare the baseline, one changed variable, expected observation, falsifying
observation, completion rule, and interpretation boundary before launch.
Retain failed and interrupted runs under their exact status. Do not tune
physical parameters during a reproduction run.

### Gate 7 — convergence and method-gap study

Only after setup identity and one completed trajectory should spatial,
temporal, element-order, and method comparisons be interpreted. A pair is a
sensitivity check; convergence requires a sequence and a declared norm.
Predeclare the multi-metric comparison suite, including vector RMSE, relative
L2, RED with its near-zero-reference warning, onset, recovery, and envelope
coverage. Simpler patch or Land-type tests may diagnose a known operator gap,
but they do not replace the setup-matched Case B comparison.

### Gate 8 — evidence and release

Retain the exact source/Core revisions and tree states, command, environment,
input identities, all-step solver diagnostics, output hashes, comparison code,
and limitations. A dirty-tree run remains development evidence even if the
source is committed later. Rerun from the clean revision rather than relabeling
metadata.

## Design changes resulting from the postmortem

The public application now encodes several of these lessons:

- historical truncated and benchmark-target closed topologies have different
  names and metadata;
- the closed mesh uses a noncollapsed five-block construction and the
  benchmark toolkit's straight Cartesian wall ruling;
- a pre-solve geometry audit checks face ownership, labels, extended
  Jacobians, measures, and extents;
- Case B pressure audit checks signed resultant, transverse fraction, and
  moment, including a reversed-load broken control;
- Robin operators are checked before integration;
- consistent Q1 Hex8 mass is the default, while historical commands request
  lumped mass explicitly;
- mass, viscosity, parameter variant, fiber policy, and portable transmural
  field identity are retained in result metadata;
- the closed-mesh structural frame is reconstructed from physical coordinates
  by one shared serial/MPI implementation, and its reconstruction identity is
  retained in result metadata;
- the pre-solve Robin audit reports rigid translation and rotation quadratic
  forms so mesh-dependent boundary restraint is visible without changing the
  benchmark law;
- Q1/P0 invalid-deformation trials are guarded and failed states cannot reach
  the completed-result writer; the distributed pointwise-`std-kappa` batch
  still needs an equivalent explicit trial-state `det(F)>0` guard;
- the existing FEniCS stress arrays are quarantined from quantitative claims;
  and
- public retention requires a clean-source run rather than a metadata edit.

These changes do not prove that the remaining trajectory differences are
resolved. They make future differences interpretable.

## Current result and remaining work

The current accepted numerical evidence now includes eight-rank
`2x20x17` and `4x20x17 tip_refine=6.0` full cycles (1,000/1,000 increments
each). Both use 1 ms, source-matched
generalized alpha, consistent mass, condensed Q1/P0 mean `log(J)`, the
canonical 1 s load horizon, and unchanged paper physical parameters.

The current interpretation is numerical sensitivity near snap-through, not a
physical branch. Relative to the two-layer full cycle, the four-layer control
has pre-snap RMSE 0.009831/0.002622 mm, reaches 1.633/1.406 mm maximum
separation in the snap window, has full-cycle RMSE 0.407604/0.361832 mm, and
returns to 0.034598/0.051776 mm separation at cycle end (`p0`/`p1`). Its
FEniCS full-shared-history RMSE improves from 1.093/1.169 to 1.017/1.073 mm,
and its
maximum snap-window gap improves from 5.552/5.024 to 4.188/3.753 mm, while
relaxation RMSE worsens 16.2%/14.9%. This phase-mixed outcome supports
transient numerical/timing sensitivity without establishing convergence. The discrete Robin
association and physical-frame null result remain useful bounded diagnostics,
not exclusive explanations.

The first four-layer full-cycle candidate is an independently audited,
non-retained dirty-tree provenance diagnostic. A clean isolated replay is
complete under Core `454f73c` and runtime-source SHA-256
`f8d9469a101709d11460a0803b6c031001192ccce61c472d032355b01070da05`;
its states match to roundoff. The retained archive records **application
`2458e7c`, NPZ
`1e333b29b05f01dedce9272b32b82ce6ccfda56036c1c8f57eb395b9b4494800`,
stdout `0810a9d0c944345464f48917c707a7d3df7187b543c74efd4b0384542b95aed7`,
elapsed `1778.4 s (29.6 min)`**.

The historical full fine archive with SHA-256
`63d41a1f69dceaa8c1fe7f3c7d46a6de4e40c270a35977214de741e06fc580a3`
remains a valid record of its pre-straight-wall, pre-frame execution. Its
localized pressure-peak `det(F)` comparison is likewise historical evidence;
it cannot substitute for the current full trajectory. The current full-cycle
run is a comparison record, not a validation or convergence claim.

The next research actions, in priority order, are:

1. keep the checked-in compact
   [`Step 0B prefix report`](../examples/cardiac_benchmark/results/step0b_case_b_clean_frame_0p32.report.json)
   as the reviewed binding for the external prefix archives, comparison
   figure and metrics, Robin operator/state audit, exact revisions, commands,
   and hashes;
2. retain the completed `4x20x17 tip_refine=6.0` prefix as the through-wall
   continuity record and do not repeat it without a source or method trigger;
3. preserve the completed clean four-layer archive and do not promote the
   dirty candidate;
4. repair and independently verify FEniCS stress reconstruction before any
   element-interior or quadrature-point stress comparison; and
5. reserve convergence and physical-validation language for proper sequences
   and experiment-matched or real-device evidence, respectively.

The exact run status, hashes, metrics, and release boundary are maintained in
[`CASE_B_STATUS.md`](CASE_B_STATUS.md). The broader application rules are in
[`lessons_learned.md`](lessons_learned.md), and the runnable interface is in
[`API.md`](API.md).
