# CoupFE-Cardiac benchmark test matrix

This document organizes the evidence required for a new CoupFE-Cardiac
implementation. The benchmark paper, participant-team curves, toolkit data,
retained local FEniCS reference, and Simula calculation are reference
evidence. They do not need to be rerun merely to exercise established external
software. The validation burden is on CoupFE: identify external inputs by hash,
then test the new material, elements, operators, distributed assembly, and
complete trajectory progressively.

This is application-owned verification. CoupFE Core remains a small generic
finite-element/code-generation library; cardiac constitutive choices,
benchmark geometries, reference-data contracts, solver profiles, and retained
scientific evidence live in CoupFE-Cardiac.

## Test layers

| Layer | Purpose | Required evidence | When to rerun |
|---|---|---|---|
| 0. Source and data identity | Prove what code and reference data are being used | Clean application/Core revisions; paper, toolkit, curve, Laplace-field, canonical node-coordinate, and element-connectivity hashes; field-to-mesh binding; structural-frame reconstruction identity | Any source or input identity change |
| 1. Material point | Check the Holzapfel--Ogden, active, volumetric, and viscous laws without FE assembly | Stress-free state; energy/stress consistency; complex-step versus finite-difference tangent; rotated fibers; nonzero active and viscous states | Material or code-generation change |
| 2. Element | Check Hex8 interpolation, quadrature, state layout, and formulation-specific tangent | Affine patch; distorted-cell tangent; positive reference Jacobian; accepted-state commit; separate F-bar, `std-kappa`, and local-pressure identities | Element, quadrature, state, or formulation change |
| 3. Operators | Check mass, Robin support, follower pressure, Laplace field, point sampling, and condensed pressure independently | Analytic resultants; symmetry; row sums; finite differences; exact affine sampling; boundary audit; rigid-translation and long-axis-rotation `q^T K q` forms, with base full-vector and epicardial normal-only Robin contributions separated | Operator or geometry-boundary code change |
| 4. Simple CoupFE application | Expose coupling mistakes before the ventricle | Uniform-fiber block with the intended formulation; zero-load preservation; active contraction; positive `det(F)`; a mesh sequence only when the simple result is not already known | Coupling, solver, or time-method change |
| 5. Ventricular setup | Prove that the intended continuum problem was discretized | Closed geometry; exactly endocardial, epicardial, and base roles; exact mesh and Laplace-field identities; declared physical-coordinate frame identity; physical landmarks; material/fiber/mass/time manifest | Mesh, fields, boundaries, or configuration change |
| 6. MPI equivalence | Prove distributed execution does not change the selected discrete problem | 1/2/4-rank short-history comparison; complete owned-row mass; identical loads, coordinate/connectivity hashes, mesh-bound field, and frame identity; accepted-state parity | MPI partition, assembly, solver, or state-commit change |
| 7. Full trajectory | Measure agreement, not just completion | Complete 1 s current-source archive; all-step diagnostics; `p0`/`p1` vectors; formulation state; common-time comparison with the local FEniCS reference and ten-team data; `det(F)` distribution, low-`J` reference-volume fractions, global wall volume, and declared sampling | Production source or configuration change |

Passing one layer does not imply the next. In particular, nonlinear
convergence is not trajectory validation, equal DOF counts do not make Q1
Hex8 and P2 tetrahedral spaces equivalent, and a stable volumetric
formulation is not automatically the paper's volumetric law.

## Package commands

Run the inexpensive package gates after every code change:

```bash
python -m pytest -q
```

Run compilation and reduced-application gates after material, element, state,
or solver changes:

```bash
python -m pytest -q -m slow
```

Run MPI gates in the matched PETSc environment after distributed assembly,
preconditioner, partition, or commit changes:

```bash
python -m pytest -q -m mpi
```

The closed local-pressure parser and reporting contract is covered by:

```bash
python -m pytest -q \
  tests/test_cardiac_mpi_closed_configuration.py \
  tests/test_cardiac_reporting.py
```

## Configuration registry

Every result-producing command must select one row explicitly. A new
combination receives a new identity; it must not borrow another row's label.

| Identity | Geometry and fields | Volumetric formulation | Mass and time | Evidence role |
|---|---|---|---|---|
| `cardiac-owned-distributed-q1p0-step0` | Historical polar ring/open apex; CG1 fibers | Condensed Q1/P0, `p=kappa*log(J)` | Lumped; backward Euler | Historical MPI compatibility only |
| `cardiac-owned-distributed-closed-std-kappa-step0` | Closed five-block; Laplace field; GP-direct fibers | Pointwise paper penalty, `p=kappa*(J**2-1)/2` | Consistent; backward Euler | Closed paper-law spatial/time diagnostic |
| `cardiac-owned-distributed-closed-local-pressure-step0` | Closed five-block; Laplace field; GP-direct fibers | Condensed Q1/P0, `p=kappa*log(J)` | Consistent; backward Euler | Selected robust CoupFE Case A formulation |
| `cardiac-owned-distributed-closed-local-pressure-mean-logj-paper-j2-step0` | Closed five-block; Laplace field; GP-direct fibers | Condensed Q1/P0, `m=<log(J)>`, `p=kappa*(exp(2*m)-1)/2` | Consistent; backward Euler | Explicit paper-volume-law local-pressure experiment |
| `cardiac-owned-distributed-closed-std-kappa-generalized-alpha-step0` | Closed five-block; Laplace field; GP-direct fibers | Pointwise paper penalty, `p=kappa*(J**2-1)/2` | Consistent; source-matched generalized alpha | Case A temporal/spatial diagnostic; historical pre-correction full Step 0B trajectory, but no current-source full trajectory retained |
| `cardiac-owned-distributed-closed-local-pressure-generalized-alpha-step0` | Closed five-block; mesh-bound Laplace field; declared GP-direct frame reconstruction | Condensed Q1/P0, `p=kappa*log(J)` at the alpha-f stage | Consistent; source-matched generalized alpha | Historical fine 8-rank Case A 1 s trajectory predating the straight-wall/frame corrections; current-source Step 0B 8-rank 0.32 s gates and completed clean `2x20x17`/`4x20x17`, `tip_refine=6.0` full cycles |
| `cardiac-owned-distributed-closed-local-pressure-mean-logj-paper-j2-generalized-alpha-step0` | Closed five-block; Laplace field; GP-direct fibers | Condensed Q1/P0 paper law at `exp(<log(J)>)` and the alpha-f stage | Consistent; source-matched generalized alpha | Retained fine 8-rank 0.70 s volume-law null trial |
| `hex8_fbar` archives | Configuration recorded per archive | Historical F-bar | Recorded per archive | Legacy/application comparison; never relabel as local pressure |

The selected historical fine Case A local-pressure mesh is
`(n_t,n_core,n_radial)=(4,36,32)`: 23,616 Hex8 cells, 29,885 nodes, and
89,655 displacement DOFs. It keeps the closed geometry, Case A zero cavity
pressure, epicardial normal Robin support, full-vector base Robin support,
Laplace transmural field, GP-direct fibers, consistent mass, `eta=100 Pa s`,
`dt=0.001 s`, and the canonical 1 s load horizon. Historical/current BE
archives must state `integrator=be` and `viscous_rate=backward_difference`.
Its retained 1 s archive predates the toolkit-matched straight-wall mapping
and physical-coordinate frame reconstruction and remains evidence only for
that recorded source. The separate source-matched generalized-alpha Case A
path passes focused kinematic, load-stage, material-point, compiled-element,
condensed-pressure, and analytic PETSc gates, plus a loaded coarse ventricular
1/2/4-rank gate to roundoff. The selected log-law formulation also has a
retained fine 8-rank 1 s trajectory; the paper-law variant has a controlled
fine 8-rank trajectory through 0.70 s. The coarse rank gate does not establish
spatial convergence or equate its direct solver profile with the fine
iterative profile, so every archive keeps its exact implementation and solver
identities and no BE archive is relabeled.

The controlled paper-law trial used the same 4x36x32 mesh, geometry, boundary
conditions, fields, material parameters, consistent mass, generalized-alpha
parameters, time step, and 8-rank solver profile as the log-law baseline. It
completed 700/700 steps with zero domain rejections. Against official Simula
on the common 0--0.70 s grid, the p0 relative L2 changed from 8.273% to 8.348%
and p1 from 11.480% to 11.578%; maximum vector gaps were unchanged to less
than 0.001 mm. The candidate archive SHA-256 is
`a3c19a6bfa040a6780466eae20a86ee653982015dc7d1fd049baf9ca1a2d21a9`.
This is a null result for the scalar volume-law hypothesis, not a 1 s
benchmark claim.

The Case A campaign was stopped by user direction on 2026-08-04. The retained
fine log-law generalized-alpha result completed 1,000/1,000 steps on 8 ranks
and has full-history relative L2 of 8.566% at `p0` and 12.263% at `p1` against
official Simula. Its archive SHA-256 is
`ba9b31ec533398be1f39fc9a898e72f77d9587c90f9b7d9e00ce91e4d2ae6a6c`.
The checked-in compact report
`case_a_local_pressure_4x36x32_dt0p001.report.json` binds that archive and
retains the ten-team mean comparison; benchmark-paper RED is
0.3337402/0.5024615 at `p0`/`p1`. Because the archive predates explicit
benchmark identity fields, the report labels Step 0A `legacy-inferred` rather
than recorded.
See the [Case A stopping record](CASE_A_STATUS.md) for the complete
configuration, metrics, and claim boundary.

The retained pre-straight-wall, pre-physical-frame Case A local-pressure
generalized-alpha mesh evidence does not show monotonic improvement toward
Simula on the shared 0--0.20 s prefix. From
`2x20x17` to `4x36x32`, `p0` Simula relative L2 changes from 4.386% to 7.103%
and `p1` from 6.956% to 10.033%. On the native 0.001 s grid, the two mesh
trajectories differ by 5.331% and 2.821%, respectively, normalized by the finer
histories. No coarser local-pressure archive reaches the 0.65--0.67 s maximum-
gap interval. This supports neither a convergence rate nor a claim that
refinement would remove the gap. No further Case A refinement is planned under
the closeout decision.

## 2026-08-07 Step 0B current Layer 7 status

The closed `2x20x17` and `4x20x17`, `tip_refine=6.0` Q1/P0 runs completed all
1,000/1,000 increments on eight MPI ranks with consistent mass, source-matched
generalized-alpha, and `dt=0.001 s`. Their checked-in dual-run
[`comparison report`](../examples/cardiac_benchmark/results/step0b_tip6p0_full_cycle_comparison.report.json)
and [`figure`](figures/step0b_tip_refine_full_cycle.svg) bind the external
archives and exact local FEniCS/ten-team inputs. Relative to two layers, the
four-layer trajectory has pre-snap/full-cycle/late pairwise RMSE
0.009831/0.407604/0.051556 mm at `p0` and
0.002622/0.361832/0.062702 mm at `p1`; maximum separation is
1.633083/1.405866 mm in the snap window and cycle-end separation is
0.034598/0.051776 mm. FEniCS full-shared-history RMSE improves from
1.092675/1.168502 to 1.016944/1.072951 mm, while relaxation RMSE worsens
16.2%/14.9%. This satisfies the full-history comparison layer for a controlled
pair, but not spatial convergence, rank equivalence, a physical-branch claim,
or a validation threshold.

The first four-layer full-cycle candidate is a non-retained provenance
diagnostic. Its clean isolated replay completed, matches it to roundoff, and
uses Core `454f73c` plus runtime-source SHA-256
`f8d9469a101709d11460a0803b6c031001192ccce61c472d032355b01070da05`.
Retained provenance is **application `2458e7c`, NPZ
`1e333b29b05f01dedce9272b32b82ce6ccfda56036c1c8f57eb395b9b4494800`,
stdout `0810a9d0c944345464f48917c707a7d3df7187b543c74efd4b0384542b95aed7`,
elapsed `1778.4 s (29.6 min)`**.

## Historical Step 0B prefix evidence and former Layer 7 boundary

> Superseded status: this section preserves the evidence and stop decision as
> recorded on 2026-08-05. The separate `2x20x17`/`4x20x17`,
> `tip_refine=6.0` full-cycle pair later completed as described above.

Step 0B generalized-alpha does have completed trajectories; the earlier
statement that none had been run is obsolete. The available evidence consists
of two explicitly separated 0.32 s generations, not a current-source full
cycle.

First, four clean `a2006b7`/`e2f42ed`, eight-rank prefixes completed all
320/320 increments with the corrected straight-wall geometry, condensed Q1/P0
log-`J` pressure, consistent mass, source-matched generalized alpha, and paper
parameters. The controlled meshes were coarse `2x20x17`, wall-only
`4x20x17`, surface-only `2x36x32`, and fine `4x36x32`. At 0.32 s, adding only
wall layers on the coarse surface changed the `p0`/`p1` endpoint vectors by
0.669/0.635 mm, whereas surface refinement with two wall layers changed them
by 9.994/8.772 mm. Coarse and wall-only form one response pair; surface-only
and fine form the other. This isolates the surface discretization as the
dominant variable in that split, but it is not a convergence rate or proof
that one mesh approaches a continuum solution.

A hash-bound reconstruction of the same reference Robin operator gives the
corresponding rigid-mode split. For
`q_x(X)=e_x cross X=[0,-X_z,X_y]`, coarse and wall-only both have total
`q_x^T K q_x=1.44455 N m/rad`, including `1.28510 N m/rad` from the
epicardial normal-only term. Surface-only and fine both give
`0.555828 N m/rad`, including `0.396153 N m/rad` from the epicardium; the
base full-vector term stays near `0.1595 N m/rad`. The actual 0.32 s nonlinear
states show the same grouping: coarse/wall-only epicardial spring energies are
0.308407/0.305358 J and long-axis moments are -0.771763/-0.774477 N m, while
surface-only/fine values are 0.116757/0.117304 J and
-0.248906/-0.250147 N m. This is strong mechanism evidence, not proof that
Robin faceting is the sole source of the CoupFE/reference difference. The
prescribed full-vector base and normal-only epicardial Robin law and its
coefficients are unchanged.

Those four trajectories used the old stored-parametric Gauss-point frame. The
subsequent physical-frame audit uses straight-wall meshes and their
mesh-bound Laplace fields, whose SHA-256 values are
`de24749a85b458c039a16a1e4b24422cf35d54cf853a40daa763e3137cb930a4`
for `2x20x17` and
`7fbdbcf9a6b6c5135ef87fe998ec23a3dbf44957ccd3184584e8f6d60768c2b6`
for `4x36x32`. The current identity is
`toolkit-physical-coordinate-u-v-v1`: it reconstructs the complete fiber,
sheet, and normal frame from each physical quadrature point and its Laplace
value. Compared with the old rule, fiber-axis p99/max changes are
0.386/2.823 degrees on the coarse mesh and 0.322/2.692 degrees on the fine
mesh. The near reversal of the old directed sheet/normal vectors is
material-invariant for the present Holzapfel--Ogden law; the small local axis
rotations are the physical material change. This field audit establishes the
frame formula and identity, not a mechanical trajectory by itself.

Second, clean current-source `056c02d`/`e2f42ed` eight-rank coarse
`2x20x17` and wall-only `4x20x17` runs use that physical-frame identity and
completed 320/320 increments with zero domain rejections. Their archive
SHA-256 values are
`efc7e42a60218ab275df2250ae4383f081554e25c46beee1c25cff3248c47785`
and `dc79faa158f04007c2430592024c69030a0894bde02967caac6dc338194b706f`.
The wall-only endpoint differs from coarse by 0.670 mm at `p0` and 0.621 mm
at `p1`, preserving the through-wall pairing. Against the hash-gated local
FEniCS history over 0--0.32 s, vector RMSE is 2.5889/2.5304 mm for coarse and
2.3499/2.2157 mm for wall-only at `p0`/`p1`. Relative to the matching old-frame
prefixes, the maximum landmark-vector changes are at most 0.0344/0.0307 mm
for coarse and 0.0348/0.0357 mm for wall-only. Thus the frame correction is
required for source fidelity but does not explain the larger response gap on
these two meshes.

At the time of this historical decision, these runs were completed short gates
and actual response evidence; they were not checked-in full-cycle benchmark
results, a rank-equivalence proof, a spatial-convergence study, or validation.
The full 1 s current-source Step 0B simulation was then intentionally paused.
The later full-cycle pair supersedes that execution boundary without changing
the evidentiary limits of these prefix artifacts.

The checked-in compact record
[`step0b_case_b_clean_frame_0p32.report.json`](../examples/cardiac_benchmark/results/step0b_case_b_clean_frame_0p32.report.json)
hash-binds the two current archives, their logs, the reference inputs, and the
mesh/Robin diagnostics while preserving this prefix-only claim boundary.

## Established findings: do not repeat without a trigger

The following findings are already supported by retained artifacts and the
Case B postmortem. They should be cited when planning work, not automatically
rerun:

- The historical polar-ring mesh is not the closed benchmark domain; reducing
  its apex offset does not test the closed topology.
- The five-block closed Hex8 geometry and its endocardial, epicardial, and base
  partition have passed extended Jacobian and boundary audits at the declared
  refinement levels.
- Consistent versus lumped mass changed the named historical truncated-polar
  snap-window trajectory by at most about 0.03 mm and did not change onset.
  Consistent mass remains the reference-aligned choice.
- Halving the historical truncated-polar time step from 2 ms to 1 ms did not
  remove its trajectory gap. This is a configuration-specific null result,
  not proof of temporal convergence.
- The serial constant-average-acceleration Newmark experiment is not the local
  FEniCS reference's generalized-alpha method. The historical MPI identities
  remain backward Euler; the source-matched generalized-alpha identity is
  separate.
  Case A has a loaded coarse 1/2/4-rank gate and retained fine 8-rank full
  trajectory, both from sources predating the straight-wall and physical-frame
  corrections. Step 0B has four completed straight-wall/old-frame mesh-split
  prefixes, two completed clean current-source physical-frame prefixes, and a
  separate completed current-source `2x20x17`, `tip_refine=6.0` full-cycle
  trajectory. The full cycle is not a rank or spatial-convergence study.
- The completed local FEniCS reference uses quadratic P2 tetrahedral
  displacement, consistent mass, generalized alpha with
  `alpha_m=0.2`, `alpha_f=0.4`, `gamma=0.7`, `beta=0.36`, staged loads, and a
  velocity-consistent viscous rate.
- The pre-straight-wall, pre-physical-frame clean fine closed Case B
  `(4,36,32)` run completed 1,000/1,000 steps on four ranks with 89,655 Q1
  displacement DOFs. Its pre-snap FEniCS RMSE is
  only 0.1689/0.0551 mm, but it enters snap 14.293/14.506 ms early and its
  full-history relative L2 error is 56.54%/52.50%. The earlier dirty coarse
  result is not a controlled refinement partner, so the available evidence is
  nonmonotonic but does not isolate mesh as the cause.
- At the shared 0.482 s label, that pre-straight-wall, pre-physical-frame fine
  CoupFE result and the local post-hoc FEniCS reconstruction have minimum
  native-quadrature `det(F)` of 0.327287 and
  0.332702, with wall-volume ratios 0.991334 and 0.992210. The low values are
  localized near the epicardial apex; CoupFE has a broader low-`J` tail. Do
  not use one native extremum to diagnose Q1 failure or explain the landmark
  gap. The FEniCS reconstruction is not yet a packaged regression gate.
- The same fine Q1 Hex8 scale has a retained pre-straight-wall,
  pre-physical-frame Case A relative L2 error of only 8.566%/12.263%. This
  does not make Case A and B equivalent, but it rules out “Q1 Hex8 alone” as
  an adequate explanation for the much larger old Case B gap.
- CG1 versus GP-direct fibers produced a measurable field change but did not
  remove the named historical truncated-polar onset gap. GP-direct plus the
  checked Laplace field remains the closed CoupFE choice; it is still not
  identical to the local FEniCS reference's P2 structural-field
  representation.
- `p0` and `p1` are physical points evaluated by checked Hex8 inverse mapping,
  not assumed mesh nodes. Shared-face candidates agree to roundoff in the
  controlled closed meshes.
- Coordinate reconstruction and displacement accuracy are separate. The
  checked inverse map rules out a wrong landmark, but the returned value is the
  Q1 trilinear field. Shape-function evaluation adds no identified sampler
  error beyond roundoff, but the Q1 spatial approximation can differ from
  Simula's P2 tetrahedral field. The retained evidence does not quantify that
  contribution separately.
- F-bar, pointwise `std-kappa`, and condensed local pressure are different
  formulations. Local pressure has shown useful robustness, but its `log(J)`
  law must remain explicit in every comparison. The optional paper-law local
  pressure evaluates the scalar paper response at the reference-volume-
  weighted geometric mean `J=exp(<log(J)>)`; it is neither the pointwise
  `std-kappa` response nor an average of pointwise pressures.
- Trial-state admissibility instrumentation is formulation-specific. The Q1/P0
  path explicitly guards nonpositive `det(F)` trials; the distributed
  pointwise-`std-kappa` batch currently does not. A positive retained peak
  state and zero recorded domain rejections do not prove that all discarded
  line-search trials were positive.
- On the controlled fine Case A generalized-alpha trajectory through 0.70 s,
  replacing the log law by that paper-law scalar response slightly worsened
  both p0 and p1 Simula errors and left the 0.65--0.67 s gap effectively
  unchanged. Do not repeat this trial as a generic explanation for the
  remaining Case A discrepancy unless another interacting input changes.
- Existing local FEniCS point-stress postprocessing is quarantined;
  displacement histories are usable reference data, while stress must be
  reconstructed from accepted `u/v/a` with a corrected definition before comparison.
- The pre-straight-wall, pre-physical-frame clean fine Case B 0.32 s prefix
  and 1 s archive have different hashes and claim boundaries. Only that old
  full archive and its full-comparison report support full-cycle metrics for
  their recorded source, even though their shared prefix agrees to roundoff.
  That old archive does not substitute for the separate current-source
  tip-refined Layer 7 record described above.

Rerun an established test only when its named configuration no longer matches
the code under review, an input hash changes, the previous evidence was
incomplete, or the new decision depends on an interaction that the old
one-variable test did not cover.

## Run and artifact organization

Use one immutable directory per configuration and attempt:

```text
<case>-<formulation>-<mesh>-<date>/
  source/       clean reviewed worktree and revision
  input/        mesh-field files, coordinate/connectivity hashes, and metadata
  build/        rank-isolated generated kernels
  smoke/        short acceptance archive and log
  run/          full archive and log
  comparison/   local-FEniCS/team reports and figures
  manifest.json exact command, environment, hashes, status, and claim boundary
```

Statuses are `planned`, `running`, `completed`, `failed-before-result`,
`interrupted`, `rejected-comparison`, and `quarantined`. A partial log is not a
completed trajectory. A result produced with the wrong formulation remains a
diagnostic under its original identity; it is never renamed into the selected
configuration.

## Production acceptance rule

The full run may be considered only after the exact production source is clean
and a short command with the same geometry, fields, formulation, mass, time
step, MPI implementation, and solver profile has:

- retained canonical node-coordinate and element-connectivity hashes, a
  matching mesh-bound Laplace-field sidecar, and the declared structural-frame
  reconstruction identity;
- passed the pre-solve geometry, pressure, and boundary audits;
- recorded finite Robin rigid-translation diagonal forms and the long-axis
  `q^T K q` form split into base full-vector, epicardial normal-only, and total
  contributions, cross-checked against the assembled spring matrix;
- accepted a nonzero load for the selected benchmark mode without committing
  a rejected state;
- retained the formulation-specific implementation identity;
- retained complete consistent-mass ownership metadata; and
- produced finite displacement, positive Gauss-point `det(F)`, and positive
  SNES/KSP acceptance diagnostics.

The rigid-mode values diagnose the selected discrete boundary geometry. They
do not impose a mesh-dependent validation target, and they do not modify or
tune the prescribed Robin law. Passing these prerequisites does not compel a
production launch. The current-source `2x20x17`, `tip_refine=6.0` full cycle
is complete, but the rank-equivalence and controlled mesh-ladder parts of the
broader Layer 7 evidence remain incomplete.

The full archive is compared at common physical times using all three
components at both physical points. Figures accompany numeric RMSE, relative
L2, RED, and componentwise errors; visual similarity alone is not a gate.
