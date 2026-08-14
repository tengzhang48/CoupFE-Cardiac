# Cardiac application design lessons

These notes capture choices that affect how the public application should be
extended and how its results should be interpreted.

The multi-round Case B investigation has its own detailed
[debugging postmortem](CASE_B_DEBUGGING_POSTMORTEM.md). That record preserves
the experiment chronology, explains why geometry and consistent mass were
audited late, checks the remaining differences against the local FEniCS run,
and defines the setup-first benchmark protocol. This page retains the broader
application-design rules.

## Keep cardiac policy in the application

The repository began as a development fork that included a copy of CoupFE. The
public form is an application over an exact Core dependency:

- Core owns generic code generation, compiled elements, operator composition,
  assembly, linear-solver routing, and distributed infrastructure.
- CoupFE-Cardiac owns ventricular mesh and facet semantics, fiber attachment,
  benchmark parameters, cardiac boundary operators, volumetric experiments,
  output-point sampling, time-history choices, diagnostics, and external-data
  comparison.
- Serial cardiac operators consume neutral NumPy arrays. Distributed examples
  use Core's `KernelMeshView` at the partitioning boundary.

This keeps domain assumptions visible without adding cardiac mesh adapters or
benchmark-specific policy to the general FE layer.

## Geometry choices need explicit names

A collapsed structured apex closes the topology but creates degenerate Hex8
cells and boundary facets. Moving the last ring away from the apex removes the
degeneracy but creates an open boundary. The historical
`topology=polar_ring, apex_offset=0.2` records use that second construction.
They are completed truncated-domain runs, not executions of the closed
Benchmark 1 domain.

This distinction must be checked before a nonlinear solve, not inferred after
one. At 2×36×48, the historical truncation removes 1.9335 mm of tip length and
adds a 2.672330 cm² traction-free annular boundary. Relative to the retained
closed reference mesh, wall volume is 0.6816% low and endocardial,
epicardial, and base areas are 0.4983%, 0.7179%, and 0.1926% low. More
importantly, the undeformed unit-pressure resultant is 4.5937% low (about
1.324 N at the retained peak pressure). Small global measure errors therefore
do not make a topology and boundary-condition change negligible.

The application now has a five-block square-to-disk Hex8 implementation that
puts one ordinary vertex at the apex, without collapsing cells or introducing
a fourth boundary. Its current mapping follows the toolkit meridional
construction: corresponding endocardial and epicardial points are joined by
straight Cartesian segments. For `n_t=2`, `n_core=20`, `n_radial=17`, and
`core_half_width=0.36`, it has 5,403 nodes and 3,520 elements. All 3,680
exterior faces are classified exactly once. Extended vertex/edge/face/center
and Gauss sampling gives minimum `det(J)=1.23646e-9 m^3`, maximum Jacobian
condition 7.87649, and minimum scaled Jacobian 0.258135. Relative to the
retained reference, wall volume and endocardial, epicardial, and base areas
differ by -0.1229%, +0.0695%, +0.0090%, and -0.0565%; the unit-pressure
resultant differs from the analytic base projection by 0.1028%. These are
pre-solve geometry and boundary qualifications; they do not substitute for the
trajectory comparison.

The first full coarse run with these mesh counts completed 1,000/1,000 steps
and reduced the local-FEniCS full-shared-history vector RMSE to
1.5803/1.7903 mm at
`p0`/`p1`, but it came from a dirty development tree. A later clean fine
`(4,36,32)` MPI-4 run also completed 1,000/1,000 steps; its
full-shared-history RMSE
was 7.5955/7.2363 mm and onset was about 14.3/14.5 ms early. Both runs predate
the straight-wall mapping and physical-coordinate frame reconstruction. They
also differ in source and execution provenance, so they are historical
development evidence rather than a controlled refinement pair. Successful
dynamics did not retroactively validate either curve, and clean execution does
not establish convergence.

Benchmark geometry and boundary-condition identity is a mandatory pre-solve
gate. Check units and extents, positive reference Jacobians beyond the Gauss
points, ownership and exactly-one labeling of every exterior face, measures by
marker, outward pressure orientation, undeformed pressure resultant and
moment, and Robin matrix properties before spending time on dynamics. A mesh
that fails this gate can remain a named demonstration or historical record,
but it must not be called the paper benchmark geometry.

An open endocardial surface also lacks a unique closed-cavity volume. Any
volume calculation needs a documented and tested cap construction; otherwise
report endocardial displacement and deformation Jacobians directly.

Closed-mesh construction and any future imported-mesh adapter belong in this
cardiac application rather than Core.

## 2026-08-05: boundary faceting can dominate a mesh split

Four clean `a2006b7`/`e2f42ed`, eight-rank Step 0B prefixes completed through
0.32 s with the corrected straight-wall geometry, condensed Q1/P0 log-`J`
local pressure, consistent mass, source-matched generalized alpha, and the
paper physical parameters. The controlled meshes were coarse `2x20x17`,
through-wall-only `4x20x17`, surface-only `2x36x32`, and fine `4x36x32`.
At 0.32 s, adding wall layers on the coarse surface changed the `p0`/`p1`
endpoint vectors by only 0.669/0.635 mm; refining the surface with two wall
layers changed them by 9.994/8.772 mm. Coarse and wall-only form one response
pair, while surface-only and fine form the other. This isolates surface
discretization as the dominant mesh-split variable; it is not a spatial-
convergence proof.

The normal-only epicardial Robin term provides a quantified discrete mechanism
consistent with the observed mesh pairing.
Rigid rotation about the long axis is tangent to the smooth axisymmetric
epicardium, so the continuum epicardial spring has zero stiffness in that mode.
Piecewise-planar facet normals are not exactly orthogonal to the same rotation
field and create a nonzero discrete stiffness. The audited total long-axis
stiffness is 1.44455 N m/rad on both coarse-surface meshes and 0.555828 N m/rad
on both fine-surface meshes. The epicardial contribution falls from 1.28510 to
0.396153 N m/rad, while the full-vector base contribution stays near
0.1595 N m/rad.

The actual nonlinear states at 0.32 s show the same pairings:

| Mesh | Epicardial/total spring energy (J) | Epicardial/total long-axis moment (N m) |
|---|---:|---:|
| coarse `2x20x17` | 0.308407 / 0.376541 | -0.771763 / -0.853824 |
| wall-only `4x20x17` | 0.305358 / 0.375250 | -0.774477 / -0.860245 |
| surface-only `2x36x32` | 0.116757 / 0.235046 | -0.248906 / -0.332284 |
| fine `4x36x32` | 0.117304 / 0.235185 | -0.250147 / -0.338022 |

These energies and moments describe materially non-rigid states; their pairing
supports the mechanism but does not prove that Robin faceting is the sole cause
of the trajectory gap.

Do not change the prescribed Robin law or tune its coefficient to hide this
effect. The benchmark and the audited nearby FEniCS source use the same
normal-only continuum form with discrete facet normals. Instead, state the
surface discretization and distinguish a finite-mesh comparison from a
continuum target. The paper requests nominal `h` approximately 5 mm; the mean
epicardial edge lengths in this audit are 3.996 mm for the coarse CoupFE
surface, 3.880 mm for the retained FEniCS mesh, and 2.173 mm for the fine
CoupFE surface. The retained FEniCS mesh has 1.85705 N m/rad total discrete
long-axis stiffness, including 1.69804 N m/rad from its epicardium and
0.159012 N m/rad from its base. Thus the participating curves are finite-mesh
evidence, not a continuum oracle, and a finer faceted surface need not move
toward their particular discrete response.

These four prefixes used the old stored-parametric Gauss-point frame. The
subsequent physical-coordinate correction reconstructs the complete toolkit
frame from each physical point and Laplace value. Its fiber change is small:
median about 0.02 degrees, p99 0.386/0.322 degrees, and maximum 2.823/2.692
degrees on the coarse/fine meshes, localized near the apex. The old sheet and
normal directed vectors were nearly reversed, but that near-global sign
reversal is material-invariant for the current Holzapfel--Ogden law: `ss` is
unchanged and the sign changes in `I8fs` and its fiber-sheet tensor factor
cancel. The small axis rotations are the actual material change. The subsequent
clean `056c02d` coarse and wall-only eight-rank gates both completed through
0.32 s; their maximum landmark trajectory change relative to the old frame is
only 0.03566 mm. The frame correction is therefore required for source
fidelity but is a controlled null result for the frame hypothesis on these
surfaces; it does not explain the larger Step 0B gap. The compact
[`Step 0B prefix report`](../examples/cardiac_benchmark/results/step0b_case_b_clean_frame_0p32.report.json)
retains that result and the full-1-s pause decision.

## Material directions need a declared reconstruction policy

Analytic nodal fiber and sheet frames are orthonormal, but component-wise CG1
interpolation does not preserve that property at Gauss points. Historical
polar-ring paths apply Gram--Schmidt after interpolation and record
`fiber_sampling=cg1_gram_schmidt`. Closed direct-Gauss-point paths now use one
shared physical-coordinate reconstruction of the complete fiber, sheet, and
normal frame and record its reconstruction identity. Each is an explicit
CoupFE-Cardiac discretization policy; neither should be silently relabeled as
the other.

The MIT-licensed cross-check and runtime adaptation use the same pinned toolkit
formula, including its apex branch. Verify the complete frame rather than only
`f0` dot products: a fiber-only check can hide an opposite but materially
equivalent sheet/normal sign convention, while still missing a genuine local
axis rotation.

## Sample the field that was solved

The earlier result path tetrahedralized all mesh nodes with global Delaunay and
sampled `p0`/`p1` with tetrahedral barycentric weights. That produced a usable
historical execution record, but it was not the Q1 Hex8 field interpolation.

The current sampler inverts the reference trilinear Hex8 map, validates the
candidate element and reconstruction, and uses the eight Hex8 shape weights.
It records the element, natural coordinates, weights, and reconstruction error.
Shared-boundary selection is deterministic, while outside, degenerate, and
orientation-reversing candidates fail closed.

Historical Delaunay-sampled numbers remain labeled as such. Changing the
sampler requires a rerun rather than a metadata-only relabeling.

## Separate point location from interpolation accuracy

Reconstructing a benchmark coordinate to roundoff proves that the requested
physical point was found. It does not prove that the finite-element solution
at that point equals the continuum solution. For the retained pre-straight-
wall fine Case A mesh, the reconstruction errors are `4.13e-20 m` at `p0` and
`1.04e-17 m` at `p1`, and shared-face candidates agree to roundoff. This rules
out a wrong point, nearest-node substitution, and an element-choice
discontinuity for that recorded mesh.

The sampled displacement is the trilinear Q1 Hex8 field assembled from the
eight element nodes. Shape-function evaluation returns that discrete field
value up to roundoff; it does not add a separately observed sampler error.
Simula's displacement is a quadratic P2 tetrahedral field, so the cross-code
gap still mixes global discretization error with approximation-space and
element differences. A checked sampler is necessary, but it cannot isolate
those terms. Describe the Q1 spatial approximation as a plausible contributor,
not the sampler as a bug or interpolation as the sole cause, unless a
controlled same-field or same-space study measures it.

## Interpret deformation Jacobians as fields, not verdicts

A minimum `det(F)` must be reported with its time, sampling rule, location,
affected measure, distribution, and global volume. At the pre-straight-wall,
pre-physical-frame clean fine Case B pressure peak, the CoupFE minimum is
0.327287, but only 48 of 188,928 native Hex8 Gauss samples are below 0.4. They
occupy 24 epicardial-apex elements and about 0.01003% of reference volume; the
global wall-volume ratio is 0.991334. The worst reference element is well
shaped, and its principal stretches show a thin through-wall compression
rather than inversion.

A local post-hoc reconstruction of the accepted FEniCS P2 state at its native
degree-6 tetrahedral quadrature gives a minimum of 0.332702 and wall-volume
ratio 0.992210 at the same 0.482 s label. Native extrema are not the same
physical point, and the reconstruction is not yet a packaged regression gate.
CoupFE also has a broader low-`J` tail: 0.01003% versus 0.000247% of reference
volume is below 0.4. This may be an indicator of different post-snap fields;
it is not isolated causality.

The old coarse minimum near 0.66 is consistent with underresolving or smearing
the localized layer, but those pre-straight-wall dirty-coarse and clean-fine
runs do not prove that interpretation. Do not diagnose inversion, locking,
element failure, or trajectory error from one extremum. Compare a common
physical time, declared sampling, low-`J` volume fractions, global volume,
and—when needed—a shared-point field reconstruction.

## Close a benchmark campaign without overstating convergence

The retained Case A log-law generalized-alpha mesh sequence reaches only
0.20 s on coarser levels. The controlled `2x20x17` to `4x36x32` refinement
changes the two landmark trajectories, but both Simula prefix errors worsen.
No coarser local-pressure archive reaches the 0.65--0.67 s maximum-gap
interval. The retained levels therefore show spatial sensitivity without
monotonic improvement toward Simula; they are not a convergence study or
evidence that one more refinement must improve the benchmark curve.

On 2026-08-04, the user stopped Case A after retaining the complete fine
generalized-alpha trajectory and the controlled paper-volume-law null trial.
The closeout preserves the observed 8.566%/12.263% full-history relative-L2
comparison at `p0`/`p1`, while explicitly leaving spatial and temporal
convergence, along with the allocation of the remaining gap, unresolved. A
closed investigation should state both the useful result and the unanswered
question. Reopen it only for a defined trigger such as a changed approximation
method, a same-space cross-code test, or a formal convergence requirement.

## Compare volumetric formulations without rewriting history

The 2026-06-27 Case B campaign used the F-bar formulation. Its table can be
preserved as project history, but the absent raw archives mean it cannot serve
as a source-verified regression oracle. The later source-identified open-tip
records are archived separately and remain exact-configuration history, not
current Benchmark 1 validation evidence.

The application-owned Q1/P0 path composes a standard Q1 cardiac material
kernel, with its bulk penalty disabled, and one algebraically eliminated
pressure per element. Its condensed tangent includes the pressure sensitivity
to displacement. Component tests cover affine and isochoric response, symmetry,
finite differences, and invalid reference or deformed cells.

Local pressure is a discretization pattern, not a fixed scalar volume law.
When comparing with the benchmark paper, keep the element constraint
`m=<log(J)>` but evaluate the selected constitutive response explicitly. The
paper-law variant uses `Jbar=exp(m)` and
`p=kappa*(Jbar**2-1)/2`, with exact derivative
`dp/dm=kappa*exp(2*m)`. It is not valid to change only the residual pressure:
the rank-one part of the condensed tangent must use that same derivative.
Likewise, this operation is not `p(<J>)` and is not the Gauss-point average of
the pointwise paper pressure. Every variant needs a separate CLI, archive, and
MPI implementation identity so old results remain reproducible.

The pre-straight-wall, pre-physical-frame fine 4x36x32, 8-rank
generalized-alpha Case A control tested this exact change through 0.70 s while
holding all other recorded inputs fixed. It was
stable (700/700 steps, zero domain rejections), but it did not improve the
benchmark curves: Simula relative L2 changed from 8.273% to 8.348% at p0 and
from 11.480% to 11.578% at p1. The largest candidate-to-baseline landmark
difference was only 0.077 mm, and the 0.65--0.67 s maximum-gap region was
unchanged to less than 0.001 mm. A plausible constitutive correction can be a
well-posed null result; retain it as evidence instead of tuning further or
silently adopting it as the default.

F-bar and Q1/P0 should be recorded and compared as distinct discretizations.
The new operator does not retroactively turn the legacy runs into local-pressure
results, and component agreement does not replace a ventricular resolution
study.

The pre-straight-wall, pre-physical-frame clean fine Case B run did not use
local pressure: it used the pointwise paper energy
`kappa*(J**2-1-2*log(J))/4` with `kappa=1.0e6 Pa`. The audited nearby FEniCS
source implements the same law and coefficient, although its exact producer
revision was not retained with the output. CoupFE's corresponding stress term
is the exact derivative of that energy; there is no doubled bulk modulus in
this run. A finite penalty strongly resists volume change but is not an exact
pointwise `J=1` constraint. Keep the constitutive law, its scalar parameter,
and the volumetric discretization as three separately auditable facts.

Do not tune the local-pressure bulk modulus silently to rescue a run. An
external controlled 2×12×16 diagnostic at `K=0.1`, `0.3`, and `1.0 MPa`
completed 300/300 increments in all three cases with zero invalid trials and
613, 614, and 645 total nonlinear iterations. Peak reference-volume-weighted
`V/V0` changed much more, from 0.889889 through 0.958346 to 0.987143. This does
not identify `K` as the cause of the historical failed three-layer Q1/P0
attempt and supports retaining the benchmark `K=1.0 MPa`; the changed-`K`
runs are external diagnostics, not retained benchmark results. Their summary
SHA-256 is
`13cec75a0134a8599197ed1f4c371fb32509c36817d4bb3d3b911663ad1108d3`.

The same rule applies to viscosity. The paper value is `eta=100 Pa s`. An
`eta=0` diagnostic failed closed before state commit at 0.026 s. A subsequent
`eta=50 Pa s` sensitivity run was stopped by user direction after its last
printed accepted state at 0.256 s and produced no completed archive. Neither
is a benchmark result, and no parameter tuning was accepted. Reproduction runs
retain `eta=100 Pa s` unless a separately named sensitivity study says
otherwise.

The checkpoint-`62ad760` Q1/P0 and F-bar 2×24×32, `dt=0.002 s` runs hold mesh,
time step, geometry, sampling, source, and nonlinear solver fixed. Their peak
vectors and RED values differ, and the direction of the RED difference is not
the same at `p0` and `p1`. This is a controlled estimate of the formulation
effect at one discretization, not evidence that either path is more accurate.
The Q1/P0 2×12×16 time-step pair
holds the spatial problem fixed and quantifies small common-grid history
differences after halving `dt`; two time steps still do not establish temporal
convergence. The archived F-bar 2×24×32 and 2×36×48 histories differ materially,
which exposes spatial sensitivity rather than a converged mesh sequence. A
broader formulation study should repeat matched comparisons across multiple
spatial and temporal resolutions.

The quantitative comparison has two simultaneous messages. On the 2×12×16
Q1/P0 time-step pair, maximum/RMS vector-history differences are
0.144513/0.024943 mm at `p0` and 0.115627/0.020255 mm at `p1`, which is bounded
coarse-grid temporal repeatability. Yet all retained Case B RED values lie
above the participant-team ranges retained in the fine report: the team maxima
are 0.2875376036 at `p0` and 0.3701682113 at `p1`, while coarse Q1/P0 is about
0.88/0.87, 2×24×32 Q1/P0 is about 1.07/1.16, and 2×36×48 Q1/P0 is
0.551/0.654. RED against the all-team mean combines amplitude and
mesh-dependent signed point-component differences. The 2×36×48 response has
benchmark-consistent component signs and
the retained F-bar vector is close to its legacy 2×36×48 counterpart, whereas
the 2×24×32 response has the opposite `p0.u_z` sign from the recorded benchmark
comparison row. Report that spatial and signed point-component sensitivity explicitly;
do not turn it into either a blanket inconsistency claim or evidence of
quantitative reproduction, convergence, bifurcation, a global twist direction,
or a solution branch.
The retained fine-report arrays sharpen the distinction: Q1/P0 2×36×48
component-history correlations with the all-team mean are
(0.939462, 0.837752, 0.969152) at `p0` and
(0.944736, 0.825150, 0.964841) at `p1`, with component RMSE
(5.432790, 0.358965, 2.161431) mm and
(5.613262, 0.244329, 2.245171) mm. High correlations
and matching signs can coexist with amplitude error and an above-team-range
aggregate RED.

Keep source identity in the comparison boundary. The additional Q1/P0
2×24×32, `dt=0.004 s` run was produced at application checkpoint `e07993b`,
while the retained `dt=0.002 s` run was produced at `62ad760`. Their grids are
nested and their histories can be compared, but both time step and application
checkpoint differ. Treat that result as cross-checkpoint sensitivity, not as a
controlled one-factor time study or temporal convergence evidence.

The archived 2×36×48 F-bar peak is within 1.11593 mm, or 2.89%, of the
rounded surviving legacy peak vector. It is not an exact-reproduction claim:
the legacy raw result, transcript, and complete configuration are absent, while the archived run
records a different apex and output-sampling policy. Quantify both the
agreement and the provenance gap.

## Retain successful output with stated boundaries

Retain a completed, source-identified run when its exact input choices,
accepted-step diagnostics, deformation checks, and output history are
available. State separately what that record shows and which broader questions
remain open, including mesh convergence and real-device comparison.

For the archived truncated-polar Case B records, every requested increment met the recorded
nonlinear acceptance rule and the retained peak-load Gauss-point deformation
Jacobians are positive. Those are meaningful facts about the named runs. They
do not establish spatial or temporal convergence or validate the open-apex
model against a real device or clinical measurement. The matched
checkpoint-`62ad760` 2×24×32, `dt=0.002 s` pair does isolate the volumetric
formulation at that one configuration. RED should be published as calculated
rather than hidden or converted into an unjustified pass/fail label. No
repository-defined pass/fail threshold is assigned.

This separation keeps code, tests, runnable examples, benchmark comparisons,
and retained output visible without asking any one artifact to prove more than
it does.

## Time integration must cover every rate term

Backward Euler remains the coherent first-order path because its inertia
predictor, Robin dashpot, and material strain difference follow one policy. The
optional Newmark operator shares its velocity map with the Robin dashpot, but
material viscosity still uses a backward strain difference; it remains a
kinematic probe rather than source-matched generalized alpha. The separately
identified closed generalized-alpha path stages inertia, material, Robin,
loads, and the viscous rate consistently for supported Step 0 Cases A/B and
Step 2 Case B. Names and archive metadata must keep all three paths distinct.

## Solver success must be independently checked

Core's compact Newton routine returns an iterate and iteration count; the count
is not a convergence flag. The application independently checks the final free
residual before every physical-operator commit.

The PETSc path preserves the recovered 2026-06-27 parameter values but uses a
new persistent application-owned context and stronger acceptance record. A
positive SNES reason, valid KSP status, finite iterate, and independently
reassembled residual within `max(atol, rtol*|R_initial|)` are all required.
Each accepted step records reasons, iterations, norms, threshold, history, and
timings. A failed step cannot reach the completed-result writer.

Element work policy is an execution choice, not a change in formulation or
tolerance. The drivers therefore keep the established joint R/K cache as the
default and expose `split` only as an explicit opt-in. They do not infer a mode
from problem size or backend. Every completed archive records the requested
mode and whether the generated native kernel supplied a residual-only entry,
so later timing or equivalence work can reconstruct the choice.

Invalid trial deformations now have a narrow, explicit recovery path. The
Q1/P0 operator raises `InvalidDeformationError` for a trial outside its
finite-strain domain and provides `max_step`, which Core Newton uses to bound a
proposed correction before residual backtracking. The PETSc residual callback
catches only that exception, returns an IEEE positive-infinity residual, and
lets SNES `bt` shorten the trial. It records the number of such evaluations and
the last detail for each step. Jacobian exceptions and unrelated residual
exceptions remain fatal; invalid initial, accepted, or final states still fail
before commit and before completed-result writing. This added recovery changed
no nonlinear tolerances or final residual acceptance rule.

Here, an invalid `det(F)` trial means the proposed nonlinear correction
overshot the admissible deformation domain. It is not a perturbed input mesh;
the documented runs use `perturb=0`. PETSc backtracking reduces the
within-increment correction, which is also distinct from changing the physical
time step `dt`.

The distinction is important: rejecting and shortening a bad trial is not the
same as accepting a distorted state, and a recovered nonlinear step is not a
mesh-quality or convergence result. If backtracking finds no valid accepted
state, reduce or revise the discretization or increment rather than weakening
the final residual rule.

Admissibility guards are formulation-specific. The Q1/P0 path has the explicit
trial-state determinant guard described above, but the distributed pointwise-
`std-kappa` batch used by the pre-straight-wall clean fine Case B run currently
does not. Its retained pressure-peak state has positive Gauss-point `det(F)`,
and the exact kinematic reconstruction agrees with the archive, but a zero
recorded domain-rejection count does not prove that every discarded line-search
trial was positive. This is an instrumentation and robustness gap to close; it
is not evidence that an invalid accepted state caused the Case B trajectory.

The `e07993b` Q1/P0 2×24×32, `dt=0.004 s` execution provides a retained example
of that distinction. PETSc rejected 46 trial residual evaluations at steps
132–133 before backtracking to valid trials. Every accepted/final state met the
unchanged domain and residual checks, and no rejected trial was committed.
This demonstrates that the recovery mechanism operated in one named run; it
does not validate the mesh, tune a tolerance, establish convergence or
accuracy, or support a bifurcation or clinical claim.

The Q1/P0 2×36×48, `dt=0.002 s` execution at the same checkpoint rejected 168
trial residual evaluations at steps 277 and 279, then completed 500/500 valid
accepted/final states under the unchanged rule. The largest final-residual-to-
threshold ratio was 0.973921. Retaining that near-bound acceptance ratio makes
the evidence more auditable; it does not relax the rule or establish accuracy.
Its F-bar counterpart was produced at `62ad760`, so their side-by-side
difference cannot be called a controlled one-factor formulation comparison or
a reproduction of the legacy F-bar observation.

The fail-closed path also matters. For a Q1/P0 3×36×48, `dt=0.002 s` attempt,
backtracking did not yield a converged increment; it ended with SNES reason
`-5` before state commit and produced no completed NPZ or report, so it is not
retained as result evidence. The surviving historical row specifies the
smaller physical time step `dt=0.001 s`. A detached clean-checkpoint-`e07993b`
Q1/P0 rerun at that value was stopped by user direction on 2026-08-01, after it
printed accepted steps 1 and 2, to prioritize an MPI companion. It produced no
completed NPZ or report; its external partial log has SHA-256
`d47ee0ce04a312db3808a2a3b373bb02a0c7f71ff5192f9dd88999e53b50554a`.
This is an interruption, not a numerical failure or result. See the
[numerical execution record](CASE_B_STATUS.md#numerical-execution-record) for
the explicit completed, failed, and interrupted boundaries.

The distributed Core API exposes a different set of diagnostics. Distributed
scripts therefore make claims only about their exact serial/rank comparisons
and reported solver fields; they do not borrow the serial application's
all-step SNES record. Keep the existing serial driver/example and
retained results unchanged as the reference. Add MPI as a companion path, and
require rank-equivalent outputs against that serial reference before treating a
distributed run as benchmark evidence.

## Residual callbacks should not construct an unused pressure tangent

The first MPI Q1/P0 adapter fused its material and condensed-pressure R/K
blocks in every callback. That preserved same-iterate reuse, but it also built
the full 24×24 pressure tangent during line-search probes and the independent
initial/final residual checks. The retained serial pressure operator already
had separate residual and tangent methods; the MPI adapter had regressed that
design.

The corrected adapter exposes explicit `split` and `joint` modes. Split caches
shared kinematics, evaluates pressure residual and tangent only on their
corresponding callbacks, and uses Core's residual-only compiled material entry
when present. Joint remains available for paired-callback comparison. Invalid
pressure kinematics are checked before invoking the stateful material kernel,
and the accepted final residual still refreshes trial state before one commit.
Parity, call-count, exact-equilibrium, and serial-versus-rank tests gate the
change. Performance must still be reported from a resource-matched retained
run rather than inferred from the design alone.

## A zero-residual increment must not require a linear solve

Case A contains a normal initial interval with zero active tension and zero
cavity pressure. On those increments the predictor can already satisfy the
nonlinear equations: SNES reports convergence without invoking KSP, leaving the
KSP reason at `0` and the linear-iteration count at zero.

The first distributed wrapper incorrectly required a positive KSP reason on
every accepted increment. It therefore rejected this normal Case A state even
though SNES and the independently assembled residual both reported
convergence. The active Case A serial-versus-MPI parity test found the defect
before release. The corrected acceptance rule is explicit:

- a negative KSP reason is a failure;
- KSP reason `0` is valid only when no linear iteration was requested; and
- KSP reason `0` after one or more linear iterations is a failure.

SNES convergence, finite diagnostics, and the independent global residual rule
remain mandatory in all three cases. This is a solver-state distinction, not a
special numerical tolerance and not an exception for Case A.

## A regression oracle is valid only for the law that produced it

The Case A slow gate initially appeared to reveal a 0.266 mm displacement
regression. It compared the current driver against a clean `62ad760` report,
so the first suspicion fell on later Core changes to fused element evaluation
and state commit. That diagnosis was wrong because the test crossed a
constitutive-law boundary before it compared any numbers.

A controlled source matrix resolved the ambiguity:

- application `62ad760` with Core `454f73c` reproduced every retained array
  exactly;
- application `62ad760` with Core `e2f42ed` also reproduced those arrays to
  roundoff, including nonlinear iterations and peak `det(F)`; and
- application `6839c13` with Core `e2f42ed` reproduced the observed drift.

Commit `6839c13` added the missing derivative of the smooth compression switch
to the fiber and sheet stresses. That change makes the stress the derivative
of the declared complete energy. The earlier output is therefore evidence for
an earlier constitutive expression, not an oracle that the corrected law
should be forced to match. Relaxing the displacement tolerance would have
hidden the semantic mismatch; reverting the derivative would have made the
implementation internally inconsistent again.

The release keeps both records. The `62ad760` report remains byte-identical and
historical. A separately named clean `6839c13` report is the numerical oracle
for the corrected law. Its predecessor NPZ predates the material-identity
field, so the report reconstructs that identity only under the immutable
reviewed-checkpoint contract and records
`method_metadata_origin=reviewed-predecessor-source-checkpoint`. Runs made from
the current checkpoint write
`material_model_id=holzapfel-ogden-smooth-switch-complete-energy-derivative-v1`
directly into the NPZ. Reporting code validates either qualified origin before
comparing curves. This makes a physical-law change visible even when peak
values happen to remain close.

The same investigation exposed a second reproducibility trap: the current
serial driver defaults to consistent mass, whereas the historical Case A run
used row-summed mass. Commands intended to reproduce an old record now state
`--mass lumped` explicitly. Defaults describe the current interface; retained
evidence commands must pin every result-defining choice.

## Recover the established formulation before launching a benchmark campaign

On 2026-08-03, a refinement campaign on the then-current topologically closed,
pre-straight-wall Case A geometry was launched with the pointwise `std-kappa`
option before the retained CoupFE configuration history and earlier
investigation were checked. The choice was physically motivated by the paper,
but it was not the established application baseline: the old driver default
and retained Case A demonstration used F-bar, while later retained Case B runs
used the condensed Q1/P0 local-pressure path. The resulting `std-kappa`
calculations remain useful diagnostics, but they must not be presented as if
they continued the established local-pressure campaign.

F-bar, pointwise `std-kappa`, and condensed Q1/P0 local pressure are three
different finite-element formulations. Q1/P0 local pressure names the
condensation pattern, not one scalar volume law: the `local-pressure` CLI
variant uses `p = kappa*log(J)`, while the separately identified paper-law
variant uses the mean-log-volume response documented above. The pointwise
benchmark penalty gives `p = kappa*(J**2 - 1)/2`. Greater nonlinear robustness is a
valid reason to select local pressure for a CoupFE study, but it is not evidence
that the two volumetric laws are identical. Time integration is independent:
historical distributed runs remain backward Euler, while a source-matched
generalized-alpha path requires its own implementation identity and evidence.

Before a production benchmark launch:

- inspect the newest clean source, retained archive metadata, historical run
  commands, and existing collaborator notes;
- write a configuration manifest that pins geometry, boundary conditions,
  formulation, mass, fibers, viscosity, time integration, time step, sampling,
  ranks, and solver profile;
- assign a formulation-specific implementation identity and fail closed when
  the requested combination has no validated identity;
- make a short acceptance run from the exact production command before
  submitting the full horizon; and
- preserve an incorrectly configured completed run as a clearly labeled
  diagnostic rather than renaming it or rewriting its provenance.

The closed MPI driver therefore gives closed local pressure its own explicit
implementation identity instead of reusing either the historical open-mesh
Q1/P0 label or the closed `std-kappa` label.  Parser and reporting gates check
that identity against the archived formulation and consistent-mass contract.

## Match every generalized-alpha stage, not only Newmark kinematics

A generalized-alpha comparison cannot be created by changing only the inertia
predictor. The audited Simula source at peeled commit
`325d17d850c2e2032abb85a4191a5795d3008ab7` evaluates acceleration at the
`1-alpha_m` stage and displacement, velocity, material, Robin terms, and loads
at the `1-alpha_f` stage. Its viscous Green--Lagrange rate is reconstructed
from stage velocity, not from `(E[n+1]-E[n])/dt`.

The application-owned closed MPI implementation therefore stores accepted
Gauss-point gradients of displacement, velocity, and acceleration in a
dedicated generated material, stages the condensed Q1/P0 pressure outside the
kernel with the endpoint chain rule, stages Robin spring and dashpot terms,
applies consistent mass to staged acceleration, and integrates the prescribed
load from `t=0` over the fixed horizon while sampling the shifted evaluation
grid. The old BE material and solver branches remain separate. A parser/
implementation mismatch or a BE material batch carrying a generalized-alpha
label fails before PETSc setup. The code contract covers closed Step 0 Cases
A/B and Step 2 Case B; that support is not itself a retained-result claim.

Focused kinematic, load, material-point, generated-element tangent/commit,
local-pressure, and analytic PETSc gates remain necessary but are not a
ventricular validation by themselves. In the Case A campaign, the subsequent
loaded coarse 1/2/4-rank gate passed to roundoff and the fine 8-rank 1 s
trajectory completed. Those records predate the straight-wall and physical-
frame corrections. That ordered evidence validates the named historical
execution layers, not spatial convergence. Repeat only the affected layer when
a future implementation, solver profile, or configuration changes.

## Retained comparisons need a complete provenance chain

Reference curves remain in the separately downloaded, CC BY 4.0 Zenodo
archive. A reviewable comparison should retain the exact application and Core
revisions, command, environment, formulation, mesh, time discretization,
nonlinear configuration, all accepted-step diagnostics, result hash, external
archive and curve-file identities, and full-precision comparison output.

Binary NPZ archives and trusted external pickles can remain outside the public
repository while their hashes and derived text/JSON records are retained. The
archived reports include the CC BY 4.0 derived all-team mean curves and RED
values with citation while leaving the raw pickle inputs external. A
legacy table whose raw artifacts are absent should be labeled as historical
observation, not silently promoted to current-source validation.

Prefixes and full trajectories are distinct evidence objects even when their
shared arrays agree to roundoff. For the pre-straight-wall, pre-physical-frame
clean fine Case B run, the 0.32 s prefix archive has SHA-256
`a7de58b1acf507fbdb059f87e14db040e514711736381a06e649a02b529d6a0b`, while
the completed 1 s archive has SHA-256
`63d41a1f69dceaa8c1fe7f3c7d46a6de4e40c270a35977214de741e06fc580a3`.
Only the latter supports the full comparison report with SHA-256
`602a1e904973d9700b6cd99f1b922fa81651815c7722c157c9a571bd9d7640f5`.
Never let a shortened diagnostic inherit a full-cycle claim because its
filename, revision, or initial history resembles the completed archive.

## Audit a reference postprocessor before treating its output as an oracle

The locally retained FEniCS `result.h5` stores accepted displacement,
velocity, and acceleration fields. Its supplied point-stress postprocessor,
however, reloads only displacement for each timestamp into a newly constructed
problem whose old displacement, velocity, and acceleration fields remain zero.
For the recorded generalized-alpha parameters, this reconstructs velocity as
approximately `1944.44*u` per second instead of using the stored accepted
velocity. It then projects the result to DG1 and samples that projection. The
implemented von Mises expression also adds unsquared shear components to
squared normal-stress differences, so its terms are not dimensionally
consistent.

The existing FEniCS point-stress arrays and XDMF stress field are therefore
quarantined as quantitative comparison oracles; they are not evidence that
the underlying solve is wrong. A future stress comparison must reload the
accepted `u/v/a` state, evaluate a corrected Cauchy-stress definition at a
declared element-interior or quadrature point, and record element identity,
coordinates, projection status, and the physical distance to the matched
CoupFE Gauss point.

## Signed response needs more than one component

Fiber handedness, sheet/normal convention, apex treatment, mesh and time
resolution, formulation, and nonlinear path can all affect Case B displacement
and twist. A sign change between discrete runs is a reason to investigate those
dependencies; it is not by itself evidence of bifurcation, bistability, or a
unique physical twist direction.

Report complete `p0`/`p1` histories, deformation Jacobians, the exact fiber
convention, and a defined twist measure together with any peak-component table.

## Treat graded-node studies as sensitivity until locality and convergence are shown

`--tip-refine` keeps element counts, connectivity, labels, and the base rim
fixed, but it moves almost every non-base meridional node: cells are clustered
toward the apex and coarsened toward the base. It is global r-adaptation, not
local h-refinement. On `2x20x17`, increasing its strength from
1.0/2.5/4.0/6.0 reduced the 0.32 s p0 endpoint gap against local FEniCS
6.15 -> 3.61 -> 2.59 -> 2.00 mm. That useful endpoint trend does not locate
the error at the apex or establish a convergence order.

History metrics are the necessary countercheck. Over the same prefix, p0
vector RMSE was 2.593/1.535/1.421/1.424 mm, while maximum history error for
the nonuniform runs grew 4.19 -> 4.93 -> 5.55 mm. At `4x32x48`, changing
the grading produced only a 0.23 mm p0 endpoint difference but a 2.39 mm
maximum transient difference. Endpoint improvement, history convergence, and
spatial error ownership are separate claims; report and test them separately.

The later `4x20x17 tip=6.0` full cycle makes the complementary point. Relative
to two layers, its pairwise RMSE is only 0.0098/0.0026 mm before snap, grows to
a 1.633/1.406 mm maximum separation during snap, and falls to
0.0346/0.0518 mm at cycle end (`p0`/`p1`). Against FEniCS,
full-shared-history RMSE and the maximum snap gap improve, but relaxation RMSE
worsens 16.2%/14.9%.
Endpoint agreement, transient agreement, and phase-specific agreement can rank
the same pair differently; report all three before judging a mesh change.

## A full-cycle pair can distinguish timing sensitivity from a branch claim

For the two-/four-layer pair, downward `u_z=-5 mm` snap crossings move about
1.2 ms toward FEniCS when two wall layers are added, while upward crossings
remain late and change by less than 1.2 ms. The histories nearly reunite late
in the cycle. This is the pattern of a discretization-dependent transient
timing/amplitude shift in the observed data; it is not evidence of two
physical solution branches. A branch claim needs continuation or stability
evidence, not curve clustering around a sharp transient.

The result also reinforces a provenance rule. A numerically complete run from
a dirty application tree can be a useful diagnostic, but it is not a retained
public result. Reproduce it from a
clean isolated source identity, verify prefix/state continuity to roundoff,
and bind its application, Core, runtime-source, NPZ, stdout, and elapsed-time
identities before promotion.

## Snap-window response clustering is not proof of solution branches

The coarse graded trajectory, the local FEniCS record, and published team
curves cluster more closely at the 0.32 s endpoint than the fine facet and
smooth-normal CoupFE trajectories. Their separation grows quickly in the snap
window. Snap-through can amplify small changes in discretization and boundary
representation, so descriptive "near" and "far" groups are useful shorthand
but are not evidence of distinct mathematical branches.

Branch claims require continuation, stability, or equivalent evidence that
these runs do not contain. Before reporting mesh convergence, compare complete
histories, snap timing, transient extrema, and declared phase metrics under one
fixed operator. Agreement of two endpoints is neither a proof of convergence
nor proof that another trajectory is on a different branch.

## Boundary operator representations are not corrections

The largest single lever in the Case B study is the epicardial normal
representation: facet vs smooth ellipsoid-gradient normals changed the
endpoint by 11.6 mm on one mesh. Faceting explains most of the
coarse-vs-fine *split* (the Robin-normal gate's pre-registered ratios),
but the smooth-normal mode moves the trajectory *away* from the ten-team
consensus. A representation diagnostic answers "is the trajectory
sensitive to this operator choice"; it does not automatically make the
alternative choice more faithful, especially when the reference
implementation itself uses finite facet normals.

## Round-keyed deduplication is a seam bug waiting for a fine mesh

Deduplicating shared block-seam points by rounding floating-point
coordinates is fragile: adjacent blocks compute shared corners through
different expressions whose last-bit differences can split a rounded key
at some resolutions. At `64x72` this duplicated six seam points and
created twelve overlapping full-size base quads (+3.62% base area) - a
silent geometry corruption that only the fail-closed pre-solve measures
audit caught before any solve.

The fix has two parts that must go together: tolerance-based matching
that preserves the first writer's bytes (so previously-working
resolutions reproduce bit-for-bit and archive/hash provenance survives),
and a hard structural invariant (the disk rim edge count) that turns any
future seam break into a build-time error instead of a subtle area error.
Whenever a dedup, renumbering, or sealing step merges points, verify
bit-identity against retained archives at the resolutions those archives
were produced with.

## The run environment and solver profile are part of the experiment

Three avoidable launch failures in this campaign were environment, not
physics: `petsc4py` missing from the default interpreter (the reference
runs pin petsc4py 3.18.4 in a specific env); a launch command missing one
flag (`--tip-refine`), caught by the fail-closed field/mesh identity
gate; and the driver's default linear-solver profile
(`direct-superlu-dist`) being different from the reference run's
`fgmres-gamg-rigid-rebuild` - a ~50x wall-clock difference that is
invisible in result files unless read from `solver_configuration_json`.
Before launching a control, copy the exact interpreter, profile, and
flags from the reference run's recorded configuration, and let
fail-closed gates (mesh identity, field sidecars, pre-solve audits)
reject mismatches rather than patching around them.

## Archives exist only at completion; budget runs for the machine you have

The fail-closed result writer saves an NPZ only after every requested
step completes. Two multi-hour MPI runs were lost to timeouts under
external machine load (load average ~400 on 64 cores, ~6x
oversubscribed). For expensive runs: measure the current load, disable or
generously size timeouts, launch one heavy job at a time, and treat a
killed run as a total loss (there is no partial archive). Machine
contention is experiment metadata; record it when runtimes are quoted.

## Define the level of reproduction before declaring success

"We reproduced the benchmark" is too broad to be auditable. Use separate
levels and report which ones the evidence supports:

1. **Problem identity:** geometry, boundary partition, loads, parameters,
   material points, and time history are identified.
2. **Numerical execution:** the declared discrete method completes under its
   solver, admissibility, MPI, and provenance contract.
3. **Qualitative response:** the expected contraction, snap-through, recovery,
   signs, and event order are present.
4. **Quantitative trajectory agreement:** complete vector histories agree
   under declared metrics and phase windows.
5. **Convergence or validation:** controlled refinements or independent
   physical evidence support a stronger conclusion.

Passing an earlier level does not imply a later one. For the retained evidence
as of 2026-08-07, Step 0 Case A provides historical approximate trajectory
agreement near the project's roughly 10% working level, but its pre-correction
geometry/frame means it is not a current-setup reproduction. Current Step 0
Case B passes physical-case identity, numerical execution, and qualitative
response, but only partially passes quantitative trajectory agreement. Its
Q1/P0 mean-`log(J)` volume response is a documented application variant rather
than the paper's pointwise scalar law. That distinction limits
method-equivalence claims; it does not make the volume law the diagnosed cause
of the gap. The separate paper-law change was a null only through 0.70 s on one
named fine, pre-correction Case A configuration and must not be generalized to
Step 0B. Neither case has a current formal convergence result, and neither is
experimental or clinical validation. The answer-first numbers and exact claim
wording are retained in
[`BENCHMARK_REPRODUCTION_STATUS.md`](BENCHMARK_REPRODUCTION_STATUS.md).

The Case B full-cycle pair adds four practical rules:

- compare the whole three-component history at both physical landmarks, not a
  single peak value or the largest component alone;
- report absolute error for small components because percentage error can
  become visually and numerically misleading near zero;
- separate pre-snap, snap, peak, relaxation, and late-cycle windows because a
  mesh change can improve one phase and worsen another; and
- do not call a finer mesh converged merely because its aggregate RMSE is
  lower: the current four-layer run improves full-shared-history FEniCS RMSE
  while worsening relaxation RMSE.

The benchmark's ten submitted teams do not define one exact discrete curve,
and the paper supplies no pass threshold. A named FEniCS/Simula comparison and
the official ten-team envelope answer different questions; retain both and
state the metric instead of converting either into an invented pass/fail rule.

## Aggregate metrics must carry their worst component

A full-history aggregate can look acceptable while one component fails
qualitatively. The retained Step 2 Case B development run reported a 9.80%
global relative L2 against the ten-team mean, yet its `p1`-z plateau had the
opposite sign from all ten official curves (ours −0.13…−0.67 mm; every team
+0.33…+0.79 mm), a median pointwise relative error of ~185% in that window.
Quoting the 9.80% alone inverted the meaning of the result.

Before quoting any aggregate (global relative L2, vector RMSE, RED), compute
the per-component, per-phase breakdown and report the worst entry with its
number and sign. A figure caption or verdict that states an aggregate without
its worst component is misleading even when every number in it is correct.

## Rerun with the corrected setup before diagnosing a discrepancy

The `p1`-z sign error looked like a convention or formulation defect: all ten
teams were positive, our run was negative, a classic signature of a fiber or
frame convention error. The retained run, however, predated the straight-wall
geometry and physical-coordinate-frame corrections. A later corrected-setup
Q1/P0 diagnostic reported a positive plateau against the official band and
relative L2 of 3.76%/5.94% at `p0`/`p1`. Its compact report does not bind the
complete release-grade source, execution, solver/deformation, and ten-team
input provenance, so this is a promising diagnostic observation rather than a
current reproduction result.

When a discrepancy appears in an old development run, the first experiment is
a current-setup reproduction, not a convention audit. Geometry and frame
corrections are live variables, and a qualitative failure in a stale run can
be a stale-setup artifact. Only if the discrepancy survives the corrected
stack do the sign-sensitive audits (fiber chirality, sheet frame, active-load
direction, landmark sampling) become the next step.

## Name the attribution boundary when one rerun changes several factors

The later Step 2 Case B diagnostic changed the geometry/frame corrections and
the volumetric formulation in the same execution (the old run used pointwise
`kappa`; the new one uses condensed Q1/P0 mean-`log(J)` local pressure).
Therefore "the geometry/frame correction fixed the sign" is not an isolated
attribution: the formulation change rides along. Its incomplete public
provenance also prevents promoting it as the best release trajectory. The
reproduction log records both boundaries, and a current-driver `std-kappa`
control remains available to separate the two if attribution ever matters.

When the goal is the best current answer, bundling corrections is acceptable
if recorded. When the goal is attribution to one factor, it is not. State
which goal a run serves, and keep the isolating control identified before
declaring the mechanism understood.
