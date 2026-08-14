# Cardiac benchmark application

This directory contains the cardiac-specific implementation described in the
repository [README](../../README.md). All model quantities use SI units.

See the [application interface reference](../../docs/API.md) for the complete
driver options, NPZ fields, helper signatures, and comparison-report schema.
See the [benchmark test matrix](../../docs/BENCHMARK_TEST_MATRIX.md) for the
progressive CoupFE validation layers, configuration identities, established
findings, and artifact layout.
See the answer-first
[benchmark reproduction status](../../docs/BENCHMARK_REPRODUCTION_STATUS.md)
for which Case A/Case B claims the retained comparisons do and do not support.
The [closed Case A stopping record](../../docs/CASE_A_STATUS.md) identifies the
selected fine generalized-alpha scientific result, its comparison metrics,
the interpolation and mesh-evidence limits, and the 2026-08-04 closeout
decision. The retained levels do not show monotonic improvement toward Simula.
The selected fine result now has a checked-in compact 101-point report and
figure. The eight-element F-bar reports remain archived open-tip
exact-configuration regression/history records, not current Benchmark 1
validation evidence.

## Paper benchmark comparisons

- [Case A active-contraction comparison](../../docs/figures/case_a_comparison.svg):
  selected CoupFE-Cardiac topologically closed 4×36×32 Hex8 Q1/P0 local pressure,
  consistent-mass source-matched generalized-alpha, `dt=0.001 s`, from the
  compact [`016a4f9` report](results/case_a_local_pressure_4x36x32_dt0p001.report.json),
  versus the benchmark all-team mean. Its Step 0A identity is explicitly
  `legacy-inferred` because the NPZ predates the later identity fields. Its
  producing source also predates the toolkit-matched straight-wall mapping and
  physical-coordinate structural-frame reconstruction, so the report remains
  bound to that earlier geometry/frame implementation. The
  archived [`6839c13` corrected-law report](results/archive/truncated_polar/case_a/case_a_fbar_1x2x4_dt0p002_corrected_switch.report.json)
  and [`62ad760` historical-law report](results/archive/truncated_polar/case_a/case_a_fbar_1x2x4_dt0p002.report.json)
  remain unchanged eight-element open-tip records and do not drive this
  figure.
- [Historical truncated-polar Case B pressure-loading comparison](../../docs/figures/archive/truncated_polar/case_b_comparison.svg):
  Step 0 CoupFE-Cardiac Q1/P0 local pressure, 2×36×48 Hex8, `dt=0.002 s`, backward
  Euler/PETSc SNES, from the
  [`e07993b` report](results/archive/truncated_polar/case_b/case_b_local_pressure_2x36x48_dt0p002.report.json),
  versus the benchmark all-team mean. This retained run used the historical
  truncated polar domain.
- [Current Step 0B tip-refined full-cycle comparison](../../docs/figures/step0b_tip_refine_full_cycle.svg):
  the closed `2x20x17` and `4x20x17`, `tip_refine=6.0` Q1/P0 runs completed
  1,000/1,000 increments on eight MPI ranks with consistent mass,
  source-matched generalized-alpha, and `dt=0.001 s`. The
  [dual-run comparison report](results/step0b_tip6p0_full_cycle_comparison.report.json)
  binds both external archives and the exact FEniCS/team inputs. The
  four-layer trajectory reduces full-shared-history FEniCS RMSE from
  1.0927/1.1685 to 1.0169/1.0730 mm at `p0`/`p1`, but its relaxation RMSE is
  16.2%/14.9% worse. This mixed phase response is numerical-sensitivity
  evidence, not convergence, rank equivalence, a physical-branch claim, or a
  validation pass.
- [Step 2 Case B full-cycle comparison](../../docs/figures/step2_case_b_comparison.svg):
  a completed four-rank, topologically closed 2×20×17, pointwise-`kappa`,
  source-matched generalized-alpha development run from the
  [`e9b7d90` report](results/step2_case_b_std_kappa_2x20x17_dt0p001.report.json),
  compared with the official ten-team range/mean and named Simula curve. The
  run has 9.8038% global relative L2, while its `p1-z` plateau has the opposite
  official sign. Its content-identified source predates the straight-wall and
  physical-frame corrections. The dirty-tree result is development evidence,
  not validation, convergence, or rank-independence evidence; see
  the [dated reproduction log](../../docs/STEP2_CASE_B_REPRODUCTION_LOG.md).
- [Step 2 Case B corrected-setup diagnostic](../../docs/figures/step2b_current_rerun_comparison.svg):
  a later Q1/P0 full-cycle trajectory whose compact
  [comparison report](results/step2b_current_rerun_comparison.report.json)
  lacks the complete source, execution, solver/deformation, and ten-team
  role/hash provenance required for release-grade evidence. It is retained as
  a promising diagnostic only, not a reproduction claim.

As of 2026-08-07, the [result index](results/) separates the pre-straight-wall
closed Case A and Step 2 records from the
[truncated-polar archive](results/archive/truncated_polar/),
which lists exact RED values, source identities, and diagnostics for the
historical open-tip configurations. It also indexes the current-source
dual-run `2x20x17`/`4x20x17`, `tip_refine=6.0` full-cycle Step 0B comparison
and the separate 0.32 s prefix diagnostic, plus the provenance-incomplete Step
2 corrected-setup diagnostic. The two Step 0B full-cycle trajectories are
nearly identical before snap (0.0098/0.0026 mm pairwise RMSE at `p0`/`p1`),
separate by at most 1.6331/1.4059 mm in the snap window, and recover to
0.0346/0.0518 mm separation at cycle end. Their full-cycle pairwise RMSE is
0.4076/0.3618 mm. The clean isolated four-layer replay completed and matched
the independently audited non-retained candidate to roundoff.
This is numerical sensitivity evidence, not convergence or evidence for
physically distinct solution branches.

The separate [Case B debugging
postmortem](../../docs/CASE_B_DEBUGGING_POSTMORTEM.md) records the failed and
null experiments, the late geometry and mass audits, and the discrete-method
differences still visible against the local FEniCS trajectory.

## Code map

- `geometry.py`: historical polar-ring geometry and the noncollapsed five-block
  closed Hex8 mesh with the toolkit-matched straight Cartesian wall ruling,
  boundary facets, output points, and topology-specific frame routing.
- `boundary_audit.py`: fail-fast exterior-face, reference-Jacobian, geometry
  measure, pressure-resultant/moment, Robin, and rigid-mode Robin checks before
  dynamics.
- `tbar_laplace.py`: reproducible Q1-Hex Laplace transmural-field generator
  for the closed multiblock mesh, with checked residual, exact mesh-coordinate/
  connectivity identity, and hashed metadata.
- `structural_directions.py`: shared serial/MPI reconstruction of the complete
  fiber, sheet, and sheet-normal frame from each physical point and Laplace
  coordinate using the pinned toolkit formula, including its apex branch.
- `material.py`: Holzapfel–Ogden, active, and viscous material definitions;
  compiled F-bar and standard Q1 kernel construction; and Gauss-point
  Gram–Schmidt sampling.
- `activation.py`: Benchmark 1 active-tension and pressure histories.
- `benchmark_parameters.py`: fail-closed Step 0 A/B and Step 2 Case B
  physical identities, exact parameter dictionaries, and runtime-source
  content manifests.
- `robin.py` and `pressure.py`: pericardial Robin and endocardial
  follower-pressure operators.
- `local_pressure.py`: application-owned Q1/P0 volumetric operator with one
  algebraically eliminated pressure per Hex8.
- `sampling.py`: checked reference Hex8 point location and isoparametric
  displacement interpolation for `p0` and `p1`.
- `newmark.py`: optional constant-average-acceleration Newmark kinematics.
- `generalized_alpha.py`: fixed Simula-source-matched generalized-alpha
  parameters, stages, and kinematic chain rules used by the closed Step 0 Case
  A/B and Step 2 Case B MPI companion paths.
- `solver.py`: guarded Core Newton and an optional persistent PETSc SNES path,
  both accepted before physical-state commit.
- `result_io.py`: atomic completed-result archive writing.
- `run.py`: serial Step 0 Case A/B and Step 2 Case B driver with backward Euler
  as the default and an explicit `--benchmark-step` identity.
- `run_mpi.py`, `distributed_solver.py`, `distributed_local_pressure.py`,
  `distributed_material.py`, and `distributed_mass.py`: additive PETSc
  companion with separate historical Q1/P0/lumped and closed
  `{pointwise-kappa,Q1/P0-local-pressure}`/consistent-mass contracts, including
  explicit `split` and `joint` material-evaluation modes. The source-matched
  generalized-alpha contract covers closed Step 0 Cases A/B and Step 2 Case B;
  only the pointwise-`kappa` Step 2 profile is labeled `paper-source-matched`.
- `diagnose.py`, `mesh_quality.py`, and `fiber_crosscheck.py`: deformation,
  mesh, and fiber-convention diagnostics.
- `post.py`: guarded comparison with separately downloaded Zenodo curves.
- `compare_step2b_case_b.py`: hash-pinned Step 2 Case B publisher comparator
  with exact physical/source identity and separate plumbing and numerical
  agreement results.
- `step2b_case_b_runtime_source_hashes.json`: byte-pinned dirty-tree source
  manifest for the retained Step 2 run, linked to repository snapshot
  `d06c3e9` so later compatibility edits do not relabel its provenance.
- `plot_retained_comparisons.py` and `plot_step2b_case_b.py`: distinct retained
  Step 0 and development Step 2 figure renderers from checked-in reports.
- `plot_step2b_current_rerun.py`: renderer for the separately labeled,
  provenance-incomplete corrected-setup Step 2 diagnostic; it accepts only the
  exact hash-manifested publisher files through the restricted NumPy unpickler
  and does not produce release-grade evidence.
- `compare_fenics_case_b.py`: five-input, hash-gated direct displacement
  comparison with the locally completed FEniCS P2-tetrahedron output; it
  accepts either qualified backward-Euler or source-matched generalized-alpha
  Step 0B archives; see the
  [comparison protocol](../../docs/CASE_B_FENICS_COMPARISON.md).
- `compare_mpi_rank_gate.py`: fail-closed comparison of the clean serial and
  one-, two-, and four-rank 2×20×17 trajectories through 0.32 s; see the
  [rank-gate protocol](../../docs/CASE_B_MPI_RANK_GATE.md).

See [`examples/REFERENCES.md`](../REFERENCES.md) for scientific sources,
file-level licenses, local checks, and reproducibility limits.

## Formulations

`--formulation fbar` selects the historical single-displacement-field Hex8
F-bar path.

`--formulation std-kappa` selects the standard Q1 Hex8 material kernel with
the paper's pointwise `kappa` penalty. It is the direct application path for
the paper volumetric energy; it is distinct from both F-bar and Q1/P0.

`--formulation local-pressure` builds the standard Q1 material kernel with its
bulk penalty set to zero, then composes `LocalPressureHex8Operator`. For each
element,

```text
p_e = K / V_e * integral(log(det(F))) dV,
P_vol = p_e F^(-T).
```

The MPI companion also exposes the separately identified experimental option
`--formulation local-pressure-paper`. It retains the same one-pressure-per-
element condensation but evaluates the paper scalar volume law at the
reference-volume-weighted geometric mean dilatation:

```text
m_e = 1 / V_e * integral(log(det(F))) dV,
p_e = K / 2 * (exp(2*m_e) - 1).
```

Its exact condensed tangent uses `dp_e/dm_e = K*exp(2*m_e)`. This is distinct
from both the historical linear `p_e=K*m_e` local law and the pointwise
`std-kappa` formulation; each has its own result and MPI implementation label.

A controlled fine Case A trial held the 4x36x32 closed mesh, boundary
conditions, fields, material parameters, consistent mass, generalized-alpha
scheme, 0.001 s step, and 8-rank iterative solver fixed and changed only this
scalar law and its exact tangent. Through 0.70 s it completed 700/700 steps
with zero domain rejections, but slightly worsened both Simula landmark
relative-L2 metrics (p0 8.273% to 8.348%; p1 11.480% to 11.578%) and did not
change the 0.65--0.67 s maximum-gap region materially. Treat
`local-pressure-paper` as an explicit experimental formulation, not the
default and not an established fix for the remaining Case A gap.

The element pressure is algebraic and is eliminated on every evaluation. The
returned displacement tangent includes both the derivative of `F^(-T)` and the
dependence of `p_e` on every element displacement degree of freedom. The
compiled material element continues to own cardiac fiber, viscous-history, and
active-stress state. This Q1/P0 path belongs to the cardiac application and is
being evaluated alongside F-bar; it is not a claim that the 2026-06-27 F-bar
campaign used local pressure.

The focused tests check affine volumetric response, isochoric response,
tangent symmetry and finite differences, and failures for degenerate reference
cells, inverted deformation, and non-finite inputs. These are component checks,
not a ventricular mesh-convergence result.

## Output-point sampling

The current driver samples `p0` and `p1` through the reference Hex8 field. It
first selects candidate elements by reference-coordinate bounding boxes, then
uses a checked damped Newton inversion of the trilinear isoparametric map. A
location is accepted only with a positive nonsingular Jacobian, natural
coordinates inside the reference cube, finite shape weights, and a
scale-aware reconstruction check. If a point lies on a shared boundary, the
lowest valid element index is selected deterministically.

The result archive records `point_sampling=hex8_reference_isoparametric` plus
the chosen element, natural coordinates, eight weights, and reconstruction
error for each output point. This replaces the global Delaunay-tetra policy in
the existing historical Case A record; that older record remains labeled with
the policy that produced it.

This check separates point-location correctness from field accuracy. The
recorded reconstruction error verifies that the requested physical coordinate
was located; it does not make the computed displacement exact. The returned
value is the trilinear Q1 Hex8 displacement field itself, so evaluating the
shape functions adds no identified sampler error beyond roundoff. The Q1
spatial approximation can still differ from a quadratic P2-tetrahedral result
even when the physical point is reconstructed exactly. The current Case A
evidence does not isolate that contribution.

## Nonlinear solvers and result integrity

`--nonlinear-solver core-newton` uses Core's serial Newton routine and an
application-level final residual check. `--nonlinear-solver petsc-snes` uses a
persistent serial PETSc context with parameter values recovered from the
2026-06-27 development adapter:

```text
SNES newtonls, line search bt, KSP preonly, PC lu,
rtol=1e-9, atol=1e-10, stol=1e-12, max_it=60.
```

The PETSc implementation is new even though those parameter values are
preserved. A positive SNES reason is necessary but not sufficient: the
application independently reassembles the accepted residual and requires
`|R| <= max(atol, rtol*|R_initial|)`. KSP status and displacement finiteness are
also checked before one state commit.

Both drivers default to joint element evaluation, retaining the established
paired material R/K cache. With explicit `--element-evaluation split`, residual
callbacks use the current Core-generated material residual-only entry. In the
MPI Q1/P0 companion, split residual callbacks also evaluate only the algebraic
pressure residual; Jacobian callbacks construct the condensed pressure tangent
and obtain the joint material R/K block.
Both modes use the same equations, tolerances, invalid-trial policy, and final
residual rule.

For the Q1/P0 path, Core Newton calls the operator's `max_step` method to bound
a correction to a determinant-valid trial before residual backtracking. The
PETSc residual callback catches only `InvalidDeformationError`, writes IEEE
positive infinity as that trial residual, and lets `bt` shorten the trial. Its
per-step diagnostics count these rejections and retain the last detail.
Jacobian exceptions and unrelated residual exceptions abort; invalid initial,
accepted, or final states remain fail-closed. This mechanism does not alter the
listed tolerances or the independently reassembled final residual rule.

Successful archives retain the formulation, material and local-pressure bulk
parameters, mesh and apex treatment, integrator, fiber and point-sampling
policies, source identities, package versions, solver configuration, and one
nonlinear diagnostic record per step. At peak load they also retain
Gauss-point `det(F)` for both formulations and the eliminated element pressure
for the Q1/P0 path. A failed solve exits before the completed-result writer;
the writer itself uses atomic replacement.

## Run

State the formulation, integrator, solver, mesh, and time step explicitly in a
retained command. For example:

```bash
python examples/cardiac_benchmark/tbar_laplace.py \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --out tbar_closed_nt2_core20_rad17.npy

python examples/cardiac_benchmark/run.py \
  --benchmark-step 0 --case B --formulation fbar --integrator be \
  --nonlinear-solver petsc-snes \
  --element-evaluation joint \
  --nt 2 --nmu 24 --ntheta 32 --apex-offset 0.2 \
  --mass lumped --fiber-sampling cg1 \
  --dt 0.002 --tend 1.0 --out caseB_fbar.npz

python examples/cardiac_benchmark/run.py \
  --benchmark-step 0 --case B --formulation local-pressure --integrator be \
  --nonlinear-solver petsc-snes \
  --element-evaluation joint \
  --nt 2 --nmu 24 --ntheta 32 --apex-offset 0.2 \
  --mass lumped --fiber-sampling cg1 \
  --dt 0.002 --tend 1.0 --out caseB_local_pressure.npz

python examples/cardiac_benchmark/run.py \
  --benchmark-step 0 --case B --formulation std-kappa --integrator be \
  --nonlinear-solver petsc-snes --element-evaluation joint \
  --mesh-topology closed-multiblock \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --mass consistent --material-eta 100 \
  --tbar-laplace tbar_closed_nt2_core20_rad17.npy \
  --fiber-sampling gp-direct \
  --dt 0.001 --tend 1.0 --out caseB_closed_candidate.npz

mpiexec -n 4 python examples/cardiac_benchmark/run_mpi.py \
  --benchmark-step 0 --case B --element-evaluation split \
  --nt 2 --nmu 24 --ntheta 32 --apex-offset 0.2 \
  --dt 0.002 --tend 1.0 --out caseB_local_pressure_mpi.npz

# Closed-configuration rank gate: first create the exact serial prefix.
python examples/cardiac_benchmark/run.py \
  --benchmark-step 0 --case B --mesh-topology closed-multiblock \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --formulation std-kappa --mass consistent --material-eta 100 \
  --tbar-laplace tbar_closed_nt2_core20_rad17.npy \
  --fiber-sampling gp-direct --element-evaluation joint \
  --nonlinear-solver petsc-snes \
  --dt 0.001 --tend 0.32 --load-horizon 1.0 \
  --out caseB_closed_snap_serial.npz

# Repeat this MPI command at 1, 2, and 4 ranks with distinct output names.
mpiexec -n 4 python examples/cardiac_benchmark/run_mpi.py \
  --benchmark-step 0 --case B --mesh-topology closed-multiblock \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --formulation std-kappa --mass consistent --material-eta 100 \
  --tbar-laplace tbar_closed_nt2_core20_rad17.npy \
  --fiber-sampling gp-direct --element-evaluation joint \
  --dt 0.001 --tend 0.32 --load-horizon 1.0 \
  --out caseB_closed_snap_mpi_rank4.npz

# The same closed discretization supports Step 0 Case A. Case A applies the
# benchmark active-tension history and keeps endocardial pressure identically
# zero.
mpiexec -n 4 python examples/cardiac_benchmark/run_mpi.py \
  --benchmark-step 0 --case A --mesh-topology closed-multiblock \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --formulation std-kappa --mass consistent --material-eta 100 \
  --tbar-laplace tbar_closed_nt2_core20_rad17.npy \
  --fiber-sampling gp-direct --element-evaluation joint \
  --dt 0.001 --tend 0.32 --load-horizon 1.0 \
  --out caseA_closed_snap_mpi_rank4.npz

# Source-matched Case A temporal diagnostic. This is deliberately a separate
# implementation identity; local pressure is evaluated at the same alpha_f
# stage as the generated cardiac material.
mpiexec -n 4 python examples/cardiac_benchmark/run_mpi.py \
  --benchmark-step 0 --case A --mesh-topology closed-multiblock \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --formulation local-pressure --mass consistent --material-eta 100 \
  --tbar-laplace tbar_closed_nt2_core20_rad17.npy \
  --fiber-sampling gp-direct --element-evaluation joint \
  --integrator generalized-alpha \
  --dt 0.001 --tend 0.32 --load-horizon 1.0 \
  --out caseA_closed_local_pressure_ga_mpi_rank4.npz

# Step 2 Case B: active stress plus pressure, canonical full-cycle profile.
mpiexec -n 4 python examples/cardiac_benchmark/run_mpi.py \
  --benchmark-step 2 --case B --mesh-topology closed-multiblock \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --formulation std-kappa --mass consistent --material-eta 100 \
  --tbar-laplace tbar_closed_nt2_core20_rad17.npy \
  --fiber-sampling gp-direct --element-evaluation joint \
  --integrator generalized-alpha \
  --dt 0.001 --tend 1.0 --load-horizon 1.0 \
  --out mpi4_ga_full.npz

python examples/cardiac_benchmark/compare_mpi_rank_gate.py \
  --serial caseB_closed_snap_serial.npz \
  --mpi1 caseB_closed_snap_mpi_rank1.npz \
  --mpi2 caseB_closed_snap_mpi_rank2.npz \
  --mpi4 caseB_closed_snap_mpi_rank4.npz \
  --report caseB_closed_snap_rank_gate.json
```

The field generator writes both the NPY and
`tbar_closed_nt2_core20_rad17.meta.json`. Keep them together: the serial
driver validates the native sidecar's field hash, mesh and boundary identity,
and linear-solve diagnostics before compilation, and retains portable hashes
for both files in the completed archive.

The F-bar and local-pressure serial dynamics commands intentionally exercise
historical truncated-domain code paths. The closed `std-kappa` serial command
is the current straight-wall/physical-frame analogue of the Step 0
configuration used by the completed 2026-08-02 development run and keeps the
paper physical parameters; it is not the exact historical configuration. That
dirty-tree historical run cannot serve as public retained evidence. Each MPI
run is one member of the rank-equivalence protocol, not rank-independence
evidence by itself. `--load-horizon 1.0` makes every shortened pressure history
the exact prefix of the production schedule. The comparison command accepts
only the fixed 2×20×17 setup and writes a report only after every rank passes.
That record is required before claiming rank equivalence, but
configuration-specific eight-rank diagnostic evidence does not depend on such
a claim. Retain convergence and deformation diagnostics, and record the actual command rather
than describing a nominal “fine” mesh.

The optional `--integrator newmark` path shares Newmark velocity kinematics
between inertia and Robin damping while material viscosity remains a backward
strain difference. Use `--integrator be` for the coherent first-order path.

The MPI-only `--integrator generalized-alpha` path is fail-closed to closed
Step 0 Cases A/B or Step 2 Case B with consistent mass, GP-direct fibers, a
validated Laplace field, and either an explicitly named local-pressure variant
or the `std-kappa` formulation. These are supported code contracts; support
does not by itself claim a retained trajectory for every case/formulation.
It matches the source defaults at peeled upstream commit
`325d17d850c2e2032abb85a4191a5795d3008ab7`:
`alpha_m=0.2`, `alpha_f=0.4`, `gamma=0.7`, and `beta=0.36`. Consistent inertia
is evaluated at `a[n+1-alpha_m]`; material, condensed local pressure, Robin
spring/dashpot, and active tension are evaluated at the `1-alpha_f` stage;
the load ODE is integrated from `t=0` over the fixed load horizon and sampled
at `t[n+1]-alpha_f*dt`; and viscous stress uses
`sym(F_stage.T @ grad(v_stage))`. The generated material stores accepted
Gauss-point displacement, velocity, and acceleration gradients so this rate is
not approximated by the older Green--Lagrange backward difference. Full and
short-prefix runs both use the source's canonical 1 s load-integration horizon;
a different generalized-alpha horizon is rejected.

Focused material-point, compiled-element tangent/state-commit, local-pressure
chain-rule, and analytic PETSc gates cover this path. The retained loaded
coarse ventricular 1/2/4-rank gate and selected fine 8-rank, 1 s log-law
trajectory are Step 0A records whose producing source predates the straight-
wall and physical-frame corrections. They validate only their named execution
layers; they do not establish spatial convergence or equivalence of Q1 Hex8
and P2 tetrahedral approximation spaces. Step 0B is supported by the current
staged implementation. Its two clean current-source 0.32 s gates are bound by
the compact prefix diagnostic below. Separate closed `2x20x17` and `4x20x17`,
`tip_refine=6.0` runs have completed the full 1 s cycle on eight ranks and are
included in the dual-run full-cycle comparison. The clean isolated replay of
the four-layer archive is retained; the first dirty-tree full-cycle candidate
remains only a provenance diagnostic.

For clean source archives or other non-Git execution contexts, source identity
can be asserted explicitly:

```bash
export COUPFE_CARDIAC_APP_REVISION=<40-hex-public-app-commit>
export COUPFE_CARDIAC_TREE_STATE=clean
export COUPFE_CORE_REVISION=<40-hex-public-core-commit>
export COUPFE_CORE_TREE_STATE=clean
```

These variables are assertions, not independent discovery or reachability
checks. Do not label an edited tree `clean`.

## Result records

The [`results/`](results/) directory separates the code's generated NPZ files
from the records suitable for source control. Its root retains the selected
compact topologically closed Case A report and the topologically closed Step 2
development report and log. Both producing sources predate the straight-wall
geometry mapping and physical-coordinate frame reconstruction; they are not
relabeled as current-source output. A separate compact Step 0B prefix report
binds the clean current-source coarse and wall-only gates and its historical
decision to pause the 1 s extension. That decision is superseded for the
separately configured `tip_refine=6.0` runs, whose dual-run full-cycle
comparison report and figure are also retained. The earlier `4x20x17`,
`tip_refine=6.0` prefix remains an external forensic archive; its shared
history agrees with the full result to roundoff.
Historical Case A and Step 0 Case B records made with the open-tip
`polar_ring`, `apex_offset=0.2` geometry live in the
[`truncated-polar archive`](results/archive/truncated_polar/); they are not
current Benchmark 1 validation evidence. The compact fine Case A report binds
the external NPZ/manifest/stdout identities without distributing the
machine-local transcript. Generated NPZ archives and external reference
pickles remain unbundled.

| App | Case and formulation | Mesh, `dt` | Accepted steps | Retained files |
|---|---|---|---:|---|
| **`2458e7c`** | Step 0 B, Q1/P0 local pressure, generalized-alpha | current closed 4×20×17, `tip_refine=6.0`, 0.001 s | 1,000/1,000; clean isolated replay complete | Dual-run report/figure; **clean NPZ SHA-256 `1e333b29b05f01dedce9272b32b82ce6ccfda56036c1c8f57eb395b9b4494800`, stdout SHA-256 `0810a9d0c944345464f48917c707a7d3df7187b543c74efd4b0384542b95aed7`, elapsed `1778.4 s (29.6 min)`** |
| `ae2c2eb` | Step 0 B, Q1/P0 local pressure, generalized-alpha | current closed 2×20×17, `tip_refine=6.0`, 0.001 s | 1,000/1,000 | `step0b_tip6p0_full_cycle_comparison.report.json` and `docs/figures/step0b_tip_refine_full_cycle.svg` (external archive and references hash-bound) |
| `a5824d3` | Step 0 B, Q1/P0 local pressure, generalized-alpha | current closed 4×20×17, `tip_refine=6.0`, 0.001 s | 320/320; 0.32 s controlled prefix | External NPZ SHA-256 `774a7dc5fc970bb744ff0188f0f428cff54c532b401a365643f4b626584d7acf`; prefix-continuity record |
| `056c02d` | Step 0 B, Q1/P0 local pressure, generalized-alpha | current straight-wall/physical-frame 2×20×17 and wall-only 4×20×17, 0.001 s | 320/320 each; diagnostic prefix | `step0b_case_b_clean_frame_0p32.report.json` (historical stop decision; external archives/logs/metrics hash-bound) |
| `016a4f9` | Step 0 A, Q1/P0 local pressure, generalized-alpha | pre-straight-wall closed 4×36×32, 0.001 s | 1,000/1,000 | `case_a_local_pressure_4x36x32_dt0p001.report.json` (compact derived report; external NPZ/log hash-bound) |
| `e9b7d90` dirty/content-identified | Step 2 B, pointwise `kappa`, generalized-alpha | pre-straight-wall closed 2×20×17, 0.001 s | 1,000/1,000 | `step2_case_b_std_kappa_2x20x17_dt0p001.report.json` and `.raw.stdout.txt` |
| `62ad760` | A, F-bar, historical material law | 1×2×4, 0.002 s | 500/500 | `archive/truncated_polar/case_a/case_a_fbar_1x2x4_dt0p002.report.json` and `.raw.stdout.txt` |
| `6839c13` | A, F-bar, corrected material law | 1×2×4, 0.002 s | 500/500 | `archive/truncated_polar/case_a/case_a_fbar_1x2x4_dt0p002_corrected_switch.report.json` and `.raw.stdout.txt` |
| `62ad760` | B, Q1/P0 local pressure | 2×12×16, 0.004 s | 250/250 | `archive/truncated_polar/case_b/case_b_local_pressure_2x12x16_dt0p004.report.json` and `.raw.stdout.txt` |
| `62ad760` | B, Q1/P0 local pressure | 2×12×16, 0.002 s | 500/500 | `archive/truncated_polar/case_b/case_b_local_pressure_2x12x16_dt0p002.report.json` and `.raw.stdout.txt` |
| `62ad760` | B, Q1/P0 local pressure | 2×24×32, 0.002 s | 500/500 | `archive/truncated_polar/case_b/case_b_local_pressure_2x24x32_dt0p002.report.json` and `.raw.stdout.txt` |
| `e07993b` | B, Q1/P0 local pressure | 2×24×32, 0.004 s | 250/250 | `archive/truncated_polar/case_b/case_b_local_pressure_2x24x32_dt0p004.report.json` and `.raw.stdout.txt` |
| `e07993b` | B, Q1/P0 local pressure | 2×36×48, 0.002 s | 500/500 | `archive/truncated_polar/case_b/case_b_local_pressure_2x36x48_dt0p002.report.json` and `.raw.stdout.txt` |
| `62ad760` | B, F-bar | 2×24×32, 0.002 s | 500/500 | `archive/truncated_polar/case_b/case_b_fbar_2x24x32_dt0p002.report.json` and `.raw.stdout.txt` |
| `62ad760` | B, F-bar | 2×36×48, 0.002 s | 500/500 | `archive/truncated_polar/case_b/case_b_fbar_2x36x48_dt0p002.report.json` and `.raw.stdout.txt` |

Six reports identify clean application checkpoint `62ad760`; the newer
2×24×32, `dt=0.004 s` and 2×36×48, `dt=0.002 s` reports identify clean
checkpoint `e07993b`; those two groups identify Core checkpoint `454f73`.
The compact selected fine Case A report identifies clean application/Core
checkpoints `016a4f9`/`e2f42ed`, retains the exact 101-point CoupFE and
ten-team-mean curves, and drives the public Case A figure. It labels Step 0A
`legacy-inferred` because the external 16.2 MB archive predates explicit
benchmark fields. The corrected-law `6839c13` report remains the regression
record for its eight-element F-bar configuration, and the `62ad760` Case A
report remains historical-law evidence. See the
[Case A stopping record](../../docs/CASE_A_STATUS.md). The Step 0 v2 reports use
reference-Hex8 output sampling, retain the full `p0`/`p1` histories and solver
diagnostics, and identify the exact external curves used for RED. The distinct
Step 2 report binds its external full NPZ and runtime-source contents while
retaining the derived 101-point CoupFE/publisher curves and full-precision
comparison metrics needed for review and figure regeneration. Its frozen
runtime-source file map remains separately reviewable after current-source
compatibility changes. The
checkpoint-`62ad760` 2×24×32,
`dt=0.002 s` Case B rows hold the mesh, time step, geometry, sampler, source,
and solver fixed while changing the volumetric formulation. That controlled
pair reports a formulation-dependent output difference at one configuration;
it is not an accuracy ranking or mesh/time convergence study. The two 2×12×16
Q1/P0 runs provide a bounded two-step time-sensitivity check. The additional
2×36×48 F-bar output is a second retained spatial resolution, not a
spatial-convergence claim.

All Step 0 Case B rows in this table use `apex_offset=0.2`. They are preserved,
source-identified truncated-domain results; none is silently relabeled as
output from the new closed mesh. The separately labeled Step 2 row uses the
closed multiblock mesh.

In the `e07993b` run, PETSc rejected 46 invalid Q1/P0 trial residuals during
steps 132–133 and backtracked to valid trials. All 250 accepted/final states
met the unchanged domain and residual checks; no rejected trial was committed.
This is evidence that the recovery path operated for that named run, not a
tolerance change or an accuracy, convergence, validation, or bifurcation
claim. The 2×24×32 Q1/P0 `dt=0.004 s` and `dt=0.002 s` records have different
application checkpoints and therefore form a cross-checkpoint sensitivity
record, not a controlled time-step study.

The `e07993b` 2×36×48 Q1/P0 run rejected 168 invalid trial residuals at steps
277 and 279 before backtracking to valid trials. All 500 accepted/final states
met the unchanged checks, with largest final-residual/threshold ratio 0.973921.
The F-bar record at the same nominal mesh and time step comes from `62ad760`,
so the two outputs can be inspected side by side but do not form a controlled
one-factor formulation comparison or an accuracy, convergence, reproduction,
validation, or bifurcation result.

The directory also preserves the reduced historical Case A result tied to its
own recorded revisions. It uses 24 nodes, eight open-apex Hex8 elements, 500
backward-Euler steps, and the older Delaunay-tetra output sampler. It is not a
paper-curve comparison and is not output from the current Hex8 sampler.

The legacy 2026-06-27 Case B F-bar observations and retained Step 0 Case B
comparisons are documented together in
[`docs/CASE_B_STATUS.md`](../../docs/CASE_B_STATUS.md). The old table remains
legacy-reported history because its raw archives and logs are absent. Step 2
has a separate [reproduction log](../../docs/STEP2_CASE_B_REPRODUCTION_LOG.md).

## External comparison

The reference archive is about 23.2 GB and is not vendored:

```bash
curl -L \
  -o benchmark_article_data.zip \
  https://zenodo.org/api/records/14260459/files/benchmark_article_data.zip/content
echo "75602be4777c4ca2262c2bcfd2134b15  benchmark_article_data.zip" | md5sum -c -
unzip benchmark_article_data.zip
export CARDIAC_BENCHMARK_DATA_DIR="$PWD/benchmark_article_data"
```

Only load pickle files from this verified, trusted Zenodo archive. It is
licensed CC BY 4.0; cite the dataset and record any exact curve files used.
Install Matplotlib only when requesting a plot:

```bash
python -m pip install -e ".[dev,reference]"
python examples/cardiac_benchmark/post.py \
  caseB_local_pressure.npz --case step_0B --plot caseB.png
```

That `post.py` command is the Step 0 path. Step 2 Case B uses its exact
publisher manifest and distinct comparator:

```bash
python examples/cardiac_benchmark/compare_step2b_case_b.py \
  --publisher-data-dir /path/to/results_time_curves/data \
  --coupfe-run mpi4_ga_full.npz \
  --output step2_case_b_std_kappa_2x20x17_dt0p001.report.json

python examples/cardiac_benchmark/plot_step2b_case_b.py --canonical
```

The checked-in Step 2 report embeds the derived 101-point CoupFE, official
ten-team range/mean, and named Simula curves, so the plot command needs neither
the external NPZ nor the publisher pickles. Recomputing the report still needs
both. The comparator quantifies agreement without inventing a validation
threshold.

Retain the comparison command, result checksum, complete RED output, external
archive and input-file identities, source revisions, and environment. RED is
reported without a repository-defined pass threshold. See the [benchmark data
and comparison guide](../../docs/BENCHMARK_COMPARISON.md).

## Geometry and interpretation limits

The historical `polar-ring` topology with `apex_offset=0.2` avoids degenerate
elements by truncating the domain. At 2×36×48 it removes 1.9335 mm of the tip,
adds a 2.672330 cm² traction-free annular surface, and makes the undeformed
pressure resultant 4.5937% low. Setting the offset to zero instead creates
degenerate collapsed cells. These are retained historical/demonstration
choices, not the closed paper domain.

The pre-solve pressure audit keeps that distinction explicit. For the
historical topology it independently closes the discrete endocardial surface
with its polygonal base and terminal rings, then checks the signed resultant
and moment against that declared truncated surface. The analytic closed-base
projection remains recorded only as a comparator. For the closed topology it
remains the pressure-audit target, including the signed-load reversal check.

`--mesh-topology closed-multiblock` selects the five-block square-to-disk
implementation. The current mapping follows the toolkit meridional geometry:
corresponding endocardial and epicardial points are joined by straight
Cartesian segments. Its default 2×20×17 construction has 5,403 nodes and 3,520
Hex8 elements, no additional tip boundary, and all exterior faces classified
exactly once. Extended sampling gives minimum `det(J)=1.23646e-9 m³`, maximum
Jacobian condition `7.87649`, and minimum scaled Jacobian `0.258135`. Relative
to the retained reference, wall volume and endocardial, epicardial, and base
areas differ by `-0.1229%`, `+0.0695%`, `+0.0090%`, and `-0.0565%`; the unit-
pressure resultant differs from the analytic base projection by `0.1028%`.
The driver audits those properties, the pressure resultant and moment, and
Robin matrices before integration and stores the record as
`pre_solve_audit_json`. Passing the gate establishes the checked setup; it does
not by itself establish convergence or agreement with the result curves.

The Robin record also reports the spring stiffness seen by rigid translations
and by rotation about the long axis, with the latter split into base
full-vector and epicardial normal-only contributions. These values diagnose
the discrete boundary geometry: in particular, a faceted epicardium can resist
a rotation that is tangential on the smooth surface. They do not add a
mesh-dependent pass/fail or validation target and do not alter the boundary
condition.

The completed 2026-08-02 development run used the same default mesh counts but
predates both the straight-wall mapping and physical-coordinate frame
reconstruction. It used `dt=0.001 s`, consistent mass, Laplace tbar, GP-direct
fibers, and the paper `kappa` and `eta`. It finished 1,000/1,000 increments.
Vector RMSE against the
local FEniCS trajectory is 1.5803/1.7903 mm at `p0`/`p1`, with snap onset
about 5.3/5.6 ms early. The dominant x/z paths capture the qualitative snap
and recovery in this development comparison, while the y histories and parts
of unloading remain outside the ten-team spread. Because
the application tree was dirty, this is development evidence awaiting an exact
clean-source rerun; it is not placed in `results/`.

The separate 2026-08-05 Step 2 Case B run used the same topologically closed,
pre-straight-wall 2×20×17 geometry and pre-physical-frame implementation with
pointwise `kappa`, consistent mass, GP-direct Laplace fibers, and source-matched
generalized-alpha, applying active stress and pressure together.
It completed 1,000/1,000 increments on four ranks and is represented in
`results/` by a content-identified dirty-tree comparison report and normalized
stdout; the external NPZ is not bundled. Its bounded global errors coexist with
an opposite-sign `p1-z` plateau, so it remains development evidence rather than
a validation or convergence result.

Paper-parameter reproduction keeps `kappa=1.0e6 Pa` and `eta=100 Pa s`.
Earlier changed-viscosity and isotropic runs are historical forensic records;
the eta=50 attempt stopped without a completed archive and no tuning was
accepted. The current public driver rejects nonpaper `eta` and isotropic
requests before setup.
For backward-Euler archives, the trilinear Q1 Hex8/backward-Euler versus
quadratic P2-tetrahedral/generalized-alpha gap must remain stated. For the new
source-matched temporal path, the Q1 Hex8 versus quadratic P2-tetrahedral
spatial/formulation gap remains and must still be stated with any result.

The supplied FEniCS point-stress archive is quarantined as a quantitative
oracle because its postprocessor does not reload the accepted velocity and
acceleration state, projects to DG1, and uses a dimensionally inconsistent von
Mises expression. A stress comparison requires corrected accepted-state,
element-interior or quadrature-point postprocessing.

An open-apex mesh has no unique closed-cavity volume unless an explicit cap
policy is defined and checked.

The local-pressure operator, follower-load controls, and retained runs show
what these code paths compute at named configurations. Whole-case
interpretation still depends on
spatial and temporal sensitivity, apex treatment, deformation Jacobians,
complete `p0`/`p1` histories, fiber convention, and comparison with the trusted
`step_0B` curves.
