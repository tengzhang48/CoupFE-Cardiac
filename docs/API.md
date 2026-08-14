# Application interface reference

This page documents the executable and Python interfaces used by the public
examples. CoupFE-Cardiac is a source application, not a separately versioned
Python library. Names and signatures under `examples/` may change as the
benchmark application develops; downstream code should pin a repository commit
and treat the result schema and recorded source identity as part of its input.

## Serial driver CLI

Run the driver from a repository clone:

```bash
python examples/cardiac_benchmark/run.py [OPTIONS]
```

The driver has no positional arguments.

| Option | Default | Meaning |
|---|---:|---|
| `--case {A,B}` | `A` | Case identity within the selected benchmark step: Step 0 A is active-only, Step 0 B is pressure-only, and Step 2 B combines active stress and ventricular pressure. |
| `--benchmark-step {0,2}` | `0` | Select the Benchmark 1 protocol identity. Step 0 supports Cases A/B; Step 2 currently supports Case B only and records the distinct material/load contract. |
| `--formulation {fbar,local-pressure,std-kappa}` | `fbar` | Select historical F-bar, application-owned Q1/P0 condensed local pressure, or the standard Q1 kernel with the paper's pointwise `kappa` penalty. |
| `--integrator {be,newmark}` | `be` | Backward Euler is the coherent first-order path. Newmark applies constant-average-acceleration kinematics to inertia and Robin damping while material viscosity remains a backward strain difference. |
| `--nonlinear-solver {core-newton,petsc-snes}` | `core-newton` | Use guarded Core Newton or the application-owned serial PETSc SNES path. PETSc requires the `mpi` optional dependencies even though this driver uses one rank. |
| `--element-evaluation {joint,split}` | `joint` | Keep one cached paired material R/K evaluation, or explicitly use Core's native residual-only material entry until a tangent is needed. No automatic selection is performed. |
| `--nt INT` | `2` | Number of elements through the wall. |
| `--mesh-topology {polar-ring,closed-multiblock}` | `closed-multiblock` | Select the noncollapsed five-block benchmark domain. `polar-ring` is an explicit historical open-tip option retained only for archived-result reproduction. |
| `--nmu INT`, `--ntheta INT` | `12`, `16` | Longitudinal and circumferential counts for `polar-ring`. |
| `--ncore INT`, `--nradial INT` | `20`, `17` | Even central-square count per side and outer radial layers for `closed-multiblock`. |
| `--core-half-width FLOAT` | `0.36` | Central-square half-width in unit-disk coordinates for `closed-multiblock`; must lie in `(0.1, 1/sqrt(2))`. |
| `--dt FLOAT` | `0.002` | Time increment in seconds; must be finite and positive. |
| `--tend FLOAT` | `1.0` | End time in seconds; must be a positive integer multiple of `dt`. |
| `--load-horizon FLOAT` | `tend` | Integrate the benchmark load once through this horizon, then use the exact prefix through `tend`. It must be at least `tend` and share the `dt` grid. Use `1.0` for shortened rank gates so their load is byte-identical to the production schedule. |
| `--apex-offset FLOAT` | `0.2` | Historical `polar-ring` tip truncation in radians. Zero collapses that topology and creates degenerate cells; the option does not control `closed-multiblock`. |
| `--mass {consistent,lumped}` | `consistent` | Use the assembled Q1 Hex8 consistent mass or the historical row-summed diagonal mass. |
| `--material-eta FLOAT` | `100.0` | Material viscosity in Pa s. The public benchmark driver requires the paper value `100` and stops before setup for any other value. |
| `--isotropic` | off | Retired forensic-control spelling. The public benchmark driver stops before setup if requested; the paper anisotropic coefficients remain fixed. |
| `--tbar-laplace PATH.npy` | none | Use a checked per-node Laplace transmural field instead of the analytic layer coordinate. The sibling `PATH.meta.json` emitted by `tbar_laplace.py` is mandatory and is validated before compilation. |
| `--fiber-sampling {cg1,gp-direct}` | `cg1` | Interpolate nodal directions then re-orthonormalize, or evaluate the analytic rule directly at each Hex8 Gauss point. |
| `--viscous-evidence-out PATH.npz` | none | Retired forensic-control spelling. The public benchmark driver stops before setup if it or its former window options are requested; no eta split or physical-parameter sweep is run. The historical implementation remains documented and tested separately. |
| `--raw-helix` | off | Use `asign=+1` instead of the default pinned-formula-matched helix convention. This changes prescribed helix handedness. |
| `--perturb FLOAT` | `0.0` | For a positive value, add deterministic Gaussian node noise with this standard deviation in metres. Intended for diagnostics, not an accuracy claim. |
| `--out PATH` | case/step-derived name | Completed NPZ result destination. Step 0 defaults to `caseA_full.npz` or `caseB_full.npz`; Step 2 uses `benchmark1_step2_caseB_full.npz`. |
| `--build-dir PATH` | temporary directory | Persistent directory for generated/compiled kernel files. Generated kernels are build artifacts, not public source files. |

Example with all result-defining choices explicit:

```bash
python examples/cardiac_benchmark/run.py \
  --benchmark-step 0 --case B --formulation std-kappa --integrator be \
  --nonlinear-solver petsc-snes \
  --element-evaluation joint \
  --mesh-topology closed-multiblock \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --mass consistent --material-eta 100 \
  --tbar-laplace tbar_closed_nt2_core20_rad17.npy \
  --fiber-sampling gp-direct \
  --dt 0.001 --tend 1.0 --load-horizon 1.0 \
  --build-dir build --out caseB_closed_candidate.npz
```

This is the Step 0 configuration used by the completed 2026-08-02
closed-domain development comparison. It remains a Q1 Hex8/backward-Euler discretization,
whereas the local reference uses a quadratic P2 displacement element on
tetrahedral cells and generalized alpha. The default closed candidate has
16,209 vector displacement DOFs and two Q1 layers through the wall. A spatial
study should include a wall-only refinement and at least three fully qualified
levels; approaching the reference's 89,700 vector DOFs is a resolution
comparison, not proof that the two approximation spaces are equivalent. The
development run came from a dirty application tree; a public retained report
requires an exact clean-source rerun.

The driver writes an archive only after every requested increment completes.
Solver, model, sampling, or serialization failures propagate as a nonzero
process exit. The completed-result writer uses atomic replacement so a failed
write does not replace an existing result.

### Closed-mesh Laplace transmural field

Generate the nodal field consumed by `--tbar-laplace` with the same closed-mesh
parameters as the dynamics run:

```bash
python examples/cardiac_benchmark/tbar_laplace.py \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --out tbar_closed_nt2_core20_rad17.npy
```

The utility solves Q1-Hex Laplace with `tbar=0` on the endocardium, `tbar=1`
on the epicardium, and the natural zero-flux condition on the base. It fails
on a nonpositive Gauss-point Jacobian, invalid bounds or boundary values, or a
linear residual outside its declared tolerance. It atomically writes the NPY
field and adjacent `.meta.json` record; the record includes the NPY SHA-256 and
uses only the output filename, not a machine-local absolute path. The dynamics
driver requires that native sidecar and checks its schema, field name and hash,
mesh identity, boundary conditions, and passing solver diagnostics before it
compiles an element or starts time integration.

## Distributed companion CLI

Launch the additive MPI companion with a matched MPI/PETSc environment:

```bash
mpiexec -n 4 python examples/cardiac_benchmark/run_mpi.py [OPTIONS]
```

The driver accepts explicit, non-interchangeable contracts:

- historical `polar-ring` + `local-pressure` + lumped mass + CG1 fibers; or
- closed Case A/Case B `closed-multiblock` +
  `{std-kappa,local-pressure,local-pressure-paper}` +
  consistent mass + GP-direct fibers + a validated `--tbar-laplace` field.

Mixed combinations stop before MPI setup and are not labeled as either
implementation. The closed contract also fixes `eta=100 Pa s`, the benchmark
helix convention, and zero perturbation. Backward Euler remains available for
the supported Step 0 and Step 2 identities. The distinct
`--integrator generalized-alpha` path accepts closed Step 0 Cases A/B and
Step 2 Case B with fixed source parameters `alpha_m=0.2`, `alpha_f=0.4`,
`gamma=0.7`, and `beta=0.36`. Complete consistent-mass rows are assembled
from every Hex8 that touches a rank-owned row.

In the source-matched path, inertia is staged at `1-alpha_m`; material,
condensed local pressure, Robin spring/dashpot, and active tension are staged
at `1-alpha_f`; load histories are evaluated at
`t[n+1]-alpha_f*dt`; and the viscous Green--Lagrange rate is computed from the
stage velocity. These semantics are application-owned and do not enlarge
CoupFE Core's time-integrator API.

In addition to the shared mesh/time/output options,
`--element-evaluation {joint,split}` selects callback work:

- ``joint`` (default) constructs and caches the paired element blocks;
- ``split`` evaluates material residual and tangent callbacks separately and
  uses Core's compiled-material residual-only entry; on the historical Q1/P0
  contract it also defers the condensed-pressure tangent until requested.

The archive and `solver_configuration_json` record the selected mode and
whether the compiled material actually exposed the residual-only entry. Split
and joint are numerical execution modes, not different formulations or
tolerances. They also record the MPI implementation ID, rank count, balanced
cell partition, SuperLU_DIST factor solver, and—on the closed contract—the
owned-row ranges, local consistent-mass nonzero counts, and touching-element
counts.

A shortened closed rank-gate command is:

```bash
mpiexec -n 4 python examples/cardiac_benchmark/run_mpi.py \
  --benchmark-step 0 --case B --mesh-topology closed-multiblock \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --formulation std-kappa --mass consistent --material-eta 100 \
  --tbar-laplace tbar_closed_nt2_core20_rad17.npy \
  --fiber-sampling gp-direct --element-evaluation joint \
  --integrator be \
  --dt 0.001 --tend 0.32 --load-horizon 1.0 \
  --out caseB_closed_snap_mpi_rank4.npz
```

Run the configuration serially and at one, two, and four MPI ranks and retain
the equality record across the declared snap window before using MPI for wall
or in-surface refinement:

```bash
python examples/cardiac_benchmark/compare_mpi_rank_gate.py \
  --serial caseB_closed_snap_serial.npz \
  --mpi1 caseB_closed_snap_mpi_rank1.npz \
  --mpi2 caseB_closed_snap_mpi_rank2.npz \
  --mpi4 caseB_closed_snap_mpi_rank4.npz \
  --report caseB_closed_snap_rank_gate.json
```

The comparison utility is hard-pinned to the 2×20×17 closed configuration,
the paper parameters, the 0.001 s step, the exact 1 s load prefix through
0.32 s, and fixed machine-level tolerances. See
[`CASE_B_MPI_RANK_GATE.md`](CASE_B_MPI_RANK_GATE.md). Existing automated
1/2/4-rank gates cover the historical reduced application; focused tests
cover the new closed operators and one-step integration. No completed closed
cross-snap gate or full-resolution scaling claim is bundled yet.

The distinct full-cycle Step 2 Case B development configuration is:

```bash
mpiexec -n 4 python examples/cardiac_benchmark/run_mpi.py \
  --benchmark-step 2 --case B \
  --mesh-topology closed-multiblock \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --formulation std-kappa --mass consistent --material-eta 100 \
  --tbar-laplace tbar_closed_nt2_core20_rad17.npy \
  --fiber-sampling gp-direct --element-evaluation joint \
  --integrator generalized-alpha \
  --dt 0.001 --tend 1.0 --load-horizon 1.0 \
  --out mpi4_ga_full.npz
```

This profile applies active stress and ventricular pressure together and
records `paper-source-matched-full-cycle` only when the canonical Step 2
identity, pointwise `std-kappa` formulation, and time/load horizon are present.
The condensed local-pressure variants remain `diagnostic-noncanonical`. The
bundled comparison report was
derived from a completed four-rank run of this profile, but its application
tree was dirty and the result is development evidence—not validation,
convergence, or rank-independence evidence. The generated NPZ remains
external.

### Source-identity assertions

In a Git checkout, the driver discovers the application and Core revisions and
tree states. A clean source archive or another non-Git context may assert them:

```bash
export COUPFE_CARDIAC_APP_REVISION=<40-hex-public-app-commit>
export COUPFE_CARDIAC_TREE_STATE=clean
export COUPFE_CORE_REVISION=<40-hex-public-core-commit>
export COUPFE_CORE_TREE_STATE=clean
```

These variables are assertions, not remote-reachability checks. Do not mark an
edited source tree clean.

The current dependency/runtime guard requires Core
`e2f42ed5772850a0a23a2ce434f430c287eae5c8`. Historical retained JSON/SVG
evidence continues to identify Core
`454f73ce2de284262b214a2b37bd676c6aca3c0a`; it is not metadata-relabeled, and
the compiled rerun gates check current-runtime compatibility separately.

## Result archive fields

Current runs use `result_schema=coupfe-cardiac-result-v1`. NPZ arrays follow
NumPy conventions. A historical v1 archive may omit fields introduced after it
was produced; the field values and source revision determine the actual
contract for that file.

### Completion, histories, and peak field

| Field | Shape or type | Meaning |
|---|---|---|
| `result_schema`, `driver` | scalar text | Schema name and source driver path. |
| `converged` | scalar bool | Written as true only by the completed-result path. |
| `completed_steps`, `expected_steps` | scalar integer | Must agree; each also equals `len(times)-1`. |
| `times`, `tau`, `pres` | `(n_time,)` float | Time, active tension, and pressure histories in seconds and pascals. |
| `load_horizon` | scalar float | Horizon used to generate the load before taking the exact prefix through `t_end`; current shortened gates use `1.0 s`. Historical archives without this field are reported as implicitly using `t_end`. |
| `u0`, `u1` | `(n_time,3)` float | Displacement histories at benchmark material points in metres. |
| `p0`, `p1` | `(3,)` float | Reference coordinates of the output points in metres. |
| `n_peak` | scalar integer | Index of the largest-magnitude case load. |
| `U_peak` | `(ndof,)` float | Global displacement vector at `n_peak`. |
| `nodes`, `elems` | `(n_node,3)`, `(n_elem,8)` | Reference coordinates and Hex8 connectivity. |
| `fiber` | `(n_elem,3)` | Element-level fiber direction used by the diagnostic output. |
| `facets_endo` | `(n_facet,4)` | Endocardial quadrilateral facets. |

### Model and discretization

| Field | Meaning |
|---|---|
| `case`, `dt`, `t_end`, `load_horizon`, `integrator` | Selected benchmark case, time discretization, and load-generation horizon. |
| `benchmark_step`, `benchmark_configuration_id`, `benchmark_load_contract`, `benchmark_active_stress_enabled`, `benchmark_pressure_enabled` | Fail-closed physical identity. In particular, Step 2 Case B is distinct from the historical Step 0 pressure-only Case B label. |
| `benchmark_material_parameters_json`, `benchmark_activation_parameters_json`, `benchmark_pressure_parameters_json` | Exact parameter dictionaries selected by the benchmark-step/case contract. |
| `generalized_alpha_alpha_m`, `generalized_alpha_alpha_f`, `generalized_alpha_gamma`, `generalized_alpha_beta`, `generalized_alpha_stage_contract` | Fixed source parameters and stage-contract identity for a generalized-alpha archive; zero/not-applicable for BE. |
| `load_evaluation_times_s` | Per-record time at which the applied load ODE was sampled. It equals endpoint time for BE. Generalized alpha retains the source ODE initial condition at `t=0` and fixed load horizon, then samples at `t-alpha_f*dt` after the initial record. |
| `mesh_topology`, `n_t`, `n_mu`, `n_theta`, `n_core`, `n_radial`, `core_half_width`, `apex_offset`, `perturb` | Mesh construction identity. `apex_offset` is zero in a closed-multiblock archive and is not a claim that the polar mesh was collapsed. |
| `flip_helix`, `fiber_sampling`, `fiber_sampling_option`, `fiber_direction_reconstruction`, `viscous_rate` | Material-direction and rate policies, including both the stored sampling method and its CLI selection. Closed meshes record `toolkit-physical-coordinate-u-v-v1`; the historical polar topology records `historical-parametric-mu-theta-v1`. |
| `density`, `a_top`, `b_top`, `a_epi`, `b_epi` | Density and Robin parameters in SI units. |
| `mass_representation`, `material_model_id`, `material_eta_pa_s`, `viscous_term_active`, `parameter_variant`, `isotropic` | Inertia choice, constitutive-law identity, viscosity, named parameter variant, and anisotropy-control choices. The model ID distinguishes the corrected complete smooth-switch-energy derivative from historical records made before that correction. |
| `method_metadata_origin` | Whether result-defining method metadata was recorded by the archive or narrowly reconstructed from an immutable reviewed source checkpoint. |
| `tbar_definition`, `tbar_source_filename`, `tbar_source_sha256`, `tbar_metadata_filename`, `tbar_metadata_sha256`, `tbar_metadata_schema` | Transmural-coordinate method and, for a presolved Laplace field, portable identities for both the NPY and its validated native sidecar. No machine-local source path is retained. |
| `formulation` | `hex8_fbar`, `hex8_local_pressure_p0_condensed_logj`, `hex8_local_pressure_p0_condensed_mean_logj_paper_j2`, or the standard pointwise-`kappa` Hex8 identifier recorded by the driver. |
| `material_kernel_formulation` | Generated material-kernel choice: `fbar_mechanics` or `standard`. |
| `material_kappa_pa` | Bulk penalty retained in the material kernel; zero for the Q1/P0 path. |
| `local_pressure_bulk_modulus_pa` | Bulk modulus used by the local-pressure operator; zero for F-bar and std-kappa. |
| `local_pressure_volume_law` | Exact local-pressure response identity. The paper-law variant records that it evaluates paper J² response at the reference-volume-weighted geometric mean J; nonlocal formulations record `not-applicable`. |
| `det_f_gauss_peak` | `(n_elem,8)` deformation Jacobians at the eight Hex8 Gauss points at peak load. |
| `element_pressure_peak_pa` | `(n_elem,)` eliminated P0 pressure for local pressure; an empty array for F-bar and std-kappa. This is not a penalty-pressure reconstruction. |
| `element_evaluation_mode` | `joint` or `split` compiled-material callback policy. |
| `compiled_material_residual_only_available` | Boolean reporting whether the generated Core material kernel exposed its R-only batch entry. |
| `pre_solve_audit_json` | JSON record of geometry/facet identity, extended reference-Jacobian checks, Robin symmetry, and, for Case B, the undeformed pressure resultant and moment. A failed audit stops before time integration and no completed archive is written. |

### Output sampling

Current results record `point_sampling=hex8_reference_isoparametric`. For each
of `p0` and `p1`, the archive includes:

- `<point>_sampling_element`: selected zero-based element index;
- `<point>_sampling_natural`: three natural coordinates;
- `<point>_sampling_weights`: eight Hex8 shape weights; and
- `<point>_sampling_reconstruction_error_m`: reference-map reconstruction
  error in metres.

The existing historical Case A record instead identifies its older
`global_delaunay_tetra` sampler. A sampler label describes how that result was
actually produced; it must not be changed without rerunning the calculation.

### Solver and source provenance

| Field | Meaning |
|---|---|
| `nonlinear_solver` | `core-newton`, serial `petsc-snes`, or the guarded companion value `petsc-snes-mpi`. |
| `solver_configuration_json` | Scalar JSON text containing the exact solver policy and available PETSc/version information, including time integrator, acceleration/force stages, material-batch time identity, and source-matched generalized-alpha parameters when selected. |
| `nonlinear_step_diagnostics_json` | Scalar JSON text containing one diagnostic object per completed increment. |
| `mpi_enabled`, `mpi_ranks`, `mpi_world_size`, `mpi_implementation` | Present for the MPI companion; identify distributed execution, rank count, and the historical or closed contract. |
| `mpi_local_element_counts`, `mpi_partition`, `mpi_build_layout`, `mpi_factor_solver_type` | Balanced cell counts and declared partition/build/factor-solver policies. |
| `mpi_mass_partition`, `mpi_mass_owned_row_ranges`, `mpi_mass_local_nnz`, `mpi_mass_touching_element_counts` | Closed consistent-mass ownership evidence; postprocessing requires contiguous complete global row coverage and nonempty rank contributions. |
| `app_revision`, `app_tree_state`, `app_source_kind` | Application source identity. |
| `core_revision`, `core_tree_state`, `core_source_kind`, `core_source_url` | Core source identity. |
| `benchmark_runtime_source_manifest_json`, `benchmark_runtime_source_sha256` | Canonical path-to-SHA-256 manifest and aggregate digest for every application-owned result-producing source, required to content-identify an intentional dirty-tree Step 2 run. |
| `python_version`, `numpy_version`, `scipy_version`, `coupfe_version` | Runtime package/version context. |

The NPZ is a working result artifact, not the public retained-record format.
Public examples keep normalized text/JSON evidence and hashes while leaving
large generated NPZ files unbundled.

## Mesh construction and pre-solve identity gate

`geometry.py` keeps the historical and benchmark-target topologies separate:

```python
historical = build_mesh(
    n_t=2, n_mu=36, n_theta=48, apex_offset=0.2
)
closed = build_closed_mesh(
    n_t=2, n_core=20, n_radial=17, core_half_width=0.36
)
```

`historical.topology == "polar_ring"`; a positive `apex_offset` means a
truncated domain with an additional free tip boundary. It is retained for
history and demonstrations. `closed.topology == "closed_multiblock_disk"`;
the five-block disk puts the apex at one ordinary surface vertex and has only
the endocardial, epicardial, and base boundary classes.

`run.py` calls `audit_geometry` before kernel compilation, then calls
`audit_robin` and, for Case B, `audit_pressure` before time integration. The gate
checks every exterior-face owner and label, Jacobians at Gauss and extended
natural-coordinate samples, reference measures and extents, Robin symmetry,
and the undeformed unit-pressure resultant and moment. The closed topology is
also compared with retained reference-mesh volume and surface measures under
the declared meshing tolerance. Any failure raises and prevents completed
result writing. Passing establishes setup identity at the recorded resolution;
it does not establish nonlinear convergence or agreement with the benchmark
curves.

The pressure target is topology-aware. `closed_multiblock_disk` uses the
signed analytic projected-base resultant. The historical truncated `polar_ring`
mesh uses an independent vector-area closure of its actual polygonal base and
traction-free terminal rings; the closed analytic projection is retained in
the audit record as a comparator, not used as the truncated mesh's pass/fail
target. Both policies reject a reversed load orientation.

## `build_group` example API

`run.py` exposes this helper for tests and application experiments:

```python
grp, elem, ta_idx, local_pressure = build_group(
    mesh,
    dt,
    build_dir=None,
    formulation="fbar",
    evaluation_mode="joint",
    fiber_sampling="cg1",
    fiber_tbar=None,
    fiber_asign=-1.0,
)
```

The return values are:

- `grp`: Core `ElementGroup` for the compiled cardiac material element;
- `elem`: the stateful `CompiledElement` instance;
- `ta_idx`: integer index of active tension in `elem.props`; and
- `local_pressure`: a `LocalPressureHex8Operator` only for the Q1/P0 path;
  `None` for F-bar and standard pointwise-`kappa`.

The helper initializes fiber-related state at all eight Gauss points and caches
the generated module by kernel formulation for the current Python process.
`fiber_sampling` selects nodal-CG1 interpolation or direct rule evaluation;
`fiber_tbar` supplies the per-node Laplace coordinate used by the direct rule;
and `fiber_asign` records the helix convention. For local pressure the helper
copies the material properties, sets material `kappa` to zero, and returns the
separate volumetric operator so the caller can compose it exactly once.
`evaluation_mode` is passed directly to Core's native
`ElementGroup`; unsupported names and a requested `split` mode without a
residual-only entry fail explicitly.

This is an internal example helper, not a stable package API. Its tuple layout,
cache, and property indexing may change. Prefer the CLI for reproducible runs.

## Hex8 sampling functions

The source module `sampling.py` provides application-owned NumPy helpers.

### Shape functions

```python
weights = hex8_shape(natural_coordinates)             # (8,)
derivatives = hex8_shape_derivatives(natural_coordinates)  # (8,3)
```

Natural coordinates must be a finite length-three vector. Node ordering is the
standard CoupFE/Abaqus Hex8 order recorded by
`HEX8_NATURAL_COORDINATES`.

### Candidate search and checked location

```python
indices = candidate_hex8_elements(
    nodes, elements, point, bbox_tolerance=None
)

location = locate_hex8_point(
    nodes,
    elements,
    point,
    bbox_tolerance=None,
    inside_tolerance=1e-9,
    reconstruction_atol=0.0,
    reconstruction_rtol=1e-11,
    jacobian_rtol=1e-12,
    max_iterations=30,
)
```

`candidate_hex8_elements` returns sorted element indices whose reference
axis-aligned bounding boxes include the point. It is only a search filter.
`locate_hex8_point` applies a damped inverse-isoparametric Newton solve and
validates the reference Jacobian, natural-coordinate bounds, shape weights,
and reconstruction. If more than one candidate is valid, it returns the lowest
element index deterministically.

The immutable `Hex8PointLocation` result contains:

```text
element_index
node_ids                  # 8 integers
natural_coordinates       # 3 floats
weights                   # 8 floats
reconstruction_error
iterations
```

Bad array shapes or tolerances raise `ValueError`. A point outside the mesh,
degenerate/orientation-reversing candidates, or inverse-map failure raises
`RuntimeError`; the function does not silently extrapolate.

### Displacement interpolation

```python
u = interpolate_displacement(displacement, location, dof_per_node=3)
```

`displacement` may be a flat global vector or an `(n_node,dof_per_node)` array.
The function interpolates the first three Cartesian components at the eight
selected nodes and returns a finite `(3,)` vector. `location` must be a checked
`Hex8PointLocation`.

## `LocalPressureHex8Operator`

Construction:

```python
operator = LocalPressureHex8Operator(
    nodes,
    elements,
    ndof,
    bulk_modulus=K,
    dof_per_node=3,
)
```

Inputs must describe at least one finite, orientation-preserving Hex8 with
eight-node integer connectivity, exactly three displacement DOFs per node,
`ndof == 3*len(nodes)`, and a finite positive bulk modulus. Reference
Jacobians must be positive at all eight 2x2x2 Gauss points. The collapsed-apex
mesh therefore is not accepted for this operator.

For each element,

```text
p_e = K / V_e * integral(log(det(F))) dV
P_vol = p_e F^(-T)
```

The main methods are:

| Method | Return or effect |
|---|---|
| `element_pressure(U)` | `(n_elem,)` copy of the algebraically eliminated pressures. |
| `deformation_jacobians(U)` | `(n_elem,8)` copy of Gauss-point `det(F)`. |
| `max_step(U,dU,t,dt=None)` | Conservative factor in `(0,1]` for a determinant-valid Core Newton trial; the implementation repeatedly halves the proposed correction. |
| `residual(U,state,t,dt)` | Core `Residual` contribution for the volumetric stress. |
| `tangent(U,state,t,dt)` | Core `Tangent` containing the condensed displacement derivative, including `dp_e/dU`. |
| `commit(U,state,t,dt)` | Returns `state` unchanged; eliminated pressure has no committed history. |

An inverted, singular, or non-finite trial deformation raises the dedicated
`InvalidDeformationError`. Core consumes `max_step` before its residual
backtracking, while the application PETSc adapter handles that exact exception
as described below. Bad inputs and other numerical or programming exceptions
are not converted into recoverable trials. Non-finite accepted residual or
tangent data fail closed. The operator supplies only the volumetric
contribution. Compose it with the standard cardiac material kernel with
material `kappa=0`; do not compose it with an active bulk penalty or count the
volumetric stress twice.

This class is an application experiment with focused component checks. It is
not a general stable Core interface or, by itself, evidence of ventricular
mesh convergence.

## Nonlinear solver interfaces

### Guarded Core Newton

```python
U, committed_states, iterations = checked_newton_solve(
    operators,
    U0,
    state,
    ndof,
    dirichlet,
    t=1.0,
    dt=1.0,
    rtol=1e-9,
    maxit=60,
)
```

The wrapper computes the initial free-residual norm and installs a commit guard
ahead of the physical operators. It rejects a non-finite or insufficiently
reduced final residual before physical state commit. The CLI currently calls
this path with `rtol=1e-8` and `maxit=40` and records those values. When the
Q1/P0 operator is present, Core's existing optional-operator hook calls its
`max_step` method to bound each proposed correction to a determinant-valid
trial before Core's ordinary residual backtracking.

### PETSc SNES

Default `PetscSnesSettings`:

| Field | Value |
|---|---:|
| `snes_type` | `newtonls` |
| `line_search_type` | `bt` |
| `ksp_type` | `preonly` |
| `pc_type` | `lu` |
| `rtol` | `1e-9` |
| `atol` | `1e-10` |
| `stol` | `1e-12` |
| `max_it` | `60` |

Usage:

```python
solver = PetscSnesSolver()
try:
    U, committed_states, diagnostics = solver.solve(
        operators, U0, state, ndof, {}, t=time, dt=dt
    )
    configuration = solver.configuration()
finally:
    solver.close()
```

The solver is serial (`PETSc.COMM_SELF`), accepts no nonempty Dirichlet map,
creates PETSc objects lazily, and reuses them only for the same `ndof`. Use one
instance per driver run. `close()` is safe before or after initialization.

Acceptance requires a finite displacement, positive SNES reason, acceptable
KSP reason, finite diagnostic values, and independently reassembled
`final_residual_norm <= max(atol, rtol*initial_residual_norm)` before state
commit. `SnesSolveError` carries `diagnostics` when available.

In the residual callback, and only there, `PetscSnesSolver` catches
`InvalidDeformationError`, writes IEEE positive infinity to every entry of the
trial residual, and lets SNES `newtonls`/`bt` shorten that trial. This is the
portable application policy used with supported petsc4py versions; it does not
claim that petsc4py exposed a native function-domain-error call. Each solve
counts these trial rejections and keeps the last exception detail. The
Jacobian callback does not catch the exception, and unrelated exceptions in
either callback propagate and abort before state commit. An invalid initial,
accepted, or independently reassembled final state also fails closed. If
backtracking cannot find an acceptable valid state, no state is committed and
no completed archive is written. None of these controls changes `rtol`,
`atol`, `stol`, `max_it`, or the final residual rule.

`configuration()` is JSON-compatible. Before initialization it contains the
declared settings, `function_domain_rejection_api`, and available `petsc4py`
version; after initialization it can also report PETSc version, the line-search
configuration API, and factor solver type.

`SnesStepDiagnostics.as_dict()` contains:

```text
time, dt
initial_residual_norm, final_residual_norm
residual_acceptance_threshold, petsc_function_norm
snes_converged_reason, ksp_converged_reason
nonlinear_iterations, linear_iterations
residual_history
assembly_seconds, solve_seconds
function_domain_rejections, last_function_domain_error
```

`function_domain_rejections` is zero when SNES did not evaluate an invalid
local-pressure residual trial; otherwise `last_function_domain_error` contains
the final recorded detail. These fields describe rejected trials, not accepted
invalid states. Historical v1 results may omit fields introduced after their
source revision. Timing fields describe that process execution and are not a
scaling claim.

## External comparison CLIs and JSON

For Step 0, `post.py` consumes a completed NPZ plus trusted external curves:

```bash
python examples/cardiac_benchmark/post.py RESULT.npz \
  [--case {step_0A,step_0B}] \
  [--reference-dir DIR] \
  [--json REPORT.json] \
  [--run-log STDOUT.txt] \
  [--supersedes-report-sha256 SHA256] \
  [--plot FIGURE.png]
```

- `--case` is an optional cross-check; the result case otherwise selects the
  reference step.
- `--reference-dir` accepts the curve-data directory or an extracted archive
  root. `CARDIAC_BENCHMARK_DATA_DIR` is the fallback.
- `--json` writes an atomic, strict-finite JSON report.
- `--run-log` normalizes a UTF-8 transcript to LF with a final LF when nonempty,
  then records its basename, size, hash, and normalization policy. It does not
  embed or redistribute the log.
- `--supersedes-report-sha256` records correction lineage to an exact predecessor
  report when regenerating corrected evidence.
- `--plot` writes an optional PNG and imports Matplotlib only for that request.

The utility validates completion, time grids and histories, formulation and
sampling labels, source identity, current Hex8 sampling metadata, solver
configuration and per-step diagnostics, positive peak Gauss-point `det(F)`,
and local-pressure metadata before reading reference pickles. It accepts the
older Delaunay label for backward compatibility with historical archives; the
current driver writes the Hex8 label and metadata.

The loader requires the exact 10 Case A/Case B files selected by upstream
`results_time_curves/figures.py`; missing selected files and unexpected matching
files fail closed. The archive's unselected base-name SimVascular file is
excluded only after a byte-for-byte comparison with selected SimVascular P2.
Each selected file must parse. Reference histories and the result history are
resampled to the canonical 101-point grid from 0 to 1 s. The report records the
selection policy, upstream `figures.py` identity, selected-file identities, and
the excluded alias identity and reason. A caller must still verify its own 23.2
GB download before treating the trusted-pickle precondition as satisfied.

The JSON report uses
`schema=coupfe-cardiac-reference-comparison-v2` and has these top-level fields:

- `bounded_claim`: the scope statement for the generated report;
- `result`: NPZ identity, configuration, source identity, solver record,
  histories, peak data, deformation-Jacobian summary, optional element-pressure
  summary, Case B peak circumferential-ring rotation profile when complete
  structured rings are present, and optional normalized-log identity;
- `reference`: DOI, license, canonical grid, exact 10-file selection provenance,
  excluded-alias provenance, and all-team mean curves;
- `comparison`: RED values for this result and each team plus the result curves
  on the canonical grid; and
- optional `correction`: exact predecessor report SHA-256, repository revision,
  and correction reason.

The console rounds RED for readability; use the JSON values for retained
full-precision evidence. RED is a measurement, not an automatic validation
threshold. `main(argv)` also returns the report dictionary when called from
Python, but that convenience is not promised as a stable package API.

Concrete Step 0 v2 reports under `examples/cardiac_benchmark/results/` retain the full
result histories and accepted-step diagnostics, the derived all-team mean
curves, exact reference-file identities, and team/result RED values. The raw
CC BY 4.0 team pickle files and generated CoupFE NPZ archives remain external.
The report is therefore sufficient to review the recorded comparison, while a
new comparison run still requires a regenerated NPZ and a caller-verified copy
of the external dataset. The separately retained console file is bound to the
report by `result.normalized_run_log` size, SHA-256, and normalization policy.

### Step 2 Case B publisher comparator and renderer

Step 2 Case B has a separate fail-closed comparator because its active-stress-
plus-pressure identity, source-matched generalized-alpha contract, publisher
selection, and phase/branch diagnostics differ from the Step 0 report schema:

```bash
python examples/cardiac_benchmark/compare_step2b_case_b.py \
  --publisher-data-dir /path/to/results_time_curves/data \
  --coupfe-run mpi4_ga_full.npz \
  --output step2_case_b_std_kappa_2x20x17_dt0p001.report.json
```

The default `--hash-manifest` is the reviewed
`examples/cardiac_benchmark/step2b_case_b_reference_hashes.json`. The
comparator accepts only a completed 1,000-step, 1 ms Step 2 Case B archive
with the exact physical/load identity, closed-mesh audits, source-matched
generalized-alpha metadata, Core revision, runtime-source content manifest,
and one accepted nonlinear diagnostic per increment. It verifies the exact
ten publisher files before restricted pickle decoding, then writes schema
`coupfe-cardiac-step2b-publisher-comparison-v1` with full-history/phase errors,
paper Eq. 21 RED, event/branch diagnostics, and derived 101-point curves.
For the retained dirty-tree run, the exact source file map is separately
byte-pinned in
`examples/cardiac_benchmark/step2b_case_b_runtime_source_hashes.json` and
points to repository snapshot `d06c3e9`. The comparator accepts either that
exact reviewed identity or the exact current runtime-source identity for a new
run; it never substitutes current hashes into the retained result.

The checked-in
[`full-cycle development report`](../examples/cardiac_benchmark/results/step2_case_b_std_kappa_2x20x17_dt0p001.report.json)
and its [normalized stdout](../examples/cardiac_benchmark/results/step2_case_b_std_kappa_2x20x17_dt0p001.raw.stdout.txt)
bind an external four-rank NPZ by name, size, and SHA-256. The report records a
completed full cycle, 9.8038% global relative L2 and 28.3004%/35.2774% paper
RED at `p0`/`p1`, while retaining the wrong-sign `p1-z` plateau as an explicit
discrepancy. Its application revision reports a dirty tree whose exact runtime
source contents are hash-identified. These are quantified development results;
the comparator defines no acceptance threshold and makes no validation,
convergence, or rank-independence claim.

Because the report embeds the reviewed CoupFE, ten-team range/mean, and named
Simula curves, the public figure does not require the external NPZ or pickles:

```bash
python examples/cardiac_benchmark/plot_step2b_case_b.py --canonical
```

The plotter validates the report byte identity and run contract before
rendering `docs/figures/step2_case_b_comparison.svg`. `--canonical` additionally
requires the exact renderer stack used for the checked-in SVG; without it,
supported Matplotlib stacks can change SVG text coordinates or IDs without
changing the report data. Recomputing the comparison itself still requires
the external NPZ and verified publisher inputs.

## Peak-field diagnostic

Run:

```bash
python examples/cardiac_benchmark/diagnose.py RESULT.npz
```

The importable interface separates calculation from printing:

```python
diagnostics = analyze_result(path)
diagnostics = main(path)  # prints the report and returns the same mapping
```

With no CLI argument the script looks for `caseB_full.npz` in the current
directory. It validates the completed history, connectivity, field shapes, and
finite values. `analyze_result` returns the peak index and time, element fiber
stretch, the selected `det(F)` field and its sampling label, reconstructed
centroid `det(F)`, endocardial radial displacement, and local-pressure values
when present.

Current archives retain the eight-Gauss-point `det_f_gauss_peak` array, which
the diagnostic prefers when auditing element validity because a centroid
summary cannot bound all Gauss points. Older archives fall back to a labeled
centroid reconstruction. For the Q1/P0 path the printed report also summarizes
the retained condensed element pressure.

`diagnose.py` is a reporting convenience. Its finite values and displacement
signs are not a paper-curve oracle, a cavity-volume calculation, or a clinical
interpretation.
