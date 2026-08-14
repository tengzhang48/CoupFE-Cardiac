# Case B landmark comparison with the local FEniCS output

`examples/cardiac_benchmark/compare_fenics_case_b.py` creates a public JSON
record and a six-panel displacement overlay for one completed CoupFE Case B
run and the retained local FEniCS p0/p1 trajectories. The paths are always
supplied by the caller; the utility has no machine-specific input defaults.

This is a landmark-trajectory comparison between different spatial
discretizations. CoupFE may use the retained backward-Euler contract or the
source-matched generalized-alpha MPI contract; FEniCS uses P2 tetrahedra with
generalized alpha. It is not a stress comparison, a mesh/time-convergence
result, evidence that the discretizations are equivalent, clinical validation,
or broad solver validation.

## Inputs and provenance boundary

The utility reads exactly five files:

- a completed CoupFE Case B NPZ archive;
- the FEniCS `parameters.json`;
- the FEniCS `time_stamps.npy`;
- `componentwise_displacement_up0.npy`; and
- `componentwise_displacement_up1.npy`.

All NumPy files are loaded with `allow_pickle=False`. The comparison does not
read `pressure_model.npy`, pickle, HDF5, retained stress fields, or a solver
executable. Each input is recorded by role, basename, byte count, and SHA-256.
Resolved paths and the `outdir`/`outpath` strings inside `parameters.json` are
never copied into the report or figure metadata. A geometry path is reduced to
its basename for both POSIX and Windows inputs, while every actual input name
must already be a portable basename. The required Core source URL
must be the public `https://github.com/tengzhang48/CoupFE.git` URL; a missing,
credential-bearing, or different URL fails closed.

The local FEniCS files and parameters have exact hashes, but the producing
source revision was not retained contemporaneously. A nearby clean source
clone is useful context; it is not claimed as the producer of these output
files. The report states this limitation directly.

The checked example manifest
`examples/cardiac_benchmark/fenics_case_b_reference_hashes.example.json`
contains the identities of the known local displacement inputs:

| Role | Basename | SHA-256 |
|---|---|---|
| `fenics-parameters` | `parameters.json` | `c1cd4c8d2521fd6c28774975843740a8af12568edd1240f5daa133d469e6fb76` |
| `fenics-times` | `time_stamps.npy` | `ddba330b1c8f8c1bb61282e187047f3aa99d0df37b2c4ed2139ea1b0e0ff0f0c` |
| `fenics-p0` | `componentwise_displacement_up0.npy` | `4344a4f599a6eabb16159682339a735bff572eaa18eedd1fe2a97ebd3ee7f4a0` |
| `fenics-p1` | `componentwise_displacement_up1.npy` | `88a679de2189bc137de5d64186c698f1702e9df20333849377cdfd01aac8bf1e` |

The manifest verifies file identity; it does not supply or infer a filesystem
location. A final retained comparison must use `--retained` and add the
completed CoupFE run hash with `--expect-sha256 coupfe-run=...`. Retained mode
requires expected hashes for all five roles and records every verified role,
digest, basename, and byte count in the report. Without `--retained`, the
output is explicitly labeled as a development comparison and is not a public
retained-result candidate, even if some expected hashes were supplied.
Retained mode establishes the integrity of a named comparison; it does not by
itself choose the repository's canonical public mesh. That selection is bound
separately by the release artifact specification and its exact result/report/
figure identities.

## Run the comparison

Install the optional plotting dependency with `pip install -e ".[reference]"`,
then provide every input and output explicitly:

```bash
FENICS_OUTPUT_DIR=/path/to/fenics-case-b-output
COUPFE_RESULT=/path/to/completed-clean-case-b.npz
python examples/cardiac_benchmark/compare_fenics_case_b.py \
  --retained \
  --coupfe-run "$COUPFE_RESULT" \
  --fenics-parameters "$FENICS_OUTPUT_DIR/parameters.json" \
  --fenics-times "$FENICS_OUTPUT_DIR/time_stamps.npy" \
  --fenics-p0 "$FENICS_OUTPUT_DIR/componentwise_displacement_up0.npy" \
  --fenics-p1 "$FENICS_OUTPUT_DIR/componentwise_displacement_up1.npy" \
  --expected-hashes examples/cardiac_benchmark/fenics_case_b_reference_hashes.example.json \
  --expect-sha256 coupfe-run=REPLACE_WITH_64_HEX_DIGITS \
  --report /path/to/output/case_b_fenics_comparison.json \
  --figure /path/to/output/case_b_fenics_displacements.png
```

The CoupFE archive must identify clean application and Core Git checkouts and
must contain a completed 1000/1000-step fixed-parameter Case B run. The loader
checks the pointwise-kappa material, `eta=100 Pa s`, `kappa=1 MPa`, zero
activation, the exact benchmark pressure history, consistent mass,
GP-direct/Laplace fibers, closed geometry,
pre-solve boundary audits, accepted nonlinear residuals, positive retained
`det(F)`, the benchmark p0/p1 coordinates, and material-model ID
`holzapfel-ogden-smooth-switch-complete-energy-derivative-v1`. A serial
`petsc-snes` result is accepted. A `petsc-snes-mpi` result is accepted only
when its closed-Case-B implementation, rank counts, balanced element
partition, owned-row consistent-mass metadata, PETSc configuration, and
per-step rank records agree. Missing or mixed serial/MPI provenance fails
closed. These are identity and evidence checks, not adjustable acceptance
parameters.

For a generalized-alpha archive, the loader additionally requires explicit
Benchmark 1 Step 0 Case B identity, the fixed `alpha_m=0.2`, `alpha_f=0.4`,
`gamma=0.7`, `beta=0.36` parameters, the dedicated MPI implementation and
material-rate contract, and the pressure history sampled at
`t_np1-alpha_f*dt`. Step 2 Case B remains routed to its separate blinded
publisher comparator.

## Common grid and reported quantities

The common grid is the retained FEniCS output convention: 999 samples from
`0.001 s` through `0.999 s` at `0.001 s`. A required CoupFE run contains the
`0.000 s` and `1.000 s` endpoints as well, so all common samples coincide
exactly and the report records `identity_sampling: true`. The mapping function
supports checked interpolation for unit tests and future compatible inputs,
but endpoint extrapolation is forbidden.

For p0 and p1, the report gives:

- vector RMSE, `sqrt(mean_t(||u_CoupFE-u_FEniCS||²))`, in millimetres;
- relative L2, `||u_CoupFE-u_FEniCS||_F / ||u_FEniCS||_F`;
- component RMSE and maximum vector error;
- full, pre-snap, snap, and post-snap records; and
- snap onset for both trajectories and their time difference.

The snap window is fixed at `0.20–0.32 s`. Onset is fixed as the first downward
crossing of landmark `u_z=-5 mm`, linearly interpolated between adjacent
samples. Neither definition is exposed as a tuning option.

The figure shows FEniCS as a solid blue curve, CoupFE as a dashed orange curve,
and the declared snap window as gold shading. It contains displacement units,
method labels, and the landmark-only/non-equivalence boundary.
