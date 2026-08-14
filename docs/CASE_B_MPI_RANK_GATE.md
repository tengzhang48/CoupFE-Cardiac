# Closed Case B MPI rank-equivalence gate

This gate checks that the closed Case B trajectory produced by the validated
MPI companion agrees with the serial PETSc-SNES trajectory at 1, 2, and 4 MPI
ranks. It is a numerical-equivalence check for one fixed method and solve
prefix. It is not a performance benchmark, mesh-convergence study, or claim of
agreement with the published reference solution.

The retained gate is intentionally not configurable. It accepts only:

- Case B with exactly zero activation and zero mesh perturbation;
- the closed five-block mesh at `n_t=2`, `n_core=20`, `n_radial=17`, and
  `core_half_width=0.36` (5,403 nodes, 3,520 Hex8 elements, 16,209 displacement
  degrees of freedom);
- standard Hex8 pointwise `kappa=1e6 Pa`, the corrected smooth-switch material
  derivative, and `eta=100 Pa s`;
- assembled consistent Q1-Hex8 mass;
- GP-direct fibers using a pre-solved Laplace transmural field;
- backward Euler with `dt=0.001 s`, `t_end=0.32 s`, and a fixed `1.0 s` load
  integration horizon;
- serial `petsc-snes` and the closed Case B MPI implementation
  `cardiac-owned-distributed-closed-std-kappa-step0` at exactly 1, 2, and 4
  ranks.

The 0.32 s endpoint covers the pressure-driven snap interval while keeping the
rank gate distinct from the later full-duration production comparison.

## Generate the four inputs

Run all four archives from the same clean public CoupFE-Cardiac revision and
the same clean public CoupFE Core revision. The Laplace `tbar` field and its
sibling metadata must also be identical. Use separate output and build
directories for every run.

The serial command has this configuration:

```bash
python examples/cardiac_benchmark/run.py \
  --case B \
  --mesh-topology closed-multiblock \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --formulation std-kappa \
  --mass consistent \
  --material-eta 100 \
  --tbar-laplace tbar_laplace_closed_multiblock_nt2_core20_rad17.npy \
  --fiber-sampling gp-direct \
  --integrator be \
  --nonlinear-solver petsc-snes \
  --element-evaluation joint \
  --apex-offset 0 --perturb 0 \
  --dt 0.001 --tend 0.32 --load-horizon 1.0 \
  --build-dir build/serial \
  --out run/serial.npz
```

Use the same physical and discretization arguments for each MPI run:

```bash
mpiexec -n 1 python examples/cardiac_benchmark/run_mpi.py \
  --case B --mesh-topology closed-multiblock \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --formulation std-kappa --mass consistent --material-eta 100 \
  --tbar-laplace tbar_laplace_closed_multiblock_nt2_core20_rad17.npy \
  --fiber-sampling gp-direct --element-evaluation joint \
  --apex-offset 0 --perturb 0 \
  --dt 0.001 --tend 0.32 --load-horizon 1.0 \
  --build-dir build/mpi1 --out run/mpi1.npz
```

Repeat with `mpiexec -n 2` and `mpiexec -n 4`, distinct build directories,
and distinct output archives. Set `PYTHONPATH` or use an installed Core build
as required by the local environment; those machine paths are execution
details and do not belong in the retained record.

## Create the retained gate record

```bash
python examples/cardiac_benchmark/compare_mpi_rank_gate.py \
  --serial run/serial.npz \
  --mpi1 run/mpi1.npz \
  --mpi2 run/mpi2.npz \
  --mpi4 run/mpi4.npz \
  --report results/case_b_closed_rank_gate.report.json
```

There are no mesh or tolerance options. A smaller smoke run cannot be labeled
as this retained gate, and a mismatch cannot be made to pass by changing a
command-line threshold.

## What is checked

Before comparing displacements, the tool requires every archive to be complete
and converged for all 320 steps. It then checks:

- clean, 40-character Git revisions for the application and Core, identical
  across all four inputs;
- the corrected material model identifier and every fixed physical and
  discretization setting listed above;
- the exact mesh counts, closed-geometry identity, and benchmark landmarks;
- identical SHA-256 identities for the Laplace field and its metadata;
- identical passed geometry, pressure, and Robin pre-solve audits;
- the serial PETSc-SNES execution identity;
- MPI ranks of exactly 1, 2, and 4, the validated implementation identifier,
  balanced element partitions, and complete owned-row consistent-mass
  provenance;
- a positive SNES and KSP convergence reason and independent residual
  acceptance at every step;
- an empty `element_pressure_peak_pa` array, because pointwise `kappa` rather
  than the condensed local-pressure formulation is used.

The following arrays must be bit-for-bit identical:

- `times`, `tau`, `pres`, `nodes`, `elems`, `p0`, and `p1`;
- retained fiber, endocardial-facet, and landmark-sampling arrays;
- all other common integer, boolean, and text fields that do not identify the
  intentionally different serial/MPI execution paths.

The serial result is the reference for `u0`, `u1`, `U_peak`, and
`det_f_gauss_peak`. Each value must satisfy

```text
abs(mpi - serial) <= 2e-13 + 2e-11 * abs(serial)
```

Those tolerances come from the historical CoupFE serial/MPI equivalence gate.
The report retains, for every field and rank, the maximum absolute difference,
the maximum difference relative to the maximum serial magnitude, the maximum
fraction of the elementwise tolerance, and the pass result.

## Failure and retention behavior

Any missing field, unsafe object array, dirty or mismatched source identity,
configuration mismatch, incomplete mass provenance, failed audit, unaccepted
step, exact-invariant difference, or numerical-tolerance violation stops the
tool before it writes or replaces the JSON report.

A passing report contains only portable input basenames, byte counts, SHA-256
identities, exact source revisions, checked setup, completion and solver
summaries, comparison definitions, and numerical differences. Caller paths,
build directories, and raw diagnostic records are not copied into it.

The focused broken-control and path-leak checks are in
`tests/test_cardiac_mpi_rank_gate.py`.
