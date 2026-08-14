# Controlled benchmark runs

How to launch a controlled Step 0 Case B run that the fail-closed
provenance and audit machinery will accept, and that stays a one-variable
comparison against the retained clean gate. This document is operational:
it records the runtime, the exact command shape, and the failure modes
that have already cost runs. It does not reinterpret results; see
`CASE_B_MESH_ERROR_LAYERS.md` and `BENCHMARK_REPRODUCTION_STATUS.md` for
evidence and claims.

## Runtime environment

- MPI runs require an interpreter providing `mpi4py` and `petsc4py`
  **3.18.4**, the version used by the retained clean gate and gate runs.
  The dependency solver may install a different petsc4py; check
  `petsc4py.PETSc.Sys.getVersion()` before launching. A version mismatch
  is recorded in the archive's solver configuration and silently changes
  linear-solver behavior.
- The CoupFE Core revision is pinned in `pyproject.toml`; verify it with
  `python .github/scripts/check_runtime_core.py`. The resolved Core is
  recorded in every result archive.
- Runs must start from a **clean, committed source tree**. The result
  archive records the application/Core revisions and tree state, and the
  release guard rejects uncommitted trees.

## Command template

The retained clean Step 0B gate used exactly this shape (8 ranks):

```bash
mpiexec -n 8 python examples/cardiac_benchmark/run_mpi.py \
  --benchmark-step 0 --case B --mesh-topology closed-multiblock \
  --nt 2 --ncore 20 --nradial 17 --core-half-width 0.36 \
  --formulation local-pressure --mass consistent --material-eta 100 \
  --tbar-laplace FIELD.npy \
  --fiber-sampling gp-direct --element-evaluation joint \
  --integrator generalized-alpha \
  --linear-solver-profile fgmres-gamg-rigid-rebuild \
  --dt 0.001 --tend 0.32 --load-horizon 1.0 \
  --out RESULT.npz
```

Hard requirements that differ from driver defaults:

- `--linear-solver-profile fgmres-gamg-rigid-rebuild` — the driver default
  (`direct-superlu-dist`) is roughly 50x slower on this model and is not
  the reference configuration. Copy the profile from the reference run's
  recorded `solver_configuration_json`, not from the parser default.
- `--integrator generalized-alpha --dt 0.001 --load-horizon 1.0` — the
  source-matched scheme. `--tend 0.32` covers the snap window;
  `--tend 1.0` is the full cycle.
- Mesh changes are parameterized by `--nt`, `--ncore`, `--nradial`,
  `--core-half-width`, and the tip-grading control `--tip-refine F`
  (node relocation toward the apex at fixed element count; `F=1.0` is
  the uniform benchmark mesh bit-for-bit). See
  `MESH_REFINEMENT_GUIDE.md` for the three mesh axes.

## Laplace transmural field

The closed benchmark configuration uses a precomputed Laplace transmural
field (`--tbar-laplace`). The field is bound to the exact mesh by content
hash and mesh parameters, and the driver rejects a mismatch fail-closed:

- regenerate the field **after any geometry-affecting change** (mesh
  counts, `core_half_width`, `tip_refine`, geometry code) with
  `python examples/cardiac_benchmark/tbar_laplace.py --nt ... --ncore ...
  --nradial 17 --tip-refine F --out FIELD.npy`, using the same mesh
  parameters as the run;
- never reuse a field across mesh variants: the gate is designed to
  reject that, and bypassing it silently changes the fiber field;
- a `tip_refine` mismatch between field and mesh is a launch error, not
  a warning.

## Acceptance and budgeting

- Pre-solve audits (geometry/topology/labels/Jacobians/measures, pressure
  resultant, Robin symmetry) are fail-closed on rank 0. A launch that
  fails an audit has never produced a result archive; treat the audit
  message as the defect report.
- Result archives are written **only after every requested step
  completes**. A killed or timed-out run leaves no archive. Size wall
  clocks for the machine as loaded, not as idle: check `uptime`/load
  first, do not rely on default timeouts for multi-hour MPI runs, and
  avoid running several heavy jobs concurrently on shared hosts.
- Record the runtime (interpreter, MPI, PETSc) and machine load with the
  run; they are experiment metadata.

## One-variable discipline

A controlled run changes exactly one input from its reference (mesh axis,
`tip_refine`, or operator option) and keeps the command template, load,
time integration, solver profile, rank count, and field policy identical.
If two things differ, the result is not assignable to either.

## Observed launch-failure modes (all fail-closed, all cheap)

1. `petsc4py` missing or wrong version in the invoking interpreter.
2. Mesh/field mismatch from a stale or mis-graded `--tbar-laplace` file.
3. Driver-default solver profile (`direct-superlu-dist`) instead of the
   reference profile.
4. Pre-solve audit rejection of a genuinely defective mesh (this caught
   a seam deduplication defect at 64x72 before any solve; see
   `lessons_learned.md`).
5. Wall-clock timeouts on an oversubscribed host; archives exist only on
   completion.
