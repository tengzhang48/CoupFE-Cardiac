# Distributed cardiac implementation checks

These scripts exercise distributed element assembly, transient Robin terms,
per-element viscous state, follower pressure, source-matched generalized-alpha
staging, and a timing-oriented rank comparison. They are executable
implementation checks, not whole-case cardiac validation or general
performance evidence.

All five retained checks use:

- application revision `44cbfed9e09d4150203faae3087f2e4617d1fc47`;
- Core revision `454f73ce2de284262b214a2b37bd676c6aca3c0a`;
- the noncollapsed open-apex mesh (`apex_offset=0.2`);
- softened volumetric penalty `kappa=1e3`; and
- the scripts' `1e-7` agreement gate.

## Retained release-configuration output

The per-script numbers below were independently reproduced from the exact
release revisions after the aggregate 5-test release gate passed.

| Script | Problem size / steps | Reference | Maximum observed difference at ranks 1/2/4 | Largest final reported residual at ranks 1/2/4 | Status |
|---|---|---|---:|---:|---|
| `distributed_cardiac_passive.py` | 1,872 DOF, 384 elements, 3 increments | serial | `4.869e-17` | `5.44e-14` | passed |
| `distributed_cardiac_dynamics.py` | 1,872 DOF, 384 elements, 4 steps | serial | `9.758e-19` | `2.76e-13` | passed |
| `distributed_cardiac_viscous.py` | 1,872 DOF, 384 elements, 5 steps | serial | `8.222e-14` | `7.70e-13` | passed |
| `distributed_cardiac_pressure.py` | 1,872 DOF, 384 elements, 192 pressure facets, 4 steps | serial | `1.375e-13` | `5.81e-15` | passed |
| `distributed_cardiac_scaling.py` | 648 DOF, 96 elements, 4 steps | distributed rank 1 | `6.47e-19` (ranks 2/4) | `4.43e-14` | passed |

The release gate ran every script at 1, 2, and 4 ranks:

```text
.....                                                                    [100%]
5 passed, 31 deselected in 110.32s (0:01:50)
```

The scripts auto-selected SuperLU_DIST in the retained environment. Wall times
are intentionally not reported because this configuration is not a scaling
benchmark. Approved Core exposes cumulative linear-solver divergence and the
final reported nonlinear residual, but not an all-step nonlinear status; the
claim remains limited to the exact scripts and ranks above.

## Source-matched generalized-alpha implementation gate

`distributed_cardiac_generalized_alpha.py` is a newer, application-owned gate
and is not part of the five-result historical table above. Its two-step linear
reference isolates the exact Simula/FEniCS parameters and stages: consistent
mass at `alpha_m`, material/Robin terms at `alpha_f`, and committed `u/v/a`
history. It also partitions one load element per scalar DOF independently from
PETSc row ownership, so remote rank-local insertion is exercised. The current
1/2/4-rank check produced maximum endpoint/history errors of `3.886e-15`,
`2.220e-15`, and `1.998e-15`, respectively. This is an assembly and state gate;
it does not establish ventricular agreement with the FEniCS benchmark.

## Run

With mutually compatible MPI, PETSc, `mpi4py`, and `petsc4py` installations:

```bash
python -m pip install -e ".[dev,mpi]"
OMP_NUM_THREADS=1 python -m pytest -q -m mpi
```

See [`../REFERENCES.md`](../REFERENCES.md) for provenance and the boundary
between implementation agreement and scientific validation.
