"""Timing-oriented distributed cardiac rank comparison on a finer mesh.

Times ONLY the distributed solve (passive HO bulk + lumped-mass inertia + Robin) at the given
rank count on a finer LV mesh, and saves the gathered U so a separate pass can confirm the
solution agrees with the rank-1 result within 1e-7 for this script and
configuration (the correctness gate must hold before interpreting timing). The
auto policy prefers reproducible SuperLU_DIST, then MUMPS, then GMRES/GAMG;
MUMPS's parallel pivoting is not run-to-run reproducible. Timings are specific to the selected
mesh, solver, machine, and MPI build. A multi-rank invocation fails unless a
rank-1 result already exists in ``CARDIAC_SCALING_DIR`` and the gathered
solution agrees; it does not establish a general scaling threshold or
performance claim.

    OMP_NUM_THREADS=1 <mpiexec> -n {1,2,4} python examples/mpi_smoke/distributed_cardiac_scaling.py [nt nmu nth]
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

import numpy as np
from petsc4py import PETSc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "cardiac_benchmark"))
import geometry as geom                                    # noqa: E402
from solver import require_distributed_success             # noqa: E402

from coupfe.mesh import KernelMeshView                      # noqa: E402
from coupfe.runtime.compiled_element import CompiledElement  # noqa: E402
from coupfe.assembly.distributed import (                   # noqa: E402
    element_partition, partition_elements, solve_dynamics_distributed,
)
from distributed_cardiac_passive import (                   # noqa: E402
    _passive_props, _init_fiber_svars_subset, _kernel, _solver_options,
)
from distributed_cardiac_dynamics import _lumped_mass, _robin, DT

N_STEPS = 4
V0 = 0.01
NT, NMU, NTH = (int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 \
    else (3, 32, 48)
SAVE_DIR = Path(os.environ.get("CARDIAC_SCALING_DIR", ".")).expanduser()


def _atomic_save(path, values):
    temporary = path.with_name(f".{path.name}.tmp")
    with temporary.open("wb") as stream:
        np.save(stream, values)
    os.replace(temporary, path)


def main():
    comm = PETSc.COMM_WORLD
    rank, size = comm.getRank(), comm.getSize()

    mesh = geom.build_mesh(NT, NMU, NTH, flip_helix=True, apex_offset=0.2)
    view = KernelMeshView(mesh.nodes, mesh.elems, dof_per_node=3)
    ndof = view.ndof
    M = _lumped_mass(mesh, ndof)
    v0 = np.zeros(ndof); v0[2::3] = V0

    parts = partition_elements(view, size)
    mine = np.where(parts == rank)[0]
    my_gm, my_coords, _ = element_partition(view, rank, size)
    props, schema, per_gp = _passive_props()
    elem = CompiledElement(_kernel(), props=props, dof_per_node=3, n_svars=per_gp * 8,
                           mcrd=3, n_elem=len(mine), dt=DT, backend="native",
                           state_schema=schema)
    _init_fiber_svars_subset(elem, mesh, schema, per_gp, mine)
    robin = _robin(mesh, ndof)

    comm.barrier()
    t0 = time.time()
    U_par, info = solve_dynamics_distributed(
        ndof, my_gm, my_coords, 3, elem.element_rk_batch, {}, M,
        dt=DT, n_steps=N_STEPS, v0=v0, robin=robin,
        tol=1e-9, max_newton=60, **_solver_options())
    U_par = require_distributed_success(
        U_par, info, tol=1e-9, context="distributed cardiac scaling smoke"
    )
    comm.barrier()
    wall = time.time() - t0

    failed = False
    if rank == 0:
        SAVE_DIR.mkdir(parents=True, exist_ok=True)
        save_path = SAVE_DIR / f"cardiac_scale_U_{size}.npy"
        ref_path = SAVE_DIR / "cardiac_scale_U_1.npy"
        drift = "  rank-1 reference saved"
        error = float("nan")
        if size != 1:
            if not ref_path.is_file():
                failed = True
                drift = "  rank-1 reference missing -> FAIL"
            else:
                U1 = np.load(ref_path)
                if U1.shape != U_par.shape:
                    failed = True
                    drift = (
                        f"  rank-1 shape {U1.shape} != {U_par.shape} -> FAIL"
                    )
                else:
                    error = float(np.max(np.abs(U_par - U1)))
                    failed = error >= 1e-7
                    drift = (
                        f"  max|U_{size} - U_1|={error:.2e}"
                        f" -> {'OK' if not failed else 'FAIL'}"
                    )
        if not failed:
            _atomic_save(save_path, U_par)
        print(f"[size={size}] ndof={ndof} elems={mesh.n_elem} wall={wall:.2f}s "
              f"newton={info['n_newton']} rnorm={info['rnorm']:.2e}{drift}", flush=True)
    comm.barrier()
    if failed:
        if np.isfinite(error):
            raise RuntimeError(f"distributed/rank-1 mismatch: {error:.3e}")
        raise RuntimeError("distributed scaling run lacks a usable rank-1 reference")


if __name__ == "__main__":
    main()
