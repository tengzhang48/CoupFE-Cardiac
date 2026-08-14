"""Distributed dynamics + Robin + follower-pressure serial-reference check.

The endocardial follower pressure (deformation-dependent, paper Eq. 1 Neumann term) is a surface
load: each rank computes the (cheap, O(surface)) full residual/tangent on the gathered U and adds
only its OWNED rows to the distributed system — the surface-replication pattern. Constant cavity
pressure inflates the LV against the Robin springs; a short backward-Euler transient.

Gate: gathered U(t) agrees with serial
`solve_dynamics([grp, inertia, robin, pressure])` within 1e-7 for this script
and configuration at 1/2/4 ranks.

    OMP_NUM_THREADS=1 <mpiexec> -n {1,2,4} python examples/mpi_smoke/distributed_cardiac_pressure.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
from petsc4py import PETSc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "cardiac_benchmark"))
import geometry as geom                                    # noqa: E402
from pressure import FollowerPressureOperator              # noqa: E402
from solver import require_distributed_success, solve_dynamics_checked  # noqa: E402

from coupfe.mesh import KernelMeshView                      # noqa: E402
from coupfe.operators.element_group import ElementGroup     # noqa: E402
from coupfe.operators.inertia import InertiaOperator        # noqa: E402
from coupfe.runtime.compiled_element import CompiledElement  # noqa: E402
from coupfe.assembly.distributed import (                   # noqa: E402
    element_partition, partition_elements, solve_dynamics_distributed,
)
from distributed_cardiac_passive import (                   # noqa: E402
    NT, NMU, NTHETA, _passive_props, _init_fiber_svars_subset, _kernel,
    _solver_options,
)
from distributed_cardiac_dynamics import _lumped_mass, _robin, DT

N_STEPS = 4
PRES = 500.0                    # constant cavity pressure (Pa) — inflates the LV


def _build_element(mesh, elem_ids):
    props, schema, per_gp = _passive_props()
    elem = CompiledElement(_kernel(), props=props, dof_per_node=3, n_svars=per_gp * 8,
                           mcrd=3, n_elem=len(elem_ids), dt=DT, backend="native",
                           state_schema=schema)
    _init_fiber_svars_subset(elem, mesh, schema, per_gp, elem_ids)
    return elem


def _pressure(mesh, ndof):
    interior = mesh.nodes[mesh.elems[mesh.facets_endo_elem]].mean(axis=1)
    return FollowerPressureOperator(mesh.nodes, ndof, mesh.facets_endo, p=PRES,
                                    dof_per_node=3, interior=interior)


def serial_reference(mesh, ndof, M):
    elem = _build_element(mesh, np.arange(mesh.n_elem))
    grp = ElementGroup(elem, mesh.nodes, mesh.elems, dof_per_node=3, comps=(0, 1, 2))
    inertia = InertiaOperator(M, ndof, damping=0.0)
    robin = _robin(mesh, ndof)
    pressure = _pressure(mesh, ndof)
    U, _ = solve_dynamics_checked(
        [grp, inertia, robin, pressure], np.zeros(ndof), ndof, {},
        dt=DT, n_steps=N_STEPS, rtol=1e-10, maxit=60,
    )
    return U


def main():
    comm = PETSc.COMM_WORLD
    rank, size = comm.getRank(), comm.getSize()

    mesh = geom.build_mesh(
        NT, NMU, NTHETA, flip_helix=True, apex_offset=0.2
    )
    view = KernelMeshView(mesh.nodes, mesh.elems, dof_per_node=3)
    ndof = view.ndof
    M = _lumped_mass(mesh, ndof)

    parts = partition_elements(view, size)
    mine = np.where(parts == rank)[0]
    my_gm, my_coords, _ = element_partition(view, rank, size)
    elem = _build_element(mesh, mine)
    robin = _robin(mesh, ndof)
    pressure = _pressure(mesh, ndof)

    U_par, info = solve_dynamics_distributed(
        ndof, my_gm, my_coords, 3, elem.element_rk_batch, {}, M,
        dt=DT, n_steps=N_STEPS, robin=robin, pressure=pressure,
        tol=1e-9, max_newton=60, **_solver_options())
    U_par = require_distributed_success(
        U_par, info, tol=1e-9, context="distributed cardiac pressure solve"
    )

    failed = False
    if rank == 0:
        U_ser = serial_reference(mesh, ndof, M)
        err = float(np.max(np.abs(U_par - U_ser)))
        umax = float(np.max(np.abs(U_ser)))
        print(f"[size={size}] ndof={ndof} my_ne={info['my_ne']} endo_facets={len(mesh.facets_endo)} "
              f"newton={info['n_newton']} rnorm={info['rnorm']:.2e}")
        print(f"[size={size}] max|U|={umax:.6e}  max|U_par - U_ser|={err:.3e}  "
              f"-> {'OK' if err < 1e-7 else 'FAIL'}", flush=True)
        failed = err >= 1e-7
    comm.barrier()
    if failed:
        raise RuntimeError(f"distributed/serial mismatch: {err:.3e}")


if __name__ == "__main__":
    main()
