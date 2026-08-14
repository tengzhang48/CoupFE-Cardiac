"""Distributed cardiac dynamics + Robin serial-reference check. Run under mpiexec.

Backward-Euler implicit dynamics: cardiac HO bulk (passive, no evolving state) + lumped-mass
inertia + the pericardial Robin spring-dashpot (base full + epi normal), distributed through
`solve_dynamics_distributed(..., robin=...)`. Driven by a small initial velocity so the LV
oscillates on its Robin springs and the dashpot damps it — a clean transient with no body force.

The Robin K/C are constant reference-config sparse matrices; the distributed solve adds each
rank's OWNED rows of (K + C/dt) to the tangent and the owned slice of K·u + C·(u−u_prev)/dt to
the residual (full u gathered — Robin couples boundary nodes that can split across ranks).

Gate: gathered U(t) agrees with serial `solve_dynamics` within 1e-7 for this
script and configuration at 1/2/4 ranks.

    OMP_NUM_THREADS=1 <mpiexec> -n {1,2,4} python examples/mpi_smoke/distributed_cardiac_dynamics.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
from petsc4py import PETSc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "cardiac_benchmark"))
import geometry as geom                                    # noqa: E402
from robin import RobinOperator                            # noqa: E402
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

DENSITY = 1.0e3
A_TOP, B_TOP = 1.0e5, 5.0e3
A_EPI, B_EPI = 1.0e8, 5.0e3
DT, N_STEPS = 1.0e-3, 4
V0 = 0.01                       # uniform initial z-velocity (m/s)


def _lumped_mass(mesh, ndof):
    nodal = np.zeros(mesh.n_node)
    for e in mesh.elems:
        np.add.at(nodal, e, DENSITY * geom._hex_volume(mesh.nodes[e]) / 8.0)
    M = np.zeros(ndof)
    for c in range(3):
        M[c::3] = nodal
    return M


def _robin(mesh, ndof):
    return RobinOperator(
        mesh.nodes, ndof,
        [(mesh.facets_base, A_TOP, B_TOP, "full"),
         (mesh.facets_epi, A_EPI, B_EPI, "normal")],
        dof_per_node=3)


def _build_element(mesh, elem_ids):
    props, schema, per_gp = _passive_props()
    elem = CompiledElement(_kernel(), props=props, dof_per_node=3, n_svars=per_gp * 8,
                           mcrd=3, n_elem=len(elem_ids), dt=DT, backend="native",
                           state_schema=schema)
    _init_fiber_svars_subset(elem, mesh, schema, per_gp, elem_ids)
    return elem


def serial_reference(mesh, ndof, M, v0):
    elem = _build_element(mesh, np.arange(mesh.n_elem))
    grp = ElementGroup(elem, mesh.nodes, mesh.elems, dof_per_node=3, comps=(0, 1, 2))
    inertia = InertiaOperator(M, ndof, v0=v0, damping=0.0)
    robin = _robin(mesh, ndof)
    U, _ = solve_dynamics_checked(
        [grp, inertia, robin], np.zeros(ndof), ndof, {},
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
    v0 = np.zeros(ndof); v0[2::3] = V0

    parts = partition_elements(view, size)
    mine = np.where(parts == rank)[0]
    my_gm, my_coords, _ = element_partition(view, rank, size)
    elem = _build_element(mesh, mine)
    robin = _robin(mesh, ndof)

    U_par, info = solve_dynamics_distributed(
        ndof, my_gm, my_coords, 3, elem.element_rk_batch, {}, M,
        dt=DT, n_steps=N_STEPS, v0=v0, robin=robin,
        tol=1e-9, max_newton=60, **_solver_options())
    U_par = require_distributed_success(
        U_par, info, tol=1e-9, context="distributed cardiac dynamics"
    )

    failed = False
    if rank == 0:
        U_ser = serial_reference(mesh, ndof, M, v0)
        err = float(np.max(np.abs(U_par - U_ser)))
        umax = float(np.max(np.abs(U_ser)))
        print(f"[size={size}] ndof={ndof} my_ne={info['my_ne']} owned={info['n_owned']} "
              f"newton={info['n_newton']} rnorm={info['rnorm']:.2e}")
        print(f"[size={size}] max|U|={umax:.6e}  max|U_par - U_ser|={err:.3e}  "
              f"-> {'OK' if err < 1e-7 else 'FAIL'}", flush=True)
        failed = err >= 1e-7
    comm.barrier()
    if failed:
        raise RuntimeError(f"distributed/serial mismatch: {err:.3e}")


if __name__ == "__main__":
    main()
