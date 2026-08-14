"""Distributed cardiac viscous-state serial-reference check.

This stateful case carries per-Gauss-point viscous history E_prev
(Green-Lagrange strain at the previous step) as committed svars: S_visco = (eta/dt)(E − E_prev).
Each rank owns its elements' svars; the state commit at step end is RANK-LOCAL (elements don't
cross ranks) — done via `solve_dynamics_distributed(step_callback=…)` calling `elem.commit()`
(svars_trial -> svars) exactly as the serial `ElementGroup.commit` does.

Gate: gathered U(t) agrees with the serial cardiac run within 1e-7 over
multiple steps for this script and configuration at 1/2/4 ranks. A
state-scatter or commit bug can appear as drift after step 1.

    OMP_NUM_THREADS=1 <mpiexec> -n {1,2,4} python examples/mpi_smoke/distributed_cardiac_viscous.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
from petsc4py import PETSc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "cardiac_benchmark"))
import geometry as geom                                    # noqa: E402
from material import CardiacHex8                           # noqa: E402
from solver import require_distributed_success, solve_dynamics_checked  # noqa: E402

from coupfe.mesh import KernelMeshView                      # noqa: E402
from coupfe.operators.element_group import ElementGroup     # noqa: E402
from coupfe.operators.inertia import InertiaOperator        # noqa: E402
from coupfe.runtime.compiled_element import CompiledElement  # noqa: E402
from coupfe.assembly.distributed import (                   # noqa: E402
    element_partition, partition_elements, solve_dynamics_distributed,
)
from distributed_cardiac_passive import (                   # noqa: E402
    NT, NMU, NTHETA, _init_fiber_svars_subset, _kernel, _solver_options,
)
from distributed_cardiac_dynamics import _lumped_mass, _robin, DT

N_STEPS = 5
V0 = 0.02                       # initial z-velocity (m/s) — drives viscous state evolution


def _viscous_props():
    """Cardiac props with viscosity ON (state evolves) but softened kappa for convergence."""
    prob = CardiacHex8()
    keys = list(prob._mat.props.keys())
    props = np.asarray(prob._mat.props_array, float).copy()
    props[keys.index("eta")] = 100.0                        # viscous history ON
    props[keys.index("Ta")] = 0.0
    props[keys.index("kappa")] = 1.0e3
    schema = prob._mat.state_schema
    per_gp = sum(v["size"] for v in schema.values())
    return props, schema, per_gp


def _build_element(mesh, elem_ids):
    props, schema, per_gp = _viscous_props()
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

    # RANK-LOCAL per-element state commit at each step end (svars_trial -> svars)
    def commit_state(step, U_full, rk):
        elem.commit()

    U_par, info = solve_dynamics_distributed(
        ndof, my_gm, my_coords, 3, elem.element_rk_batch, {}, M,
        dt=DT, n_steps=N_STEPS, v0=v0, robin=robin, step_callback=commit_state,
        tol=1e-9, max_newton=60, **_solver_options())
    U_par = require_distributed_success(
        U_par, info, tol=1e-9, context="distributed cardiac viscous solve"
    )

    failed = False
    if rank == 0:
        U_ser = serial_reference(mesh, ndof, M, v0)
        err = float(np.max(np.abs(U_par - U_ser)))
        umax = float(np.max(np.abs(U_ser)))
        print(f"[size={size}] ndof={ndof} my_ne={info['my_ne']} steps={N_STEPS} "
              f"newton={info['n_newton']} rnorm={info['rnorm']:.2e}")
        print(f"[size={size}] max|U|={umax:.6e}  max|U_par - U_ser|={err:.3e}  "
              f"-> {'OK' if err < 1e-7 else 'FAIL'}", flush=True)
        failed = err >= 1e-7
    comm.barrier()
    if failed:
        raise RuntimeError(f"distributed/serial mismatch: {err:.3e}")


if __name__ == "__main__":
    main()
