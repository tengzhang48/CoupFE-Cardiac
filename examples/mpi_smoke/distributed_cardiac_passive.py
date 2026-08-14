"""Distributed passive cardiac serial-reference check. Run under mpiexec.

The cardiac Holzapfel-Ogden Hex8 bulk kernel (the same `element_rk_batch` the serial
`ElementGroup` uses), distributed by element partition through `solve_distributed`, is compared
with the serial load-stepped solve of the identical problem.

"No state" here means path-INDEPENDENT: viscosity eta=0 and active Ta=0, so the material is pure
anisotropic hyperelastic (HO + volumetric). The gate requires serial-versus-rank
agreement within 1e-7 for this script and configuration at 1/2/4 ranks. The
per-GP fiber structural tensors (ff/ss/fssym) are set
(the anisotropy) but do not evolve. The viscous-state script checks evolving history.

BC (well-posed without Robin): base facets fully fixed; endocardial nodes squeezed radially
inward (toward the z-axis) by a small fraction. Both serial and distributed use the same BC.

    OMP_NUM_THREADS=1 <mpiexec> -n {1,2,4} python examples/mpi_smoke/distributed_cardiac_passive.py
"""
from __future__ import annotations

import atexit
import os
import sys
import tempfile

import numpy as np
from petsc4py import PETSc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "cardiac_benchmark"))
import geometry as geom                                  # noqa: E402
from material import CardiacHex8, struct_tensors, build_kernel  # noqa: E402
from solver import require_distributed_success, solve_increments_checked  # noqa: E402

from coupfe.mesh import KernelMeshView                    # noqa: E402
from coupfe.operators.element_group import ElementGroup   # noqa: E402
from coupfe.runtime.compiled_element import CompiledElement  # noqa: E402
from coupfe.assembly.distributed import (                 # noqa: E402
    element_partition, partition_elements, solve_distributed,
)

NT, NMU, NTHETA = 2, 12, 16          # open-apex smoke mesh (~1872 dof)
SQUEEZE = 0.001                        # endo radial-inward fraction
N_STEPS = 3
_NAT = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
                 [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)


def _hex8_gp_shape():
    g1 = 1.0 / np.sqrt(3.0)
    GP = np.empty((8, 3))
    for g in range(8):
        GP[g, 0] = g1 if g % 2 == 0 else -g1
        GP[g, 1] = g1 if (g // 2) % 2 == 0 else -g1
        GP[g, 2] = g1 if g // 4 == 0 else -g1
    N = np.empty((8, 8))
    for g in range(8):
        for a in range(8):
            N[g, a] = 0.125 * np.prod(1.0 + _NAT[a] * GP[g])
    return N


def _init_fiber_svars_subset(elem, mesh, schema, per_gp, elem_ids):
    """Init per-GP fiber structural tensors for a SUBSET of global elements (in order).
    Local element i corresponds to global element elem_ids[i]."""
    n_gp = elem.svars.shape[1] // per_gp
    o = {k: schema[k]["offset"] for k in schema}
    N = _hex8_gp_shape()
    for li, e in enumerate(elem_ids):
        fn = mesh.fiber_node[mesh.elems[e]]
        sn = mesh.sheet_node[mesh.elems[e]]
        for g in range(n_gp):
            f0 = N[g] @ fn
            s0 = N[g] @ sn
            ff, ss, fssym = struct_tensors(f0, s0)
            base = g * per_gp
            elem.svars[li, base + o["ff"]:base + o["ff"] + 9] = ff.ravel()
            elem.svars[li, base + o["ss"]:base + o["ss"] + 9] = ss.ravel()
            elem.svars[li, base + o["fssym"]:base + o["fssym"] + 9] = fssym.ravel()
    elem.svars_trial = elem.svars.copy()


def _passive_props():
    """Cardiac material props with viscosity + active stress OFF (path-independent)."""
    prob = CardiacHex8()
    keys = list(prob._mat.props.keys())
    props = np.asarray(prob._mat.props_array, float).copy()
    props[keys.index("eta")] = 0.0
    props[keys.index("Ta")] = 0.0
    props[keys.index("kappa")] = 1.0e3   # soften near-incompressibility for the 1-vs-N gate
    schema = prob._mat.state_schema
    per_gp = sum(v["size"] for v in schema.values())
    return props, schema, per_gp


def _dirichlet_fn(mesh):
    base_nodes = np.unique(mesh.facets_base.ravel())
    endo_nodes = np.unique(mesh.facets_endo.ravel())
    xy = mesh.nodes[endo_nodes, :2]

    def fn(frac):
        d = {}
        for i, n in enumerate(endo_nodes):
            d[int(n) * 3 + 0] = -SQUEEZE * frac * xy[i, 0]
            d[int(n) * 3 + 1] = -SQUEEZE * frac * xy[i, 1]
            d[int(n) * 3 + 2] = 0.0
        for n in base_nodes:                              # base wins on shared nodes
            d[int(n) * 3 + 0] = 0.0
            d[int(n) * 3 + 1] = 0.0
            d[int(n) * 3 + 2] = 0.0
        return d
    return fn


_KMOD = {}
_BUILD_TMP = None


def _solver_options():
    """Select an available PETSc solver; prefer reproducible SuperLU_DIST."""
    requested = os.environ.get("COUPFE_CARDIAC_LINEAR_SOLVER", "auto").lower()
    if requested == "auto":
        for candidate in ("superlu_dist", "mumps"):
            if PETSc.Sys.hasExternalPackage(candidate):
                return {"pc": "lu", "ksp_type": "preonly", "solver": candidate}
        return {"pc": "gamg", "ksp_type": "gmres", "solver": None}
    if requested == "gamg":
        return {"pc": "gamg", "ksp_type": "gmres", "solver": None}
    return {"pc": "lu", "ksp_type": "preonly", "solver": requested}


def _kernel():
    global _BUILD_TMP
    if "mod" not in _KMOD:
        # Each MPI process gets an isolated temporary build. This avoids both
        # rank races and generated kernels in the source checkout.
        _BUILD_TMP = tempfile.TemporaryDirectory(prefix="coupfe-cardiac-mpi-build-")
        atexit.register(_BUILD_TMP.cleanup)
        _, _KMOD["mod"] = build_kernel(
            tmpdir=_BUILD_TMP.name, module_name="cardiac_hex8_passive_dist"
        )
    return _KMOD["mod"]


def serial_reference(mesh, ndof, dirichlet_fn):
    props, schema, per_gp = _passive_props()
    elem = CompiledElement(_kernel(), props=props, dof_per_node=3, n_svars=per_gp * 8,
                           mcrd=3, n_elem=mesh.n_elem, dt=1.0, backend="native",
                           state_schema=schema)
    _init_fiber_svars_subset(elem, mesh, schema, per_gp, np.arange(mesh.n_elem))
    grp = ElementGroup(elem, mesh.nodes, mesh.elems, dof_per_node=3, comps=(0, 1, 2))
    full_bc = dirichlet_fn(1.0)
    U, tot = solve_increments_checked(
        [grp], np.zeros(ndof), ndof, full_bc,
        n_steps=N_STEPS, rtol=1e-9, maxit=40,
    )
    print(f'[serial] total newton iters={tot} (n_steps={N_STEPS})', flush=True)
    return U


def main():
    comm = PETSc.COMM_WORLD
    rank, size = comm.getRank(), comm.getSize()

    mesh = geom.build_mesh(
        NT, NMU, NTHETA, flip_helix=True, apex_offset=0.2
    )
    view = KernelMeshView(mesh.nodes, mesh.elems, dof_per_node=3)
    ndof = view.ndof
    dirichlet_fn = _dirichlet_fn(mesh)

    parts = partition_elements(view, size)
    mine = np.where(parts == rank)[0]
    my_gm, my_coords, ndof = element_partition(view, rank, size)

    props, schema, per_gp = _passive_props()
    elem = CompiledElement(_kernel(), props=props, dof_per_node=3, n_svars=per_gp * 8,
                           mcrd=3, n_elem=len(mine), dt=1.0, backend="native",
                           state_schema=schema)
    _init_fiber_svars_subset(elem, mesh, schema, per_gp, mine)

    U_par, info = solve_distributed(
        ndof, my_gm, my_coords, 3, elem.element_rk_batch,
        dirichlet_fn, n_steps=N_STEPS, max_newton=100,
        tol=1e-10, **_solver_options(),
    )
    U_par = require_distributed_success(
        U_par, info, tol=1e-10, context="distributed cardiac passive solve"
    )

    failed = False
    if rank == 0:
        U_ser = serial_reference(mesh, ndof, dirichlet_fn)
        err = float(np.max(np.abs(U_par - U_ser)))
        umax = float(np.max(np.abs(U_ser)))
        print(f"[size={size}] nodes={mesh.n_node} elems={mesh.n_elem} ndof={ndof} "
              f"my_ne={info['my_ne']} owned={info['n_owned']} ghost={info['n_ghost']}")
        print(f"[dist] n_newton={info.get('n_newton')} rnorm={info.get('rnorm'):.2e}", flush=True)
        print(f"[size={size}] max|U|={umax:.6e}  max|U_par - U_ser|={err:.3e}  "
              f"-> {'OK' if err < 1e-7 else 'FAIL'}", flush=True)
        failed = err >= 1e-7
    comm.barrier()
    if failed:
        raise RuntimeError(f"distributed/serial mismatch: {err:.3e}")


if __name__ == "__main__":
    main()
