"""Source-matched generalized-alpha cardiac assembly check under MPI.

This deliberately small linear problem isolates the application-owned transient
assembly used by the closed Case A companion.  It exercises, at 1/2/4 ranks:

* rank-local material residual insertion into remotely owned PETSc rows;
* an owned-row consistent-mass action at the ``alpha_m`` acceleration stage;
* Robin spring and dashpot terms at the ``alpha_f`` displacement/velocity stage;
* the exact source parameters ``(.2, .4, .7, .36)``; and
* transactional displacement, velocity, acceleration, and material commits.

Each scalar DOF is independent, so the two accepted endpoints have a closed
form reference.  This is an MPI implementation gate, not ventricular or
FEniCS scientific validation.

    OMP_NUM_THREADS=1 <mpiexec> -n {1,2,4} \
        python examples/mpi_smoke/distributed_cardiac_generalized_alpha.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import scipy.sparse as sp
from mpi4py import MPI
from petsc4py import PETSc


sys.path.insert(
    0, os.path.join(os.path.dirname(__file__), "..", "cardiac_benchmark")
)
from distributed_mass import OwnedRowMassMatrix  # noqa: E402
from distributed_solver import (  # noqa: E402
    DistributedPetscSnesSolver,
    FGMRES_ASM_ILU1_PROFILE,
)
from generalized_alpha import SOURCE_MATCHED_GENERALIZED_ALPHA  # noqa: E402


NDOF = 12
DT = 0.05


class _PartitionedLoadBatch:
    """One zero-tangent load element for each rank-local global DOF."""

    evaluation_mode = "split"
    material_residual_only_available = True
    time_integrator = "generalized-alpha-source-matched"
    pressure_law = "log"
    pressure_stage_alpha_f = 0.4

    def __init__(self, global_dofs):
        self.global_dofs = np.asarray(global_dofs, dtype=np.int64)
        self.force = np.zeros(len(self.global_dofs))
        self.material_dt = DT
        self.commits = 0

    def clear_cache(self):
        pass

    def element_residual_batch(self, coordinates, displacement, increment):
        del coordinates, displacement, increment
        return -self.force[:, None]

    def element_tangent_batch(self, coordinates, displacement, increment):
        del coordinates, displacement, increment
        return np.zeros((len(self.global_dofs), 1, 1))

    def commit(self):
        self.commits += 1


class _DiagonalRobin:
    def __init__(self, stiffness, damping, initial_displacement):
        self.Kmat = sp.diags(stiffness, format="csr")
        self.Cmat = sp.diags(damping, format="csr")
        self.dofs = np.arange(len(stiffness), dtype=np.int64)
        self.u_prev = np.asarray(initial_displacement, dtype=float).copy()
        self.commits = 0

    def commit(self, displacement, state, t, dt):
        del t, dt
        self.u_prev = np.asarray(displacement, dtype=float).copy()
        self.commits += 1
        return state


def _owned_diagonal_mass(comm, diagonal):
    ownership_probe = PETSc.Vec().createMPI(len(diagonal), bsize=1, comm=comm)
    row_start, row_end = ownership_probe.getOwnershipRange()
    ownership_probe.destroy()
    local_rows = np.arange(row_end - row_start, dtype=np.int64)
    global_columns = np.arange(row_start, row_end, dtype=np.int64)
    matrix = sp.csr_matrix(
        (diagonal[row_start:row_end], (local_rows, global_columns)),
        shape=(row_end - row_start, len(diagonal)),
    )
    return OwnedRowMassMatrix(
        matrix,
        row_start,
        row_end,
        len(diagonal),
        np.arange(row_start, row_end, dtype=np.int64),
    )


def _closed_form_endpoint(
    displacement,
    velocity,
    acceleration,
    force,
    mass,
    stiffness,
    damping,
):
    parameters = SOURCE_MATCHED_GENERALIZED_ALPHA

    def residual(endpoint):
        endpoint_acceleration = parameters.acceleration(
            endpoint, displacement, velocity, acceleration, DT
        )
        endpoint_velocity = parameters.velocity(
            endpoint_acceleration, velocity, acceleration, DT
        )
        return (
            mass
            * parameters.acceleration_stage(
                acceleration, endpoint_acceleration
            )
            + stiffness
            * parameters.force_stage(displacement, endpoint)
            + damping
            * parameters.force_stage(velocity, endpoint_velocity)
            - force
        )

    zero = residual(np.zeros_like(displacement))
    slope = residual(np.ones_like(displacement)) - zero
    return -zero / slope


def main():
    comm = PETSc.COMM_WORLD
    rank = int(comm.getRank())
    size = int(comm.getSize())

    all_dofs = np.arange(NDOF, dtype=np.int64)
    local_dofs = all_dofs[all_dofs % size == rank]
    my_gm = local_dofs[:, None]
    my_coordinates = local_dofs.astype(float)[:, None, None]

    mass_values = np.linspace(1.2, 2.3, NDOF)
    stiffness = np.linspace(4.0, 7.0, NDOF)
    damping = np.linspace(0.2, 0.65, NDOF)
    u_expected = np.linspace(-0.006, 0.005, NDOF)
    v_expected = np.linspace(0.018, -0.012, NDOF)
    a_expected = np.linspace(-0.11, 0.08, NDOF)

    mass = _owned_diagonal_mass(comm, mass_values)
    batch = _PartitionedLoadBatch(local_dofs)
    robin = _DiagonalRobin(stiffness, damping, u_expected)
    solver = DistributedPetscSnesSolver(
        NDOF,
        my_gm,
        my_coordinates,
        batch,
        mass,
        dof_per_node=1,
        robin=robin,
        integrator="generalized-alpha",
        implementation=(
            DistributedPetscSnesSolver.
            CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION
        ),
        linear_solver_profile=FGMRES_ASM_ILU1_PROFILE,
        u0=u_expected,
        v0=v_expected,
        a0=a_expected,
    )

    maximum_error = 0.0
    try:
        for step, force in enumerate(
            (
                np.linspace(2.0, 3.1, NDOF),
                np.linspace(-0.8, 0.9, NDOF),
            ),
            start=1,
        ):
            batch.force = force[local_dofs]
            endpoint = _closed_form_endpoint(
                u_expected,
                v_expected,
                a_expected,
                force,
                mass_values,
                stiffness,
                damping,
            )
            new_acceleration = SOURCE_MATCHED_GENERALIZED_ALPHA.acceleration(
                endpoint, u_expected, v_expected, a_expected, DT
            )
            new_velocity = SOURCE_MATCHED_GENERALIZED_ALPHA.velocity(
                new_acceleration, v_expected, a_expected, DT
            )

            displacement, diagnostics = solver.solve_step(t=step * DT, dt=DT)
            maximum_error = max(
                maximum_error,
                float(np.max(np.abs(displacement - endpoint))),
                float(np.max(np.abs(solver._v_prev - new_velocity))),
                float(np.max(np.abs(solver._a_prev - new_acceleration))),
            )
            if diagnostics.snes_converged_reason <= 0:
                raise RuntimeError("generalized-alpha SNES did not converge")
            u_expected = endpoint
            v_expected = new_velocity
            a_expected = new_acceleration

        configuration = solver.configuration()
        metadata = configuration["generalized_alpha"]
        if tuple(
            metadata[key] for key in ("alpha_m", "alpha_f", "gamma", "beta")
        ) != (0.2, 0.4, 0.7, 0.36):
            raise RuntimeError("generalized-alpha configuration metadata changed")
        if batch.commits != 2 or robin.commits != 2:
            raise RuntimeError("accepted generalized-alpha histories were not committed")
    finally:
        solver.close()

    maximum_error = float(MPI.COMM_WORLD.allreduce(maximum_error, op=MPI.MAX))
    tolerance = 2.0e-11
    if rank == 0:
        status = "OK" if maximum_error < tolerance else "FAIL"
        print(
            f"[size={size}] source-matched generalized-alpha two-step "
            f"max error={maximum_error:.3e} -> {status}",
            flush=True,
        )
    MPI.COMM_WORLD.Barrier()
    if maximum_error >= tolerance:
        raise RuntimeError(
            f"generalized-alpha distributed/reference mismatch: {maximum_error:.3e}"
        )


if __name__ == "__main__":
    main()
