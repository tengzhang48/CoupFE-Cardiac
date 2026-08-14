"""Cardiac-owned distributed PETSc SNES for the Q1/P0 Step 0 companion.

This is deliberately separate from both the retained serial SNES adapter and
CoupFE Core's generic distributed smoke solver.  The benchmark requires a
persistent ``PETSc.COMM_WORLD`` SNES, collective rejection of invalid
finite-strain line-search trials, and a second global residual check before any
history is committed.  Those experiment-specific rules live here.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import importlib.metadata
import time

import numpy as np


PETSC_FUNCTION_DOMAIN_REJECTION_API = "nonfinite residual for PETSc BT"

DIRECT_SUPERLU_DIST_PROFILE = "direct-superlu-dist"
FGMRES_GAMG_RIGID_PROFILE = "fgmres-gamg-rigid"
FGMRES_GAMG_RIGID_REBUILD_PROFILE = "fgmres-gamg-rigid-rebuild"
FGMRES_ASM_LU_PROFILE = "fgmres-asm-lu"
FGMRES_ASM_ILU1_PROFILE = "fgmres-asm-ilu1"
LINEAR_SOLVER_PROFILES = frozenset(
    {
        DIRECT_SUPERLU_DIST_PROFILE,
        FGMRES_GAMG_RIGID_PROFILE,
        FGMRES_GAMG_RIGID_REBUILD_PROFILE,
        FGMRES_ASM_LU_PROFILE,
        FGMRES_ASM_ILU1_PROFILE,
    }
)

from distributed_mass import OwnedRowMassMatrix
from generalized_alpha import SOURCE_MATCHED_GENERALIZED_ALPHA
from local_pressure import InvalidDeformationError
from solver import PetscSnesSettings, RECOVERED_PETSC_SNES_SETTINGS


def linear_solver_requires_node_aligned_ownership(profile):
    """Return whether a profile requires complete 3-DOF node ownership."""
    profile = str(profile)
    if profile not in LINEAR_SOLVER_PROFILES:
        raise ValueError(f"unsupported distributed linear solver profile {profile!r}")
    return profile in {
        FGMRES_GAMG_RIGID_PROFILE,
        FGMRES_GAMG_RIGID_REBUILD_PROFILE,
    }


def _load_parallel_runtime():
    try:
        import petsc4py
    except ImportError as error:
        raise RuntimeError(
            "the distributed cardiac companion requires petsc4py "
            "(install the application with the 'mpi' extra)"
        ) from error
    petsc4py.init([])
    from petsc4py import PETSc

    try:
        from mpi4py import MPI
    except ImportError as error:
        raise RuntimeError(
            "the distributed cardiac companion requires mpi4py "
            "(install the application with the 'mpi' extra)"
        ) from error
    return PETSc, MPI


@dataclass(frozen=True)
class DistributedSnesStepDiagnostics:
    """JSON-compatible diagnostics for one accepted distributed step."""

    time: float
    dt: float
    ranks: int
    initial_residual_norm: float
    final_residual_norm: float
    residual_acceptance_threshold: float
    petsc_function_norm: float
    snes_converged_reason: int
    ksp_converged_reason: int
    nonlinear_iterations: int
    linear_iterations: int
    residual_history: tuple
    assembly_seconds: float
    solve_seconds: float
    function_domain_rejections: int = 0
    last_function_domain_error: str | None = None
    ksp_iterations_per_solve: tuple = ()
    ksp_reasons_per_solve: tuple = ()
    ksp_residual_histories: tuple = ()

    def as_dict(self):
        return asdict(self)


class DistributedSnesSolveError(RuntimeError):
    """Fail-closed distributed solve error with optional step diagnostics."""

    def __init__(self, message, diagnostics=None):
        super().__init__(message)
        self.diagnostics = diagnostics


class DistributedAssemblyError(RuntimeError):
    """An unrelated rank-local callback error made the collective assembly fatal."""


class DistributedPetscSnesSolver:
    """Persistent ``COMM_WORLD`` SNES for partitioned cardiac Hex8 batches.

    Parameters are global except ``my_gm``, ``my_coordinates``, and
    ``element_batch``, which describe only the elements assigned to this rank.
    Element blocks are inserted with global indices and ``ADD_VALUES``.  Inertia
    and replicated boundary operators contribute only rows owned by the rank.

    The accepted displacement and velocity are retained on every rank because
    the benchmark already replicates the global mesh and because Robin and
    follower-pressure evaluation require a full displacement.  The PETSc
    solution, residual, Jacobian, KSP, and SNES remain distributed and persist
    for the lifetime of the object.
    """

    HISTORICAL_Q1P0_IMPLEMENTATION = "cardiac-owned-distributed-q1p0-step0"
    CLOSED_STD_KAPPA_IMPLEMENTATION = (
        "cardiac-owned-distributed-closed-std-kappa-step0"
    )
    CLOSED_LOCAL_PRESSURE_IMPLEMENTATION = (
        "cardiac-owned-distributed-closed-local-pressure-step0"
    )
    CLOSED_LOCAL_PRESSURE_PAPER_IMPLEMENTATION = (
        "cardiac-owned-distributed-closed-local-pressure-mean-logj-"
        "paper-j2-step0"
    )
    CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION = (
        "cardiac-owned-distributed-closed-std-kappa-generalized-alpha-step0"
    )
    CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION = (
        "cardiac-owned-distributed-closed-local-pressure-generalized-alpha-step0"
    )
    CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION = (
        "cardiac-owned-distributed-closed-local-pressure-mean-logj-"
        "paper-j2-generalized-alpha-step0"
    )
    SUPPORTED_IMPLEMENTATIONS = frozenset(
        {
            HISTORICAL_Q1P0_IMPLEMENTATION,
            CLOSED_STD_KAPPA_IMPLEMENTATION,
            CLOSED_LOCAL_PRESSURE_IMPLEMENTATION,
            CLOSED_LOCAL_PRESSURE_PAPER_IMPLEMENTATION,
            CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION,
            CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION,
            CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION,
        }
    )
    # Backward-compatible class spelling used by historical callers/tests.
    IMPLEMENTATION = HISTORICAL_Q1P0_IMPLEMENTATION

    def __init__(
        self,
        ndof,
        my_gm,
        my_coordinates,
        element_batch,
        mass,
        *,
        dof_per_node=3,
        robin=None,
        pressure=None,
        settings=RECOVERED_PETSC_SNES_SETTINGS,
        factor_solver_type=None,
        linear_solver_profile=DIRECT_SUPERLU_DIST_PROFILE,
        node_coordinates=None,
        implementation=HISTORICAL_Q1P0_IMPLEMENTATION,
        integrator="be",
        u0=None,
        v0=None,
        a0=None,
    ):
        if not isinstance(settings, PetscSnesSettings):
            raise TypeError("settings must be a PetscSnesSettings instance")
        if settings != RECOVERED_PETSC_SNES_SETTINGS:
            raise ValueError(
                "the distributed Step 0 companion requires the recovered exact "
                "PETSc SNES settings"
            )
        self.settings = settings
        self.linear_solver_profile = str(linear_solver_profile)
        if self.linear_solver_profile not in LINEAR_SOLVER_PROFILES:
            raise ValueError(
                "unsupported distributed cardiac linear solver profile"
            )
        self.implementation = str(implementation)
        if self.implementation not in self.SUPPORTED_IMPLEMENTATIONS:
            raise ValueError(
                "unsupported distributed cardiac implementation identity"
            )
        self.integrator = str(integrator)
        if self.integrator not in {"be", "generalized-alpha"}:
            raise ValueError(
                "distributed cardiac integrator must be 'be' or "
                "'generalized-alpha'"
            )
        generalized_alpha_implementations = {
            self.CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION,
            self.CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION,
            self.CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION,
        }
        if (self.integrator == "generalized-alpha") != (
            self.implementation in generalized_alpha_implementations
        ):
            raise ValueError(
                "distributed cardiac integrator and implementation identity "
                "must describe the same time-discretization contract"
            )
        expected_pressure_law = {
            self.HISTORICAL_Q1P0_IMPLEMENTATION: "log",
            self.CLOSED_STD_KAPPA_IMPLEMENTATION: "not-applicable",
            self.CLOSED_LOCAL_PRESSURE_IMPLEMENTATION: "log",
            self.CLOSED_LOCAL_PRESSURE_PAPER_IMPLEMENTATION: "paper",
            self.CLOSED_STD_KAPPA_GENERALIZED_ALPHA_IMPLEMENTATION: (
                "not-applicable"
            ),
            self.CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION: "log",
            self.CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION: (
                "paper"
            ),
        }[self.implementation]
        pressure_law_default = (
            "log"
            if self.implementation == self.HISTORICAL_Q1P0_IMPLEMENTATION
            else "not-applicable"
        )
        actual_pressure_law = getattr(
            element_batch, "pressure_law", pressure_law_default
        )
        if actual_pressure_law != expected_pressure_law:
            raise ValueError(
                "distributed cardiac implementation identity and local-pressure "
                "law disagree"
            )
        self.generalized_alpha = (
            SOURCE_MATCHED_GENERALIZED_ALPHA
            if self.integrator == "generalized-alpha"
            else None
        )
        batch_integrator = getattr(element_batch, "time_integrator", "be")
        expected_batch_integrator = (
            "generalized-alpha-source-matched"
            if self.generalized_alpha is not None
            else "be"
        )
        if batch_integrator != expected_batch_integrator:
            raise ValueError(
                "distributed solver and material batch time integrators disagree"
            )
        material_dt = getattr(element_batch, "material_dt", None)
        self.material_dt = None if material_dt is None else float(material_dt)
        if self.generalized_alpha is not None and (
            self.material_dt is None
            or not np.isfinite(self.material_dt)
            or self.material_dt <= 0.0
        ):
            raise ValueError(
                "source-matched generalized-alpha requires the compiled "
                "material time step"
            )
        local_pressure_generalized_alpha_implementations = {
            self.CLOSED_LOCAL_PRESSURE_GENERALIZED_ALPHA_IMPLEMENTATION,
            self.CLOSED_LOCAL_PRESSURE_PAPER_GENERALIZED_ALPHA_IMPLEMENTATION,
        }
        if self.implementation in (
            local_pressure_generalized_alpha_implementations
        ) and not np.isclose(
            getattr(element_batch, "pressure_stage_alpha_f", np.nan),
            SOURCE_MATCHED_GENERALIZED_ALPHA.alpha_f,
            rtol=0.0,
            atol=1.0e-15,
        ):
            raise ValueError(
                "source-matched local pressure must be staged with alpha_f=0.4"
            )
        expected_factor = (
            "superlu_dist"
            if self.linear_solver_profile == DIRECT_SUPERLU_DIST_PROFILE
            else None
        )
        requested_factor = (
            expected_factor if factor_solver_type is None else str(factor_solver_type)
        )
        if requested_factor != expected_factor:
            raise ValueError(
                "factor_solver_type does not match the selected distributed "
                "linear solver profile"
            )
        self.factor_solver_type = expected_factor or "none"
        if not isinstance(ndof, int) or ndof < 1:
            raise ValueError("ndof must be a positive integer")
        if not isinstance(dof_per_node, int) or dof_per_node < 1:
            raise ValueError("dof_per_node must be a positive integer")
        self.ndof = ndof
        self.dof_per_node = dof_per_node
        self.node_aligned_ownership = linear_solver_requires_node_aligned_ownership(
            self.linear_solver_profile
        )
        self.node_coordinates = None
        if self.node_aligned_ownership:
            if dof_per_node != 3 or ndof % dof_per_node:
                raise ValueError(
                    "the GAMG rigid-body profile requires three displacement "
                    "DOFs per complete node"
                )
            self.node_coordinates = np.asarray(node_coordinates, dtype=float)
            expected_coordinates = (ndof // dof_per_node, 3)
            if (
                self.node_coordinates.shape != expected_coordinates
                or not np.all(np.isfinite(self.node_coordinates))
            ):
                raise ValueError(
                    "the GAMG rigid-body profile requires finite global node "
                    f"coordinates with shape {expected_coordinates}"
                )
        self.my_gm = np.asarray(my_gm, dtype=np.int64)
        if self.my_gm.ndim != 2:
            raise ValueError("my_gm must have shape (n_local_element, ndof_element)")
        if self.my_gm.size and (
            np.min(self.my_gm) < 0 or np.max(self.my_gm) >= ndof
        ):
            raise ValueError("my_gm contains a DOF outside the global system")
        self.my_coordinates = np.asarray(my_coordinates, dtype=float)
        if self.my_coordinates.ndim < 2 or len(self.my_coordinates) != len(
            self.my_gm
        ):
            raise ValueError(
                "my_coordinates must contain one coordinate block per local element"
            )
        if not np.all(np.isfinite(self.my_coordinates)):
            raise ValueError("my_coordinates must be finite")
        required_batch_methods = (
            "element_residual_batch",
            "element_tangent_batch",
            "commit",
        )
        if not all(hasattr(element_batch, method) for method in required_batch_methods):
            raise TypeError(
                "element_batch must supply element_residual_batch(...), "
                "element_tangent_batch(...), and commit()"
            )
        self.element_batch = element_batch

        self._consistent_mass = isinstance(mass, OwnedRowMassMatrix)
        if self._consistent_mass:
            if mass.ndof != ndof:
                raise ValueError("owned consistent mass does not match ndof")
            self.mass = mass
            self.mass_representation = mass.representation
            self.mass_partition = mass.assembly_policy
        else:
            self.mass = np.asarray(mass, dtype=float)
            if self.mass.shape != (ndof,):
                raise ValueError(
                    "mass must be an owned-row consistent matrix or have "
                    f"diagonal shape {(ndof,)}"
                )
            if not np.all(np.isfinite(self.mass)) or np.any(self.mass < 0.0):
                raise ValueError("mass must be finite and nonnegative")
            self.mass_representation = "lumped_row_sum"
            self.mass_partition = "replicated-diagonal-owned-entry-insertion"
        self.robin = robin
        self.pressure = pressure
        if robin is not None:
            for attribute in ("Kmat", "Cmat", "dofs", "u_prev"):
                if not hasattr(robin, attribute):
                    raise TypeError(
                        f"robin operator is missing required attribute {attribute!r}"
                    )
            if robin.Kmat.shape != (ndof, ndof) or robin.Cmat.shape != (
                ndof,
                ndof,
            ):
                raise ValueError("Robin matrices must match the global system")
        if pressure is not None and not (
            hasattr(pressure, "residual") and hasattr(pressure, "tangent")
        ):
            raise TypeError("pressure must implement residual(...) and tangent(...)")

        self._PETSc, self._MPI = _load_parallel_runtime()
        PETSc = self._PETSc
        self.comm = PETSc.COMM_WORLD
        self.mpi_comm = self._MPI.COMM_WORLD
        self.rank = int(self.comm.getRank())
        self.size = int(self.comm.getSize())
        if self.mpi_comm.Get_rank() != self.rank or self.mpi_comm.Get_size() != self.size:
            raise RuntimeError("petsc4py and mpi4py COMM_WORLD communicators disagree")
        if (
            self.linear_solver_profile == DIRECT_SUPERLU_DIST_PROFILE
            and not PETSc.Sys.hasExternalPackage(self.factor_solver_type)
        ):
            raise RuntimeError(
                "PETSc does not provide the required distributed factor solver "
                f"{self.factor_solver_type!r}"
            )

        self._u_prev = self._initial_state(u0, "u0")
        self._v_prev = self._initial_state(v0, "v0")
        if self.integrator == "generalized-alpha":
            self._a_prev = self._initial_state(a0, "a0")
        else:
            if a0 is not None:
                raise ValueError("a0 is only valid for generalized-alpha")
            self._a_prev = None
        if self.robin is not None:
            robin_initial = np.asarray(self.robin.u_prev, dtype=float)
            if robin_initial.shape != (ndof,) or not np.all(
                np.isfinite(robin_initial)
            ):
                raise ValueError("Robin u_prev must be a finite global vector")
            if not np.array_equal(robin_initial, self._u_prev):
                raise ValueError("Robin u_prev must equal the distributed initial u0")

        vector_block_size = dof_per_node if self.node_aligned_ownership else 1
        self._x = PETSc.Vec().createMPI(
            ndof, bsize=vector_block_size, comm=self.comm
        )
        self._x.set(0.0)
        self._f = self._x.duplicate()
        self._rs, self._re = self._x.getOwnershipRange()
        if self._consistent_mass and (
            self.mass.row_start != self._rs or self.mass.row_end != self._re
        ):
            raise ValueError(
                "owned consistent-mass rows do not match PETSc vector ownership: "
                f"mass=[{self.mass.row_start}, {self.mass.row_end}), "
                f"PETSc=[{self._rs}, {self._re})"
            )
        self._owned = np.arange(self._rs, self._re, dtype=PETSc.IntType)
        self._scatter_all, self._x_all = PETSc.Scatter.toAll(self._x)
        self._rows_all = self.my_gm.astype(PETSc.IntType)

        self._jacobian = PETSc.Mat().createAIJ(
            [ndof, ndof], bsize=vector_block_size, comm=self.comm
        )
        # A cardiac Hex8 row normally couples to 27 nodes x 3 components;
        # follower-pressure rows add a small boundary stencil.  Dynamic
        # allocation remains enabled so this estimate is never correctness-
        # critical.
        self._jacobian.setPreallocationNNZ(max(1, dof_per_node * 30))
        self._jacobian.setOption(
            PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False
        )

        self._snes = PETSc.SNES().create(comm=self.comm)
        self._snes.setFunction(self._form_function, self._f)
        self._snes.setJacobian(self._form_jacobian, self._jacobian)
        self._snes.setMonitor(self._monitor)
        self._snes.setType(settings.snes_type)
        self._line_search_configuration = self._configure_line_search()
        self._petsc_option_prefix = f"coupfe_cardiac_{id(self):x}_"
        self._petsc_options = {}
        self._coordinate_vector = None
        self._near_nullspace = None
        self._near_nullspace_modes = 0
        self._linear_solver_configuration = self._configure_linear_solver()
        self._snes.setTolerances(
            rtol=settings.rtol,
            atol=settings.atol,
            stol=settings.stol,
            max_it=settings.max_it,
        )

        self._active = None
        self._residual_history = []
        self._ksp_iterations_per_solve = []
        self._ksp_reasons_per_solve = []
        self._ksp_residual_histories = []
        self._step_assembly_seconds = 0.0
        self._function_domain_rejections = 0
        self._last_function_domain_error = None
        self._callback_fatal_error = None
        self._accepted_steps = 0
        self.last_diagnostics = None
        self._closed = False

    def _initial_state(self, values, name):
        if values is None:
            return np.zeros(self.ndof, dtype=float)
        result = np.asarray(values, dtype=float)
        if result.shape != (self.ndof,) or not np.all(np.isfinite(result)):
            raise ValueError(f"{name} must be a finite vector of shape {(self.ndof,)}")
        return result.copy()

    def _configure_line_search(self):
        PETSc = self._PETSc
        if hasattr(self._snes, "getLineSearch"):
            self._snes.getLineSearch().setType(self.settings.line_search_type)
            return "SNES.getLineSearch"
        prefix = f"coupfe_cardiac_mpi_{id(self):x}_"
        self._snes.setOptionsPrefix(prefix)
        options = PETSc.Options()
        option_name = f"{prefix}snes_linesearch_type"
        options[option_name] = self.settings.line_search_type
        try:
            self._snes.setFromOptions()
        finally:
            del options[option_name]
        return "namespaced PETSc option"

    def _install_petsc_options(self, values):
        """Install namespaced options retained until PC teardown.

        PETSc creates GAMG levels and ASM sub-KSPs lazily.  Their options must
        therefore remain in the database until :meth:`close`, rather than only
        around the initial ``setFromOptions`` call.
        """
        options = self._PETSc.Options()
        for name, value in values.items():
            full_name = f"{self._petsc_option_prefix}{name}"
            options[full_name] = value
            self._petsc_options[full_name] = str(value)

    def _attach_rigid_body_near_nullspace(self, pc):
        PETSc = self._PETSc
        coordinates = PETSc.Vec().createMPI(
            self.ndof, bsize=self.dof_per_node, comm=self.comm
        )
        row_start, row_end = coordinates.getOwnershipRange()
        if row_start % self.dof_per_node or row_end % self.dof_per_node:
            coordinates.destroy()
            raise RuntimeError(
                "GAMG coordinate ownership split a displacement node"
            )
        local_coordinates = self.node_coordinates[
            row_start // self.dof_per_node : row_end // self.dof_per_node
        ]
        coordinates.getArray(readonly=False)[:] = local_coordinates.reshape(-1)
        coordinates.assemble()
        near_nullspace = PETSc.NullSpace().createRigidBody(coordinates)
        modes = near_nullspace.getVecs()
        if (
            len(modes) == 2
            and isinstance(modes[0], (bool, np.bool_))
            and isinstance(modes[1], (tuple, list))
        ):
            modes = modes[1]
        if len(modes) != 6:
            near_nullspace.destroy()
            coordinates.destroy()
            raise RuntimeError(
                "PETSc did not construct the six 3-D rigid-body modes"
            )
        self._jacobian.setNearNullSpace(near_nullspace)
        pc.setCoordinates(local_coordinates)
        self._coordinate_vector = coordinates
        self._near_nullspace = near_nullspace
        self._near_nullspace_modes = len(modes)

    def _configure_linear_solver(self):
        """Configure one explicit direct or iterative linear-solver profile."""
        PETSc = self._PETSc
        ksp = self._snes.getKSP()
        pc = ksp.getPC()

        if self.linear_solver_profile == DIRECT_SUPERLU_DIST_PROFILE:
            ksp.setType(self.settings.ksp_type)
            pc.setType(self.settings.pc_type)
            pc.setFactorSolverType(self.factor_solver_type)
            return {
                "ksp_type": self.settings.ksp_type,
                "pc_type": self.settings.pc_type,
                "ksp_rtol": None,
                "ksp_atol": None,
                "ksp_divtol": None,
                "ksp_max_it": None,
                "gmres_restart": None,
                "ksp_norm_type": "default",
                "pc_side": "default",
                "preconditioner": "global LU with SuperLU_DIST",
            }

        robust_gamg = (
            self.linear_solver_profile == FGMRES_GAMG_RIGID_REBUILD_PROFILE
        )
        ksp_max_it = 400 if robust_gamg else 200
        gmres_restart = 100 if robust_gamg else 50

        ksp.setOptionsPrefix(self._petsc_option_prefix)
        ksp.setType("fgmres")
        ksp.setTolerances(
            rtol=1.0e-8,
            atol=1.0e-12,
            divtol=1.0e4,
            max_it=ksp_max_it,
        )
        ksp.setGMRESRestart(gmres_restart)
        ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
        ksp.setPCSide(PETSc.PC.Side.RIGHT)
        ksp.setInitialGuessNonzero(False)
        # Let SNES return its negative convergence reason so the application can
        # attach the final KSP history to a fail-closed, pre-commit rejection.
        ksp.setErrorIfNotConverged(False)
        ksp.setConvergenceHistory(ksp_max_it + 1, reset=True)
        ksp.setMonitor(self._ksp_monitor)

        option_values = {}
        if self.linear_solver_profile in {
            FGMRES_GAMG_RIGID_PROFILE,
            FGMRES_GAMG_RIGID_REBUILD_PROFILE,
        }:
            pc.setType("gamg")
            pc.setGAMGType("agg")
            self._attach_rigid_body_near_nullspace(pc)
            option_values = {
                "pc_gamg_type": "agg",
                "pc_gamg_agg_nsmooths": 1,
                "pc_gamg_threshold": 0.01,
                # PT-Scotch repartitioning in this PETSc build requires an
                # unavailable MPI_THREAD_MULTIPLE execution mode.
                "pc_gamg_repartition": "false",
            }
            if robust_gamg:
                # The tangent changes substantially through peak pressure.  A
                # hierarchy frozen near the undeformed state can otherwise
                # become ineffective late in the trajectory.
                option_values["pc_gamg_reuse_interpolation"] = "false"
            preconditioner = (
                "PETSc GAMG aggregation with six rigid-body near-null modes"
                + (
                    " and interpolation rebuilt for changed matrices"
                    if robust_gamg
                    else ""
                )
            )
        else:
            pc.setType("asm")
            pc.setASMOverlap(1)
            pc.setASMType(PETSc.PC.ASMType.RESTRICT)
            if self.linear_solver_profile == FGMRES_ASM_LU_PROFILE:
                option_values = {
                    "sub_ksp_type": "preonly",
                    "sub_pc_type": "lu",
                    "sub_pc_factor_mat_solver_type": "superlu",
                }
                preconditioner = (
                    "restricted additive Schwarz overlap 1 with local SuperLU"
                )
            else:
                option_values = {
                    "sub_ksp_type": "preonly",
                    "sub_pc_type": "ilu",
                    "sub_pc_factor_levels": 1,
                    "sub_pc_factor_shift_type": "nonzero",
                }
                preconditioner = (
                    "restricted additive Schwarz overlap 1 with local ILU(1)"
                )

        self._install_petsc_options(option_values)
        ksp.setFromOptions()
        # Reassert the outer convergence contract after reading only the
        # namespaced preconditioner/subsolver options above.
        ksp.setType("fgmres")
        ksp.setTolerances(
            rtol=1.0e-8,
            atol=1.0e-12,
            divtol=1.0e4,
            max_it=ksp_max_it,
        )
        ksp.setGMRESRestart(gmres_restart)
        ksp.setNormType(PETSc.KSP.NormType.UNPRECONDITIONED)
        ksp.setPCSide(PETSc.PC.Side.RIGHT)
        ksp.setErrorIfNotConverged(False)
        return {
            "ksp_type": "fgmres",
            "pc_type": pc.getType(),
            "ksp_rtol": 1.0e-8,
            "ksp_atol": 1.0e-12,
            "ksp_divtol": 1.0e4,
            "ksp_max_it": ksp_max_it,
            "gmres_restart": gmres_restart,
            "ksp_norm_type": "unpreconditioned",
            "pc_side": "right",
            "ksp_error_if_not_converged": False,
            "preconditioner": preconditioner,
            "petsc_options": dict(sorted(self._petsc_options.items())),
        }

    @property
    def accepted_steps(self):
        return self._accepted_steps

    def configuration(self):
        """Return JSON-compatible solver and MPI provenance."""
        configuration = {
            "name": "petsc-snes-mpi",
            "implementation": self.implementation,
            "settings_source": "recovered 2026-06-27 Case B development adapter",
            "communicator": "PETSc.COMM_WORLD",
            "ranks": self.size,
            "matrix_scope": "persistent one solver instance per MPI run",
            "dirichlet_support": "none",
            "rank_local_element_assembly": True,
            "global_mesh_replicated": True,
            "collective_invalid_trial_policy": (
                "all-owned-residual-entries-inf-for-bt"
            ),
            "function_domain_rejection_api": (
                PETSC_FUNCTION_DOMAIN_REJECTION_API
            ),
            "commit_policy": "independent-global-residual-check-before-commit",
            "element_evaluation_mode": getattr(
                self.element_batch, "evaluation_mode", "unspecified"
            ),
            "material_batch_time_integrator": getattr(
                self.element_batch, "time_integrator", "unspecified"
            ),
            "local_pressure_law": getattr(
                self.element_batch,
                "pressure_law",
                "log"
                if self.implementation == self.HISTORICAL_Q1P0_IMPLEMENTATION
                else "not-applicable",
            ),
            "compiled_material_dt": self.material_dt,
            "compiled_material_residual_only_available": bool(
                getattr(
                    self.element_batch,
                    "material_residual_only_available",
                    False,
                )
            ),
            "mass_representation": self.mass_representation,
            "mass_partition": self.mass_partition,
            "mass_owned_row_range": [int(self._rs), int(self._re)],
            "mass_local_nnz": int(
                self.mass.matrix.nnz
                if self._consistent_mass
                else np.count_nonzero(self.mass[self._rs : self._re])
            ),
            "factor_solver_type": self.factor_solver_type,
            "linear_solver_profile": self.linear_solver_profile,
            "node_aligned_ownership": self.node_aligned_ownership,
            "vector_block_size": int(self._x.getBlockSize()),
            "matrix_block_size": int(self._jacobian.getBlockSize()),
            "near_nullspace_kind": (
                "six-rigid-body-modes"
                if self._near_nullspace_modes
                else "none"
            ),
            "near_nullspace_mode_count": int(self._near_nullspace_modes),
            "line_search_configuration_api": self._line_search_configuration,
            "time_integrator": self.integrator,
            **asdict(self.settings),
        }
        if self.generalized_alpha is None:
            configuration.update(
                {
                    "acceleration_stage": "backward-euler-endpoint",
                    "force_stage": "backward-euler-endpoint",
                    "material_viscous_rate": "green-lagrange-backward-difference",
                    "nonlinear_initial_guess": "backward-euler-kinematic-predictor",
                }
            )
        else:
            configuration["generalized_alpha"] = (
                self.generalized_alpha.metadata()
            )
            configuration.update(
                {
                    "acceleration_stage": "1-alpha_m",
                    "force_stage": "1-alpha_f",
                    "material_viscous_rate": (
                        "sym(F_stage^T*grad(v_stage))"
                    ),
                    "nonlinear_initial_guess": "accepted-u_n-like-simula",
                    "local_pressure_stage_alpha_f": getattr(
                        self.element_batch, "pressure_stage_alpha_f", None
                    ),
                }
            )
        configuration.update(self._linear_solver_configuration)
        try:
            configuration["petsc4py_version"] = importlib.metadata.version(
                "petsc4py"
            )
        except importlib.metadata.PackageNotFoundError:
            configuration["petsc4py_version"] = "unavailable"
        version = self._PETSc.Sys.getVersion()
        configuration["petsc_version"] = ".".join(
            str(component) for component in version
        )
        if self.linear_solver_profile == DIRECT_SUPERLU_DIST_PROFILE:
            try:
                configured_factor = (
                    self._snes.getKSP().getPC().getFactorSolverType()
                )
            except (AttributeError, TypeError, RuntimeError):
                configured_factor = None
            configuration["configured_factor_solver_type"] = (
                str(configured_factor) if configured_factor else "unknown"
            )
        else:
            configuration["configured_factor_solver_type"] = "not-applicable"
        return configuration

    def predictor(self, dt):
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt must be finite and positive")
        if self.generalized_alpha is not None:
            # The audited Simula NonlinearSolver leaves its unknown at the
            # accepted endpoint between calls.  Match that initial guess; the
            # Newmark acceleration base remains inside acceleration(...).
            return self._u_prev.copy()
        return self._u_prev + float(dt) * self._v_prev

    def _generalized_alpha_trial(self, displacement):
        if self.generalized_alpha is None:
            raise RuntimeError(
                "generalized-alpha trial requested for backward Euler"
            )
        active = self._active
        if active is None:
            raise RuntimeError("distributed callback has no active step")
        parameters = self.generalized_alpha
        acceleration = parameters.acceleration(
            displacement,
            self._u_prev,
            self._v_prev,
            self._a_prev,
            active["dt"],
        )
        velocity = parameters.velocity(
            acceleration, self._v_prev, self._a_prev, active["dt"]
        )
        return {
            "acceleration": acceleration,
            "velocity": velocity,
            "displacement_stage": parameters.force_stage(
                self._u_prev, displacement
            ),
            "velocity_stage": parameters.force_stage(
                self._v_prev, velocity
            ),
            "acceleration_stage": parameters.acceleration_stage(
                self._a_prev, acceleration
            ),
        }

    def _gather_full(self, vector):
        self._scatter_all.scatter(
            vector,
            self._x_all,
            addv=self._PETSc.InsertMode.INSERT_VALUES,
            mode=self._PETSc.ScatterMode.FORWARD,
        )
        return np.asarray(self._x_all.getArray(readonly=True), dtype=float).copy()

    def _monitor(self, _snes, _iteration, residual_norm):
        self._residual_history.append(float(residual_norm))
        if int(_iteration) > 0:
            ksp = self._snes.getKSP()
            self._ksp_iterations_per_solve.append(
                self._integer_query(ksp, "getIterationNumber")
            )
            self._ksp_reasons_per_solve.append(
                self._integer_query(ksp, "getConvergedReason")
            )

    def _ksp_monitor(self, _ksp, iteration, residual_norm):
        iteration = int(iteration)
        if iteration == 0 or not self._ksp_residual_histories:
            self._ksp_residual_histories.append([])
        self._ksp_residual_histories[-1].append(float(residual_norm))

    def _record_pending_ksp_outcome(self, ksp):
        """Pair a final monitored KSP history with its convergence outcome.

        A successful KSP is recorded by the following SNES monitor call.  When
        the final KSP diverges, SNES returns without that monitor call, leaving
        exactly one residual history pending.  Preserve that outcome so the
        application can report the decisive failure before rejecting the step.
        """
        if self.linear_solver_profile == DIRECT_SUPERLU_DIST_PROFILE:
            return
        histories = len(self._ksp_residual_histories)
        outcomes = len(self._ksp_reasons_per_solve)
        if histories == outcomes:
            return
        if histories != outcomes + 1 or len(
            self._ksp_iterations_per_solve
        ) != outcomes:
            raise DistributedSnesSolveError(
                "distributed cardiac KSP diagnostics became inconsistent "
                "before state commit"
            )
        self._ksp_iterations_per_solve.append(
            self._integer_query(ksp, "getIterationNumber")
        )
        self._ksp_reasons_per_solve.append(
            self._integer_query(ksp, "getConvergedReason")
        )

    @staticmethod
    def _integer_query(owner, method, default=0):
        try:
            return int(getattr(owner, method)())
        except (AttributeError, TypeError, ValueError):
            return int(default)

    @staticmethod
    def _float_query(owner, method, default=float("nan")):
        try:
            return float(getattr(owner, method)())
        except (AttributeError, TypeError, ValueError):
            return float(default)

    def _trial_error_record(self, error):
        if error is None:
            return None
        global_element = getattr(error, "global_element_id", None)
        if global_element is None:
            return {
                "kind": "fatal",
                "rank": self.rank,
                "message": (
                    "recoverable deformation error lacks a global element ID: "
                    f"{error}"
                ),
            }
        return {
            "kind": "invalid",
            "rank": self.rank,
            "global_element": int(global_element),
            "message": str(error),
        }

    def _collective_callback_outcome(self, invalid_error, fatal_error, *, jacobian):
        record = self._trial_error_record(invalid_error)
        if fatal_error is not None:
            record = {
                "kind": "fatal",
                "rank": self.rank,
                "message": f"{type(fatal_error).__name__}: {fatal_error}",
            }
        records = self.mpi_comm.allgather(record)
        fatals = [entry for entry in records if entry and entry["kind"] == "fatal"]
        invalids = [
            entry for entry in records if entry and entry["kind"] == "invalid"
        ]
        if jacobian and invalids:
            first = min(
                invalids,
                key=lambda entry: (entry["global_element"], entry["rank"]),
            )
            fatals.append(
                {
                    "kind": "fatal",
                    "rank": first["rank"],
                    "message": (
                        "invalid deformation reached the Jacobian callback: "
                        f"global_element={first['global_element']}"
                    ),
                }
            )
        if fatals:
            first = min(fatals, key=lambda entry: entry["rank"])
            message = (
                "distributed cardiac callback failed on rank "
                f"{first['rank']}: {first['message']}"
            )
            self._callback_fatal_error = message
            raise DistributedAssemblyError(message)
        if invalids:
            first = min(
                invalids,
                key=lambda entry: (entry["global_element"], entry["rank"]),
            )
            message = (
                "distributed cardiac trial rejected: "
                f"global_element={first['global_element']}"
            )
            return message
        return None

    def _element_residual_blocks(self, displacement):
        if len(self.my_gm) == 0:
            ndofel = self.my_gm.shape[1]
            return np.empty((0, ndofel), dtype=float)
        element_u = displacement[self.my_gm]
        element_du = element_u - self._u_prev[self.my_gm]
        residual = self.element_batch.element_residual_batch(
            self.my_coordinates, element_u, element_du
        )
        residual = np.asarray(residual, dtype=float)
        ndofel = self.my_gm.shape[1]
        if residual.shape != (len(self.my_gm), ndofel):
            raise RuntimeError("rank-local element batch returned malformed residuals")
        if not np.all(np.isfinite(residual)):
            raise RuntimeError("rank-local element batch returned non-finite residuals")
        return residual

    def _element_tangent_blocks(self, displacement):
        if len(self.my_gm) == 0:
            ndofel = self.my_gm.shape[1]
            return np.empty((0, ndofel, ndofel), dtype=float)
        element_u = displacement[self.my_gm]
        element_du = element_u - self._u_prev[self.my_gm]
        tangent = self.element_batch.element_tangent_batch(
            self.my_coordinates, element_u, element_du
        )
        tangent = np.asarray(tangent, dtype=float)
        ndofel = self.my_gm.shape[1]
        if tangent.shape != (len(self.my_gm), ndofel, ndofel):
            raise RuntimeError("rank-local element batch returned malformed tangents")
        if not np.all(np.isfinite(tangent)):
            raise RuntimeError("rank-local element batch returned non-finite tangents")
        return tangent

    def _boundary_residuals(self, displacement):
        active = self._active
        if active is None:
            raise RuntimeError("distributed residual callback has no active step")
        if self.generalized_alpha is None:
            boundary_displacement = displacement
            boundary_velocity = None
            boundary_time = active["t"]
        else:
            trial = self._generalized_alpha_trial(displacement)
            boundary_displacement = trial["displacement_stage"]
            boundary_velocity = trial["velocity_stage"]
            boundary_time = active["stage_time"]
        robin_values = None
        if self.robin is not None:
            if self.generalized_alpha is None:
                # Retain the historical BE arithmetic and Robin-owned history.
                robin_values = self.robin.Kmat @ displacement + self.robin.Cmat @ (
                    displacement - np.asarray(self.robin.u_prev, dtype=float)
                ) / active["dt"]
            else:
                robin_values = (
                    self.robin.Kmat @ boundary_displacement
                    + self.robin.Cmat @ boundary_velocity
                )
            if not np.all(np.isfinite(robin_values)):
                raise RuntimeError("distributed Robin residual is non-finite")
        pressure_residual = None
        if self.pressure is not None:
            pressure_residual = self.pressure.residual(
                boundary_displacement, None, boundary_time, active["dt"]
            )
            if not np.all(np.isfinite(pressure_residual.values)):
                raise RuntimeError("distributed follower-pressure residual is non-finite")
        return robin_values, pressure_residual

    def _boundary_tangents(self, displacement):
        active = self._active
        if active is None:
            raise RuntimeError("distributed Jacobian callback has no active step")
        robin_tangent = None
        if self.robin is not None:
            if self.generalized_alpha is None:
                robin_tangent = (
                    self.robin.Kmat + self.robin.Cmat / active["dt"]
                ).tocsr()
            else:
                displacement_scale = (
                    self.generalized_alpha.force_displacement_tangent()
                )
                velocity_scale = (
                    self.generalized_alpha.force_velocity_tangent(active["dt"])
                )
                robin_tangent = (
                    displacement_scale * self.robin.Kmat
                    + velocity_scale * self.robin.Cmat
                ).tocsr()
            if not np.all(np.isfinite(robin_tangent.data)):
                raise RuntimeError("distributed Robin tangent is non-finite")
        pressure_tangent = None
        if self.pressure is not None:
            if self.generalized_alpha is None:
                pressure_displacement = displacement
                pressure_scale = 1.0
                boundary_time = active["t"]
            else:
                trial = self._generalized_alpha_trial(displacement)
                pressure_displacement = trial["displacement_stage"]
                pressure_scale = (
                    self.generalized_alpha.force_displacement_tangent()
                )
                boundary_time = active["stage_time"]
            pressure_tangent = self.pressure.tangent(
                pressure_displacement, None, boundary_time, active["dt"]
            )
            if pressure_scale != 1.0:
                pressure_tangent = type(pressure_tangent)(
                    pressure_tangent.rows,
                    pressure_tangent.cols,
                    pressure_scale * pressure_tangent.values,
                )
            if not np.all(np.isfinite(pressure_tangent.values)):
                raise RuntimeError("distributed follower-pressure tangent is non-finite")
        return robin_tangent, pressure_tangent

    def _form_function(self, _snes, X, F):
        started = time.perf_counter()
        try:
            active = self._active
            if active is None:
                raise RuntimeError("distributed residual callback has no active step")
            displacement = self._gather_full(X)
            invalid_error = None
            fatal_error = None
            element_residual = None
            try:
                element_residual = self._element_residual_blocks(displacement)
            except InvalidDeformationError as error:
                invalid_error = error
            except Exception as error:  # synchronized below before propagation
                fatal_error = error
            robin_values = pressure_residual = None
            try:
                robin_values, pressure_residual = self._boundary_residuals(
                    displacement
                )
            except Exception as error:  # synchronized below before propagation
                if fatal_error is None:
                    fatal_error = error
            invalid_message = self._collective_callback_outcome(
                invalid_error, fatal_error, jacobian=False
            )
            if invalid_message is not None:
                F.set(np.inf)
                F.assemble()
                self._function_domain_rejections += 1
                self._last_function_domain_error = invalid_message
                return

            PETSc = self._PETSc
            F.set(0.0)
            for local_index in range(len(self.my_gm)):
                F.setValues(
                    self._rows_all[local_index],
                    element_residual[local_index],
                    addv=PETSc.InsertMode.ADD_VALUES,
                )
            if self.generalized_alpha is None:
                if self._consistent_mass:
                    inertia = self.mass.action(
                        displacement - active["predictor"]
                    ) / active["dt"] ** 2
                else:
                    owned_mass = self.mass[self._rs : self._re]
                    inertia = (owned_mass / active["dt"] ** 2) * (
                        displacement[self._rs : self._re]
                        - active["predictor"][self._rs : self._re]
                    )
            else:
                inertial_vector = self._generalized_alpha_trial(
                    displacement
                )["acceleration_stage"]
                if self._consistent_mass:
                    inertia = self.mass.action(inertial_vector)
                else:
                    owned_mass = self.mass[self._rs : self._re]
                    inertia = owned_mass * inertial_vector[self._rs : self._re]
            if len(self._owned):
                F.setValues(
                    self._owned, inertia, addv=PETSc.InsertMode.ADD_VALUES
                )
            if robin_values is not None:
                robin_dofs = np.asarray(self.robin.dofs, dtype=np.int64)
                owned_robin = robin_dofs[
                    (robin_dofs >= self._rs) & (robin_dofs < self._re)
                ].astype(PETSc.IntType)
                if len(owned_robin):
                    F.setValues(
                        owned_robin,
                        robin_values[owned_robin],
                        addv=PETSc.InsertMode.ADD_VALUES,
                    )
            if pressure_residual is not None:
                gdofs = np.asarray(pressure_residual.gdofs, dtype=np.int64)
                values = np.asarray(pressure_residual.values, dtype=float)
                owned_pressure = (gdofs >= self._rs) & (gdofs < self._re)
                if np.any(owned_pressure):
                    F.setValues(
                        gdofs[owned_pressure].astype(PETSc.IntType),
                        values[owned_pressure],
                        addv=PETSc.InsertMode.ADD_VALUES,
                    )
            F.assemble()
        finally:
            self._step_assembly_seconds += time.perf_counter() - started

    def _form_jacobian(self, _snes, X, J, _P):
        started = time.perf_counter()
        try:
            if self._active is None:
                raise RuntimeError("distributed Jacobian callback has no active step")
            displacement = self._gather_full(X)
            invalid_error = None
            fatal_error = None
            element_tangent = None
            try:
                element_tangent = self._element_tangent_blocks(displacement)
            except InvalidDeformationError as error:
                invalid_error = error
            except Exception as error:  # synchronized below before propagation
                fatal_error = error
            robin_tangent = pressure_tangent = None
            try:
                robin_tangent, pressure_tangent = self._boundary_tangents(
                    displacement
                )
            except Exception as error:  # synchronized below before propagation
                if fatal_error is None:
                    fatal_error = error
            self._collective_callback_outcome(
                invalid_error, fatal_error, jacobian=True
            )

            PETSc = self._PETSc
            J.zeroEntries()
            for local_index in range(len(self.my_gm)):
                J.setValues(
                    self._rows_all[local_index],
                    self._rows_all[local_index],
                    element_tangent[local_index],
                    addv=PETSc.InsertMode.ADD_VALUES,
                )
            if self._consistent_mass:
                mass_matrix = self.mass.matrix
                if self.generalized_alpha is None:
                    scale = 1.0 / self._active["dt"] ** 2
                else:
                    scale = self.generalized_alpha.inertia_tangent(
                        self._active["dt"]
                    )
                for local_row, global_row in enumerate(self._owned):
                    start = mass_matrix.indptr[local_row]
                    stop = mass_matrix.indptr[local_row + 1]
                    columns = mass_matrix.indices[start:stop].astype(
                        PETSc.IntType, copy=False
                    )
                    values = mass_matrix.data[start:stop] * scale
                    if len(columns):
                        J.setValues(
                            [int(global_row)],
                            columns,
                            values.reshape(1, -1),
                            addv=PETSc.InsertMode.ADD_VALUES,
                        )
            else:
                if self.generalized_alpha is None:
                    diagonal = (
                        self.mass[self._rs : self._re]
                        / self._active["dt"] ** 2
                    )
                else:
                    mass_scale = self.generalized_alpha.inertia_tangent(
                        self._active["dt"]
                    )
                    diagonal = self.mass[self._rs : self._re] * mass_scale
                for global_dof, value in zip(self._owned, diagonal):
                    if value != 0.0:
                        J.setValue(
                            int(global_dof),
                            int(global_dof),
                            float(value),
                            addv=PETSc.InsertMode.ADD_VALUES,
                        )
            if robin_tangent is not None:
                robin_dofs = np.asarray(self.robin.dofs, dtype=np.int64)
                owned_robin = robin_dofs[
                    (robin_dofs >= self._rs) & (robin_dofs < self._re)
                ]
                for row in owned_robin:
                    start = robin_tangent.indptr[row]
                    stop = robin_tangent.indptr[row + 1]
                    columns = robin_tangent.indices[start:stop].astype(PETSc.IntType)
                    values = robin_tangent.data[start:stop]
                    if len(columns):
                        J.setValues(
                            [int(row)],
                            columns,
                            values.reshape(1, -1),
                            addv=PETSc.InsertMode.ADD_VALUES,
                        )
            if pressure_tangent is not None:
                rows = np.asarray(pressure_tangent.rows, dtype=np.int64)
                columns = np.asarray(pressure_tangent.cols, dtype=np.int64)
                values = np.asarray(pressure_tangent.values, dtype=float)
                owned_pressure = (rows >= self._rs) & (rows < self._re)
                for row, column, value in zip(
                    rows[owned_pressure],
                    columns[owned_pressure],
                    values[owned_pressure],
                ):
                    J.setValue(
                        int(row),
                        int(column),
                        float(value),
                        addv=PETSc.InsertMode.ADD_VALUES,
                    )
            J.assemble()
        finally:
            self._step_assembly_seconds += time.perf_counter() - started

    def _set_distributed_x(self, values):
        values = np.asarray(values, dtype=float)
        if values.shape != (self.ndof,) or not np.all(np.isfinite(values)):
            raise DistributedSnesSolveError(
                "distributed cardiac SNES has an invalid initial displacement"
            )
        local = self._x.getArray(readonly=False)
        local[:] = values[self._rs : self._re]
        self._x.assemble()

    def _independent_residual_norm(self):
        if hasattr(self.element_batch, "clear_cache"):
            self.element_batch.clear_cache()
        self._form_function(self._snes, self._x, self._f)
        norm = float(self._f.norm())
        return norm

    def solve_step(self, *, t, dt, dirichlet=None):
        """Solve, independently accept, and transactionally commit one step."""
        if self._closed:
            raise RuntimeError("distributed cardiac solver is closed")
        if dirichlet:
            raise NotImplementedError(
                "the distributed Step 0 companion supports only empty Dirichlet data"
            )
        if not np.isfinite(t):
            raise ValueError("time must be finite")
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt must be finite and positive")
        t = float(t)
        dt = float(dt)
        if self.generalized_alpha is not None and not np.isclose(
            dt,
            self.material_dt,
            rtol=0.0,
            atol=1.0e-15 * max(1.0, abs(self.material_dt)),
        ):
            raise ValueError(
                "generalized-alpha solve dt differs from compiled material dt"
            )
        predictor = self.predictor(dt)
        stage_time = (
            t
            if self.generalized_alpha is None
            else self.generalized_alpha.load_time(t, dt)
        )
        self._active = {
            "t": t,
            "stage_time": stage_time,
            "dt": dt,
            "predictor": predictor,
        }
        self._residual_history = []
        self._ksp_iterations_per_solve = []
        self._ksp_reasons_per_solve = []
        self._ksp_residual_histories = []
        self._step_assembly_seconds = 0.0
        self._function_domain_rejections = 0
        self._last_function_domain_error = None
        self._callback_fatal_error = None
        self.last_diagnostics = None
        self._set_distributed_x(predictor)

        try:
            initial_norm = self._independent_residual_norm()
        except Exception as error:
            raise DistributedSnesSolveError(
                "distributed cardiac SNES initial residual assembly failed before "
                "state commit"
            ) from error
        if not np.isfinite(initial_norm):
            raise DistributedSnesSolveError(
                "distributed cardiac SNES has a non-finite initial residual"
            )

        started = time.perf_counter()
        try:
            self._snes.solve(None, self._x)
        except Exception as error:
            message = "distributed cardiac SNES failed before state commit"
            if self._callback_fatal_error:
                message += f": {self._callback_fatal_error}"
            raise DistributedSnesSolveError(message) from error
        local_solve_seconds = time.perf_counter() - started

        snes_reason = self._integer_query(self._snes, "getConvergedReason")
        nonlinear_iterations = self._integer_query(
            self._snes, "getIterationNumber"
        )
        linear_iterations = self._integer_query(
            self._snes, "getLinearSolveIterations"
        )
        ksp = self._snes.getKSP()
        ksp_reason = self._integer_query(ksp, "getConvergedReason")
        self._record_pending_ksp_outcome(ksp)
        petsc_function_norm = self._float_query(self._snes, "getFunctionNorm")
        displacement = self._gather_full(self._x)

        try:
            final_norm = self._independent_residual_norm()
        except Exception as error:
            raise DistributedSnesSolveError(
                "distributed cardiac independent final residual assembly failed "
                "before state commit"
            ) from error
        threshold = max(self.settings.atol, self.settings.rtol * initial_norm)
        assembly_seconds = float(
            self.mpi_comm.allreduce(
                self._step_assembly_seconds, op=self._MPI.MAX
            )
        )
        solve_seconds = float(
            self.mpi_comm.allreduce(local_solve_seconds, op=self._MPI.MAX)
        )
        diagnostics = DistributedSnesStepDiagnostics(
            time=t,
            dt=dt,
            ranks=self.size,
            initial_residual_norm=float(initial_norm),
            final_residual_norm=float(final_norm),
            residual_acceptance_threshold=float(threshold),
            petsc_function_norm=float(petsc_function_norm),
            snes_converged_reason=snes_reason,
            ksp_converged_reason=ksp_reason,
            nonlinear_iterations=nonlinear_iterations,
            linear_iterations=linear_iterations,
            residual_history=tuple(self._residual_history),
            assembly_seconds=assembly_seconds,
            solve_seconds=solve_seconds,
            function_domain_rejections=int(self._function_domain_rejections),
            last_function_domain_error=self._last_function_domain_error,
            ksp_iterations_per_solve=tuple(self._ksp_iterations_per_solve),
            ksp_reasons_per_solve=tuple(self._ksp_reasons_per_solve),
            ksp_residual_histories=tuple(
                tuple(history) for history in self._ksp_residual_histories
            ),
        )
        self.last_diagnostics = diagnostics

        rejection = None
        if displacement.shape != (self.ndof,) or not np.all(
            np.isfinite(displacement)
        ):
            rejection = "distributed cardiac SNES returned an invalid displacement"
        elif snes_reason <= 0:
            rejection = (
                "distributed cardiac SNES did not report convergence before state "
                f"commit: reason={snes_reason}"
            )
        elif ksp_reason < 0 or (ksp_reason == 0 and linear_iterations > 0):
            rejection = (
                "distributed cardiac KSP did not report convergence before state "
                f"commit: reason={ksp_reason}, linear_iterations={linear_iterations}"
            )
        elif linear_iterations > 0 and (
            not diagnostics.ksp_reasons_per_solve
            or any(reason <= 0 for reason in diagnostics.ksp_reasons_per_solve)
        ):
            rejection = (
                "at least one distributed cardiac KSP solve did not report "
                "positive convergence before state commit"
            )
        finite_diagnostics = (
            diagnostics.time,
            diagnostics.dt,
            diagnostics.initial_residual_norm,
            diagnostics.final_residual_norm,
            diagnostics.residual_acceptance_threshold,
            diagnostics.petsc_function_norm,
            diagnostics.assembly_seconds,
            diagnostics.solve_seconds,
            *diagnostics.residual_history,
            *(
                value
                for history in diagnostics.ksp_residual_histories
                for value in history
            ),
        )
        if rejection is None and not np.all(
            np.isfinite(np.asarray(finite_diagnostics, dtype=float))
        ):
            rejection = (
                "distributed cardiac SNES produced non-finite diagnostics before "
                "state commit"
            )
        if rejection is None and not (
            np.isfinite(final_norm) and final_norm < threshold
        ):
            rejection = (
                "distributed cardiac SNES did not meet the independent residual "
                f"rule before state commit: |R|={final_norm:.6e}, "
                f"required < {threshold:.6e}"
            )
        if diagnostics.function_domain_rejections == 0:
            domain_record_ok = diagnostics.last_function_domain_error is None
        else:
            domain_record_ok = isinstance(
                diagnostics.last_function_domain_error, str
            ) and bool(diagnostics.last_function_domain_error)
        if rejection is None and not domain_record_ok:
            rejection = (
                "distributed cardiac SNES produced inconsistent invalid-trial "
                "diagnostics before state commit"
            )

        # Every rank must reach the same commit decision.  Carry the first
        # rank's rejection text so a future rank-local check cannot permit a
        # partial material/history commit.
        decisions = self.mpi_comm.allgather(rejection)
        collective_rejections = [entry for entry in decisions if entry is not None]
        if collective_rejections:
            raise DistributedSnesSolveError(
                collective_rejections[0], diagnostics
            )

        # The independent final residual evaluation above refreshed the
        # compiled material's trial state at the accepted displacement.  The
        # pressure tangent is intentionally unnecessary for this acceptance
        # check.  Only now may any history move forward.
        self.element_batch.commit()
        if self.generalized_alpha is None:
            old_u = self._u_prev.copy()
            self._v_prev = (displacement - old_u) / dt
            self._u_prev = displacement.copy()
        else:
            accepted = self._generalized_alpha_trial(displacement)
            self._u_prev = displacement.copy()
            self._v_prev = accepted["velocity"].copy()
            self._a_prev = accepted["acceleration"].copy()
        if self.robin is not None:
            self.robin.commit(displacement, None, t, dt)
        self._accepted_steps += 1
        return displacement, diagnostics

    def local_element_displacement(self, displacement=None):
        """Gather the accepted or supplied global displacement on local elements."""
        values = self._u_prev if displacement is None else np.asarray(
            displacement, dtype=float
        )
        if values.shape != (self.ndof,):
            raise ValueError("displacement must match ndof")
        return values[self.my_gm]

    def close(self):
        """Destroy all PETSc resources; idempotent before or after solves."""
        if self._closed:
            return
        for resource in (
            self._snes,
            self._near_nullspace,
            self._coordinate_vector,
            self._jacobian,
            self._f,
            self._scatter_all,
            self._x_all,
            self._x,
        ):
            if resource is not None:
                resource.destroy()
        options = self._PETSc.Options()
        for name in self._petsc_options:
            try:
                del options[name]
            except KeyError:
                pass
        self._petsc_options.clear()
        self._active = None
        self._closed = True
