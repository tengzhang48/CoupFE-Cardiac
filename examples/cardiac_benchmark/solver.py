"""Application-level convergence checks around CoupFE's low-level solvers.

CoupFE's compact ``newton_solve`` returns an iterate and iteration count; the
count is explicitly not a convergence flag.  These helpers independently
check the unconstrained residual before allowing any cardiac operator to
commit state.  The optional PETSc SNES path is also owned here, rather than in
Core, because its solver configuration and acceptance policy are part of the
cardiac experiment.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import importlib.metadata
import time

import numpy as np

from coupfe.assembly.assemble import assemble_residual, assemble_tangent, newton_solve
from coupfe.operators.base import Residual, Tangent
from local_pressure import InvalidDeformationError


# Core-Newton final-residual oracle: intentionally tighter than the PETSc
# guard's atol=1e-10 (adapter-matched acceptance in PetscSnesSettings).
# The two guards serve different roles: this is a strict final-residual
# check for the Core path, not a shared tolerance contract with PETSc.
_ABSOLUTE_RESIDUAL_TOL = 1.0e-14


@dataclass(frozen=True)
class PetscSnesSettings:
    """Settings recovered from the 2026-06-27 Case B development adapter."""

    snes_type: str = "newtonls"
    line_search_type: str = "bt"
    ksp_type: str = "preonly"
    pc_type: str = "lu"
    rtol: float = 1.0e-9
    atol: float = 1.0e-10
    stol: float = 1.0e-12
    max_it: int = 60


RECOVERED_PETSC_SNES_SETTINGS = PetscSnesSettings()


@dataclass(frozen=True)
class SnesStepDiagnostics:
    """Serializable nonlinear diagnostics for one accepted time step."""

    time: float
    dt: float
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

    def as_dict(self):
        return asdict(self)


class SnesSolveError(RuntimeError):
    """Fail-closed SNES error carrying diagnostics when they are available."""

    def __init__(self, message, diagnostics=None):
        super().__init__(message)
        self.diagnostics = diagnostics


def _residual_norm(operators, U, state, ndof, dirichlet, *, t, dt):
    """Return the finite unconstrained residual norm for one trial state."""
    displacement = np.asarray(U, dtype=float).copy()
    if displacement.shape != (ndof,) or not np.all(np.isfinite(displacement)):
        return float("inf")
    constrained = np.fromiter(sorted(dirichlet), dtype=int) if dirichlet else np.empty(0, int)
    if len(constrained):
        if constrained[0] < 0 or constrained[-1] >= ndof:
            raise ValueError("Dirichlet degree of freedom is outside the global system")
        values = np.array([dirichlet[int(dof)] for dof in constrained], dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError("Dirichlet values must be finite")
        displacement[constrained] = values
    residual, _ = assemble_residual(operators, displacement, state, t, dt, ndof)
    residual = np.asarray(residual, dtype=float)
    if residual.shape != (ndof,):
        raise RuntimeError(
            f"cardiac residual has shape {residual.shape}; expected {(ndof,)}"
        )
    if len(constrained):
        residual[constrained] = 0.0
    if not np.all(np.isfinite(residual)):
        return float("inf")
    return float(np.linalg.norm(residual))


class _CommitGuard:
    """Zero operator whose commit runs before every physical-state commit."""

    def __init__(
        self,
        operators,
        state,
        ndof,
        dirichlet,
        *,
        t,
        dt,
        initial_norm,
        rtol,
    ):
        self.operators = tuple(operators)
        self.state = state
        self.ndof = int(ndof)
        self.dirichlet = dict(dirichlet)
        self.t = float(t)
        self.dt = float(dt)
        self.initial_norm = float(initial_norm)
        self.rtol = float(rtol)

    def residual(self, U, state, t, dt):
        return Residual(np.empty(0, dtype=int), np.empty(0, dtype=float))

    def tangent(self, U, state, t, dt):
        empty = np.empty(0, dtype=int)
        return Tangent(empty, empty, np.empty(0, dtype=float))

    def commit(self, U, state, t, dt):
        final_norm = _residual_norm(
            self.operators,
            U,
            self.state,
            self.ndof,
            self.dirichlet,
            t=self.t,
            dt=self.dt,
        )
        threshold = max(
            _ABSOLUTE_RESIDUAL_TOL,
            self.rtol * max(self.initial_norm, 1.0e-300),
        )
        if not np.isfinite(final_norm) or final_norm > threshold:
            raise RuntimeError(
                "cardiac Newton solve did not converge before state commit: "
                f"|R_free|={final_norm:.6e}, required <= {threshold:.6e}"
            )
        return state


def checked_newton_solve(
    operators,
    U0,
    state,
    ndof,
    dirichlet,
    *,
    t=1.0,
    dt=1.0,
    rtol=1.0e-9,
    maxit=60,
):
    """Run Core Newton and reject a nonconverged iterate before state commit."""
    if not np.isfinite(rtol) or not 0.0 < rtol < 1.0:
        raise ValueError("rtol must be finite and satisfy 0 < rtol < 1")
    if not isinstance(maxit, int) or maxit < 1:
        raise ValueError("maxit must be a positive integer")
    physical_operators = tuple(operators)
    initial_norm = _residual_norm(
        physical_operators, U0, state, ndof, dirichlet, t=t, dt=dt
    )
    if not np.isfinite(initial_norm):
        raise RuntimeError("cardiac Newton solve has a non-finite initial residual")
    guard = _CommitGuard(
        physical_operators,
        state,
        ndof,
        dirichlet,
        t=t,
        dt=dt,
        initial_norm=initial_norm,
        rtol=rtol,
    )
    U, committed, nit = newton_solve(
        (guard, *physical_operators),
        U0,
        state,
        ndof,
        dirichlet,
        t=t,
        dt=dt,
        rtol=rtol,
        maxit=maxit,
    )
    return U, committed[1:], nit


def _load_petsc():
    """Import and initialize petsc4py only when the SNES path is selected."""
    try:
        import petsc4py
    except ImportError as error:
        raise RuntimeError(
            "--nonlinear-solver petsc-snes requires the optional petsc4py "
            "dependency (install this application with the 'mpi' extra)"
        ) from error
    petsc4py.init([])
    from petsc4py import PETSc

    return PETSc


class PetscSnesSolver:
    """Persistent, application-owned PETSc SNES solver for serial cardiac runs.

    The matrix, vectors, and SNES object are created lazily on the first step
    and reused for the rest of this solver instance.  They are never cached
    globally.  The initial implementation intentionally supports only the
    unconstrained cardiac benchmark: nonempty Dirichlet data fail before any
    PETSc callback or physical-state commit.
    """

    def __init__(self, settings=RECOVERED_PETSC_SNES_SETTINGS):
        if not isinstance(settings, PetscSnesSettings):
            raise TypeError("settings must be a PetscSnesSettings instance")
        self.settings = settings
        self._petsc = None
        self._ndof = None
        self._snes = None
        self._x = None
        self._f = None
        self._jacobian = None
        self._active = None
        self._residual_history = []
        self._step_assembly_seconds = 0.0
        self._line_search_configuration = "not-initialized"
        self._function_domain_rejection_api = "nonfinite residual for PETSc BT"
        self._function_domain_rejections = 0
        self._last_function_domain_error = None
        self.last_diagnostics = None

    def configuration(self):
        """Return JSON-compatible provenance for the exact solver settings."""
        configuration = {
            "name": "petsc-snes",
            "settings_source": "recovered 2026-06-27 Case B development adapter",
            "matrix_scope": "one solver instance per application run",
            "dirichlet_support": "none",
            "function_domain_rejection_api": self._function_domain_rejection_api,
            **asdict(self.settings),
        }
        try:
            configuration["petsc4py_version"] = importlib.metadata.version("petsc4py")
        except importlib.metadata.PackageNotFoundError:
            configuration["petsc4py_version"] = "unavailable"
        if self._petsc is not None:
            try:
                version = self._petsc.Sys.getVersion()
            except (AttributeError, TypeError):
                version = None
            if version is not None:
                configuration["petsc_version"] = ".".join(
                    str(component) for component in version
                )
            configuration["line_search_configuration_api"] = (
                self._line_search_configuration
            )
            try:
                factor_solver = (
                    self._snes.getKSP().getPC().getFactorSolverType()
                )
            except (AttributeError, TypeError, RuntimeError):
                factor_solver = None
            configuration["factor_solver_type"] = (
                str(factor_solver) if factor_solver else "unknown"
            )
        return configuration

    @staticmethod
    def _validate_tangent(matrix, ndof):
        matrix = matrix.tocsr()
        if matrix.shape != (ndof, ndof):
            raise RuntimeError(
                f"cardiac tangent has shape {matrix.shape}; expected {(ndof, ndof)}"
            )
        if not np.all(np.isfinite(matrix.data)):
            raise RuntimeError("cardiac tangent contains non-finite values")
        return matrix

    def _form_function(self, _snes, X, F):
        started = time.perf_counter()
        try:
            active = self._active
            if active is None:
                raise RuntimeError("PETSc residual callback has no active time step")
            displacement = np.asarray(X.getArray(readonly=True), dtype=float)
            residual, _ = assemble_residual(
                active["operators"],
                displacement,
                active["state"],
                active["t"],
                active["dt"],
                active["ndof"],
            )
            residual = np.asarray(residual, dtype=float)
            if residual.shape != (active["ndof"],):
                raise RuntimeError(
                    f"cardiac residual has shape {residual.shape}; "
                    f"expected {(active['ndof'],)}"
                )
            target = F.getArray(readonly=False)
            target[:] = residual
        except InvalidDeformationError as error:
            # A finite-strain trial outside the element domain is not a Python
            # callback failure. A nonfinite trial residual makes PETSc's bt
            # line search reject and shorten the step, including on petsc4py
            # 3.18 where no public function-domain-error method is exposed.
            target = F.getArray(readonly=False)
            target[:] = np.inf
            self._function_domain_rejections += 1
            self._last_function_domain_error = str(error)
        finally:
            self._step_assembly_seconds += time.perf_counter() - started

    def _form_jacobian(self, _snes, X, J, _P):
        started = time.perf_counter()
        try:
            active = self._active
            if active is None:
                raise RuntimeError("PETSc tangent callback has no active time step")
            displacement = np.asarray(X.getArray(readonly=True), dtype=float)
            matrix = self._validate_tangent(
                assemble_tangent(
                    active["operators"],
                    displacement,
                    active["state"],
                    active["t"],
                    active["dt"],
                    active["ndof"],
                ),
                active["ndof"],
            )
            int_type = self._petsc.IntType
            J.zeroEntries()
            J.setValuesCSR(
                np.asarray(matrix.indptr, dtype=int_type),
                np.asarray(matrix.indices, dtype=int_type),
                matrix.data,
            )
            J.assemble()
        finally:
            self._step_assembly_seconds += time.perf_counter() - started

    def _monitor(self, _snes, _iteration, residual_norm):
        self._residual_history.append(float(residual_norm))

    def _ensure_context(self, operators, U0, state, ndof, *, t, dt):
        if self._snes is not None:
            if ndof != self._ndof:
                raise RuntimeError(
                    "one PetscSnesSolver instance cannot be reused with a different ndof"
                )
            return

        self._petsc = _load_petsc()
        PETSc = self._petsc
        try:
            world_size = PETSc.COMM_WORLD.getSize()
        except AttributeError:
            world_size = 1
        if world_size != 1:
            raise RuntimeError(
                "the cardiac petsc-snes path is serial; run it with one MPI rank"
            )

        initial_tangent = self._validate_tangent(
            assemble_tangent(operators, U0, state, t, dt, ndof), ndof
        )
        int_type = PETSc.IntType
        csr = (
            np.asarray(initial_tangent.indptr, dtype=int_type),
            np.asarray(initial_tangent.indices, dtype=int_type),
            initial_tangent.data,
        )
        self._jacobian = PETSc.Mat().createAIJ(
            (ndof, ndof), csr=csr, comm=PETSc.COMM_SELF
        )
        self._jacobian.setOption(
            PETSc.Mat.Option.NEW_NONZERO_ALLOCATION_ERR, False
        )
        self._x = PETSc.Vec().createSeq(ndof, comm=PETSc.COMM_SELF)
        self._f = self._x.duplicate()
        self._snes = PETSc.SNES().create(comm=PETSc.COMM_SELF)
        self._snes.setFunction(self._form_function, self._f)
        self._snes.setJacobian(self._form_jacobian, self._jacobian)
        self._snes.setMonitor(self._monitor)
        self._snes.setType(self.settings.snes_type)
        if hasattr(self._snes, "getLineSearch"):
            self._snes.getLineSearch().setType(self.settings.line_search_type)
            self._line_search_configuration = "SNES.getLineSearch"
        else:
            # petsc4py 3.18 does not expose SNES.getLineSearch. A temporary,
            # namespaced PETSc option configures the identical bt line search;
            # remove it immediately so no process-global option leaks into a
            # later solver instance.
            prefix = f"coupfe_cardiac_{id(self):x}_"
            self._snes.setOptionsPrefix(prefix)
            options = PETSc.Options()
            option_name = f"{prefix}snes_linesearch_type"
            options[option_name] = self.settings.line_search_type
            try:
                self._snes.setFromOptions()
            finally:
                del options[option_name]
            self._line_search_configuration = "namespaced PETSc option"
        ksp = self._snes.getKSP()
        ksp.setType(self.settings.ksp_type)
        ksp.getPC().setType(self.settings.pc_type)
        self._snes.setTolerances(
            rtol=self.settings.rtol,
            atol=self.settings.atol,
            stol=self.settings.stol,
            max_it=self.settings.max_it,
        )
        self._ndof = int(ndof)

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

    def solve(
        self,
        operators,
        U0,
        state,
        ndof,
        dirichlet,
        *,
        t=1.0,
        dt=1.0,
    ):
        """Solve, independently accept, then commit one cardiac time step."""
        if dirichlet:
            raise NotImplementedError(
                "cardiac petsc-snes currently supports only empty Dirichlet data"
            )
        if not isinstance(ndof, int) or ndof < 1:
            raise ValueError("ndof must be a positive integer")
        if not np.isfinite(t):
            raise ValueError("time must be finite")
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("time increment must be finite and positive")
        displacement0 = np.asarray(U0, dtype=float)
        if displacement0.shape != (ndof,) or not np.all(np.isfinite(displacement0)):
            raise SnesSolveError(
                "cardiac PETSc SNES solve has an invalid initial displacement"
            )
        physical_operators = tuple(operators)
        initial_norm = _residual_norm(
            physical_operators, displacement0, state, ndof, {}, t=t, dt=dt
        )
        if not np.isfinite(initial_norm):
            raise SnesSolveError(
                "cardiac PETSc SNES solve has a non-finite initial residual"
            )

        self.last_diagnostics = None
        self._function_domain_rejections = 0
        self._last_function_domain_error = None
        self._active = {
            "operators": physical_operators,
            "state": state,
            "t": float(t),
            "dt": float(dt),
            "ndof": int(ndof),
        }
        self._residual_history = []
        self._step_assembly_seconds = 0.0
        self._ensure_context(
            physical_operators, displacement0, state, ndof, t=t, dt=dt
        )
        x_values = self._x.getArray(readonly=False)
        x_values[:] = displacement0

        started = time.perf_counter()
        try:
            self._snes.solve(None, self._x)
        except Exception as error:
            raise SnesSolveError(
                "cardiac PETSc SNES solve failed before state commit"
            ) from error
        solve_seconds = time.perf_counter() - started

        displacement = np.asarray(
            self._x.getArray(readonly=True), dtype=float
        ).copy()
        snes_reason = self._integer_query(self._snes, "getConvergedReason")
        nonlinear_iterations = self._integer_query(
            self._snes, "getIterationNumber"
        )
        linear_iterations = self._integer_query(
            self._snes, "getLinearSolveIterations"
        )
        ksp = self._snes.getKSP()
        ksp_reason = self._integer_query(ksp, "getConvergedReason")
        petsc_function_norm = self._float_query(self._snes, "getFunctionNorm")
        final_norm = _residual_norm(
            physical_operators, displacement, state, ndof, {}, t=t, dt=dt
        )
        threshold = max(
            self.settings.atol, self.settings.rtol * initial_norm
        )
        diagnostics = SnesStepDiagnostics(
            time=float(t),
            dt=float(dt),
            initial_residual_norm=float(initial_norm),
            final_residual_norm=float(final_norm),
            residual_acceptance_threshold=float(threshold),
            petsc_function_norm=float(petsc_function_norm),
            snes_converged_reason=snes_reason,
            ksp_converged_reason=ksp_reason,
            nonlinear_iterations=nonlinear_iterations,
            linear_iterations=linear_iterations,
            residual_history=tuple(self._residual_history),
            assembly_seconds=float(self._step_assembly_seconds),
            solve_seconds=float(solve_seconds),
            function_domain_rejections=int(self._function_domain_rejections),
            last_function_domain_error=self._last_function_domain_error,
        )
        self.last_diagnostics = diagnostics

        if displacement.shape != (ndof,) or not np.all(np.isfinite(displacement)):
            raise SnesSolveError(
                "cardiac PETSc SNES returned a non-finite or malformed displacement",
                diagnostics,
            )
        if snes_reason <= 0:
            raise SnesSolveError(
                "cardiac PETSc SNES did not report convergence before state commit: "
                f"reason={snes_reason}",
                diagnostics,
            )
        if ksp_reason < 0 or (linear_iterations > 0 and ksp_reason == 0):
            raise SnesSolveError(
                "cardiac PETSc linear solve did not report convergence before "
                "state commit: "
                f"reason={ksp_reason}",
                diagnostics,
            )
        finite_petsc_diagnostics = (
            diagnostics.time,
            diagnostics.dt,
            diagnostics.initial_residual_norm,
            diagnostics.residual_acceptance_threshold,
            diagnostics.petsc_function_norm,
            diagnostics.assembly_seconds,
            diagnostics.solve_seconds,
            *diagnostics.residual_history,
        )
        if diagnostics.function_domain_rejections < 0:
            raise SnesSolveError(
                "cardiac PETSc SNES produced invalid domain-rejection diagnostics "
                "before state commit",
                diagnostics,
            )
        if diagnostics.function_domain_rejections == 0:
            domain_message_is_consistent = (
                diagnostics.last_function_domain_error is None
            )
        else:
            domain_message_is_consistent = isinstance(
                diagnostics.last_function_domain_error, str
            ) and bool(diagnostics.last_function_domain_error)
        if not domain_message_is_consistent:
            raise SnesSolveError(
                "cardiac PETSc SNES produced inconsistent domain-rejection "
                "diagnostics before state commit",
                diagnostics,
            )
        if not np.all(
            np.isfinite(np.asarray(finite_petsc_diagnostics, dtype=float))
        ):
            raise SnesSolveError(
                "cardiac PETSc SNES produced non-finite solver diagnostics before "
                "state commit",
                diagnostics,
            )
        if not np.isfinite(final_norm) or final_norm > threshold:
            raise SnesSolveError(
                "cardiac PETSc SNES did not meet the recovered residual rule before "
                f"state commit: |R|={final_norm:.6e}, required <= {threshold:.6e}",
                diagnostics,
            )

        committed = [
            operator.commit(displacement, state, t, dt)
            for operator in physical_operators
        ]
        return displacement, committed, diagnostics

    def close(self):
        """Destroy owned PETSc objects; safe to call before or after first use."""
        for name in ("_snes", "_f", "_x", "_jacobian"):
            resource = getattr(self, name)
            if resource is not None:
                resource.destroy()
                setattr(self, name, None)
        self._ndof = None
        self._active = None


def solve_dynamics_checked(
    operators,
    U0,
    ndof,
    dirichlet,
    *,
    dt,
    n_steps,
    rtol=1.0e-9,
    maxit=60,
):
    """Small cardiac backward-Euler driver using checked Newton steps."""
    if not isinstance(n_steps, int) or n_steps < 1:
        raise ValueError("n_steps must be a positive integer")
    if not np.isfinite(dt) or dt <= 0.0:
        raise ValueError("dt must be finite and positive")
    U = np.asarray(U0, dtype=float).copy()
    total = 0
    for step in range(1, n_steps + 1):
        t = step * dt
        bc = dirichlet(t) if callable(dirichlet) else dict(dirichlet)
        predictor = U
        for operator in operators:
            if hasattr(operator, "predictor"):
                predictor = np.asarray(operator.predictor(dt), dtype=float)
                break
        U, _, nit = checked_newton_solve(
            operators,
            predictor,
            None,
            ndof,
            bc,
            t=t,
            dt=dt,
            rtol=rtol,
            maxit=maxit,
        )
        total += nit
    return U, {
        "steps": n_steps,
        "t": n_steps * dt,
        "total_newton": total,
    }


def solve_increments_checked(
    operators,
    U0,
    ndof,
    dirichlet,
    *,
    n_steps=4,
    rtol=1.0e-9,
    maxit=60,
):
    """History-free proportional load stepping with checked Newton increments."""
    if not isinstance(n_steps, int) or n_steps < 1:
        raise ValueError("n_steps must be a positive integer")
    target = dict(dirichlet)
    U = np.asarray(U0, dtype=float).copy()
    total = 0
    for step in range(1, n_steps + 1):
        fraction = step / n_steps
        bc = {dof: fraction * value for dof, value in target.items()}
        U, _, nit = checked_newton_solve(
            operators,
            U,
            None,
            ndof,
            bc,
            t=fraction,
            dt=1.0 / n_steps,
            rtol=rtol,
            maxit=maxit,
        )
        total += nit
    return U, total


def require_distributed_success(U, info, *, tol, context):
    """Reject a distributed smoke result whose reported final solve failed."""
    if not np.isfinite(tol) or tol <= 0.0:
        raise ValueError("tol must be finite and positive")
    displacement = np.asarray(U, dtype=float)
    residual_norm = float(info.get("rnorm", float("inf")))
    if not np.all(np.isfinite(displacement)):
        raise RuntimeError(f"{context} returned a non-finite displacement")
    if bool(info.get("ksp_diverged", False)):
        raise RuntimeError(f"{context} reported a diverged PETSc linear solve")
    if not np.isfinite(residual_norm) or residual_norm >= tol:
        raise RuntimeError(
            f"{context} did not meet its final residual tolerance: "
            f"|R|={residual_norm:.6e}, required < {tol:.6e}"
        )
    return displacement
