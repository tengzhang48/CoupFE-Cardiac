"""Constant-average-acceleration Newmark inertia and velocity kinematics.

CoupFE's compact ``InertiaOperator`` uses backward-Euler kinematics. This
application also provides the Newmark-β scheme with β=0.25 and γ=0.5. The
same velocity map is exposed to the Robin dashpot so inertia and boundary
damping use one time discretization.

Backward-Euler residual:  M/dt²·(u − û)
Newmark residual:         M·a_{n+1},  a_{n+1} = (u − u_pred)/(β dt²)
  u_pred = u_n + dt·v_n + (0.5−β) dt²·a_n
  v_{n+1} = v_n + dt[(1−γ) a_n + γ a_{n+1}]

Internal forces (element, Robin, pressure) are evaluated at u_{n+1} (α_f=0), so
this composes with the existing operators unchanged.
"""
from __future__ import annotations

import numpy as np

from coupfe.operators.base import Residual, Tangent


class NewmarkInertia:
    def __init__(self, mass, ndof, *, beta=0.25, gamma=0.5,
                 u0=None, v0=None, a0=None):
        self.M = np.asarray(mass, float)
        self.ndof = int(ndof)
        self.beta = float(beta)
        self.gamma = float(gamma)
        if self.M.shape != (self.ndof,) or not np.all(np.isfinite(self.M)):
            raise ValueError("mass must be a finite vector with shape (ndof,)")
        if np.any(self.M < 0.0):
            raise ValueError("mass entries must be nonnegative")
        if not np.isfinite(self.beta) or self.beta <= 0.0:
            raise ValueError("beta must be finite and positive")
        if not np.isfinite(self.gamma) or self.gamma <= 0.0:
            raise ValueError("gamma must be finite and positive")
        self.u_prev = np.zeros(ndof) if u0 is None else np.asarray(u0, float).copy()
        self.v_prev = np.zeros(ndof) if v0 is None else np.asarray(v0, float).copy()
        self.a_prev = np.zeros(ndof) if a0 is None else np.asarray(a0, float).copy()
        for name, value in (
            ("u0", self.u_prev),
            ("v0", self.v_prev),
            ("a0", self.a_prev),
        ):
            if value.shape != (self.ndof,) or not np.all(np.isfinite(value)):
                raise ValueError(f"{name} must be a finite vector with shape (ndof,)")
        self._idx = np.nonzero(self.M != 0.0)[0]
        self._m = self.M[self._idx]

    def _u_pred(self, dt):
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt must be finite and positive")
        return (self.u_prev + dt * self.v_prev
                + (0.5 - self.beta) * dt * dt * self.a_prev)

    def predictor(self, dt):
        """Warm-start displacement predictor (constant-average-acceleration)."""
        return self.u_prev + dt * self.v_prev + 0.5 * dt * dt * self.a_prev

    def _accel(self, U, dt):
        return (U - self._u_pred(dt)) / (self.beta * dt * dt)

    def velocity(self, U, dt):
        """Velocity at the trial displacement under the Newmark update."""
        acceleration = self._accel(np.asarray(U, float), dt)
        return self.v_prev + dt * (
            (1.0 - self.gamma) * self.a_prev + self.gamma * acceleration
        )

    def velocity_tangent(self, dt):
        """Scalar ``dv_{n+1}/du_{n+1}`` for boundary dashpot tangents."""
        if not np.isfinite(dt) or dt <= 0.0:
            raise ValueError("dt must be finite and positive")
        return self.gamma / (self.beta * dt)

    def residual(self, U, state, t, dt) -> Residual:
        U = np.asarray(U, float)
        a = self._accel(U, dt)
        return Residual(self._idx, self._m * a[self._idx])

    def tangent(self, U, state, t, dt) -> Tangent:
        d = self._m / (self.beta * dt * dt)
        return Tangent(self._idx, self._idx, d)

    def commit(self, U, state, t, dt):
        U = np.asarray(U, float)
        a_new = self._accel(U, dt)
        self.v_prev = self.velocity(U, dt)
        self.a_prev = a_new
        self.u_prev = U.copy()
        return state
