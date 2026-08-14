"""Source-matched generalized-alpha kinematics for the cardiac benchmark.

The local Simula/FEniCS benchmark uses ``alpha_m=0.2`` and ``alpha_f=0.4``.
Its definitions are

``a[n+1] = (u[n+1] - u[n] - dt*v[n] - (0.5-beta)*dt**2*a[n]) /
            (beta*dt**2)``

``v[n+1] = v[n] + dt*((1-gamma)*a[n] + gamma*a[n+1])``

with ``gamma=0.5+alpha_f-alpha_m=0.7`` and
``beta=(gamma+0.5)**2/4=0.36``.  Inertia is evaluated at the
``1-alpha_m`` stage; displacement, velocity, constitutive terms, Robin terms,
and prescribed loads are evaluated at the ``1-alpha_f`` stage.

This application-owned module intentionally does not add generalized-alpha to
CoupFE Core.  The benchmark parameters are a closed contract rather than user
tuning knobs: a mismatched value is rejected instead of being labeled as the
source-matched path.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np


@dataclass(frozen=True)
class GeneralizedAlphaParameters:
    """The fixed Simula/FEniCS generalized-alpha parameter set."""

    alpha_m: float = 0.2
    alpha_f: float = 0.4
    gamma: float = 0.7
    beta: float = 0.36

    def __post_init__(self):
        values = np.asarray(
            [self.alpha_m, self.alpha_f, self.gamma, self.beta], dtype=float
        )
        if not np.all(np.isfinite(values)):
            raise ValueError("generalized-alpha parameters must be finite")
        expected_gamma = 0.5 + self.alpha_f - self.alpha_m
        expected_beta = 0.25 * (expected_gamma + 0.5) ** 2
        if not np.isclose(self.gamma, expected_gamma, rtol=0.0, atol=1.0e-15):
            raise ValueError(
                "generalized-alpha gamma must equal 0.5 + alpha_f - alpha_m"
            )
        if not np.isclose(self.beta, expected_beta, rtol=0.0, atol=1.0e-15):
            raise ValueError(
                "generalized-alpha beta must equal (gamma + 0.5)^2 / 4"
            )
        if not (0.0 <= self.alpha_m < 1.0 and 0.0 <= self.alpha_f < 1.0):
            raise ValueError("generalized-alpha alpha values must lie in [0, 1)")
        if self.gamma <= 0.0 or self.beta <= 0.0:
            raise ValueError("generalized-alpha gamma and beta must be positive")
        source_values = np.array([0.2, 0.4, 0.7, 0.36])
        if not np.allclose(values, source_values, rtol=0.0, atol=1.0e-15):
            raise ValueError(
                "source-matched generalized-alpha requires exactly "
                "alpha_m=0.2, alpha_f=0.4, gamma=0.7, beta=0.36"
            )

    def metadata(self) -> dict:
        return {
            **asdict(self),
            "parameter_source": (
                "finsberg/cardiac_benchmark problem.py defaults"
            ),
            "acceleration_stage": "alpha_m*a_n + (1-alpha_m)*a_np1",
            "force_stage": "alpha_f*x_n + (1-alpha_f)*x_np1",
            "load_time": "t_np1 - alpha_f*dt",
        }

    @staticmethod
    def _dt(dt) -> float:
        value = float(dt)
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError("dt must be finite and positive")
        return value

    @staticmethod
    def _state(name, value, shape):
        result = np.asarray(value, dtype=float)
        if result.shape != shape or not np.all(np.isfinite(result)):
            raise ValueError(f"{name} must be finite with shape {shape}")
        return result

    def displacement_predictor(self, u_n, v_n, a_n, dt):
        """Return a second-order Newton warm start, not a committed state."""
        dt = self._dt(dt)
        u_n = np.asarray(u_n, dtype=float)
        shape = u_n.shape
        u_n = self._state("u_n", u_n, shape)
        v_n = self._state("v_n", v_n, shape)
        a_n = self._state("a_n", a_n, shape)
        return u_n + dt * v_n + 0.5 * dt * dt * a_n

    def acceleration(self, u_np1, u_n, v_n, a_n, dt):
        dt = self._dt(dt)
        u_n = np.asarray(u_n, dtype=float)
        shape = u_n.shape
        u_np1 = self._state("u_np1", u_np1, shape)
        u_n = self._state("u_n", u_n, shape)
        v_n = self._state("v_n", v_n, shape)
        a_n = self._state("a_n", a_n, shape)
        base = u_n + dt * v_n + (0.5 - self.beta) * dt * dt * a_n
        return (u_np1 - base) / (self.beta * dt * dt)

    def velocity(self, a_np1, v_n, a_n, dt):
        dt = self._dt(dt)
        v_n = np.asarray(v_n, dtype=float)
        shape = v_n.shape
        a_np1 = self._state("a_np1", a_np1, shape)
        v_n = self._state("v_n", v_n, shape)
        a_n = self._state("a_n", a_n, shape)
        return v_n + dt * (
            (1.0 - self.gamma) * a_n + self.gamma * a_np1
        )

    def force_stage(self, old, new):
        old = np.asarray(old, dtype=float)
        shape = old.shape
        old = self._state("old force-stage state", old, shape)
        new = self._state("new force-stage state", new, shape)
        return self.alpha_f * old + (1.0 - self.alpha_f) * new

    def acceleration_stage(self, a_n, a_np1):
        a_n = np.asarray(a_n, dtype=float)
        shape = a_n.shape
        a_n = self._state("a_n", a_n, shape)
        a_np1 = self._state("a_np1", a_np1, shape)
        return self.alpha_m * a_n + (1.0 - self.alpha_m) * a_np1

    def load_time(self, endpoint_time, dt) -> float:
        dt = self._dt(dt)
        endpoint_time = float(endpoint_time)
        if not np.isfinite(endpoint_time):
            raise ValueError("endpoint time must be finite")
        return endpoint_time - self.alpha_f * dt

    def acceleration_tangent(self, dt) -> float:
        dt = self._dt(dt)
        return 1.0 / (self.beta * dt * dt)

    def velocity_tangent(self, dt) -> float:
        dt = self._dt(dt)
        return self.gamma / (self.beta * dt)

    def force_displacement_tangent(self) -> float:
        return 1.0 - self.alpha_f

    def force_velocity_tangent(self, dt) -> float:
        return (1.0 - self.alpha_f) * self.velocity_tangent(dt)

    def inertia_tangent(self, dt) -> float:
        return (1.0 - self.alpha_m) * self.acceleration_tangent(dt)


SOURCE_MATCHED_GENERALIZED_ALPHA = GeneralizedAlphaParameters()
