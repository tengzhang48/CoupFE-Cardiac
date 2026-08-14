"""Active tension τ(t) and endocardial pressure p(t) time courses.

Bestel-Clément-Sorine activation + pressure ODEs (paper Eqs. 5-8), adapted
from `activation_model.py` and `pressure_model.py` in Henrik Finsberg,
Joakim Sundnes, and Jonas van den Brink, "finsberg/cardiac_benchmark:
v1.0.0," Zenodo, https://doi.org/10.5281/zenodo.10875818. The exact archive
is `finsberg/cardiac_benchmark-v1.0.0.zip` (MD5
`be92da5dbc1fd26d424bf88ef7db13b4`; peeled Git commit
`325d17d850c2e2032abb85a4191a5795d3008ab7`).

Changes made for CoupFE-Cardiac: the two upstream modules were combined;
parameters were routed through the benchmark-configuration registry; scalar
math was vectorized with NumPy; input, integration-span, solver-completion,
shape, and finiteness checks were added; and the public functions were renamed
to `tau_of_t` and `p_of_t`. The Radau integration method is retained.
τ(t)/p(t) are precomputed in Python and fed to the solver as per-step scalars.

The cited Zenodo deposit is licensed CC BY 4.0. License text and the complete
source/attribution record are in `LICENSES/CC-BY-4.0.txt` and
`THIRD_PARTY_NOTICES.md`.

SPDX-License-Identifier: CC-BY-4.0
"""
from __future__ import annotations

import numpy as np
import scipy.integrate

try:  # package import
    from .benchmark_parameters import benchmark_configuration
except ImportError:  # direct example/script import
    from benchmark_parameters import benchmark_configuration

# --- Benchmark-1 activation parameters (Table 3, SI) ----------------------
_STEP0_CONFIGURATION = benchmark_configuration(0, "A")
ACT_PARAMS = dict(_STEP0_CONFIGURATION.activation_parameters)

# --- Benchmark-1 endocardial-pressure parameters (Table 4, SI) ------------
PRES_PARAMS = dict(_STEP0_CONFIGURATION.pressure_parameters)


def _time_samples(times):
    values = np.asarray(times, dtype=float)
    if (
        values.ndim != 1
        or len(values) == 0
        or not np.all(np.isfinite(values))
        or (len(values) > 1 and not np.all(np.diff(values) > 0.0))
    ):
        raise ValueError("times must be a nonempty, finite, strictly increasing vector")
    return values


def _integrate_history(rhs, times, label, *, t_span=None):
    times = _time_samples(times)
    if t_span is None:
        integration_span = (float(times[0]), float(times[-1]))
    else:
        span = np.asarray(t_span, dtype=float)
        if (
            span.shape != (2,)
            or not np.all(np.isfinite(span))
            or not span[0] < span[1]
            or times[0] < span[0]
            or times[-1] > span[1]
        ):
            raise ValueError(
                "t_span must be a finite increasing pair containing all samples"
            )
        integration_span = (float(span[0]), float(span[1]))
    if len(times) == 1 and integration_span[0] == times[0]:
        return np.zeros(1)
    result = scipy.integrate.solve_ivp(
        rhs, integration_span, [0.0], t_eval=times, method="Radau"
    )
    if (
        not result.success
        or result.t.shape != times.shape
        or not np.allclose(result.t, times, rtol=0.0, atol=1.0e-14)
        or result.y.shape != (1, len(times))
        or not np.all(np.isfinite(result.y))
    ):
        message = getattr(result, "message", "incomplete or non-finite output")
        raise RuntimeError(f"{label} integration failed: {message}")
    return result.y[0].copy()


def tau_of_t(times, p=ACT_PARAMS, *, t_span=None):
    r"""Active tension: τ̇ = -|a|τ + σ₀⟨a⟩₊, a = a_max f + a_min(1-f),
    f = S⁺(t-t_sys) S⁻(t-t_dias)."""
    g = p["gamma"]
    f = lambda t: 0.25 * (1 + np.tanh((t - p["t_sys"]) / g)) * (1 - np.tanh((t - p["t_dias"]) / g))
    a = lambda t: p["a_max"] * f(t) + p["a_min"] * (1 - f(t))
    rhs = lambda t, tau: -abs(a(t)) * tau + p["sigma_0"] * max(a(t), 0.0)
    return _integrate_history(rhs, times, "activation", t_span=t_span)


def p_of_t(times, p=PRES_PARAMS, *, t_span=None):
    r"""Endocardial pressure: ṗ = -|b|p + σ_mid⟨b⟩₊ + σ_pre⟨g_pre⟩₊,
    b = a(t) + α_pre g_pre + α_mid,  g_pre = S⁻(t-t_dias_pre)."""
    g = p["gamma"]
    f = lambda t: 0.25 * (1 + np.tanh((t - p["t_sys_pre"]) / g)) * (1 - np.tanh((t - p["t_dias_pre"]) / g))
    a = lambda t: p["a_max"] * f(t) + p["a_min"] * (1 - f(t))
    g_pre = lambda t: 0.5 * (1 - np.tanh((t - p["t_dias_pre"]) / g))
    b = lambda t: a(t) + p["alpha_pre"] * g_pre(t) + p["alpha_mid"]
    rhs = lambda t, P: (-abs(b(t)) * P + p["sigma_mid"] * max(b(t), 0.0)
                        + p["sigma_pre"] * max(g_pre(t), 0.0))
    return _integrate_history(rhs, times, "pressure", t_span=t_span)


if __name__ == "__main__":
    t = np.arange(0.0, 1.0 + 1e-9, 1e-3)
    tau = tau_of_t(t)
    print(f"τ(t): peak = {tau.max():.2f} Pa at t={t[tau.argmax()]:.3f}s "
          f"(paper ≈ 118817.07 Pa)")
    pr = p_of_t(t)
    print(f"p(t): peak = {pr.max():.2f} Pa at t={t[pr.argmax()]:.3f}s "
          f"(paper ≈ 16117.52 Pa)")
