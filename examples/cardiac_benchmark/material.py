"""Holzapfel-Ogden passive myocardium + active stress + viscous damping.

Benchmark 1 of Aróstica et al. (CMAME 2025). The same finite-strain material
can be generated as a Hex8 F-bar element or as a standard Hex8 whose bulk term
is replaced by the driver's condensed local-pressure operator. All quantities
are SI (Pa, m, s).

2nd Piola-Kirchhoff stress (paper Eqs. 2-4):

    S = S_passive(C; f0,s0) + Ta * (f0 ⊗ f0) + eta * (E - E_prev)/dt

with the passive Holzapfel-Ogden energy (Eq. 3), nearly-incompressible via the
isochoric first invariant Ī1 = J^(-2/3) tr C and the κ/4 (J²-1-2lnJ) penalty:

    Ψ = a/(2b) exp(b(Ī1-3))
      + Σ_{i=f,s} a_i/(2b_i) χ(I4i) (exp(b_i (I4i-1)²) - 1)
      + a_fs/(2b_fs) (exp(b_fs I8fs²) - 1)
      + κ/4 (J² - 1 - 2 ln J)

I4f = f0·C f0, I4s = s0·C s0, I8fs = f0·C s0.  χ is a smooth fiber-compression
switch ≈ H(I4-1) (sigmoid, steepness k).  We supply S analytically; the
framework complex-steps it for the consistent tangent.

State per Gauss point (36 scalars) -- stored as 3x3 structural tensors so the
codegen traces them as matrices (array-valued vector state is unsupported):
    ff    (9) : f0 ⊗ f0          -- prescribed fiber structural tensor
    ss    (9) : s0 ⊗ s0          -- prescribed sheet structural tensor
    fssym (9) : sym(f0 ⊗ s0)     -- prescribed fiber-sheet coupling tensor
    E_prev(9) : Green-Lagrange strain at previous step (viscous history)

Invariants via I4f = tr(C·ff), I4s = tr(C·ss), I8fs = tr(C·fssym).

Active tension Ta is a *property* updated each time step by the driver.

The separate generalized-alpha material has 54 scalars per Gauss point: the
three structural tensors plus accepted ``grad(u)``, ``grad(v)``, and
``grad(a)``. It reconstructs the Simula ``alpha_f`` stage and the
velocity-consistent Green--Lagrange rate inside the generated kernel; the BE
material and its ``E_prev`` update are not repurposed.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

import numpy as np

import coupfe.codegen as au
from coupfe.codegen.core import tensor
from coupfe.codegen.core.tensor import det, inv, trace, eye, dyad, exp
from coupfe.codegen.generators.uel_gen import generate_element, generate_uel
from coupfe.runtime.compiled_element import CompiledElement, build_element_kernel

try:  # package import
    from .benchmark_parameters import benchmark_configuration
except ImportError:  # direct example/script import
    from benchmark_parameters import benchmark_configuration


# --- material parameters (Tables 1-3, SI) ---------------------------------
# Compatibility copy for the retained Step 0 API. Immutable per-step physical
# definitions live in benchmark_parameters.py and are selected by the drivers.
HO_PARAMS = dict(benchmark_configuration(0, "A").material_parameters)

# Stable semantic identity for result archives.  Commit 6839c13 corrected the
# stress to differentiate the complete smooth-switch energy, so records made
# before and after that correction must never be compared as if they used one
# constitutive law.
MATERIAL_MODEL_ID = (
    "holzapfel-ogden-smooth-switch-complete-energy-derivative-v1"
)


class HolzapfelOgdenActive(au.Material):
    """Holzapfel-Ogden + active stress + viscous damping (finite strain)."""

    props = dict(HO_PARAMS)
    state_vars = dict(
        ff=np.diag([1.0, 0.0, 0.0]),    # f0 ⊗ f0 (prescribed)
        ss=np.diag([0.0, 1.0, 0.0]),    # s0 ⊗ s0 (prescribed)
        fssym=np.zeros((3, 3)),         # sym(f0 ⊗ s0) (prescribed)
        E_prev=np.zeros((3, 3)),        # previous-step Green-Lagrange strain
    )

    def _switch(self, x):
        """Smooth Heaviside H(x-1) (numerically stable sigmoid)."""
        z = self.k_sw * (x - 1.0)
        # branch on real part to avoid exp overflow (complex-step safe)
        return tensor.where(z.real > 0.0,
                            1.0 / (1.0 + exp(-z)),
                            exp(z) / (1.0 + exp(z)))

    def stress_PK1(self, F, ff_old, ss_old, fssym_old, E_prev_old, dt):
        I = eye(3)
        C = tensor.transpose(F) @ F
        J = det(F)
        Cinv = inv(C)
        I1 = trace(C)
        Jm23 = J ** (-2.0 / 3.0)
        I1bar = Jm23 * I1

        # anisotropic invariants from structural tensors
        I4f = trace(C @ ff_old)
        I4s = trace(C @ ss_old)
        I8fs = trace(C @ fssym_old)

        # --- isochoric isotropic term ---
        S_iso = (self.a * exp(self.b * (I1bar - 3.0)) * Jm23
                 * (I - (I1 / 3.0) * Cinv))

        # --- anisotropic fiber/sheet terms (with smooth compression switch) ---
        chi_f = self._switch(I4f)
        chi_s = self._switch(I4s)
        dchi_f = self.k_sw * chi_f * (1.0 - chi_f)
        dchi_s = self.k_sw * chi_s * (1.0 - chi_s)
        exp_f = exp(self.b_f * (I4f - 1.0) ** 2)
        exp_s = exp(self.b_s * (I4s - 1.0) ** 2)

        # Differentiate the complete smooth-switch energy.  The first term in
        # each coefficient comes from d(chi)/dI4; omitting it differentiates a
        # different constitutive law near I4=1.
        S_f = (
            (self.a_f / self.b_f) * dchi_f * (exp_f - 1.0)
            + 2.0 * self.a_f * chi_f * (I4f - 1.0) * exp_f
        ) * ff_old
        S_s = (
            (self.a_s / self.b_s) * dchi_s * (exp_s - 1.0)
            + 2.0 * self.a_s * chi_s * (I4s - 1.0) * exp_s
        ) * ss_old

        # --- fiber-sheet coupling term ((f⊗s + s⊗f) = 2 sym(f⊗s)) ---
        S_fs = (2.0 * self.a_fs * I8fs
                * exp(self.b_fs * I8fs ** 2)) * fssym_old

        # --- volumetric penalty ---
        S_vol = 0.5 * self.kappa * (J * J - 1.0) * Cinv

        # --- active stress along fibers ---
        S_act = self.Ta * ff_old

        # --- viscous pseudo-potential: S_visco = eta * (E - E_prev)/dt ---
        E = 0.5 * (C - I)
        S_visco = (self.eta / dt) * (E - E_prev_old)

        S = S_iso + S_f + S_s + S_fs + S_vol + S_act + S_visco
        P = F @ S
        return P, {"ff": ff_old, "ss": ss_old, "fssym": fssym_old, "E_prev": E}


class CardiacHex8(au.WeakForm):
    """Finite-strain cardiac material weak form on a Hex8 displacement field."""

    material = HolzapfelOgdenActive
    ndim = 3

    def define_fields(self):
        self.u = au.VectorField("u", degree=1)

    def momentum_equation(self, v, F):
        # dummy body for signature detection (framework calls stress_PK1)
        return self.material.stress_PK1(F)


class HolzapfelOgdenActiveGeneralizedAlpha(HolzapfelOgdenActive):
    """Source-matched generalized-alpha cardiac constitutive update.

    The Simula reference evaluates the constitutive response at
    ``u[n+1-alpha_f]`` and computes the viscous Green--Lagrange rate from the
    velocity at the same stage.  Keeping only ``E_prev`` cannot reproduce that
    rate.  This material therefore stores the accepted displacement, velocity,
    and acceleration *gradients* at each Gauss point and reconstructs the exact
    stage kinematics from the endpoint displacement supplied to the kernel.

    The returned PK1 stress is the stage stress itself.  Code generation
    differentiates it with respect to ``u[n+1]``, so all ``1-alpha_f`` and
    Newmark velocity-chain-rule factors are retained in the tangent.
    """

    state_vars = dict(
        ff=np.diag([1.0, 0.0, 0.0]),
        ss=np.diag([0.0, 1.0, 0.0]),
        fssym=np.zeros((3, 3)),
        grad_u_prev=np.zeros((3, 3)),
        grad_v_prev=np.zeros((3, 3)),
        grad_a_prev=np.zeros((3, 3)),
    )

    def stress_PK1(
        self,
        F,
        ff_old,
        ss_old,
        fssym_old,
        grad_u_prev_old,
        grad_v_prev_old,
        grad_a_prev_old,
        dt,
    ):
        # Literals are deliberate: this generated material is the closed
        # Simula-source-matched contract, not a tunable generalized-alpha API.
        alpha_f = 0.4
        gamma = 0.7
        beta = 0.36

        I = eye(3)
        grad_u_new = F - I
        grad_a_new = (
            grad_u_new
            - (
                grad_u_prev_old
                + dt * grad_v_prev_old
                + (0.5 - beta) * dt * dt * grad_a_prev_old
            )
        ) / (beta * dt * dt)
        grad_v_new = grad_v_prev_old + dt * (
            (1.0 - gamma) * grad_a_prev_old + gamma * grad_a_new
        )
        grad_u_stage = (
            alpha_f * grad_u_prev_old + (1.0 - alpha_f) * grad_u_new
        )
        grad_v_stage = (
            alpha_f * grad_v_prev_old + (1.0 - alpha_f) * grad_v_new
        )
        F_stage = I + grad_u_stage

        C = tensor.transpose(F_stage) @ F_stage
        J = det(F_stage)
        Cinv = inv(C)
        I1 = trace(C)
        Jm23 = J ** (-2.0 / 3.0)
        I1bar = Jm23 * I1

        I4f = trace(C @ ff_old)
        I4s = trace(C @ ss_old)
        I8fs = trace(C @ fssym_old)

        S_iso = (
            self.a
            * exp(self.b * (I1bar - 3.0))
            * Jm23
            * (I - (I1 / 3.0) * Cinv)
        )

        chi_f = self._switch(I4f)
        chi_s = self._switch(I4s)
        dchi_f = self.k_sw * chi_f * (1.0 - chi_f)
        dchi_s = self.k_sw * chi_s * (1.0 - chi_s)
        exp_f = exp(self.b_f * (I4f - 1.0) ** 2)
        exp_s = exp(self.b_s * (I4s - 1.0) ** 2)
        S_f = (
            (self.a_f / self.b_f) * dchi_f * (exp_f - 1.0)
            + 2.0 * self.a_f * chi_f * (I4f - 1.0) * exp_f
        ) * ff_old
        S_s = (
            (self.a_s / self.b_s) * dchi_s * (exp_s - 1.0)
            + 2.0 * self.a_s * chi_s * (I4s - 1.0) * exp_s
        ) * ss_old
        S_fs = (
            2.0 * self.a_fs * I8fs * exp(self.b_fs * I8fs**2)
        ) * fssym_old
        S_vol = 0.5 * self.kappa * (J * J - 1.0) * Cinv
        S_act = self.Ta * ff_old

        # Holzapfel Eq. 2.163 as implemented by the Simula reference:
        # E_dot = F^T sym(grad(v) F^-1) F = sym(F^T grad(v)).
        E_dot_stage = 0.5 * (
            tensor.transpose(F_stage) @ grad_v_stage
            + tensor.transpose(grad_v_stage) @ F_stage
        )
        S_visco = self.eta * E_dot_stage

        P = F_stage @ (
            S_iso + S_f + S_s + S_fs + S_vol + S_act + S_visco
        )
        return P, {
            "ff": ff_old,
            "ss": ss_old,
            "fssym": fssym_old,
            "grad_u_prev": grad_u_new,
            "grad_v_prev": grad_v_new,
            "grad_a_prev": grad_a_new,
        }


class CardiacHex8GeneralizedAlpha(CardiacHex8):
    """Hex8 weak form with the source-matched staged constitutive update."""

    material = HolzapfelOgdenActiveGeneralizedAlpha


_PER_GP = 36   # ff(9) + ss(9) + fssym(9) + E_prev(9)
_PER_GP_GENERALIZED_ALPHA = 54  # three tensors + grad(u), grad(v), grad(a)


def struct_tensors(f0, s0):
    """Return structural tensors after a local Gram--Schmidt projection.

    Interpolating nodal directions component-wise does not preserve unit length
    or orthogonality at a Gauss point.  The constitutive invariants assume a
    local orthonormal fiber/sheet pair, so enforce that application-owned
    convention before constructing the tensors.
    """
    f0 = np.asarray(f0, dtype=float)
    s0 = np.asarray(s0, dtype=float)
    if f0.shape != (3,) or s0.shape != (3,):
        raise ValueError("fiber and sheet directions must have shape (3,)")
    f_norm = float(np.linalg.norm(f0))
    if not np.isfinite(f_norm) or f_norm <= 1.0e-12:
        raise ValueError("fiber direction is degenerate")
    f0 = f0 / f_norm
    s0 = s0 - float(np.dot(s0, f0)) * f0
    s_norm = float(np.linalg.norm(s0))
    if not np.isfinite(s_norm) or s_norm <= 1.0e-12:
        raise ValueError("sheet direction is degenerate after fiber projection")
    s0 = s0 / s_norm
    ff = np.outer(f0, f0)
    ss = np.outer(s0, s0)
    fssym = 0.5 * (np.outer(f0, s0) + np.outer(s0, f0))
    return ff, ss, fssym


def verification_state():
    """A representative material-point state for CS-vs-FD tangent check."""
    F = np.array([[1.06, 0.03, -0.02],
                  [0.01, 0.97, 0.04],
                  [-0.02, 0.03, 1.02]])
    f0 = np.array([0.6, 0.8, 0.0]); f0 = f0 / np.linalg.norm(f0)
    s0 = np.array([-0.8, 0.6, 0.0]); s0 = s0 / np.linalg.norm(s0)
    ff, ss, fssym = struct_tensors(f0, s0)
    E_prev = np.array([[0.01, 0.002, 0.0],
                       [0.002, -0.005, 0.001],
                       [0.0, 0.001, 0.003]])
    return dict(F=F, ff=ff, ss=ss, fssym=fssym, E_prev=E_prev, dt=1.0e-3)


def verify_material_tangent(eps: float = 1e-6, tol: float = 1e-6,
                            verbose: bool = True) -> bool:
    """Self-contained complex-step vs finite-difference dP/dF check."""
    m = HolzapfelOgdenActive()
    st = verification_state()
    F0, ff, ss, fssym, Ev, dt = (st["F"], st["ff"], st["ss"], st["fssym"],
                                 st["E_prev"], st["dt"])

    def P_of(Fc):
        P, _ = m.stress_PK1(Fc, ff, ss, fssym, Ev, dt)
        return np.asarray(P)

    # stress-free reference: F=I, no active, no viscous history -> P≈0
    Z = np.zeros((3, 3))
    P_ref = m.stress_PK1(np.eye(3) + 0j, ff + 0j, ss + 0j, fssym + 0j, Z + 0j, dt)[0]
    free = float(np.max(np.abs(np.asarray(P_ref).real)))

    h = 1e-20
    cs = np.zeros((3, 3, 3, 3))
    fd = np.zeros((3, 3, 3, 3))
    for a in range(3):
        for b in range(3):
            Fc = F0.astype(complex).copy(); Fc[a, b] += 1j * h
            cs[:, :, a, b] = P_of(Fc).imag / h
            Fp = F0.copy(); Fp[a, b] += eps
            Fm = F0.copy(); Fm[a, b] -= eps
            fd[:, :, a, b] = (P_of(Fp).real - P_of(Fm).real) / (2 * eps)
    norm = float(np.max(np.abs(cs)))
    rel = float(np.max(np.abs(cs - fd))) / (norm if norm > 0 else 1.0)
    ok = (free < 1e-6) and (rel < tol)
    if verbose:
        print(f"  stress-free |P(F=I)| = {free:.3e}")
        print(f"  dP/dF CS-vs-FD rel err = {rel:.3e}  "
              f"[{'PASS' if ok else 'FAIL'}]")
    return ok


def verify_generalized_alpha_material_tangent(
    eps: float = 1.0e-6, tol: float = 2.0e-6, verbose: bool = True
) -> bool:
    """Check the staged PK1 derivative with old kinematic state held fixed."""
    model = HolzapfelOgdenActiveGeneralizedAlpha()
    st = verification_state()
    F0 = st["F"]
    ff, ss, fssym = st["ff"], st["ss"], st["fssym"]
    dt = st["dt"]
    grad_u_old = np.array(
        [[0.012, -0.004, 0.001], [0.003, -0.008, 0.002], [0.0, 0.001, 0.006]]
    )
    grad_v_old = np.array(
        [[0.08, -0.03, 0.01], [0.02, -0.05, 0.015], [0.0, 0.01, 0.03]]
    )
    grad_a_old = np.array(
        [[0.4, -0.2, 0.1], [0.1, -0.3, 0.05], [0.0, 0.04, 0.2]]
    )

    def P_of(Fc):
        P, _ = model.stress_PK1(
            Fc,
            ff,
            ss,
            fssym,
            grad_u_old,
            grad_v_old,
            grad_a_old,
            dt,
        )
        return np.asarray(P)

    h = 1.0e-20
    cs = np.zeros((3, 3, 3, 3))
    fd = np.zeros((3, 3, 3, 3))
    for a in range(3):
        for b in range(3):
            Fc = F0.astype(complex).copy()
            Fc[a, b] += 1j * h
            cs[:, :, a, b] = P_of(Fc).imag / h
            Fp = F0.copy()
            Fm = F0.copy()
            Fp[a, b] += eps
            Fm[a, b] -= eps
            fd[:, :, a, b] = (P_of(Fp).real - P_of(Fm).real) / (2.0 * eps)
    scale = max(1.0, float(np.max(np.abs(cs))))
    relative = float(np.max(np.abs(cs - fd))) / scale
    ok = relative < tol
    if verbose:
        print(
            "  generalized-alpha staged dP/dF CS-vs-FD rel err = "
            f"{relative:.3e}  [{'PASS' if ok else 'FAIL'}]"
        )
    return ok


def build_kernel(tmpdir: Optional[str] = None, backend: str = "native",
                 formulation: str = "fbar_mechanics",
                 module_name: Optional[str] = None) -> Tuple[Path, object]:
    if formulation not in {"fbar_mechanics", "standard"}:
        raise ValueError(
            "cardiac kernel formulation must be 'fbar_mechanics' or 'standard'"
        )
    problem = CardiacHex8()
    if not verify_material_tangent(verbose=False):
        raise RuntimeError("HolzapfelOgdenActive tangent verification failed")
    here = Path(__file__).resolve().parent
    out = Path(tmpdir) if tmpdir else here
    out.mkdir(parents=True, exist_ok=True)
    tag = "fbar" if formulation == "fbar_mechanics" else "std"
    for_path = out / f"cardiac_hex8_{tag}.for"
    if module_name is None:
        module_name = f"cardiac_hex8_{tag}_module"
    if backend == "native":
        generate_element(problem, str(for_path), element="Hex8",
                         formulation=formulation, backend="native")
    else:
        generate_uel(problem, str(for_path), element="Hex8",
                     formulation=formulation)
    mod = build_element_kernel(str(for_path), module_name, workdir=str(out))
    return for_path, mod


def build_generalized_alpha_kernel(
    tmpdir: Optional[str] = None,
    backend: str = "native",
    formulation: str = "standard",
    module_name: Optional[str] = None,
) -> Tuple[Path, object]:
    """Build the endpoint-unknown/source-staged generalized-alpha kernel."""
    if formulation != "standard":
        raise ValueError(
            "generalized-alpha cardiac kernel requires the standard Hex8 "
            "kernel; local pressure is fused by the application batch"
        )
    problem = CardiacHex8GeneralizedAlpha()
    if not verify_generalized_alpha_material_tangent(verbose=False):
        raise RuntimeError(
            "HolzapfelOgdenActiveGeneralizedAlpha tangent verification failed"
        )
    here = Path(__file__).resolve().parent
    out = Path(tmpdir) if tmpdir else here
    out.mkdir(parents=True, exist_ok=True)
    for_path = out / "cardiac_hex8_generalized_alpha_std.for"
    if module_name is None:
        module_name = "cardiac_hex8_generalized_alpha_std_module"
    if backend == "native":
        generate_element(
            problem,
            str(for_path),
            element="Hex8",
            formulation="standard",
            backend="native",
        )
    else:
        generate_uel(
            problem, str(for_path), element="Hex8", formulation="standard"
        )
    module = build_element_kernel(
        str(for_path), module_name, workdir=str(out)
    )
    return for_path, module


if __name__ == "__main__":
    print("Verifying Holzapfel-Ogden active material tangent...")
    ok = verify_material_tangent(verbose=True)
    if not ok:
        raise RuntimeError("tangent verification failed")
    print("Generating + compiling Hex8 F-bar kernel...")
    path, _ = build_kernel()
    print(f"OK -> {path}")
