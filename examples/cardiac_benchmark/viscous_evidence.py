"""Experiment-only accepted-state split of cardiac material forces.

The decorator probes the production residual-only kernel immediately before
the ordinary material-state commit.  It changes only ``eta`` while holding the
accepted displacement and committed state fixed, then restores every mutable
input and delegates to the normal ``ElementGroup.commit``.  This is campaign
instrumentation, not a replacement constitutive law.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np


_NAT = np.array(
    [
        [-1, -1, -1],
        [1, -1, -1],
        [1, 1, -1],
        [-1, 1, -1],
        [-1, -1, 1],
        [1, -1, 1],
        [1, 1, 1],
        [-1, 1, 1],
    ],
    dtype=float,
)


def _shape_derivatives():
    """Return dN/dxi in the generated kernel's eight-point order."""
    g1 = 1.0 / np.sqrt(3.0)
    points = np.empty((8, 3), dtype=float)
    for gp in range(8):
        points[gp, 0] = g1 if gp % 2 == 0 else -g1
        points[gp, 1] = g1 if (gp // 2) % 2 == 0 else -g1
        points[gp, 2] = g1 if gp // 4 == 0 else -g1
    derivatives = np.empty((8, 8, 3), dtype=float)
    for gp, point in enumerate(points):
        for node, natural in enumerate(_NAT):
            for direction in range(3):
                others = [axis for axis in range(3) if axis != direction]
                derivatives[gp, node, direction] = (
                    0.125
                    * natural[direction]
                    * (1.0 + natural[others[0]] * point[others[0]])
                    * (1.0 + natural[others[1]] * point[others[1]])
                )
    return derivatives


def _norm(vector):
    return float(np.linalg.norm(np.asarray(vector, dtype=float)))


class ViscousEvidenceGroup:
    """Decorate one material ElementGroup with fail-closed eta probes."""

    requires_converged_commit = True

    def __init__(
        self,
        group,
        *,
        ndof,
        eta_index,
        ta_index,
        e_prev_offset,
        per_gp,
        output_path,
        probe_start,
        probe_end,
        case,
        formulation,
        checkpoint_interval=10,
    ):
        if case != "B" or formulation != "std-kappa":
            raise ValueError(
                "viscous evidence is scoped only to passive Case B std-kappa"
            )
        if not bool(getattr(group.element, "has_element_r_batch", False)):
            raise RuntimeError("viscous evidence requires the residual-only ABI")
        if not np.isfinite(probe_start) or not np.isfinite(probe_end):
            raise ValueError("viscous probe bounds must be finite")
        if probe_start < 0.0 or probe_end < probe_start:
            raise ValueError("viscous probe bounds are invalid")
        self.group = group
        self.element = group.element
        self.ndof = int(ndof)
        self.eta_index = int(eta_index)
        self.ta_index = int(ta_index)
        self.e_prev_offset = int(e_prev_offset)
        self.per_gp = int(per_gp)
        self.output_path = Path(output_path).expanduser().resolve()
        if self.output_path.suffix != ".npz":
            raise ValueError("viscous evidence output must end in .npz")
        if self.output_path.exists():
            raise FileExistsError(
                f"refusing to overwrite viscous evidence {self.output_path}"
            )
        self.output_path.parent.mkdir(parents=True, exist_ok=True)
        self.probe_start = float(probe_start)
        self.probe_end = float(probe_end)
        self.checkpoint_interval = int(checkpoint_interval)
        if self.checkpoint_interval < 1:
            raise ValueError("checkpoint interval must be positive")
        self.previous_u = np.zeros(self.ndof, dtype=float)
        self.records = []
        self.other_operators = ()
        self.solver = None
        self._finalized = False
        self._prepare_reference_kinematics()

    def _prepare_reference_kinematics(self):
        derivatives = _shape_derivatives()
        coordinates = np.asarray(self.group._coords(), dtype=float)
        if coordinates.shape[1:] != (8, 3):
            raise ValueError("viscous evidence requires three-dimensional Hex8")
        jacobian = np.einsum("gai,eaj->egij", derivatives, coordinates)
        determinants = np.linalg.det(jacobian)
        if not np.all(np.isfinite(determinants)) or np.any(determinants <= 0.0):
            raise RuntimeError("reference Hex8 Jacobians must be finite and positive")
        inverse = np.linalg.inv(jacobian)
        # Generated convention: Jac(i,j)=dX_j/dxi_i and
        # dN/dX_j=sum_i Jinv(j,i)*dN/dxi_i.
        gradients = np.einsum("gai,egji->egaj", derivatives, inverse)
        self._gradients = gradients
        self._determinants = determinants
        self._reference_volume = float(np.sum(determinants))

    def configure_closure(self, *, other_operators, solver):
        self.other_operators = tuple(other_operators)
        self.solver = solver

    def residual(self, U, state, t, dt):
        return self.group.residual(U, state, t, dt)

    def tangent(self, U, state, t, dt):
        return self.group.tangent(U, state, t, dt)

    def __getattr__(self, name):
        # Preserve optional Operator capabilities without changing the wrapped
        # ElementGroup's normal residual/tangent/commit path.
        return getattr(self.group, name)

    def _scatter(self, element_residual):
        vector = np.zeros(self.ndof, dtype=float)
        np.add.at(
            vector,
            self.group.gm.ravel(),
            np.asarray(element_residual, dtype=float).ravel(),
        )
        return vector

    def _independent_viscous_residual(self, U, committed_state, dt, eta):
        element_u = np.asarray(U, dtype=float)[self.group.gm].reshape(-1, 8, 3)
        deformation_gradient = np.eye(3)[None, None, :, :] + np.einsum(
            "eai,egaj->egij", element_u, self._gradients
        )
        strain = 0.5 * (
            np.einsum(
                "egki,egkj->egij", deformation_gradient, deformation_gradient
            )
            - np.eye(3)[None, None, :, :]
        )
        previous_strain = np.asarray(committed_state, dtype=float).reshape(
            len(element_u), 8, self.per_gp
        )[:, :, self.e_prev_offset : self.e_prev_offset + 9].reshape(-1, 8, 3, 3)
        strain_rate = (strain - previous_strain) / float(dt)
        viscous_second_piola = float(eta) * strain_rate
        viscous_first_piola = np.einsum(
            "egik,egkj->egij", deformation_gradient, viscous_second_piola
        )
        element_residual = np.einsum(
            "egij,egaj,eg->eai",
            viscous_first_piola,
            self._gradients,
            self._determinants,
        ).reshape(-1, 24)
        strain_rate_squared = np.einsum(
            "egij,egij->eg", strain_rate, strain_rate
        )
        weighted_integral = float(
            np.sum(self._determinants * strain_rate_squared)
        )
        return {
            "element_residual": element_residual,
            "strain": strain,
            "strain_rate": strain_rate,
            "strain_rate_frobenius_rms_per_s": float(
                np.sqrt(weighted_integral / self._reference_volume)
            ),
            "strain_rate_frobenius_max_per_s": float(
                np.sqrt(np.max(strain_rate_squared))
            ),
            "dissipation_rate_w": float(eta) * weighted_integral,
            "viscous_potential_w": 0.5 * float(eta) * weighted_integral,
        }

    def _global_closure(self, U, state, t, dt, material_residual):
        if self.solver is None or self.solver.last_diagnostics is None:
            raise RuntimeError("viscous evidence requires accepted PETSc diagnostics")
        total = np.asarray(material_residual, dtype=float).copy()
        component_norms = {}
        for operator in self.other_operators:
            contribution = operator.residual(U, state, t, dt)
            vector = np.zeros(self.ndof, dtype=float)
            np.add.at(vector, contribution.gdofs, contribution.values)
            name = type(operator).__name__
            component_norms[name] = component_norms.get(name, 0.0) + _norm(vector)
            total += vector
        observed = _norm(total)
        expected = float(self.solver.last_diagnostics.final_residual_norm)
        absolute_error = abs(observed - expected)
        if not np.isclose(observed, expected, rtol=1.0e-8, atol=1.0e-12):
            raise RuntimeError(
                "accepted-state global residual does not close: "
                f"recomputed={observed:.6e}, solver={expected:.6e}"
            )
        return observed, expected, absolute_error, component_norms

    def _probe(self, U, state, t, dt):
        if not np.isclose(float(self.element.dt), float(dt), rtol=0.0, atol=1.0e-15):
            raise RuntimeError("compiled element dt differs from accepted step dt")
        props0 = np.asarray(self.element.props, dtype=float).copy()
        if not np.array_equal(props0, self.element.props):
            raise RuntimeError("material properties are not a stable float array")
        eta = float(props0[self.eta_index])
        active_tension = float(props0[self.ta_index])
        if not np.isclose(eta, 100.0, rtol=0.0, atol=1.0e-12):
            raise RuntimeError(f"unexpected viscosity eta={eta}")
        if active_tension != 0.0:
            raise RuntimeError("passive Case B evidence requires Ta=0")
        committed0 = np.asarray(self.element.svars, dtype=float).copy()
        trial0 = np.asarray(self.element.svars_trial, dtype=float).copy()
        gathered_u, gathered_du = self.group._gather(U, state)
        coordinates = self.group._coords()

        def evaluate(test_eta):
            self.element.props[:] = props0
            self.element.props[self.eta_index] = float(test_eta)
            element_residual = self.element.element_r_batch(
                coordinates, gathered_u, gathered_du
            ).copy()
            trial = np.asarray(self.element.svars_trial, dtype=float).copy()
            if not np.array_equal(self.element.svars, committed0):
                raise RuntimeError("residual-only probe mutated committed material state")
            if not np.all(np.isfinite(element_residual)) or not np.all(
                np.isfinite(trial)
            ):
                raise RuntimeError("viscous probe produced non-finite data")
            return element_residual, trial

        try:
            r_eta_element, state_eta = evaluate(eta)
            r_half_element, state_half = evaluate(0.5 * eta)
            r_zero_element, state_zero = evaluate(0.0)
            r_repeat_element, state_repeat = evaluate(eta)
            if not (
                np.array_equal(state_eta, state_half)
                and np.array_equal(state_eta, state_zero)
                and np.array_equal(state_eta, state_repeat)
            ):
                raise RuntimeError("material state update depends on eta")
            if not np.array_equal(self.element.svars, committed0):
                raise RuntimeError("eta probes changed committed state")
        finally:
            self.element.props[:] = props0
            self.element.svars_trial = trial0.copy()
            self.group._rk_cache = None
        if not np.array_equal(self.element.props, props0):
            raise RuntimeError("viscous probe did not restore material properties")

        r_eta = self._scatter(r_eta_element)
        r_half = self._scatter(r_half_element)
        r_zero = self._scatter(r_zero_element)
        r_repeat = self._scatter(r_repeat_element)
        r_viscous = r_eta - r_zero
        half_error = (r_half - r_zero) - 0.5 * r_viscous
        eta_linearity_relative = _norm(half_error) / max(_norm(r_viscous), 1.0e-30)
        element_half_error = (r_half_element - r_zero_element) - 0.5 * (
            r_eta_element - r_zero_element
        )
        element_linearity_relative = _norm(element_half_error) / max(
            _norm(r_eta_element - r_zero_element), 1.0e-30
        )
        repeat_relative = _norm(r_repeat - r_eta) / max(_norm(r_eta), 1.0e-30)
        if max(eta_linearity_relative, element_linearity_relative) >= 1.0e-11:
            raise RuntimeError("material residual is not linear in eta")
        if repeat_relative >= 1.0e-12:
            raise RuntimeError("restored eta residual is not reproducible")

        independent = self._independent_viscous_residual(U, committed0, dt, eta)
        state_strain = state_eta.reshape(-1, 8, self.per_gp)[
            :, :, self.e_prev_offset : self.e_prev_offset + 9
        ].reshape(-1, 8, 3, 3)
        state_strain_error = float(
            np.max(np.abs(state_strain - independent["strain"]))
        )
        if state_strain_error >= 1.0e-11:
            raise RuntimeError("accepted E_prev state does not equal reconstructed E_n")
        symmetry_error = float(
            np.max(
                np.abs(
                    independent["strain_rate"]
                    - np.swapaxes(independent["strain_rate"], -1, -2)
                )
            )
        )
        if symmetry_error >= 1.0e-12:
            raise RuntimeError("reconstructed Green-strain rate is not symmetric")
        independent_relative = _norm(
            independent["element_residual"] - (r_eta_element - r_zero_element)
        ) / max(_norm(r_eta_element - r_zero_element), 1.0e-30)
        independent_max_abs = float(
            np.max(
                np.abs(
                    independent["element_residual"]
                    - (r_eta_element - r_zero_element)
                )
            )
        )
        if independent_relative >= 5.0e-11:
            raise RuntimeError("independent quadrature does not close viscous residual")
        if independent["dissipation_rate_w"] < -1.0e-14:
            raise RuntimeError("viscous dissipation is negative")

        closure = self._global_closure(U, state, t, dt, r_eta)
        increment = np.asarray(U, dtype=float) - self.previous_u
        velocity = increment / float(dt)
        norm_eta = _norm(r_eta)
        norm_zero = _norm(r_zero)
        norm_viscous = _norm(r_viscous)
        denominator = norm_zero * norm_viscous
        cosine = (
            float(np.dot(r_zero, r_viscous) / denominator)
            if denominator > 0.0
            else 0.0
        )
        record = {
            "time_s": float(t),
            "accepted_displacement_m": np.asarray(U, dtype=float).copy(),
            "accepted_increment_m": increment.copy(),
            "material_total_eta100_n": r_eta,
            "material_nonviscous_eta0_n": r_zero,
            "material_viscous_difference_n": r_viscous,
            "unassembled_element_viscous_l2_n": np.linalg.norm(
                r_eta_element - r_zero_element, axis=1
            ),
            "material_total_l2_n": norm_eta,
            "material_nonviscous_l2_n": norm_zero,
            "material_viscous_l2_n": norm_viscous,
            "viscous_to_nonviscous_l2_ratio": norm_viscous / max(norm_zero, 1.0e-30),
            "viscous_bounded_l2_fraction": norm_viscous
            / max(norm_zero + norm_viscous, 1.0e-30),
            "viscous_to_total_l2_ratio": norm_viscous / max(norm_eta, 1.0e-30),
            "nonviscous_viscous_cosine": cosine,
            "material_total_max_nodal_l2_n": float(
                np.max(np.linalg.norm(r_eta.reshape(-1, 3), axis=1))
            ),
            "material_nonviscous_max_nodal_l2_n": float(
                np.max(np.linalg.norm(r_zero.reshape(-1, 3), axis=1))
            ),
            "material_viscous_max_nodal_l2_n": float(
                np.max(np.linalg.norm(r_viscous.reshape(-1, 3), axis=1))
            ),
            "material_total_incremental_power_w": float(np.dot(r_eta, velocity)),
            "material_nonviscous_incremental_power_w": float(
                np.dot(r_zero, velocity)
            ),
            "material_viscous_incremental_power_w": float(
                np.dot(r_viscous, velocity)
            ),
            "eta_linearity_relative_error": eta_linearity_relative,
            "eta_linearity_element_relative_error": element_linearity_relative,
            "eta_repeat_relative_error": repeat_relative,
            "accepted_state_strain_max_abs_error": state_strain_error,
            "strain_rate_symmetry_max_abs_error_per_s": symmetry_error,
            "independent_quadrature_relative_error": independent_relative,
            "independent_quadrature_max_abs_n": independent_max_abs,
            "full_residual_recomputed_l2_n": closure[0],
            "full_residual_solver_l2_n": closure[1],
            "full_residual_norm_absolute_error_n": closure[2],
            "other_operator_residual_l2_json": json.dumps(
                closure[3], sort_keys=True, separators=(",", ":")
            ),
            "strain_rate_frobenius_rms_per_s": independent[
                "strain_rate_frobenius_rms_per_s"
            ],
            "strain_rate_frobenius_max_per_s": independent[
                "strain_rate_frobenius_max_per_s"
            ],
            "viscous_dissipation_rate_w": independent["dissipation_rate_w"],
            "viscous_potential_w": independent["viscous_potential_w"],
        }
        scalar_values = [
            value
            for value in record.values()
            if isinstance(value, (int, float, np.integer, np.floating))
        ]
        if not np.all(np.isfinite(np.asarray(scalar_values, dtype=float))):
            raise RuntimeError("viscous evidence contains non-finite scalar data")
        return record, state_eta

    def commit(self, U, state, t, dt):
        probe = self.probe_start - 1.0e-12 <= float(t) <= self.probe_end + 1.0e-12
        record = None
        expected_state = None
        if probe:
            record, expected_state = self._probe(U, state, t, dt)
        committed = self.group.commit(U, state, t, dt)
        if expected_state is not None and not np.array_equal(
            self.element.svars, expected_state
        ):
            raise RuntimeError("normal material commit differs from probed accepted state")
        self.previous_u = np.asarray(U, dtype=float).copy()
        if record is not None:
            self.records.append(record)
            if len(self.records) % self.checkpoint_interval == 0:
                checkpoint = self.output_path.with_name(
                    f"{self.output_path.stem}.checkpoint-{len(self.records):04d}.npz"
                )
                self._write(checkpoint, final=False)
        return committed

    def _payload(self, *, final):
        if not self.records:
            raise RuntimeError("viscous evidence has no accepted probe records")
        keys = tuple(self.records[0])
        payload = {}
        for key in keys:
            values = [record[key] for record in self.records]
            if isinstance(values[0], np.ndarray):
                payload[key] = np.stack(values)
            else:
                payload[key] = np.asarray(values)
        payload.update(
            {
                "evidence_schema": np.asarray("coupfe-cardiac-viscous-split-v1"),
                "complete": np.asarray(bool(final)),
                "record_count": np.asarray(len(self.records), dtype=int),
                "probe_start_s": np.asarray(self.probe_start),
                "probe_end_s": np.asarray(self.probe_end),
                "eta_pa_s": np.asarray(float(self.element.props[self.eta_index])),
                "reference_volume_m3": np.asarray(self._reference_volume),
                "force_decomposition": np.asarray(
                    "assembled material residual at fixed accepted U/state: "
                    "r_visc=r_eta100-r_eta0; r_nonvisc=r_eta0"
                ),
                "power_definition": np.asarray(
                    "assembled material residual dot accepted secant velocity"
                ),
                "independent_quadrature_definition": np.asarray(
                    "Pvis=F*eta*(E_n-E_prev)/dt integrated with 2x2x2 Hex8 quadrature"
                ),
            }
        )
        return payload

    def _write(self, path, *, final):
        path = Path(path)
        if path.exists():
            raise FileExistsError(f"refusing to overwrite viscous evidence {path}")
        np.savez_compressed(path, **self._payload(final=final))

    def finalize(self):
        if self._finalized:
            return self.output_path
        self._write(self.output_path, final=True)
        self._finalized = True
        return self.output_path
