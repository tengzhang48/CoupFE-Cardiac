"""Rank-local Q1/P0 condensation for the distributed cardiac companion.

The retained serial benchmark owns its existing
``LocalPressureHex8Operator``.  This module provides the element-batch form
needed by PETSc distributed assembly without changing that serial operator or
CoupFE Core.  A :class:`FusedQ1P0Batch` caches the standard cardiac material
kernel's joint residual/tangent result (with ``kappa=0``), while evaluating the
exactly condensed one-pressure-per-element residual and tangent only when the
corresponding PETSc callback asks for them.

Global element identifiers are carried explicitly.  Consequently an invalid
finite-strain trial can be reduced across ranks to one deterministic element
identifier rather than reporting a rank-local array offset.
"""
from __future__ import annotations

import numpy as np

from local_pressure import (
    InvalidDeformationError,
    LOCAL_PRESSURE_LAWS,
    LOG_MEAN_DILATATION_LAW,
    _DN_DXI,
    _pressure_response,
)


class RankLocalInvalidDeformationError(InvalidDeformationError):
    """Recoverable invalid trial tagged with its global element identifier."""

    def __init__(self, global_element_id, minimum_determinant):
        self.global_element_id = int(global_element_id)
        self.minimum_determinant = float(minimum_determinant)
        super().__init__(
            "distributed local-pressure Hex8 deformation is inverted or "
            "non-finite: "
            f"global_element={self.global_element_id}, "
            f"min(det(F))={self.minimum_determinant:.6e}"
        )


def _as_element_ids(element_ids, count):
    identifiers = np.asarray(element_ids, dtype=np.int64)
    if identifiers.shape != (count,):
        raise ValueError(
            f"global element IDs have shape {identifiers.shape}; expected {(count,)}"
        )
    if len(np.unique(identifiers)) != count:
        raise ValueError("global element IDs must be unique on each rank")
    if count and np.min(identifiers) < 0:
        raise ValueError("global element IDs must be nonnegative")
    return identifiers.copy()


class RankLocalQ1P0PressureBatch:
    """Exactly condensed volumetric Q1/P0 blocks for this rank's Hex8 cells."""

    def __init__(
        self,
        coordinates,
        *,
        bulk_modulus,
        global_element_ids,
        pressure_law=LOG_MEAN_DILATATION_LAW,
    ):
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.ndim != 3 or coordinates.shape[1:] != (8, 3):
            raise ValueError("coordinates must have shape (n_element, 8, 3)")
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("coordinates must be finite")
        self.coordinates = coordinates.copy()
        self.n_element = len(coordinates)
        self.global_element_ids = _as_element_ids(
            global_element_ids, self.n_element
        )
        self.bulk_modulus = float(bulk_modulus)
        self.pressure_law = str(pressure_law)
        if not np.isfinite(self.bulk_modulus) or self.bulk_modulus <= 0.0:
            raise ValueError("bulk_modulus must be finite and positive")
        if self.pressure_law not in LOCAL_PRESSURE_LAWS:
            raise ValueError(
                f"pressure_law must be one of {sorted(LOCAL_PRESSURE_LAWS)}"
            )

        jacobian = np.einsum("eai,gaj->egij", coordinates, _DN_DXI)
        with np.errstate(over="ignore", invalid="ignore"):
            determinant = np.linalg.det(jacobian)
        bad = np.any(~np.isfinite(determinant) | (determinant <= 0.0), axis=1)
        if np.any(bad):
            local = self._lowest_global_local_index(bad)
            raise ValueError(
                "distributed local-pressure reference Hex8 has a non-positive "
                "or non-finite Gauss-point Jacobian determinant in global "
                f"element {self.global_element_ids[local]}"
            )
        try:
            inverse = np.linalg.inv(jacobian)
        except np.linalg.LinAlgError as error:
            raise ValueError(
                "distributed local-pressure reference Hex8 Jacobian is singular"
            ) from error
        self._gradient = np.einsum("gaj,egjk->egak", _DN_DXI, inverse)
        self._weight = determinant
        self._volume = determinant.sum(axis=1)
        self._kinematics_cache = None
        if (
            not np.all(np.isfinite(inverse))
            or not np.all(np.isfinite(self._gradient))
            or not np.all(np.isfinite(self._volume))
            or np.any(self._volume <= 0.0)
        ):
            raise ValueError(
                "distributed local-pressure reference Hex8 produced non-finite "
                "gradients or volume"
            )

    def _lowest_global_local_index(self, mask):
        local_indices = np.flatnonzero(mask)
        if len(local_indices) == 0:
            raise RuntimeError("internal distributed validity check is inconsistent")
        return int(
            local_indices[
                np.argmin(self.global_element_ids[local_indices])
            ]
        )

    def clear_cache(self):
        """Discard the one-iterate kinematics shared by residual and tangent."""
        self._kinematics_cache = None

    def _kinematics(self, element_displacement):
        displacement = np.ascontiguousarray(
            np.asarray(element_displacement, dtype=float)
        )
        expected = (self.n_element, 24)
        if displacement.shape != expected:
            raise ValueError(
                f"element displacement has shape {displacement.shape}; "
                f"expected {expected}"
            )
        if not np.all(np.isfinite(displacement)):
            raise RuntimeError("distributed local-pressure displacement is non-finite")
        key = displacement.tobytes()
        if self._kinematics_cache is not None and self._kinematics_cache[0] == key:
            return self._kinematics_cache[1:]
        element_u = displacement.reshape(self.n_element, 8, 3)
        deformation = np.eye(3)[None, None, :, :] + np.einsum(
            "eai,egaj->egij", element_u, self._gradient
        )
        with np.errstate(over="ignore", invalid="ignore"):
            determinant = np.linalg.det(deformation)
        bad = np.any(~np.isfinite(determinant) | (determinant <= 0.0), axis=1)
        if np.any(bad):
            local = self._lowest_global_local_index(bad)
            finite = determinant[local][np.isfinite(determinant[local])]
            minimum = float(np.min(finite)) if len(finite) else float("nan")
            raise RankLocalInvalidDeformationError(
                self.global_element_ids[local], minimum
            )
        try:
            inverse_transpose = np.swapaxes(
                np.linalg.inv(deformation), -1, -2
            )
        except np.linalg.LinAlgError as error:
            # A singular cell should normally have been caught by det(F), but
            # keep the recoverable domain classification deterministic.
            local = int(np.argmin(self.global_element_ids))
            raise RankLocalInvalidDeformationError(
                self.global_element_ids[local], float("nan")
            ) from error
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            mean_log_j = (
                (self._weight * np.log(determinant)).sum(axis=1) / self._volume
            )
            pressure, pressure_derivative = _pressure_response(
                mean_log_j, self.bulk_modulus, self.pressure_law
            )
        if not np.all(np.isfinite(inverse_transpose)) or not np.all(
            np.isfinite(pressure)
        ) or not np.all(np.isfinite(pressure_derivative)):
            bad_response = (
                ~np.isfinite(pressure)
                | ~np.isfinite(pressure_derivative)
                | np.any(~np.isfinite(inverse_transpose), axis=(1, 2, 3))
            )
            local = self._lowest_global_local_index(bad_response)
            raise RankLocalInvalidDeformationError(
                self.global_element_ids[local], float("nan")
            )
        b = np.einsum(
            "egij,egaj->egai",
            inverse_transpose,
            self._gradient,
        )
        q = np.einsum("egai,eg->eai", b, self._weight)
        if not np.all(np.isfinite(b)) or not np.all(np.isfinite(q)):
            local = int(np.argmin(self.global_element_ids))
            raise RankLocalInvalidDeformationError(
                self.global_element_ids[local], float("nan")
            )
        self._kinematics_cache = (
            key,
            inverse_transpose,
            determinant,
            pressure,
            pressure_derivative,
            b,
            q,
        )
        return self._kinematics_cache[1:]

    def element_residual_batch(self, element_displacement):
        """Return condensed residual blocks with shape ``(ne, 24)``."""
        (
            _inverse_transpose,
            _determinant,
            pressure,
            _pressure_derivative,
            _b,
            q,
        ) = self._kinematics(element_displacement)
        residual = (pressure[:, None, None] * q).reshape(self.n_element, 24)
        if not np.all(np.isfinite(residual)):
            raise RuntimeError("distributed local-pressure residual is non-finite")
        return residual

    def element_tangent_batch(self, element_displacement):
        """Return condensed tangent blocks with shape ``(ne, 24, 24)``."""
        (
            _inverse_transpose,
            _determinant,
            pressure,
            pressure_derivative,
            b,
            q,
        ) = self._kinematics(element_displacement)
        dp_du = (
            pressure_derivative[:, None, None]
            * q
            / self._volume[:, None, None]
        )
        condensed = np.einsum("eai,ebk->eaibk", q, dp_du)
        geometric = -pressure[:, None, None, None, None] * np.einsum(
            "egak,egbi,eg->eaibk",
            b,
            b,
            self._weight,
        )
        tangent = (condensed + geometric).reshape(self.n_element, 24, 24)
        if not np.all(np.isfinite(tangent)):
            raise RuntimeError("distributed local-pressure tangent is non-finite")
        return tangent

    def element_rk_batch(self, element_displacement):
        """Return condensed ``(R_e, K_e)`` arrays with shape ``(ne,24[,24])``."""
        return (
            self.element_residual_batch(element_displacement),
            self.element_tangent_batch(element_displacement),
        )

    def element_pressure(self, element_displacement):
        """Return the algebraically eliminated pressure for local elements."""
        return self._kinematics(element_displacement)[2].copy()

    def deformation_jacobians(self, element_displacement):
        """Return ``det(F)`` at each local element/Gauss-point pair."""
        return self._kinematics(element_displacement)[1].copy()


class FusedQ1P0Batch:
    """Fuse a standard ``kappa=0`` material batch with condensed Q1/P0 blocks.

    ``evaluation_mode="joint"`` is the compatibility-preserving default and
    caches the paired residual/tangent blocks. ``evaluation_mode="split"``
    keeps residual and tangent callbacks separate.
    It uses Core's residual-only material entry when the compiled kernel has
    one, and always avoids constructing the Python Q1/P0 pressure tangent in a
    residual-only call.
    Both modes key their caches on material properties because active tension
    changes between accepted time steps.  :meth:`commit` invalidates the caches
    after promoting only the accepted material trial state.
    """

    def __init__(
        self,
        material_element,
        coordinates,
        *,
        bulk_modulus,
        global_element_ids,
        evaluation_mode="joint",
        pressure_stage_alpha_f=None,
        time_integrator="be",
        pressure_law=LOG_MEAN_DILATATION_LAW,
    ):
        if evaluation_mode not in ("split", "joint"):
            raise ValueError("evaluation_mode must be 'split' or 'joint'")
        self.evaluation_mode = evaluation_mode
        self.material_element = material_element
        material_dt = getattr(material_element, "dt", None)
        self.material_dt = None if material_dt is None else float(material_dt)
        if self.material_dt is not None and (
            not np.isfinite(self.material_dt) or self.material_dt <= 0.0
        ):
            raise ValueError("compiled material dt must be finite and positive")
        self.pressure_batch = RankLocalQ1P0PressureBatch(
            coordinates,
            bulk_modulus=bulk_modulus,
            global_element_ids=global_element_ids,
            pressure_law=pressure_law,
        )
        self.pressure_law = self.pressure_batch.pressure_law
        self.coordinates = self.pressure_batch.coordinates
        self.global_element_ids = self.pressure_batch.global_element_ids
        self.n_element = self.pressure_batch.n_element
        self.time_integrator = str(time_integrator)
        if self.time_integrator not in {
            "be",
            "generalized-alpha-source-matched",
        }:
            raise ValueError("unsupported fused-material time integrator")
        if pressure_stage_alpha_f is None:
            self.pressure_stage_alpha_f = None
        else:
            self.pressure_stage_alpha_f = float(pressure_stage_alpha_f)
            if not np.isfinite(self.pressure_stage_alpha_f) or not (
                0.0 <= self.pressure_stage_alpha_f < 1.0
            ):
                raise ValueError(
                    "pressure_stage_alpha_f must be finite in [0, 1)"
                )
        self._material_cache = None
        self._joint_cache = None

    @property
    def props(self):
        return self.material_element.props

    @property
    def material_residual_only_available(self):
        """Whether Core can skip the compiled material tangent as well."""
        return bool(
            getattr(
                self.material_element,
                "has_element_r_batch",
                getattr(self.material_element, "has_residual_only", False),
            )
        )

    def clear_cache(self):
        """Force the next evaluation to recompute from committed state."""
        self._material_cache = None
        self._joint_cache = None
        self.pressure_batch.clear_cache()

    def _cache_key(self, coordinates, displacement, increment):
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.shape != self.coordinates.shape or not np.array_equal(
            coordinates, self.coordinates
        ):
            raise ValueError(
                "fused local-pressure batch coordinates changed after construction"
            )
        displacement = np.asarray(displacement, dtype=float)
        increment = np.asarray(increment, dtype=float)
        expected = (self.n_element, 24)
        if displacement.shape != expected or increment.shape != expected:
            raise ValueError(
                "fused element displacement and increment must both have shape "
                f"{expected}"
            )
        return (
            displacement,
            increment,
            (
                displacement.tobytes(),
                increment.tobytes(),
                np.asarray(self.material_element.props, dtype=float).tobytes(),
            ),
        )

    def _material_cache_for(self, key):
        if self._material_cache is None or self._material_cache["key"] != key:
            self._material_cache = {
                "key": key,
                "residual": None,
                "tangent": None,
            }
        return self._material_cache

    def _material_rk(self, displacement, increment, key):
        cache = self._material_cache_for(key)
        if cache["residual"] is not None and cache["tangent"] is not None:
            return cache["residual"], cache["tangent"]
        material_residual, material_tangent = (
            self.material_element.element_rk_batch(
                self.coordinates, displacement, increment
            )
        )
        material_residual = np.asarray(material_residual, dtype=float)
        material_tangent = np.asarray(material_tangent, dtype=float)
        if (
            material_residual.shape != (self.n_element, 24)
            or material_tangent.shape != (self.n_element, 24, 24)
        ):
            raise RuntimeError(
                "compiled material returned malformed cardiac element blocks"
            )
        if not np.all(np.isfinite(material_residual)) or not np.all(
            np.isfinite(material_tangent)
        ):
            raise RuntimeError("compiled cardiac material blocks are non-finite")
        cache["residual"] = material_residual
        cache["tangent"] = material_tangent
        return material_residual, material_tangent

    def _material_residual(self, displacement, increment, key):
        cache = self._material_cache_for(key)
        if cache["residual"] is not None:
            return cache["residual"]
        if not self.material_residual_only_available:
            return self._material_rk(displacement, increment, key)[0]
        material_residual = np.asarray(
            self.material_element.element_r_batch(
                self.coordinates, displacement, increment
            ),
            dtype=float,
        )
        if material_residual.shape != (self.n_element, 24):
            raise RuntimeError(
                "compiled material returned malformed cardiac residual blocks"
            )
        if not np.all(np.isfinite(material_residual)):
            raise RuntimeError("compiled cardiac material residual is non-finite")
        cache["residual"] = material_residual
        return material_residual

    def _pressure_trial(self, displacement, increment):
        """Return pressure-stage displacement and d(stage)/d(endpoint).

        The generalized-alpha material kernel consumes the endpoint unknown and
        performs its own source-matched staging.  The application-owned Q1/P0
        pressure block is outside that kernel, so it must be evaluated at
        ``u[n+1-alpha_f] = u[n+1] - alpha_f*(u[n+1]-u[n])`` explicitly and its
        tangent must include the endpoint chain-rule factor.

        ``None`` preserves the historical backward-Euler behavior bit for bit.
        """
        if self.pressure_stage_alpha_f is None:
            return displacement, 1.0
        alpha_f = self.pressure_stage_alpha_f
        return displacement - alpha_f * increment, 1.0 - alpha_f

    def _joint_rk(self, coordinates, displacement, increment):
        displacement, increment, key = self._cache_key(
            coordinates, displacement, increment
        )
        if self._joint_cache is not None and self._joint_cache[0] == key:
            return self._joint_cache[1], self._joint_cache[2]
        # Validate the algebraic-pressure kinematics before evaluating the
        # stateful material at an invalid trial.
        pressure_displacement, pressure_tangent_scale = self._pressure_trial(
            displacement, increment
        )
        pressure_residual, pressure_tangent = self.pressure_batch.element_rk_batch(
            pressure_displacement
        )
        if pressure_tangent_scale != 1.0:
            pressure_tangent = pressure_tangent_scale * pressure_tangent
        material_residual, material_tangent = self._material_rk(
            displacement, increment, key
        )
        residual = material_residual + pressure_residual
        tangent = material_tangent + pressure_tangent
        if not np.all(np.isfinite(residual)) or not np.all(np.isfinite(tangent)):
            raise RuntimeError("joint cardiac element blocks are non-finite")
        self._joint_cache = (key, residual, tangent)
        return residual, tangent

    def element_residual_batch(self, coordinates, displacement, increment):
        """Return material-plus-pressure residual blocks for one rank."""
        if self.evaluation_mode == "joint":
            return self._joint_rk(coordinates, displacement, increment)[0]
        displacement, increment, key = self._cache_key(
            coordinates, displacement, increment
        )
        pressure_displacement, _pressure_tangent_scale = self._pressure_trial(
            displacement, increment
        )
        pressure_residual = self.pressure_batch.element_residual_batch(
            pressure_displacement
        )
        material_residual = self._material_residual(
            displacement, increment, key
        )
        residual = material_residual + pressure_residual
        if not np.all(np.isfinite(residual)):
            raise RuntimeError("fused cardiac element residual is non-finite")
        return residual

    def element_tangent_batch(self, coordinates, displacement, increment):
        """Return material-plus-pressure tangent blocks for one rank."""
        if self.evaluation_mode == "joint":
            return self._joint_rk(coordinates, displacement, increment)[1]
        displacement, increment, key = self._cache_key(
            coordinates, displacement, increment
        )
        pressure_displacement, pressure_tangent_scale = self._pressure_trial(
            displacement, increment
        )
        pressure_tangent = self.pressure_batch.element_tangent_batch(
            pressure_displacement
        )
        if pressure_tangent_scale != 1.0:
            pressure_tangent = pressure_tangent_scale * pressure_tangent
        _material_residual, material_tangent = self._material_rk(
            displacement, increment, key
        )
        tangent = material_tangent + pressure_tangent
        if not np.all(np.isfinite(tangent)):
            raise RuntimeError("fused cardiac element tangent is non-finite")
        return tangent

    def element_rk_batch(self, coordinates, displacement, increment):
        """Return both blocks using the selected split or joint mode."""
        if self.evaluation_mode == "joint":
            return self._joint_rk(coordinates, displacement, increment)
        return (
            self.element_residual_batch(coordinates, displacement, increment),
            self.element_tangent_batch(coordinates, displacement, increment),
        )

    def commit(self):
        """Promote the most recent accepted material trial state exactly once."""
        self.material_element.commit()
        self.clear_cache()

    def element_pressure(self, element_displacement):
        return self.pressure_batch.element_pressure(element_displacement)

    def deformation_jacobians(self, element_displacement):
        return self.pressure_batch.deformation_jacobians(element_displacement)
