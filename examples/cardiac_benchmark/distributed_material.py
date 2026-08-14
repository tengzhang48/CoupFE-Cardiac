"""Rank-local standard-material batch for the Cardiac MPI companion.

Unlike the historical :mod:`distributed_local_pressure` adapter, this batch
does not replace the material's volumetric term.  It evaluates the standard
Hex8 kernel with the paper's pointwise ``kappa`` contribution intact, matching
the serial ``--formulation std-kappa`` path.
"""
from __future__ import annotations

import numpy as np

from distributed_local_pressure import RankLocalQ1P0PressureBatch, _as_element_ids


class RankLocalStandardMaterialBatch:
    """Adapt a native standard Hex8 kernel to distributed element callbacks."""

    def __init__(
        self,
        material_element,
        coordinates,
        *,
        global_element_ids,
        evaluation_mode="joint",
        time_integrator="be",
    ):
        if evaluation_mode not in ("joint", "split"):
            raise ValueError("evaluation_mode must be 'joint' or 'split'")
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.ndim != 3 or coordinates.shape[1:] != (8, 3):
            raise ValueError("coordinates must have shape (n_element, 8, 3)")
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("coordinates must be finite")
        self.material_element = material_element
        material_dt = getattr(material_element, "dt", None)
        self.material_dt = None if material_dt is None else float(material_dt)
        if self.material_dt is not None and (
            not np.isfinite(self.material_dt) or self.material_dt <= 0.0
        ):
            raise ValueError("compiled material dt must be finite and positive")
        self.coordinates = coordinates.copy()
        self.n_element = len(coordinates)
        self.global_element_ids = _as_element_ids(
            global_element_ids, self.n_element
        )
        self.time_integrator = str(time_integrator)
        if self.time_integrator not in {
            "be",
            "generalized-alpha-source-matched",
        }:
            raise ValueError("unsupported standard-material time integrator")
        self.evaluation_mode = evaluation_mode
        # Reuse the independently checked finite-strain kinematics solely for
        # retained det(F) output.  Its pressure never enters this batch's R/K.
        self._kinematic_probe = RankLocalQ1P0PressureBatch(
            self.coordinates,
            bulk_modulus=1.0,
            global_element_ids=self.global_element_ids,
        )
        self._cache = None

    @property
    def props(self):
        return self.material_element.props

    @property
    def material_residual_only_available(self):
        return bool(
            getattr(
                self.material_element,
                "has_element_r_batch",
                getattr(self.material_element, "has_residual_only", False),
            )
        )

    def clear_cache(self):
        self._cache = None
        self._kinematic_probe.clear_cache()

    def _inputs(self, coordinates, displacement, increment):
        coordinates = np.asarray(coordinates, dtype=float)
        if coordinates.shape != self.coordinates.shape or not np.array_equal(
            coordinates, self.coordinates
        ):
            raise ValueError(
                "rank-local standard material coordinates changed after construction"
            )
        displacement = np.asarray(displacement, dtype=float)
        increment = np.asarray(increment, dtype=float)
        expected = (self.n_element, 24)
        if displacement.shape != expected or increment.shape != expected:
            raise ValueError(
                "rank-local standard material displacement and increment must "
                f"both have shape {expected}"
            )
        if not np.all(np.isfinite(displacement)) or not np.all(
            np.isfinite(increment)
        ):
            raise ValueError("rank-local standard material inputs must be finite")
        key = (
            displacement.tobytes(),
            increment.tobytes(),
            np.asarray(self.material_element.props, dtype=float).tobytes(),
        )
        if self._cache is None or self._cache["key"] != key:
            self._cache = {"key": key, "residual": None, "tangent": None}
        return displacement, increment, self._cache

    def _joint(self, displacement, increment, cache):
        if cache["residual"] is None or cache["tangent"] is None:
            residual, tangent = self.material_element.element_rk_batch(
                self.coordinates, displacement, increment
            )
            residual = np.asarray(residual, dtype=float)
            tangent = np.asarray(tangent, dtype=float)
            if residual.shape != (self.n_element, 24) or tangent.shape != (
                self.n_element,
                24,
                24,
            ):
                raise RuntimeError(
                    "compiled standard material returned malformed element blocks"
                )
            if not np.all(np.isfinite(residual)) or not np.all(
                np.isfinite(tangent)
            ):
                raise RuntimeError(
                    "compiled standard material returned non-finite element blocks"
                )
            cache["residual"] = residual
            cache["tangent"] = tangent
        return cache["residual"], cache["tangent"]

    def element_residual_batch(self, coordinates, displacement, increment):
        displacement, increment, cache = self._inputs(
            coordinates, displacement, increment
        )
        if self.evaluation_mode == "joint" or not (
            self.material_residual_only_available
        ):
            return self._joint(displacement, increment, cache)[0]
        if cache["residual"] is None:
            residual = np.asarray(
                self.material_element.element_r_batch(
                    self.coordinates, displacement, increment
                ),
                dtype=float,
            )
            if residual.shape != (self.n_element, 24) or not np.all(
                np.isfinite(residual)
            ):
                raise RuntimeError(
                    "compiled standard material returned an invalid residual batch"
                )
            cache["residual"] = residual
        return cache["residual"]

    def element_tangent_batch(self, coordinates, displacement, increment):
        displacement, increment, cache = self._inputs(
            coordinates, displacement, increment
        )
        return self._joint(displacement, increment, cache)[1]

    def element_rk_batch(self, coordinates, displacement, increment):
        displacement, increment, cache = self._inputs(
            coordinates, displacement, increment
        )
        if self.evaluation_mode == "joint":
            return self._joint(displacement, increment, cache)
        return (
            self.element_residual_batch(coordinates, displacement, increment),
            self.element_tangent_batch(coordinates, displacement, increment),
        )

    def commit(self):
        self.material_element.commit()
        self.clear_cache()

    def deformation_jacobians(self, element_displacement):
        return self._kinematic_probe.deformation_jacobians(element_displacement)
