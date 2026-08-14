"""Fail-closed checks for the optional accepted-state viscosity instrument."""

from __future__ import annotations

import numpy as np
import pytest

from viscous_evidence import ViscousEvidenceGroup


_CUBE = np.array(
    [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ]
)


class _Element:
    has_element_r_batch = True

    def __init__(self, *, raise_during_probe=False):
        self.props = np.array([100.0, 0.0])
        self.svars = np.zeros((1, 8 * 9))
        self.svars_trial = np.full_like(self.svars, 7.0)
        self.dt = 0.1
        self.raise_during_probe = raise_during_probe

    def element_r_batch(self, coordinates, element_u, element_du):
        self.svars_trial[:] = 99.0
        if self.raise_during_probe:
            raise RuntimeError("synthetic residual-only failure")
        return np.zeros((1, 24))


class _Group:
    def __init__(self, element):
        self.element = element
        self.gm = np.arange(24).reshape(1, 24)
        self._rk_cache = object()
        self.commits = 0

    def _coords(self):
        return _CUBE[None, :, :]

    def _gather(self, U, state):
        values = np.asarray(U, dtype=float).reshape(1, 24)
        return values, np.zeros_like(values)

    def residual(self, U, state, t, dt):
        return "residual"

    def tangent(self, U, state, t, dt):
        return "tangent"

    def commit(self, U, state, t, dt):
        self.commits += 1
        return state


def _instrument(tmp_path, *, element=None, start=0.2, end=0.32):
    element = _Element() if element is None else element
    group = _Group(element)
    instrument = ViscousEvidenceGroup(
        group,
        ndof=24,
        eta_index=0,
        ta_index=1,
        e_prev_offset=0,
        per_gp=9,
        output_path=tmp_path / "viscous.npz",
        probe_start=start,
        probe_end=end,
        case="B",
        formulation="std-kappa",
    )
    return instrument, group, element


def test_instrument_delegates_outside_window_without_writing_evidence(tmp_path):
    instrument, group, _ = _instrument(tmp_path)
    displacement = np.linspace(0.0, 0.01, 24)
    state = {"accepted": True}

    assert instrument.residual(displacement, state, 0.1, 0.1) == "residual"
    assert instrument.tangent(displacement, state, 0.1, 0.1) == "tangent"
    assert instrument.commit(displacement, state, 0.1, 0.1) is state
    assert group.commits == 1
    np.testing.assert_array_equal(instrument.previous_u, displacement)
    assert instrument.records == []
    assert not instrument.output_path.exists()
    with pytest.raises(RuntimeError, match="no accepted probe records"):
        instrument.finalize()


def test_probe_exception_restores_material_inputs_and_does_not_commit(tmp_path):
    element = _Element(raise_during_probe=True)
    instrument, group, element = _instrument(tmp_path, element=element)
    props_before = element.props.copy()
    state_before = element.svars.copy()
    trial_before = element.svars_trial.copy()

    with pytest.raises(RuntimeError, match="synthetic residual-only failure"):
        instrument.commit(np.zeros(24), None, 0.25, 0.1)

    np.testing.assert_array_equal(element.props, props_before)
    np.testing.assert_array_equal(element.svars, state_before)
    np.testing.assert_array_equal(element.svars_trial, trial_before)
    assert group._rk_cache is None
    assert group.commits == 0
    assert instrument.records == []
    assert not instrument.output_path.exists()


def test_constructor_rejects_wrong_scope_missing_abi_and_output_collision(tmp_path):
    element = _Element()
    group = _Group(element)
    arguments = dict(
        ndof=24,
        eta_index=0,
        ta_index=1,
        e_prev_offset=0,
        per_gp=9,
        output_path=tmp_path / "viscous.npz",
        probe_start=0.2,
        probe_end=0.32,
        case="B",
        formulation="std-kappa",
    )
    with pytest.raises(ValueError, match="passive Case B std-kappa"):
        ViscousEvidenceGroup(group, **{**arguments, "case": "A"})
    with pytest.raises(ValueError, match="passive Case B std-kappa"):
        ViscousEvidenceGroup(group, **{**arguments, "formulation": "fbar"})

    element.has_element_r_batch = False
    with pytest.raises(RuntimeError, match="residual-only ABI"):
        ViscousEvidenceGroup(group, **arguments)
    element.has_element_r_batch = True

    collision = tmp_path / "existing.npz"
    collision.write_bytes(b"do not overwrite")
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        ViscousEvidenceGroup(group, **{**arguments, "output_path": collision})
    assert collision.read_bytes() == b"do not overwrite"

    with pytest.raises(ValueError, match="must end in .npz"):
        ViscousEvidenceGroup(
            group, **{**arguments, "output_path": tmp_path / "viscous.txt"}
        )
    with pytest.raises(ValueError, match="bounds"):
        ViscousEvidenceGroup(
            group, **{**arguments, "probe_start": 0.4, "probe_end": 0.3}
        )
