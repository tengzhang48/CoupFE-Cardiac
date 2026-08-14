from __future__ import annotations

import numpy as np
import pytest

import geometry
import run as serial_driver
import run_mpi
from distributed_local_pressure import (
    FusedQ1P0Batch,
    RankLocalInvalidDeformationError,
    RankLocalQ1P0PressureBatch,
)
from local_pressure import LocalPressureHex8Operator
from material import CardiacHex8


def _element_dofs(elements):
    return (
        np.asarray(elements)[:, :, None] * 3 + np.arange(3)[None, None, :]
    ).reshape(len(elements), 24)


@pytest.mark.parametrize("pressure_law", ["log", "paper"])
def test_rank_local_q1p0_blocks_equal_serial_operator_for_global_cell_subset(
    pressure_law,
):
    mesh = geometry.build_mesh(
        n_t=1, n_mu=2, n_theta=4, flip_helix=True, apex_offset=0.2
    )
    global_ids = np.array([5, 1], dtype=np.int64)
    elements = mesh.elems[global_ids]
    gm = _element_dofs(elements)
    rng = np.random.default_rng(20260801)
    displacement = rng.normal(scale=2.0e-6, size=3 * mesh.n_node)

    local = RankLocalQ1P0PressureBatch(
        mesh.nodes[elements],
        bulk_modulus=1.0e6,
        global_element_ids=global_ids,
        pressure_law=pressure_law,
    )
    local_residual = local.element_residual_batch(displacement[gm])
    local_tangent = local.element_tangent_batch(displacement[gm])

    serial = LocalPressureHex8Operator(
        mesh.nodes,
        elements,
        3 * mesh.n_node,
        bulk_modulus=1.0e6,
        pressure_law=pressure_law,
    )
    serial_residual = serial.residual(
        displacement, None, 0.0, 1.0
    ).values.reshape(2, 24)
    serial_tangent = serial.tangent(
        displacement, None, 0.0, 1.0
    ).values.reshape(2, 24, 24)

    np.testing.assert_allclose(
        local_residual, serial_residual, rtol=3.0e-14, atol=3.0e-11
    )
    np.testing.assert_allclose(
        local_tangent, serial_tangent, rtol=3.0e-14, atol=3.0e-9
    )
    np.testing.assert_allclose(
        local.element_pressure(displacement[gm]),
        serial.element_pressure(displacement),
        rtol=3.0e-14,
        atol=3.0e-9,
    )
    np.testing.assert_allclose(
        local.deformation_jacobians(displacement[gm]),
        serial.deformation_jacobians(displacement),
        rtol=3.0e-14,
        atol=3.0e-14,
    )


def test_rank_local_invalid_trial_reports_lowest_global_element_id():
    mesh = geometry.build_mesh(
        n_t=1, n_mu=2, n_theta=4, flip_helix=True, apex_offset=0.2
    )
    global_ids = np.array([6, 2], dtype=np.int64)
    coordinates = mesh.nodes[mesh.elems[global_ids]]
    batch = RankLocalQ1P0PressureBatch(
        coordinates,
        bulk_modulus=1.0e6,
        global_element_ids=global_ids,
    )
    reflection = np.diag([-1.0, 1.0, 1.0])
    displacement = np.einsum(
        "eai,ji->eaj", coordinates, reflection - np.eye(3)
    ).reshape(2, 24)

    with pytest.raises(
        RankLocalInvalidDeformationError,
        match=r"global_element=2",
    ) as caught:
        batch.element_rk_batch(displacement)
    assert caught.value.global_element_id == 2


class _FakeMaterialBatch:
    def __init__(self, n_element):
        self.props = np.array([0.0, 2.0])
        self.n_element = n_element
        self.evaluations = 0
        self.commits = 0

    def element_rk_batch(self, coordinates, displacement, increment):
        self.evaluations += 1
        return (
            np.full((self.n_element, 24), self.props[1]),
            np.full((self.n_element, 24, 24), 3.0),
        )

    def commit(self):
        self.commits += 1


class _ResidualOnlyFakeMaterialBatch(_FakeMaterialBatch):
    has_residual_only = True

    def __init__(self, n_element):
        super().__init__(n_element)
        self.residual_evaluations = 0

    def element_r_batch(self, coordinates, displacement, increment):
        self.residual_evaluations += 1
        return np.full((self.n_element, 24), self.props[1])


def test_fused_q1p0_batch_evaluates_material_once_per_iterate_and_commit():
    mesh = geometry.build_mesh(
        n_t=1, n_mu=2, n_theta=4, flip_helix=True, apex_offset=0.2
    )
    global_ids = np.array([0, 7], dtype=np.int64)
    coordinates = mesh.nodes[mesh.elems[global_ids]]
    material = _FakeMaterialBatch(len(global_ids))
    fused = FusedQ1P0Batch(
        material,
        coordinates,
        bulk_modulus=1.0e6,
        global_element_ids=global_ids,
    )
    displacement = np.zeros((2, 24))
    increment = np.zeros_like(displacement)

    first_residual = fused.element_residual_batch(
        coordinates, displacement, increment
    )
    first_tangent = fused.element_tangent_batch(
        coordinates, displacement, increment
    )
    second_residual = fused.element_residual_batch(
        coordinates, displacement, increment
    )
    second_tangent = fused.element_tangent_batch(
        coordinates, displacement, increment
    )
    assert material.evaluations == 1
    np.testing.assert_array_equal(first_residual, second_residual)
    np.testing.assert_array_equal(first_tangent, second_tangent)

    material.props[1] = 4.0
    fused.element_residual_batch(coordinates, displacement, increment)
    assert material.evaluations == 2
    fused.commit()
    assert material.commits == 1
    fused.element_residual_batch(coordinates, displacement, increment)
    assert material.evaluations == 3


def test_fused_q1p0_residual_does_not_construct_pressure_tangent(monkeypatch):
    mesh = geometry.build_mesh(
        n_t=1, n_mu=2, n_theta=4, flip_helix=True, apex_offset=0.2
    )
    global_ids = np.array([0, 7], dtype=np.int64)
    coordinates = mesh.nodes[mesh.elems[global_ids]]
    fused = FusedQ1P0Batch(
        _FakeMaterialBatch(len(global_ids)),
        coordinates,
        bulk_modulus=1.0e6,
        global_element_ids=global_ids,
        evaluation_mode="split",
    )
    tangent_calls = 0
    original = fused.pressure_batch.element_tangent_batch

    def counted_tangent(displacement):
        nonlocal tangent_calls
        tangent_calls += 1
        return original(displacement)

    monkeypatch.setattr(
        fused.pressure_batch, "element_tangent_batch", counted_tangent
    )
    displacement = np.zeros((2, 24))
    increment = np.zeros_like(displacement)

    fused.element_residual_batch(coordinates, displacement, increment)
    fused.clear_cache()
    fused.element_residual_batch(coordinates, displacement, increment)
    assert tangent_calls == 0

    fused.element_tangent_batch(coordinates, displacement, increment)
    assert tangent_calls == 1


def test_split_mode_uses_core_residual_only_entry_then_joint_tangent():
    mesh = geometry.build_mesh(
        n_t=1, n_mu=2, n_theta=4, flip_helix=True, apex_offset=0.2
    )
    global_ids = np.array([0, 7], dtype=np.int64)
    coordinates = mesh.nodes[mesh.elems[global_ids]]
    material = _ResidualOnlyFakeMaterialBatch(len(global_ids))
    fused = FusedQ1P0Batch(
        material,
        coordinates,
        bulk_modulus=1.0e6,
        global_element_ids=global_ids,
        evaluation_mode="split",
    )
    displacement = np.zeros((2, 24))
    increment = np.zeros_like(displacement)

    residual = fused.element_residual_batch(
        coordinates, displacement, increment
    )
    assert material.residual_evaluations == 1
    assert material.evaluations == 0

    tangent = fused.element_tangent_batch(
        coordinates, displacement, increment
    )
    assert material.residual_evaluations == 1
    assert material.evaluations == 1
    np.testing.assert_array_equal(residual, np.full((2, 24), 2.0))
    assert tangent.shape == (2, 24, 24)
    assert np.all(np.isfinite(tangent))


def test_joint_mode_caches_pressure_and_material_pair(monkeypatch):
    mesh = geometry.build_mesh(
        n_t=1, n_mu=2, n_theta=4, flip_helix=True, apex_offset=0.2
    )
    global_ids = np.array([0, 7], dtype=np.int64)
    coordinates = mesh.nodes[mesh.elems[global_ids]]
    material = _ResidualOnlyFakeMaterialBatch(len(global_ids))
    fused = FusedQ1P0Batch(
        material,
        coordinates,
        bulk_modulus=1.0e6,
        global_element_ids=global_ids,
        evaluation_mode="joint",
    )
    tangent_calls = 0
    original = fused.pressure_batch.element_tangent_batch

    def counted_tangent(displacement):
        nonlocal tangent_calls
        tangent_calls += 1
        return original(displacement)

    monkeypatch.setattr(
        fused.pressure_batch, "element_tangent_batch", counted_tangent
    )
    displacement = np.zeros((2, 24))
    increment = np.zeros_like(displacement)

    fused.element_residual_batch(coordinates, displacement, increment)
    fused.element_tangent_batch(coordinates, displacement, increment)
    assert tangent_calls == 1
    assert material.residual_evaluations == 0
    assert material.evaluations == 1


def test_mpi_cli_exposes_split_and_joint_element_evaluation():
    parser = run_mpi._parser()
    assert parser.parse_args([]).element_evaluation == "joint"
    assert parser.parse_args(
        ["--element-evaluation", "split"]
    ).element_evaluation == "split"
    with pytest.raises(SystemExit):
        parser.parse_args(["--element-evaluation", "auto"])


class _StateBuffer:
    def __init__(self, n_element, state_size):
        self.svars = np.zeros((n_element, state_size), dtype=float)
        self.svars_trial = self.svars.copy()


def test_rank_local_fiber_initialization_uses_explicit_global_element_ids():
    mesh = geometry.build_mesh(
        n_t=1, n_mu=2, n_theta=4, flip_helix=True, apex_offset=0.2
    )
    schema = CardiacHex8()._mat.state_schema
    per_gp = sum(entry["size"] for entry in schema.values())
    full = _StateBuffer(mesh.n_elem, per_gp * 8)
    serial_driver._init_fiber_svars(full, mesh, schema, per_gp)

    global_ids = np.array([7, 1, 4], dtype=np.int64)
    local = _StateBuffer(len(global_ids), per_gp * 8)
    run_mpi._init_fiber_svars_subset(
        local, mesh, schema, per_gp, global_ids
    )
    np.testing.assert_array_equal(local.svars, full.svars[global_ids])
    np.testing.assert_array_equal(local.svars_trial, local.svars)
