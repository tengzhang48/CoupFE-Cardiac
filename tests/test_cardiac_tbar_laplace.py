"""Focused checks for the closed-mesh Laplace transmural-field utility."""

from __future__ import annotations

import hashlib
import json

import numpy as np
import pytest

import geometry
import run as cardiac_run
import tbar_laplace


@pytest.fixture(scope="module")
def small_closed_mesh():
    return geometry.build_closed_mesh(
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )


def test_q1_laplace_field_satisfies_equation_boundaries_and_bounds(
    small_closed_mesh,
):
    mesh = small_closed_mesh
    endocardial = np.unique(mesh.facets_endo)
    epicardial = np.unique(mesh.facets_epi)

    field, diagnostics = tbar_laplace.solve_q1_hex_laplace(
        mesh.nodes,
        mesh.elems,
        endocardial,
        epicardial,
    )

    assert field.shape == (mesh.n_node,)
    assert np.all(np.isfinite(field))
    np.testing.assert_array_equal(field[endocardial], 0.0)
    np.testing.assert_array_equal(field[epicardial], 1.0)
    assert field.min() >= -1.0e-9
    assert field.max() <= 1.0 + 1.0e-9
    assert diagnostics["minimum_gauss_det_j_m3"] > 0.0
    assert diagnostics["linear_residual_l2"] <= diagnostics[
        "linear_residual_limit"
    ]
    assert diagnostics["natural_boundary"] == "homogeneous Neumann on base"


def test_cli_writes_atomic_self_checking_outputs_without_absolute_paths(
    tmp_path, capsys
):
    output = tmp_path / "closed_tbar.npy"
    assert tbar_laplace.main(
        [
            "--nt", "2",
            "--ncore", "4",
            "--nradial", "1",
            "--core-half-width", "0.36",
            "--out", str(output),
        ]
    ) == 0

    metadata_path = output.with_suffix(".meta.json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    field = np.load(output, allow_pickle=False)
    digest = hashlib.sha256(output.read_bytes()).hexdigest()

    assert field.shape == (123,)
    assert metadata["schema"] == "coupfe-cardiac-laplace-tbar-v1"
    assert metadata["mesh_topology"] == "closed_multiblock_disk"
    assert metadata["output_npy"] == output.name
    assert metadata["sha256"] == digest
    generated_mesh = geometry.build_closed_mesh(
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )
    assert metadata["mesh_identity"] == tbar_laplace.closed_mesh_identity(
        generated_mesh.nodes,
        generated_mesh.elems,
    )
    assert not output.is_absolute() or not metadata["output_npy"].startswith("/")
    assert str(tmp_path) not in json.dumps(metadata)
    assert not list(tmp_path.glob(".*.tmp"))
    assert json.loads(capsys.readouterr().out)["sha256"] == digest

    provenance = cardiac_run._load_laplace_tbar_metadata(
        generated_mesh,
        output,
    )
    assert provenance == {
        "field_sha256": digest,
        "metadata_filename": metadata_path.name,
        "metadata_sha256": hashlib.sha256(metadata_path.read_bytes()).hexdigest(),
        "metadata_schema": "coupfe-cardiac-laplace-tbar-v1",
    }


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        ("field-name", "does not name this field"),
        ("field-sha", "SHA-256 disagrees"),
        ("mesh-count", "mesh parameter 'elements' disagrees"),
        ("boundary", "different boundary conditions"),
        ("residual", "failed solver check"),
    ],
)
def test_serial_tbar_input_rejects_inconsistent_sidecar(
    tmp_path, mutation, message
):
    output = tmp_path / "closed_tbar.npy"
    tbar_laplace.generate_closed_mesh_tbar(
        output,
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )
    metadata_path = output.with_suffix(".meta.json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if mutation == "field-name":
        metadata["output_npy"] = "different.npy"
    elif mutation == "field-sha":
        metadata["sha256"] = "0" * 64
    elif mutation == "mesh-count":
        metadata["mesh_parameters"]["elements"] += 1
    elif mutation == "boundary":
        metadata["boundary_conditions"]["base"] = "Dirichlet tbar=0"
    else:
        metadata["solver"]["linear_residual_l2"] = (
            2.0 * metadata["solver"]["linear_residual_limit"]
        )
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    mesh = geometry.build_closed_mesh(
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )
    with pytest.raises(ValueError, match=message):
        cardiac_run._load_laplace_tbar_metadata(mesh, output)


def test_serial_tbar_input_requires_generator_sidecar(tmp_path, small_closed_mesh):
    output = tmp_path / "orphan.npy"
    np.save(output, np.asarray(small_closed_mesh.param[:, 0], dtype=float))

    with pytest.raises(ValueError, match="requires sibling metadata"):
        cardiac_run._load_laplace_tbar_metadata(small_closed_mesh, output)


def test_mesh_identity_is_independent_of_input_byte_order(small_closed_mesh):
    expected = tbar_laplace.closed_mesh_identity(
        small_closed_mesh.nodes,
        small_closed_mesh.elems,
    )
    observed = tbar_laplace.closed_mesh_identity(
        np.asarray(small_closed_mesh.nodes, dtype=">f8"),
        np.asarray(small_closed_mesh.elems, dtype=">i8"),
    )
    assert observed == expected


def test_serial_tbar_input_rejects_legacy_geometry_without_identity(
    tmp_path, small_closed_mesh
):
    output = tmp_path / "closed_tbar.npy"
    tbar_laplace.generate_closed_mesh_tbar(
        output,
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )
    metadata_path = output.with_suffix(".meta.json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata.pop("mesh_identity")
    metadata_path.write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="not bound to mesh geometry"):
        cardiac_run._load_laplace_tbar_metadata(small_closed_mesh, output)


def test_serial_tbar_input_rejects_pre_correction_same_shape_coordinates(
    tmp_path,
):
    output = tmp_path / "closed_tbar.npy"
    tbar_laplace.generate_closed_mesh_tbar(
        output,
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )
    stale_mesh = geometry.build_closed_mesh(
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )
    interior = np.flatnonzero(stale_mesh.param[:, 0] == 0.5)
    for node in interior:
        t, mu, theta = stale_mesh.param[node]
        stale_mesh.nodes[node] = geometry.point(
            t,
            mu,
            0.0 if not np.isfinite(theta) else theta,
        )

    with pytest.raises(ValueError, match="node_coordinates_sha256.*disagrees"):
        cardiac_run._load_laplace_tbar_metadata(stale_mesh, output)


def test_serial_tbar_input_rejects_same_nodes_with_different_connectivity(
    tmp_path,
):
    output = tmp_path / "closed_tbar.npy"
    tbar_laplace.generate_closed_mesh_tbar(
        output,
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )
    changed = geometry.build_closed_mesh(
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )
    changed.elems[[0, 1]] = changed.elems[[1, 0]]

    with pytest.raises(
        ValueError,
        match="element_connectivity_sha256.*disagrees",
    ):
        cardiac_run._load_laplace_tbar_metadata(changed, output)


def test_q1_laplace_rejects_a_nonpositive_gauss_jacobian(small_closed_mesh):
    mesh = small_closed_mesh
    elems = mesh.elems.copy()
    elems[0] = elems[0, [4, 5, 6, 7, 0, 1, 2, 3]]

    with pytest.raises(ValueError, match="Gauss-point Jacobian"):
        tbar_laplace.solve_q1_hex_laplace(
            mesh.nodes,
            elems,
            np.unique(mesh.facets_endo),
            np.unique(mesh.facets_epi),
        )


def test_serial_tbar_input_rejects_a_tip_refine_mismatch(tmp_path):
    output = tmp_path / "closed_tbar.npy"
    tbar_laplace.generate_closed_mesh_tbar(
        output,
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
    )
    graded = geometry.build_closed_mesh(
        n_t=2, n_core=4, n_radial=1, core_half_width=0.36, tip_refine=2.0,
    )

    with pytest.raises(ValueError, match="tip-refine.*disagrees"):
        cardiac_run._load_laplace_tbar_metadata(graded, output)


def test_graded_mesh_tbar_field_roundtrip(tmp_path):
    output = tmp_path / "closed_tbar_graded.npy"
    graded = geometry.build_closed_mesh(
        n_t=2, n_core=4, n_radial=1, core_half_width=0.36, tip_refine=2.0,
    )
    tbar_laplace.generate_closed_mesh_tbar(
        output,
        n_t=2,
        n_core=4,
        n_radial=1,
        core_half_width=0.36,
        tip_refine=2.0,
    )

    provenance = cardiac_run._load_laplace_tbar_metadata(graded, output)
    field = cardiac_run._apply_laplace_tbar(graded, output, asign=-1.0)
    assert provenance["field_sha256"]
    assert field.shape == (graded.n_node,)
    assert np.all(np.isfinite(field))
