"""Owned-row consistent mass assembly for the Cardiac MPI companion.

The global mesh is currently replicated, but the consistent inertia matrix is
not.  Each rank assembles only its PETSc-owned rows.  An owned row receives
contributions from *every* Hex8 that touches that row, including elements whose
material block is assigned to another rank.  This distinction is essential:
partitioning mass by the material-element owner would omit shared-node terms.

The element integration is delegated to :func:`consistent_mass_coo`, the same
routine used by the serial driver.  This module only selects all touching
elements, retains the owned global rows, and sums them into a local-row/global-
column CSR matrix.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from consistent_mass import consistent_mass_coo


OWNED_ROW_ASSEMBLY_POLICY = "owned-row-csr-all-touching-elements"


def balanced_row_range(ndof: int, rank: int, size: int) -> tuple[int, int]:
    """Return a deterministic contiguous partition used by non-PETSc tests."""
    if not isinstance(ndof, int) or ndof < 0:
        raise ValueError("ndof must be a nonnegative integer")
    if not isinstance(size, int) or size < 1:
        raise ValueError("size must be a positive integer")
    if not isinstance(rank, int) or rank < 0 or rank >= size:
        raise ValueError("rank must lie in [0, size)")
    quotient, remainder = divmod(ndof, size)
    start = rank * quotient + min(rank, remainder)
    stop = start + quotient + (1 if rank < remainder else 0)
    return start, stop


@dataclass(frozen=True)
class OwnedRowMassMatrix:
    """One rank's global-column CSR rows of a consistent mass matrix."""

    matrix: sp.csr_matrix
    row_start: int
    row_end: int
    ndof: int
    touching_element_ids: np.ndarray

    def __post_init__(self):
        matrix = self.matrix
        if not sp.isspmatrix_csr(matrix):
            raise TypeError("owned mass matrix must be CSR")
        row_start = int(self.row_start)
        row_end = int(self.row_end)
        ndof = int(self.ndof)
        if ndof < 1 or row_start < 0 or row_end < row_start or row_end > ndof:
            raise ValueError("owned mass row range is outside the global system")
        if matrix.shape != (row_end - row_start, ndof):
            raise ValueError(
                "owned mass CSR must have shape "
                f"{(row_end - row_start, ndof)}"
            )
        if not np.all(np.isfinite(matrix.data)) or np.any(matrix.data < 0.0):
            raise ValueError("owned mass CSR entries must be finite and nonnegative")
        identifiers = np.asarray(self.touching_element_ids, dtype=np.int64)
        if identifiers.ndim != 1 or (
            len(identifiers) and np.min(identifiers) < 0
        ):
            raise ValueError("touching element IDs must be nonnegative and 1-D")
        if len(np.unique(identifiers)) != len(identifiers):
            raise ValueError("touching element IDs must be unique")

        canonical = matrix.copy()
        canonical.sum_duplicates()
        canonical.sort_indices()
        object.__setattr__(self, "matrix", canonical)
        object.__setattr__(self, "row_start", row_start)
        object.__setattr__(self, "row_end", row_end)
        object.__setattr__(self, "ndof", ndof)
        object.__setattr__(self, "touching_element_ids", identifiers.copy())

    @property
    def representation(self) -> str:
        return "consistent_q1_hex8"

    @property
    def assembly_policy(self) -> str:
        return OWNED_ROW_ASSEMBLY_POLICY

    def action(self, global_vector) -> np.ndarray:
        """Apply the owned rows to a finite replicated global vector."""
        vector = np.asarray(global_vector, dtype=float)
        if vector.shape != (self.ndof,) or not np.all(np.isfinite(vector)):
            raise ValueError(
                f"mass action requires a finite vector with shape {(self.ndof,)}"
            )
        result = np.asarray(self.matrix @ vector, dtype=float).reshape(-1)
        if result.shape != (self.row_end - self.row_start,) or not np.all(
            np.isfinite(result)
        ):
            raise RuntimeError("owned consistent-mass action is invalid")
        return result

    def metadata(self) -> dict:
        """Return JSON-safe rank-local assembly provenance."""
        return {
            "representation": self.representation,
            "partition": self.assembly_policy,
            "row_start": self.row_start,
            "row_end": self.row_end,
            "owned_rows": self.row_end - self.row_start,
            "global_columns": self.ndof,
            "local_nnz": int(self.matrix.nnz),
            "touching_elements": int(len(self.touching_element_ids)),
        }


def assemble_owned_consistent_mass(
    nodes,
    elements,
    density,
    ndof,
    row_start,
    row_end,
    *,
    dof_per_node=3,
) -> OwnedRowMassMatrix:
    """Assemble complete consistent-mass contributions for owned rows.

    Material cells may be partitioned independently.  This routine therefore
    selects elements by row incidence, not by material ownership.  For an
    owned row interval ``[row_start, row_end)``, every element containing at
    least one owned DOF is integrated and then filtered to those rows.
    """
    nodes = np.asarray(nodes, dtype=float)
    elements = np.asarray(elements, dtype=np.int64)
    if nodes.ndim != 2 or nodes.shape[1] != 3 or not np.all(np.isfinite(nodes)):
        raise ValueError("nodes must have finite shape (n_node, 3)")
    if elements.ndim != 2 or elements.shape[1] != 8:
        raise ValueError("elements must have shape (n_element, 8)")
    if len(elements) and (
        np.min(elements) < 0 or np.max(elements) >= len(nodes)
    ):
        raise ValueError("elements contain a node outside the mesh")
    if not isinstance(dof_per_node, int) or dof_per_node < 1:
        raise ValueError("dof_per_node must be a positive integer")
    ndof = int(ndof)
    row_start = int(row_start)
    row_end = int(row_end)
    if ndof != len(nodes) * dof_per_node:
        raise ValueError("ndof does not match nodes and dof_per_node")
    if row_start < 0 or row_end < row_start or row_end > ndof:
        raise ValueError("owned row range is outside the global system")

    element_dofs = (
        elements[:, :, None] * dof_per_node
        + np.arange(dof_per_node, dtype=np.int64)[None, None, :]
    ).reshape(len(elements), -1)
    touching = np.any(
        (element_dofs >= row_start) & (element_dofs < row_end), axis=1
    )
    touching_ids = np.flatnonzero(touching).astype(np.int64)

    local_row_count = row_end - row_start
    if local_row_count == 0 or len(touching_ids) == 0:
        matrix = sp.csr_matrix((local_row_count, ndof), dtype=float)
    else:
        rows, columns, values = consistent_mass_coo(
            nodes,
            elements[touching_ids],
            density,
            dof_per_node=dof_per_node,
        )
        owned = (rows >= row_start) & (rows < row_end)
        matrix = sp.coo_matrix(
            (
                np.asarray(values[owned], dtype=float),
                (
                    np.asarray(rows[owned] - row_start, dtype=np.int64),
                    np.asarray(columns[owned], dtype=np.int64),
                ),
            ),
            shape=(local_row_count, ndof),
        ).tocsr()

    return OwnedRowMassMatrix(
        matrix=matrix,
        row_start=row_start,
        row_end=row_end,
        ndof=ndof,
        touching_element_ids=touching_ids,
    )
