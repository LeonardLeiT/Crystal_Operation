"""Internal lattice-generation helpers shared by structure builders."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from crystalops.core import UnitCell


FloatArray = NDArray[np.float64]


def candidate_lattice(
    unit_cell: UnitCell,
    lattice_constant: float,
    box: FloatArray,
    rotation: FloatArray,
) -> FloatArray:
    """Generate a centered, rotated crystal large enough to cover ``box``."""
    primitive = unit_cell.cell * lattice_constant
    shortest = np.min(np.linalg.norm(primitive, axis=1))
    corners = np.stack(
        np.meshgrid([-0.5, 0.5], [-0.5, 0.5], [-0.5, 0.5], indexing="ij"),
        axis=-1,
    ).reshape(-1, 3)
    radius = np.max(np.linalg.norm(corners @ box, axis=1))
    extent = int(np.ceil(radius / shortest)) + 1
    indices = np.stack(
        np.meshgrid(
            np.arange(-extent, extent + 1),
            np.arange(-extent, extent + 1),
            np.arange(-extent, extent + 1),
            indexing="ij",
        ),
        axis=-1,
    ).reshape(-1, 3)
    fractional = (
        indices[:, None, :] + unit_cell.fractional_positions[None, :, :]
    ).reshape(-1, 3)
    return (fractional @ primitive) @ rotation.T


def orientation_from_rows(rows: ArrayLike) -> FloatArray:
    """Normalize an orthogonal crystallographic orientation matrix."""
    matrix = np.asarray(rows, dtype=np.float64)
    normalized = matrix / np.linalg.norm(matrix, axis=1, keepdims=True)
    if not np.allclose(normalized @ normalized.T, np.eye(3), atol=1e-10):
        raise ValueError("orientation rows must be mutually orthogonal")
    return normalized


__all__ = ["candidate_lattice", "orientation_from_rows"]
