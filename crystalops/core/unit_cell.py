"""Shared crystallographic unit-cell data model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import numpy as np
from numpy.typing import ArrayLike


@dataclass(slots=True)
class UnitCell:
    """Describe a reusable crystallographic unit-cell template.

    ``fractional_positions`` are dimensionless, while ``cell`` describes the
    cell shape relative to a lattice constant supplied by a builder.

    Examples:
        >>> from crystalops.modeling.lattices import face_centered_cubic
        >>> cell = face_centered_cubic("Cu")
        >>> cell.number_of_sites
        4
    """

    name: str
    fractional_positions: ArrayLike
    cell: ArrayLike
    atom_types: ArrayLike
    type_elements: Mapping[int, str]

    def __post_init__(self) -> None:
        """Normalize arrays and validate the unit-cell definition."""
        self.fractional_positions = np.asarray(
            self.fractional_positions,
            dtype=np.float64,
        ).copy()
        self.cell = np.asarray(self.cell, dtype=np.float64).copy()
        self.atom_types = np.asarray(self.atom_types, dtype=np.int64).copy()
        self.type_elements = {
            int(type_id): str(symbol)
            for type_id, symbol in self.type_elements.items()
        }
        self.validate()

    @property
    def number_of_sites(self) -> int:
        """Return the number of basis sites in the unit cell."""
        return int(self.fractional_positions.shape[0])

    def validate(self) -> None:
        """Validate coordinates, cell vectors, and atom type mappings."""
        if (
            self.fractional_positions.ndim != 2
            or self.fractional_positions.shape[1] != 3
        ):
            raise ValueError("fractional_positions must have shape (N, 3)")
        if self.cell.shape != (3, 3):
            raise ValueError("cell must have shape (3, 3)")
        if abs(np.linalg.det(self.cell)) <= np.finfo(float).eps:
            raise ValueError("cell vectors must define a non-zero volume")
        if self.atom_types.shape != (self.number_of_sites,):
            raise ValueError("atom_types must contain one value per site")
        if self.number_of_sites and np.any(self.atom_types < 1):
            raise ValueError("atom type identifiers must start at one")
        missing_types = set(np.unique(self.atom_types)) - set(self.type_elements)
        if missing_types:
            raise ValueError(
                f"missing element mappings for atom types {sorted(missing_types)}"
            )
