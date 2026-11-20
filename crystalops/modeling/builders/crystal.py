"""Builders for periodic crystals based on unit-cell templates."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from crystalops.modeling.builders.base import StructureBuilder
from crystalops.core import CrystalStructure, UnitCell
from crystalops.modeling.operations import repeat


@dataclass(slots=True)
class CrystalBuilder(StructureBuilder):
    """Construct a periodic supercell from a unit-cell template.

    Examples:
        >>> from crystalops.modeling.builders import CrystalBuilder
        >>> from crystalops.modeling.lattices import face_centered_cubic
        >>> structure = CrystalBuilder(
        ...     face_centered_cubic("Cu"),
        ...     repeats=(2, 2, 2),
        ...     lattice_constant=3.61,
        ... ).build()
        >>> structure.number_of_atoms
        32
    """

    unit_cell: UnitCell
    repeats: tuple[int, int, int] = (1, 1, 1)
    lattice_constant: float = 1.0
    name: str | None = None

    def build(self) -> CrystalStructure:
        """Build the requested supercell.

        Returns:
            A periodic structure with Cartesian positions.
        """
        if self.lattice_constant <= 0:
            raise ValueError("lattice_constant must be positive")
        cell = self.unit_cell.cell * self.lattice_constant
        primitive = CrystalStructure.from_fractional(
            fractional_positions=self.unit_cell.fractional_positions,
            cell=cell,
            atom_types=self.unit_cell.atom_types,
            type_elements=self.unit_cell.type_elements,
            name=self.name or self.unit_cell.name,
            metadata={"lattice_constant": float(self.lattice_constant)},
        )
        return repeat(primitive, self.repeats)
