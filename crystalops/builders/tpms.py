"""Triply periodic minimal-surface structure builder."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from crystalops.builders.base import StructureBuilder
from crystalops.builders.crystal import CrystalBuilder
from crystalops.core import CrystalStructure, UnitCell
from crystalops.operations import select


@dataclass(slots=True)
class TpmsBuilder(StructureBuilder):
    """Select a periodic crystal with a TPMS implicit equation.

    ``mode="negative"`` and ``mode="positive"`` retain either labyrinth.
    ``mode="shell"`` keeps atoms whose absolute field value is no larger
    than ``thickness``.

    Examples:
        >>> from crystalops.builders import TpmsBuilder
        >>> from crystalops.lattices import face_centered_cubic
        >>> structure = TpmsBuilder(
        ...     face_centered_cubic("Cu"),
        ...     repeats=(6, 6, 6),
        ...     lattice_constant=3.61,
        ...     surface="gyroid",
        ... ).build()
        >>> structure.number_of_atoms > 0
        True
    """

    unit_cell: UnitCell
    repeats: tuple[int, int, int]
    lattice_constant: float = 1.0
    surface: str = "d-schwarz"
    mode: str = "negative"
    thickness: float = 0.2
    periods: tuple[float, float, float] = (1.0, 1.0, 1.0)
    name: str | None = None

    def build(self) -> CrystalStructure:
        """Build a crystal and retain atoms selected by the TPMS field."""
        structure = CrystalBuilder(
            self.unit_cell,
            repeats=self.repeats,
            lattice_constant=self.lattice_constant,
            name=self.name or self.surface,
        ).build()
        coordinates = structure.fractional_positions * (
            2.0 * np.pi * np.asarray(self.periods)
        )
        x, y, z = coordinates.T
        key = self.surface.lower().replace("_", "-")
        if key in {"d", "d-schwarz", "diamond"}:
            field = (
                np.sin(x) * np.sin(y) * np.sin(z)
                + np.sin(x) * np.cos(y) * np.cos(z)
                + np.cos(x) * np.sin(y) * np.cos(z)
                + np.cos(x) * np.cos(y) * np.sin(z)
            )
        elif key in {"p", "p-schwarz", "primitive"}:
            field = np.cos(x) + np.cos(y) + np.cos(z)
        elif key == "gyroid":
            field = (
                np.sin(x) * np.cos(y)
                + np.sin(y) * np.cos(z)
                + np.sin(z) * np.cos(x)
            )
        elif key in {"iwp", "i-wp"}:
            field = 2.0 * (
                np.cos(x) * np.cos(y)
                + np.cos(y) * np.cos(z)
                + np.cos(z) * np.cos(x)
            ) - (np.cos(2.0 * x) + np.cos(2.0 * y) + np.cos(2.0 * z))
        else:
            raise ValueError("surface must be d-schwarz, p-schwarz, gyroid, or iwp")

        if self.mode == "negative":
            mask = field <= 0.0
        elif self.mode == "positive":
            mask = field >= 0.0
        elif self.mode == "shell":
            if self.thickness < 0:
                raise ValueError("thickness must be non-negative")
            mask = np.abs(field) <= self.thickness
        else:
            raise ValueError("mode must be negative, positive, or shell")
        result = select(structure, mask)
        result.metadata.update(
            {
                "builder": type(self).__name__,
                "surface": key,
                "mode": self.mode,
            }
        )
        return result
