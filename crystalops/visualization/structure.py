"""Three-dimensional atomistic structure plotting."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from crystalops.core import CrystalStructure, UnitCell
from crystalops.modeling.builders import CrystalBuilder


class StructurePlotter:
    """Render particles and simulation-cell edges in three dimensions.

    Examples:
        >>> from crystalops.modeling import CrystalBuilder, FCC
        >>> structure = CrystalBuilder(FCC("Cu"), (2, 2, 2), 3.61).build()
        >>> plotter = StructurePlotter(structure)
        >>> figure, axes = plotter.plot(show=False)
    """

    def __init__(self, structure: CrystalStructure) -> None:
        """Store the structure rendered by subsequent plot calls.

        Args:
            structure: Particle positions, types, and simulation cell.
        """
        self.structure = structure

    @classmethod
    def from_unit_cell(
        cls,
        unit_cell: UnitCell,
        *,
        lattice_constant: float = 1.0,
    ) -> StructurePlotter:
        """Create a plotter for one crystallographic unit cell."""
        structure = CrystalBuilder(
            unit_cell,
            lattice_constant=lattice_constant,
        ).build()
        return cls(structure)

    def plot(
        self,
        *,
        atom_size: float = 40.0,
        show_cell: bool = True,
        show_axes: bool = False,
        output: str | Path | None = None,
        show: bool = True,
    ) -> tuple[Figure, Axes]:
        """Create a three-dimensional scatter plot.

        Args:
            atom_size: Marker area in points squared.
            show_cell: Draw the twelve simulation-cell edges.
            show_axes: Keep Cartesian axes visible.
            output: Optional image path. Its suffix selects the format.
            show: Display the Matplotlib window.

        Returns:
            The Matplotlib figure and three-dimensional axes.
        """
        figure = plt.figure(figsize=(7, 7))
        axes = figure.add_subplot(111, projection="3d")
        palette = plt.get_cmap("tab10")
        for color_index, type_id in enumerate(np.unique(self.structure.atom_types)):
            selected = self.structure.atom_types == type_id
            points = self.structure.positions[selected]
            label = self.structure.type_elements[int(type_id)]
            axes.scatter(
                points[:, 0],
                points[:, 1],
                points[:, 2],
                s=atom_size,
                color=palette(color_index % 10),
                label=label,
                depthshade=True,
            )
        if show_cell:
            self._draw_cell(axes)
        axes.set_box_aspect(np.ptp(self._cell_corners(), axis=0))
        axes.legend()
        if not show_axes:
            axes.set_axis_off()
        if output is not None:
            figure.savefig(Path(output), bbox_inches="tight", transparent=True)
        if show:
            plt.show()
        return figure, axes

    def _cell_corners(self) -> np.ndarray:
        """Return the eight Cartesian simulation-cell corners."""
        coefficients = np.stack(
            np.meshgrid([0.0, 1.0], [0.0, 1.0], [0.0, 1.0], indexing="ij"),
            axis=-1,
        ).reshape(-1, 3)
        return coefficients @ self.structure.cell + self.structure.origin

    def _draw_cell(self, axes: Axes) -> None:
        """Draw all cell edges on an existing three-dimensional axes."""
        corners = self._cell_corners()
        coefficients = np.stack(
            np.meshgrid([0, 1], [0, 1], [0, 1], indexing="ij"),
            axis=-1,
        ).reshape(-1, 3)
        for first in range(8):
            for second in range(first + 1, 8):
                if np.count_nonzero(coefficients[first] != coefficients[second]) == 1:
                    segment = corners[[first, second]]
                    axes.plot(segment[:, 0], segment[:, 1], segment[:, 2], color="black")
