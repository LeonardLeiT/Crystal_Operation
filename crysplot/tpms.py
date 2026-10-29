"""Implicit-surface plotting for common TPMS equations."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from skimage import measure


class TpmsPlotter:
    """Render a triply periodic minimal surface with marching cubes.

    Examples:
        >>> plotter = TpmsPlotter("gyroid", resolution=40)
        >>> figure, axes = plotter.plot(show=False)
    """

    def __init__(
        self,
        surface: str = "d-schwarz",
        *,
        resolution: int = 50,
        periods: tuple[float, float, float] = (1.0, 1.0, 1.0),
    ) -> None:
        """Configure an implicit surface grid.

        Args:
            surface: ``d-schwarz``, ``p-schwarz``, ``gyroid``, or ``iwp``.
            resolution: Number of grid nodes along each dimension.
            periods: Number of periods along x, y, and z.
        """
        if resolution < 3:
            raise ValueError("resolution must be at least three")
        self.surface = surface.lower().replace("_", "-")
        self.resolution = resolution
        self.periods = periods

    def field(self) -> np.ndarray:
        """Evaluate the configured implicit field on a regular cubic grid."""
        axes = [
            np.linspace(0.0, 2.0 * np.pi * period, self.resolution)
            for period in self.periods
        ]
        x, y, z = np.meshgrid(*axes, indexing="ij")
        if self.surface in {"d", "d-schwarz", "diamond"}:
            return (
                np.sin(x) * np.sin(y) * np.sin(z)
                + np.sin(x) * np.cos(y) * np.cos(z)
                + np.cos(x) * np.sin(y) * np.cos(z)
                + np.cos(x) * np.cos(y) * np.sin(z)
            )
        if self.surface in {"p", "p-schwarz", "primitive"}:
            return np.cos(x) + np.cos(y) + np.cos(z)
        if self.surface == "gyroid":
            return (
                np.sin(x) * np.cos(y)
                + np.sin(y) * np.cos(z)
                + np.sin(z) * np.cos(x)
            )
        if self.surface in {"iwp", "i-wp"}:
            return 2.0 * (
                np.cos(x) * np.cos(y)
                + np.cos(y) * np.cos(z)
                + np.cos(z) * np.cos(x)
            ) - (np.cos(2.0 * x) + np.cos(2.0 * y) + np.cos(2.0 * z))
        raise ValueError("surface must be d-schwarz, p-schwarz, gyroid, or iwp")

    def plot(
        self,
        *,
        level: float = 0.0,
        output: str | Path | None = None,
        show: bool = True,
    ) -> tuple[Figure, Axes]:
        """Extract and render an isosurface at ``level``."""
        field = self.field()
        spacing = tuple(period / (self.resolution - 1) for period in self.periods)
        vertices, faces, _, _ = measure.marching_cubes(field, level=level, spacing=spacing)
        figure = plt.figure(figsize=(7, 7))
        axes = figure.add_subplot(111, projection="3d")
        axes.plot_trisurf(
            vertices[:, 0],
            vertices[:, 1],
            vertices[:, 2],
            triangles=faces,
            cmap="Spectral",
            linewidth=0.1,
        )
        axes.set_box_aspect(self.periods)
        axes.set_axis_off()
        if output is not None:
            figure.savefig(Path(output), bbox_inches="tight", transparent=True)
        if show:
            plt.show()
        return figure, axes
