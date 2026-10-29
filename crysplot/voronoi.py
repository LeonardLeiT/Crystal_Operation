"""Voronoi polyhedron plotting."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from numpy.typing import ArrayLike


class VoronoiPlotter:
    """Render a polyhedron from vertices and indexed polygonal faces."""

    def __init__(self, vertices: ArrayLike, faces: list[list[int]]) -> None:
        """Store polyhedron geometry for subsequent rendering."""
        self.vertices = np.asarray(vertices, dtype=np.float64)
        self.faces = [list(face) for face in faces]
        if self.vertices.ndim != 2 or self.vertices.shape[1] != 3:
            raise ValueError("vertices must have shape (N, 3)")
        if any(len(face) < 3 for face in self.faces):
            raise ValueError("each face must contain at least three vertices")

    def plot(
        self,
        *,
        face_color: str = "skyblue",
        edge_color: str = "black",
        alpha: float = 0.5,
        output: str | Path | None = None,
        show: bool = True,
    ) -> tuple[Figure, Axes]:
        """Render the stored polyhedron and optionally save it."""
        figure = plt.figure(figsize=(7, 7))
        axes = figure.add_subplot(111, projection="3d")
        polygons = [[self.vertices[index] for index in face] for face in self.faces]
        axes.add_collection3d(
            Poly3DCollection(
                polygons,
                facecolors=face_color,
                edgecolors=edge_color,
                alpha=alpha,
            )
        )
        axes.scatter(*self.vertices.T, color="black", s=15)
        center = 0.5 * (self.vertices.min(axis=0) + self.vertices.max(axis=0))
        half_range = 0.5 * np.ptp(self.vertices, axis=0).max()
        axes.set_xlim(center[0] - half_range, center[0] + half_range)
        axes.set_ylim(center[1] - half_range, center[1] + half_range)
        axes.set_zlim(center[2] - half_range, center[2] + half_range)
        axes.set_box_aspect((1.0, 1.0, 1.0))
        axes.set_axis_off()
        if output is not None:
            figure.savefig(Path(output), bbox_inches="tight", transparent=True)
        if show:
            plt.show()
        return figure, axes
