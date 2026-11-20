"""Render DXA defect meshes without periodic caps or face subsampling."""

from __future__ import annotations

from itertools import product
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from crystalops.core import CrystalStructure


def _mesh_arrays(result: Any):
    """Read the already clipped triangles; never triangulate or close a mesh."""
    vertices = np.asarray(result.clipped_vertices, dtype=float)
    faces = np.asarray(result.clipped_faces)
    cell = np.asarray(result.cell, dtype=float)
    origin = np.asarray(result.origin, dtype=float)
    pbc = np.asarray(result.pbc, dtype=bool)
    if vertices.shape == (0,):
        vertices = vertices.reshape(0, 3)
    if faces.shape == (0,):
        faces = np.empty((0, 3), dtype=np.int64)
    if vertices.ndim != 2 or vertices.shape[1:] != (3,) or not np.isfinite(vertices).all():
        raise ValueError("clipped_vertices must be finite Cartesian coordinates with shape (N, 3)")
    if faces.ndim != 2 or faces.shape[1:] != (3,) or not np.issubdtype(faces.dtype, np.integer):
        raise ValueError("clipped_faces must have shape (M, 3) and integer indices")
    if faces.size and (faces.min() < 0 or faces.max() >= len(vertices)):
        raise ValueError("clipped_faces contain an out-of-bounds vertex index")
    if cell.shape != (3, 3) or not np.isfinite(cell).all() or abs(np.linalg.det(cell)) == 0:
        raise ValueError("result.cell must contain three finite, independent row vectors")
    if origin.shape != (3,) or not np.isfinite(origin).all() or pbc.shape != (3,):
        raise ValueError("result.origin and result.pbc must each contain three values")
    return vertices, faces, cell, origin, pbc


def _cell_edges(cell: np.ndarray, origin: np.ndarray):
    """Return exact cell corners and twelve edges, including for tilted cells."""
    fractional = np.asarray(list(product((0.0, 1.0), repeat=3)))
    corners = origin + fractional @ cell
    edges = []
    for index, corner in enumerate(fractional):
        for axis in range(3):
            if corner[axis] == 0:
                other = corner.copy()
                other[axis] = 1
                target = int(np.flatnonzero(np.all(fractional == other, axis=1))[0])
                edges.append(corners[[index, target]])
    return corners, edges


def _wrapped_atoms(structure: Any, cell: np.ndarray, origin: np.ndarray, pbc: np.ndarray):
    structure_cell = np.asarray(structure.cell, dtype=float)
    structure_origin = np.asarray(structure.origin, dtype=float)
    structure_pbc = np.asarray(structure.pbc, dtype=bool)
    if (
        structure_cell.shape != (3, 3)
        or not np.isfinite(structure_cell).all()
        or not np.allclose(structure_cell, cell, rtol=1e-8, atol=1e-8)
    ):
        raise ValueError("atom overlay structure.cell does not match the defect-mesh cell")
    if (
        structure_origin.shape != (3,)
        or not np.isfinite(structure_origin).all()
        or not np.allclose(structure_origin, origin, rtol=1e-8, atol=1e-8)
    ):
        raise ValueError("atom overlay structure.origin does not match the defect-mesh origin")
    if structure_pbc.shape != (3,) or not np.array_equal(structure_pbc, pbc):
        raise ValueError("atom overlay structure.pbc does not match the defect-mesh periodic boundaries")
    positions = np.asarray(structure.positions, dtype=float)
    if positions.ndim != 2 or positions.shape[1:] != (3,) or not np.isfinite(positions).all():
        raise ValueError("structure.positions must be finite Cartesian coordinates with shape (N, 3)")
    fractional = np.linalg.solve(cell.T, (positions - origin).T).T
    fractional[:, pbc] %= 1.0
    return origin + fractional @ cell


def _smoothing_label(parameters: dict[str, Any]) -> str:
    for key in ("surface_smoothing", "defect_mesh_smoothing_level", "mesh_smoothing_level", "smoothing_level"):
        if key in parameters:
            return f"mesh smoothing: {parameters[key]}"
    return "mesh smoothing: unspecified"


def _interactive_figure(
    vertices, faces, atoms, edges, bounds, labels, title, path: Path,
):
    """Import Plotly only for callers requesting an interactive export."""
    try:
        import plotly.graph_objects as go
    except ImportError as error:
        raise ImportError(
            "Interactive defect-mesh export requires Plotly; install plotly or use interactive=False."
        ) from error

    figure = go.Figure()
    if len(faces):
        figure.add_trace(go.Mesh3d(
            x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
            i=faces[:, 0], j=faces[:, 1], k=faces[:, 2],
            color="#269FB6", opacity=1.0, flatshading=True,
            lighting={"ambient": 0.55, "diffuse": 0.7, "specular": 0.18, "roughness": 0.8},
            name="DXA defect mesh", showlegend=True,
            hovertemplate="x: %{x:.4g}<br>y: %{y:.4g}<br>z: %{z:.4g}<extra>DXA defect mesh</extra>",
        ))
    if atoms is not None and len(atoms):
        figure.add_trace(go.Scatter3d(
            x=atoms[:, 0], y=atoms[:, 1], z=atoms[:, 2], mode="markers",
            marker={"size": 2, "color": "#737A80", "opacity": 0.18},
            name="Atoms (toggle overlay)", showlegend=True,
            hovertemplate="x: %{x:.4g}<br>y: %{y:.4g}<br>z: %{z:.4g}<extra>Atom</extra>",
        ))
    edge_coordinates = np.concatenate(
        [np.vstack((edge, np.full((1, 3), np.nan))) for edge in edges], axis=0,
    )
    figure.add_trace(go.Scatter3d(
        x=edge_coordinates[:, 0], y=edge_coordinates[:, 1], z=edge_coordinates[:, 2],
        mode="lines", line={"color": "#66717B", "width": 3},
        name="Simulation cell", showlegend=True, hoverinfo="skip",
    ))
    scene = {"aspectmode": "data", "camera": {"eye": {"x": 1.45, "y": 1.45, "z": 1.0}}}
    for index, key in enumerate(("xaxis", "yaxis", "zaxis")):
        scene[key] = {
            "title": labels[index], "range": bounds[:, index].tolist(),
            "backgroundcolor": "white", "gridcolor": "#E6EAED", "showbackground": True,
            "zeroline": False,
        }
    figure.update_layout(
        title={"text": title.replace("\n", "<br>"), "x": 0.03},
        scene=scene, template="plotly_white", margin={"l": 0, "r": 0, "b": 0, "t": 85},
        legend={"x": 0.01, "y": 0.98, "bgcolor": "rgba(255,255,255,0.8)"},
        font={"family": "Arial, sans-serif", "size": 13},
        uirevision="defect-mesh",
    )
    if not len(faces):
        figure.add_annotation(
            text="No defect mesh", x=0.5, y=0.52, xref="paper", yref="paper",
            showarrow=False, font={"size": 20, "color": "#52616B"},
        )
    figure.write_html(
        str(path), include_plotlyjs=True, full_html=True, auto_open=False,
        config={"responsive": True, "displaylogo": False, "scrollZoom": True},
    )


def plot_defect_mesh(
    result: Any,
    *,
    output_prefix: str | Path,
    structure: CrystalStructure | None = None,
    show_atoms: bool = False,
    interactive: bool = True,
):
    """Save a DXA defect-mesh PNG and, optionally, a self-contained HTML view.

    Every supplied ``clipped_faces`` triangle is rendered. These triangles
    must already be clipped at periodic cell boundaries without cap faces;
    native periodic faces are deliberately not used for display. No geometric
    smoothing, decimation, surface closing, or curvature calculation occurs.

    Atom overlays require matching cell vectors, origin, and periodic flags,
    and include every supplied atom, wrapped only along periodic directions.
    Plotly is needed only when ``interactive``
    is true. The offline HTML supports rotation, zoom, and legend toggles.

    Returns:
        ``(figure, paths)`` where paths contains a PNG and optional HTML path.
        The caller owns the Matplotlib figure and may close it after viewing.
    """
    if show_atoms and structure is None:
        raise ValueError("show_atoms=True requires structure")
    vertices, faces, cell, origin, pbc = _mesh_arrays(result)
    atoms = _wrapped_atoms(structure, cell, origin, pbc) if show_atoms else None
    corners, edges = _cell_edges(cell, origin)
    extent_arrays = [corners]
    if len(vertices):
        extent_arrays.append(vertices)
    if atoms is not None and len(atoms):
        extent_arrays.append(atoms)
    points = np.concatenate(extent_arrays)
    lower, upper = points.min(axis=0), points.max(axis=0)
    span = upper - lower
    padding = span * 0.025
    bounds = np.vstack((lower - padding, upper + padding))
    parameters = dict(result.parameters)
    unit = str(parameters.get("length_unit", "input_length_unit"))
    if unit.casefold() in {"angstrom", "ångström", "å"}:
        unit = "Å"
    labels = [f"{axis} ({unit})" for axis in ("x", "y", "z")]
    title = f"DXA defect mesh | {_smoothing_label(parameters)}\n{len(faces):,} displayed triangles; periodic cuts without caps"
    prefix = Path(output_prefix)
    if prefix.suffix.lower() in {".png", ".html"}:
        prefix = prefix.with_suffix("")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    paths = {"png": Path(str(prefix) + ".png")}

    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    figure = plt.figure(figsize=(9.0, 7.5), facecolor="white")
    axes = figure.add_subplot(111, projection="3d")
    if len(faces):
        triangles = vertices[faces]
        normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
        lengths = np.linalg.norm(normals, axis=1)
        normals = np.divide(normals, lengths[:, None], out=np.zeros_like(normals), where=lengths[:, None] > 0)
        light = np.asarray([0.3, -0.5, 0.8])
        light /= np.linalg.norm(light)
        illumination = 0.52 + 0.48 * np.abs(normals @ light)
        base_color = np.asarray([0.15, 0.62, 0.71])
        face_colors = illumination[:, None] * base_color
        axes.add_collection3d(Poly3DCollection(
            triangles, facecolors=face_colors, edgecolors="none", linewidths=0,
            antialiaseds=False,
        ))
    else:
        axes.text2D(0.5, 0.52, "No defect mesh", transform=axes.transAxes,
                    ha="center", va="center", fontsize=16, color="#52616B")
    if atoms is not None and len(atoms):
        axes.scatter(
            atoms[:, 0], atoms[:, 1], atoms[:, 2], s=1.2, c="#737A80",
            alpha=0.18, linewidths=0, depthshade=False,
        )
    for edge in edges:
        axes.plot(*edge.T, color="#66717B", linewidth=0.8, alpha=0.8)
    axes.set_xlim(*bounds[:, 0])
    axes.set_ylim(*bounds[:, 1])
    axes.set_zlim(*bounds[:, 2])
    axes.set_box_aspect(bounds[1] - bounds[0])
    axes.set_xlabel(labels[0], labelpad=9)
    axes.set_ylabel(labels[1], labelpad=9)
    axes.set_zlabel(labels[2], labelpad=9)
    for axis in (axes.xaxis, axes.yaxis, axes.zaxis):
        axis.set_major_locator(MaxNLocator(nbins=4))
        axis.pane.fill = False
    axes.tick_params(labelsize=9)
    axes.view_init(elev=23, azim=-54)
    axes.set_title(title, fontsize=12, pad=20)
    handles = [Line2D([0], [0], color="#66717B", linewidth=1, label="Simulation cell")]
    if len(faces):
        handles.insert(0, Patch(facecolor="#269FB6", edgecolor="none", label="DXA defect mesh"))
    if atoms is not None and len(atoms):
        handles.append(Line2D([0], [0], marker="o", linestyle="", color="#737A80", markersize=4, label="Atoms"))
    axes.legend(handles=handles, loc="upper left", frameon=False, fontsize=9)
    figure.subplots_adjust(left=0.02, right=0.84, bottom=0.06, top=0.88)
    try:
        # Axes3D's automatic tight bounding box can omit the projected z label.
        # Include all labels explicitly and reserve room for their projection.
        figure.savefig(
            paths["png"], dpi=220, facecolor="white", bbox_inches="tight",
            bbox_extra_artists=(axes.xaxis.label, axes.yaxis.label, axes.zaxis.label),
            pad_inches=0.25,
        )
        if interactive:
            paths["html"] = Path(str(prefix) + ".html")
            _interactive_figure(vertices, faces, atoms, edges, bounds, labels, title, paths["html"])
    except Exception:
        plt.close(figure)
        raise
    return figure, paths


__all__ = ["plot_defect_mesh"]
