"""Area-aware comparison figures for reconstructed grain-boundary curvature."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np


def _surface_arrays(result) -> list[dict[str, Any]]:
    records = []
    for surface in result.surfaces:
        vertices = np.asarray(surface.vertices, dtype=float)
        faces = np.asarray(surface.faces)
        if vertices.ndim != 2 or vertices.shape[1:] != (3,) or not np.isfinite(vertices).all():
            raise ValueError("surface vertices must be finite Cartesian coordinates of shape (N, 3)")
        if faces.ndim != 2 or faces.shape[1:] != (3,) or not np.issubdtype(faces.dtype, np.integer):
            raise ValueError("surface faces must have shape (M, 3) and integer indices")
        if faces.size and (faces.min() < 0 or faces.max() >= len(vertices)):
            raise ValueError("surface face index is out of bounds")
        mean = np.asarray(surface.face_mean_curvature, dtype=float)
        valid = np.asarray(surface.face_valid, dtype=bool)
        area = np.asarray(surface.face_area, dtype=float)
        if any(value.shape != (len(faces),) for value in (mean, valid, area)):
            raise ValueError("face curvature, validity, and area arrays must match the face count")
        if not np.isfinite(area).all() or np.any(area < 0):
            raise ValueError("face areas must be finite and nonnegative")
        if np.any(valid & ~np.isfinite(mean)):
            raise ValueError("valid faces must have finite curvature")
        records.append({
            "vertices": vertices, "faces": faces, "absolute": np.abs(mean),
            "valid": valid, "area": area,
        })
    if not records or not sum(len(record["faces"]) for record in records):
        raise ValueError("a comparison result has no reconstructed faces")
    return records


def _sensitivity_series(rows, names, smoothing_length, sampling_spacing):
    selected = {name: [] for name in names}
    excluded = []
    lookup = {Path(str(name)).stem.casefold(): name for name in names}
    for index, row in enumerate(rows):
        raw_name = row.get("structure", row.get("name", row.get("sample")))
        name = lookup.get(Path(str(raw_name)).stem.casefold())
        if name is None:
            excluded.append({"row": index, "reason": "structure not in comparison"})
            continue
        scheme = row.get("grain_scheme")
        if scheme not in (None, "", "primary_mst", "primary_provided"):
            excluded.append({"row": index, "reason": "alternative grain-label scheme"})
            continue
        sigma = float(row.get("smoothing_length", row.get("sigma")))
        spacing = float(row.get(
            "sampling_spacing", row.get("requested_sampling_spacing", row.get("spacing"))
        ))
        if not np.isclose(spacing, sampling_spacing, rtol=1e-8, atol=1e-10):
            excluded.append({"row": index, "reason": "different sampling spacing"})
            continue
        value = row.get("rms_H")
        if value in (None, "") or not np.isfinite(float(value)):
            excluded.append({"row": index, "reason": "undefined valid-area RMS"})
            continue
        if not np.isfinite(sigma) or sigma <= 0:
            raise ValueError("sensitivity smoothing lengths must be positive and finite")
        selected[name].append((sigma, float(value)))
    for name, points in selected.items():
        points.sort()
        if len({point[0] for point in points}) != len(points):
            raise ValueError(f"duplicate primary sensitivity values for {name}")
    if not any(selected.values()):
        raise ValueError("no primary, constant-spacing sensitivity rows are available")
    return selected, excluded


def plot_grain_boundary_comparison(
    results: dict[str, Any],
    sensitivity_rows: list[dict[str, Any]],
    *,
    output_prefix: str | Path,
    smoothing_length: float,
    sampling_spacing: float,
):
    """Plot two full surface meshes, an area CDF, and smoothing-scale sensitivity.

    All mesh faces are rendered; invalid faces remain gray. Quantitative
    distributions use the supplied validity mask and face-area weights.
    The smoothing-sensitivity panel uses primary grain labels and holds the
    requested sampling spacing fixed. Triangles are quadrature samples from
    one configuration, not independent experimental replicates.

    Returns:
        ``(figure, paths)``. Paths include PNG (300 dpi), editable SVG/PDF,
        and a JSON audit of face counts, weights, color limits, and exclusions.
        The caller owns the returned figure and may close it after inspection.
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.cm import ScalarMappable
    from matplotlib.patches import Patch
    from matplotlib.ticker import MaxNLocator
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    if len(results) != 2:
        raise ValueError("the comparison requires exactly two named results")
    if not np.isfinite(smoothing_length) or smoothing_length <= 0:
        raise ValueError("smoothing_length must be positive and finite")
    if not np.isfinite(sampling_spacing) or sampling_spacing <= 0:
        raise ValueError("sampling_spacing must be positive and finite")
    for name, result in results.items():
        sigma = result.parameters.get("smoothing_length")
        spacing = result.parameters.get("requested_sampling_spacing")
        if sigma is None or not np.isclose(float(sigma), smoothing_length):
            raise ValueError(f"{name} result smoothing length does not match the stated baseline")
        if spacing is None or not np.isclose(float(spacing), sampling_spacing):
            raise ValueError(f"{name} result sampling spacing does not match the stated baseline")
    names = list(results)
    records = {name: _surface_arrays(result) for name, result in results.items()}
    series, excluded_rows = _sensitivity_series(
        sensitivity_rows, names, smoothing_length, sampling_spacing
    )
    valid_values = [
        record["absolute"][record["valid"]]
        for name in names for record in records[name]
        if np.any(record["valid"])
    ]
    if not valid_values:
        raise ValueError("neither structure has valid face curvature")
    maximum = float(np.max(np.concatenate(valid_values)))
    color_maximum = maximum if maximum > 0 else 1.0
    units = [str(result.parameters.get("length_unit", "input_length_unit")) for result in results.values()]
    if len(set(units)) != 1:
        raise ValueError("comparison results must use the same coordinate unit")
    unit = "Å" if units[0].casefold() in ("angstrom", "ångström", "å") else units[0]
    curvature_unit = r"$\mathrm{\AA}^{-1}$" if unit == "Å" else f"1/{unit}"
    colors = dict(zip(names, ("#3B6B8A", "#B66A3C")))
    invalid_color = "#C9CDD0"
    prefix = Path(output_prefix)
    if prefix.suffix.lower() in {".png", ".svg", ".pdf"}:
        prefix = prefix.with_suffix("")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    paths = {
        "png": Path(str(prefix) + ".png"),
        "svg": Path(str(prefix) + ".svg"),
        "pdf": Path(str(prefix) + ".pdf"),
        "qa": Path(str(prefix) + "_qa.json"),
    }
    audit = {
        "figure_width_mm": 183,
        "figure_height_mm": 154,
        "baseline_smoothing_length": float(smoothing_length),
        "baseline_requested_sampling_spacing": float(sampling_spacing),
        "length_unit": units[0],
        "color_mapping": "shared linear viridis of absolute mean curvature",
        "color_limits": [0.0, color_maximum],
        "color_clipped_faces": 0,
        "mesh_rendering": "all supplied faces; invalid gray; natural 3D occlusion",
        "quantitative_mask": "supplied face_valid; weights are supplied face_area",
        "sensitivity_filter": "primary labels (MST/provided/unspecified) and baseline requested spacing",
        "sensitivity_excluded_rows": excluded_rows,
        "replication": "one configuration per structure; no inferential statistics",
        "export_notes": "PNG preview and rasterized meshes use the requested 300 dpi; SVG/PDF text remains editable. No TIFF submission raster is requested.",
        "structures": {},
    }
    style = {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "font.size": 7,
        "axes.labelsize": 7,
        "axes.titlesize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "text.usetex": False,
    }
    with mpl.rc_context(style):
        # Fixed 183 x 154 mm canvas; explicit inch values also aid static preflight.
        fig = plt.figure(figsize=(7.2047244094, 6.0629921260), facecolor="white")
        grid = fig.add_gridspec(
            2, 2, left=0.085, right=0.955, bottom=0.095, top=0.95,
            height_ratios=(1.3, 1.0), wspace=0.32, hspace=0.40,
        )
        normalization = Normalize(vmin=0.0, vmax=color_maximum, clip=False)
        colormap = mpl.colormaps["viridis"]
        for column, name in enumerate(names):
            result = results[name]
            ax = fig.add_subplot(grid[0, column], projection="3d")
            cell = np.asarray(result.parameters.get("periodic_cell"), dtype=float)
            origin = np.asarray(result.parameters.get("periodic_origin", [0.0, 0.0, 0.0]), dtype=float)
            if cell.shape != (3, 3) or origin.shape != (3,) or not np.isfinite(cell).all():
                raise ValueError("result parameters must include a finite periodic_cell")
            if not np.allclose(cell, np.diag(np.diag(cell))) or np.any(np.diag(cell) <= 0):
                raise ValueError("this figure currently supports orthogonal simulation cells")
            lengths = np.diag(cell).copy()
            total_area = sum(float(record["area"].sum()) for record in records[name])
            valid_area = sum(float(record["area"][record["valid"]].sum()) for record in records[name])
            fraction = valid_area / total_area if total_area else 0.0
            face_count = sum(len(record["faces"]) for record in records[name])
            valid_count = sum(int(record["valid"].sum()) for record in records[name])
            for record in records[name]:
                face_colors = np.tile(mpl.colors.to_rgba(invalid_color), (len(record["faces"]), 1))
                face_colors[record["valid"]] = colormap(normalization(record["absolute"][record["valid"]]))
                collection = Poly3DCollection(
                    record["vertices"][record["faces"]], facecolors=face_colors,
                    edgecolors="none", linewidths=0, antialiaseds=False,
                )
                collection.set_rasterized(True)
                ax.add_collection3d(collection)
            for axis in range(3):
                other = [value for value in range(3) if value != axis]
                for side_one in (0, 1):
                    for side_two in (0, 1):
                        start = origin.copy()
                        start[other[0]] += side_one * lengths[other[0]]
                        start[other[1]] += side_two * lengths[other[1]]
                        end = start.copy()
                        end[axis] += lengths[axis]
                        ax.plot(*np.stack((start, end)).T, color="#969B9E", linewidth=0.35)
            ax.set_xlim(origin[0], origin[0] + lengths[0])
            ax.set_ylim(origin[1], origin[1] + lengths[1])
            ax.set_zlim(origin[2], origin[2] + lengths[2])
            ax.set_box_aspect(lengths)
            ax.view_init(elev=24, azim=-62)
            ax.grid(False)
            for dimension in (ax.xaxis, ax.yaxis, ax.zaxis):
                dimension.pane.fill = False
                dimension.pane.set_edgecolor("white")
                dimension.set_major_locator(MaxNLocator(nbins=2))
            ax.tick_params(pad=0, labelsize=6.5)
            ax.set_xlabel(f"x ({unit})", labelpad=0)
            ax.set_ylabel(f"y ({unit})", labelpad=0)
            ax.set_zlabel(f"z ({unit})", labelpad=0)
            grains = result.summary.get("grain_count", "?")
            ax.set_title(f"{name} · {grains} grains\nValid area {fraction:.1%}", pad=5, fontsize=8)
            ax.text2D(-0.08, 1.06, chr(ord("a") + column), transform=ax.transAxes, fontweight="bold", fontsize=8)
            audit["structures"][name] = {
                "grain_count": grains, "face_count": face_count,
                "valid_face_count": valid_count, "excluded_face_count": face_count - valid_count,
                "total_area": total_area, "valid_area": valid_area,
                "valid_area_fraction": fraction,
            }

        color_axis = fig.add_axes((0.29, 0.495, 0.40, 0.014))
        colorbar = fig.colorbar(ScalarMappable(norm=normalization, cmap=colormap), cax=color_axis, orientation="horizontal")
        colorbar.set_label(f"|H| ({curvature_unit})", fontsize=7, labelpad=1)
        colorbar.ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
        colorbar.ax.tick_params(labelsize=6.5, length=2, pad=1)
        fig.legend(
            handles=[Patch(facecolor=invalid_color, edgecolor="none", label="Excluded")],
            loc="center left", bbox_to_anchor=(0.715, 0.502), frameon=False, fontsize=7,
            handlelength=1.0, handletextpad=0.4,
        )
        fig.text(
            0.52, 0.438,
            f"Baseline: σ = {smoothing_length:g} {unit}; requested h = {sampling_spacing:g} {unit}",
            ha="center", va="center", fontsize=7,
        )
        cdf_axis = fig.add_subplot(grid[1, 0])
        scale_axis = fig.add_subplot(grid[1, 1])
        for name in names:
            values = np.concatenate([record["absolute"][record["valid"]] for record in records[name]])
            weights = np.concatenate([record["area"][record["valid"]] for record in records[name]])
            if weights.sum() > 0:
                order = np.argsort(values)
                cdf_axis.step(
                    np.r_[0.0, values[order]], np.r_[0.0, np.cumsum(weights[order]) / weights.sum()],
                    where="post", color=colors[name], linewidth=1.3, label=name,
                )
            points = series[name]
            if points:
                scale_axis.plot(
                    [point[0] for point in points], [point[1] for point in points],
                    color=colors[name], marker="o", markersize=3.5, linewidth=1.3, label=name,
                )
        cdf_axis.set_xlabel(f"|H| ({curvature_unit})")
        cdf_axis.set_ylabel("Fraction of valid interface area")
        cdf_axis.set_xlim(0, maximum * 1.04 if maximum > 0 else 1.0)
        cdf_axis.set_ylim(0, 1.02)
        cdf_axis.set_yticks((0, 0.25, 0.5, 0.75, 1.0))
        cdf_axis.xaxis.set_major_locator(MaxNLocator(nbins=4))
        cdf_axis.legend(loc="lower right", frameon=False)
        scale_axis.axvline(smoothing_length, color="#A5A9AB", linewidth=0.7, linestyle="--", zorder=0)
        scale_axis.set_xlabel(f"Smoothing length σ ({unit})")
        scale_axis.set_ylabel(f"RMS H ({curvature_unit})")
        scale_axis.set_ylim(bottom=0)
        scale_axis.xaxis.set_major_locator(MaxNLocator(nbins=4))
        scale_axis.yaxis.set_major_locator(MaxNLocator(nbins=4))
        scale_axis.legend(loc="best", frameon=False)
        for label, axis in (("c", cdf_axis), ("d", scale_axis)):
            axis.text(-0.18, 1.04, label, transform=axis.transAxes, fontweight="bold", fontsize=8)
            axis.tick_params(direction="out", length=3, width=0.7)
        fig.text(
            0.52, 0.015,
            "One configuration per structure; triangles are area quadrature samples.",
            ha="center", fontsize=6.5,
        )
        fig.savefig(paths["png"], dpi=300, facecolor="white")
        fig.savefig(paths["svg"], dpi=300, facecolor="white")
        fig.savefig(paths["pdf"], dpi=300, facecolor="white")
    audit["sensitivity_points"] = {
        name: [{"smoothing_length": sigma, "rms_H": value} for sigma, value in series[name]]
        for name in names
    }
    paths["qa"].write_text(json.dumps(audit, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return fig, paths


__all__ = ["plot_grain_boundary_comparison"]
