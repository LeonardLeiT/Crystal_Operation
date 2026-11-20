"""Export grain-boundary curvature results without a visualization dependency."""

from __future__ import annotations

import csv
from dataclasses import asdict, is_dataclass
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from crystalops.core import CrystalStructure


def _json_value(value: Any) -> Any:
    """Convert NumPy values and non-finite measurements to strict JSON."""
    if is_dataclass(value) and not isinstance(value, type):
        return _json_value(asdict(value))
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _array(value: Any, shape: tuple[int, ...], name: str, dtype=None) -> np.ndarray:
    result = np.asarray(value, dtype=dtype)
    if result.shape != shape:
        raise ValueError(f"{name} must have shape {shape}, got {result.shape}")
    return result


def _validated_surfaces(surfaces) -> list[dict[str, Any]]:
    records = []
    for index, surface in enumerate(surfaces):
        vertices = np.asarray(surface.vertices, dtype=float)
        faces = np.asarray(surface.faces)
        if vertices.ndim != 2 or vertices.shape[1:] != (3,):
            raise ValueError("surface vertices must have shape (N, 3)")
        if not np.isfinite(vertices).all():
            raise ValueError("surface vertex coordinates must be finite")
        if faces.ndim != 2 or faces.shape[1:] != (3,):
            raise ValueError("surface faces must have shape (M, 3)")
        if not np.issubdtype(faces.dtype, np.integer):
            raise ValueError("surface faces must contain integer vertex indices")
        if faces.size and (faces.min() < 0 or faces.max() >= len(vertices)):
            raise ValueError("surface face index is out of range")
        count = len(faces)
        pair = np.asarray(surface.grain_pair)
        if pair.shape != (2,) or not np.issubdtype(pair.dtype, np.integer):
            raise ValueError("surface grain_pair must contain two integer labels")
        mean = _array(surface.face_mean_curvature, (count,), "face_mean_curvature", float)
        valid = _array(surface.face_valid, (count,), "face_valid", bool)
        area = _array(surface.face_area, (count,), "face_area", float)
        if not np.isfinite(area).all() or np.any(area < 0):
            raise ValueError("surface face areas must be finite and nonnegative")
        if np.any(valid & ~np.isfinite(mean)):
            raise ValueError("valid surface faces must have finite mean curvature")
        centers = getattr(surface, "face_centers", None)
        centers = vertices[faces].mean(axis=1) if centers is None else _array(
            centers, (count, 3), "face_centers", float
        )
        gaussian = _array(
            getattr(surface, "face_gaussian_curvature", np.full(count, np.nan)),
            (count,), "face_gaussian_curvature", float,
        )
        junction = _array(
            getattr(surface, "face_junction", np.zeros(count, dtype=bool)),
            (count,), "face_junction", bool,
        )
        records.append({
            "index": index, "pair": pair.astype(np.int64),
            "vertices": vertices, "faces": faces.astype(np.int64),
            "mean": mean, "gaussian": gaussian, "valid": valid,
            "area": area, "centers": centers, "junction": junction,
        })
    return records


def _write_surface_vtk(path: Path, records: list[dict[str, Any]]) -> None:
    """Write standard ASCII VTK PolyData with measurements on triangle cells."""
    vertex_count = sum(len(record["vertices"]) for record in records)
    face_count = sum(len(record["faces"]) for record in records)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# vtk DataFile Version 3.0\n")
        handle.write("CrystalOps grain-boundary curvature; invalid cells carry NaN\n")
        handle.write("ASCII\nDATASET POLYDATA\n")
        handle.write(f"POINTS {vertex_count} double\n")
        for record in records:
            np.savetxt(handle, record["vertices"], fmt="%.17g")
        handle.write(f"POLYGONS {face_count} {4 * face_count}\n")
        offset = 0
        for record in records:
            faces = record["faces"] + offset
            rows = np.column_stack((np.full(len(faces), 3, dtype=int), faces))
            np.savetxt(handle, rows, fmt="%d")
            offset += len(record["vertices"])
        if not face_count:
            return
        handle.write(f"CELL_DATA {face_count}\n")
        fields = (
            ("surface_index", "int", lambda r: np.full(len(r["faces"]), r["index"])),
            ("grain_i", "long", lambda r: np.full(len(r["faces"]), r["pair"][0])),
            ("grain_j", "long", lambda r: np.full(len(r["faces"]), r["pair"][1])),
            ("valid", "unsigned_char", lambda r: r["valid"].astype(int)),
            ("junction", "unsigned_char", lambda r: r["junction"].astype(int)),
            ("area", "double", lambda r: r["area"]),
            ("mean_curvature", "double", lambda r: r["mean"]),
            ("gaussian_curvature", "double", lambda r: r["gaussian"]),
        )
        for name, kind, values in fields:
            handle.write(f"SCALARS {name} {kind} 1\nLOOKUP_TABLE default\n")
            for record in records:
                np.savetxt(handle, values(record), fmt="%.17g" if kind == "double" else "%d")


def write_grain_boundary_results(
    result,
    structure: CrystalStructure,
    output_dir: str | Path,
    *,
    grain_ids=None,
) -> dict[str, Path]:
    """Write atom/face CSV, an annotated surface VTK, and a strict JSON summary.

    ``result`` is the object returned by the grain-boundary analysis API. The
    original atom order and IDs are retained. The surface mesh is exported as
    supplied; it must already have consistent periodic geometry. No artificial
    box-face caps are added. Invalid curvature remains NaN in numerical files
    and is accompanied by an explicit validity flag; JSON uses ``null``.

    Lengths use the input structure's coordinate unit. Curvature therefore has
    units of inverse input length. This writer does not infer physical units
    from the file extension.
    """
    count = structure.number_of_atoms
    mean = _array(result.atom_mean_curvature, (count,), "atom_mean_curvature", float)
    valid = _array(result.atom_valid, (count,), "atom_valid", bool)
    pairs = _array(result.atom_grain_pairs, (count, 2), "atom_grain_pairs", np.int64)
    boundary = _array(result.atom_boundary_mask, (count,), "atom_boundary_mask", bool)
    junction = _array(result.atom_junction_mask, (count,), "atom_junction_mask", bool)
    distance = _array(result.atom_surface_distance, (count,), "atom_surface_distance", float)
    if np.any(valid & ~np.isfinite(mean)):
        raise ValueError("valid atoms must have finite mean curvature")
    atom_ids = _array(
        structure.properties.get("atom_ids", np.arange(1, count + 1)),
        (count,), "atom_ids", np.int64,
    )
    if grain_ids is None:
        grain_ids = structure.properties.get("grain_ids", getattr(result, "grain_ids", None))
    if grain_ids is None:
        raise ValueError("grain_ids must be supplied or stored on the input structure")
    labels = _array(grain_ids, (count,), "grain_ids", np.int64)
    records = _validated_surfaces(result.surfaces)
    payload = _json_value({
        "schema_version": 1,
        "structure": {
            "name": structure.name,
            "source": structure.metadata.get("source"),
            "number_of_atoms": count,
            "cell_vectors": structure.cell,
            "origin": structure.origin,
            "pbc": structure.pbc,
            "coordinate_unit": structure.metadata.get("length_unit", "input_length_unit"),
            "grain_id_source": structure.metadata.get("grain_id_source"),
            "grain_segmentation": structure.metadata.get("grain_segmentation"),
        },
        "parameters": result.parameters,
        "summary": result.summary,
        "field_interpretation": {
            "mean_curvature": "local H=(k1+k2)/2 at a surface sample; mapped to boundary atoms only",
            "mean_signed_H": "signed network area average; cancellation and grain-ID normal conventions affect it; not a minimal-surface criterion",
            "invalid_values": "NaN means unmeasured or invalid, never zero curvature; select valid=1 before interpreting values",
            "vtk_curvature_association": "CELL_DATA (triangle faces), not POINT_DATA (vertices)",
            "reconstruction_accuracy": "absolute/RMS values describe the chosen occupancy interface; validate reconstruction and smoothing against geometric controls before interpreting strain or minimality",
        },
        "files": {
            "atoms": "atoms.csv", "surface_faces": "surface_faces.csv",
            "surface_mesh": "surfaces.vtk", "summary": "summary.json",
        },
    })
    # Validate serialization before creating any output files.
    summary_text = json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    paths = {
        "atoms": output / "atoms.csv", "surface_faces": output / "surface_faces.csv",
        "surface_mesh": output / "surfaces.vtk", "summary": output / "summary.json",
    }
    with paths["atoms"].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow((
            "atom_id", "x", "y", "z", "grain_id", "grain_i", "grain_j",
            "boundary", "junction", "valid", "surface_distance", "mean_curvature",
        ))
        for index in range(count):
            writer.writerow((
                int(atom_ids[index]), *structure.positions[index], int(labels[index]),
                *pairs[index], int(boundary[index]), int(junction[index]), int(valid[index]),
                float(distance[index]), float(mean[index]),
            ))
    with paths["surface_faces"].open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow((
            "surface_index", "grain_i", "grain_j", "face_index", "x", "y", "z",
            "area", "junction", "valid", "mean_curvature", "gaussian_curvature",
        ))
        for record in records:
            for index in range(len(record["faces"])):
                writer.writerow((
                    record["index"], *record["pair"], index, *record["centers"][index],
                    record["area"][index], int(record["junction"][index]),
                    int(record["valid"][index]), record["mean"][index],
                    record["gaussian"][index],
                ))
    _write_surface_vtk(paths["surface_mesh"], records)
    paths["summary"].write_text(summary_text, encoding="utf-8")
    return paths


__all__ = ["write_grain_boundary_results"]
