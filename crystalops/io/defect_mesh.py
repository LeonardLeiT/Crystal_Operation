"""Portable DXA mesh archives and cap-free VTK previews."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from crystalops.analysis.defect_mesh import DefectMeshResult

__all__ = ["write_defect_mesh", "read_defect_mesh"]

_FORMAT = "crystalops.dxa-defect-mesh"


def write_defect_mesh(
    result: DefectMeshResult,
    output_dir: str | Path,
    *,
    prefix: str = "defect_mesh",
) -> dict[str, Path]:
    """Save native periodic geometry, a clipped VTK mesh, and JSON metadata.

    The NPZ archive is self-contained and can be read without OVITO. The VTK
    file contains the clipped surface only, without cell cap polygons, and a
    ``source_face`` index for each triangle. Its artificial domain cuts are not
    physical holes. No curvature or grain-pair interpretation is added.
    """
    result.validate()
    if not isinstance(prefix, str) or not prefix or Path(prefix).name != prefix or prefix in {".", ".."}:
        raise ValueError("prefix must be a nonempty filename stem without directories")
    destination = Path(output_dir)
    paths = {key: destination / f"{prefix}.{suffix}" for key, suffix in (
        ("mesh", "npz"), ("vtk", "vtk"), ("json", "json"),
    )}
    metadata = {
        "format": _FORMAT, "version": 1,
        "parameters": result.parameters, "diagnostics": result.diagnostics,
        "cell": result.cell.tolist(), "origin": result.origin.tolist(), "pbc": list(result.pbc),
        "array_conventions": {
            "cell": "row vectors in Cartesian coordinates",
            "vertices_faces": "native periodic mesh; use together with cell, origin and pbc",
            "clipped_vertices_faces": "Cartesian triangles cut at periodic cell boundaries, without caps",
            "clipped_face_sources": "native face index per clipped triangle",
        },
        "files": {key: path.name for key, path in paths.items()},
    }
    serialized = json.dumps(metadata, ensure_ascii=False, indent=2, allow_nan=False)
    destination.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        paths["mesh"], vertices=result.vertices, faces=result.faces,
        clipped_vertices=result.clipped_vertices, clipped_faces=result.clipped_faces,
        clipped_face_sources=result.clipped_face_sources,
        cell=result.cell, origin=result.origin, pbc=np.asarray(result.pbc, dtype=bool),
        metadata_json=np.asarray(serialized),
    )
    with paths["vtk"].open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("# vtk DataFile Version 3.0\nDXA defect mesh - clipped, no cell caps\nASCII\nDATASET POLYDATA\n")
        handle.write(f"POINTS {len(result.clipped_vertices)} double\n")
        np.savetxt(handle, result.clipped_vertices, fmt="%.17g")
        handle.write(f"POLYGONS {len(result.clipped_faces)} {4 * len(result.clipped_faces)}\n")
        if len(result.clipped_faces):
            np.savetxt(handle, np.column_stack((np.full(len(result.clipped_faces), 3), result.clipped_faces)), fmt="%d")
            handle.write(f"CELL_DATA {len(result.clipped_faces)}\nSCALARS source_face long 1\nLOOKUP_TABLE default\n")
            np.savetxt(handle, result.clipped_face_sources, fmt="%d")
    paths["json"].write_text(serialized + "\n", encoding="utf-8")
    return paths


def read_defect_mesh(filename: str | Path) -> DefectMeshResult:
    """Load a saved NPZ mesh without pickle, OVITO or the original data file."""
    from crystalops.analysis.defect_mesh import DefectMeshResult

    required = {
        "vertices", "faces", "clipped_vertices", "clipped_faces", "clipped_face_sources",
        "cell", "origin", "pbc", "metadata_json",
    }
    with np.load(filename, allow_pickle=False) as archive:
        if not required.issubset(archive.files):
            raise ValueError(f"defect mesh archive is missing arrays: {sorted(required - set(archive.files))}")
        metadata = json.loads(str(archive["metadata_json"].item()))
        if metadata.get("format") != _FORMAT or metadata.get("version") != 1:
            raise ValueError("unsupported defect mesh archive format or version")
        return DefectMeshResult(
            **{name: archive[name] for name in required - {"metadata_json"}},
            parameters=metadata["parameters"], diagnostics=metadata["diagnostics"],
        )
