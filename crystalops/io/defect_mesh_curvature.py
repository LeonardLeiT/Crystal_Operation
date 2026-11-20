"""Save native periodic curvature arrays and a clipped mesh for visualization."""

from __future__ import annotations

import json
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING

import numpy as np

from .grain_boundary_results import _json_value

if TYPE_CHECKING:
    from crystalops.analysis.defect_mesh import DefectMeshResult
    from crystalops.analysis.defect_mesh_curvature import DefectMeshCurvature

__all__ = ["write_defect_mesh_curvature"]

_FORMAT = "crystalops.dxa-defect-mesh-curvature"


def _write_vtk(
    path: Path,
    mesh: DefectMeshResult,
    clipped_mean_curvature: np.ndarray,
    clipped_valid: np.ndarray,
) -> None:
    # Binary legacy VTK preserves IEEE NaN values in readers that reject them
    # in ASCII. The legacy polygon connectivity and cell indices are int32.
    face_count = len(mesh.clipped_faces)
    limit = np.iinfo(np.int32).max
    if max(len(mesh.clipped_vertices), len(mesh.faces), 4 * face_count) > limit:
        raise ValueError("mesh exceeds legacy VTK's 32-bit indexing limits")
    with path.open("wb") as file:
        file.write(
            b"# vtk DataFile Version 3.0\n"
            b"DXA defect mesh curvature display - clipped, no cell caps\n"
            b"BINARY\nDATASET POLYDATA\n"
        )
        file.write(f"POINTS {len(mesh.clipped_vertices)} double\n".encode("ascii"))
        file.write(np.asarray(mesh.clipped_vertices, dtype=">f8").tobytes())
        file.write(f"\nPOLYGONS {face_count} {4 * face_count}\n".encode("ascii"))
        polygons = np.empty((face_count, 4), dtype=">i4")
        polygons[:, 0] = 3
        polygons[:, 1:] = mesh.clipped_faces
        file.write(polygons.tobytes())
        file.write(b"\n")
        if face_count:
            file.write(f"CELL_DATA {face_count}\n".encode("ascii"))
            for name, vtk_type, dtype, values in (
                ("mean_curvature", "double", ">f8", clipped_mean_curvature),
                ("source_face", "int", ">i4", mesh.clipped_face_sources),
                ("valid", "unsigned_char", "u1", clipped_valid),
            ):
                file.write(f"SCALARS {name} {vtk_type} 1\nLOOKUP_TABLE default\n".encode("ascii"))
                file.write(np.asarray(values, dtype=dtype).tobytes())
                file.write(b"\n")


def write_defect_mesh_curvature(
    mesh: DefectMeshResult,
    result: DefectMeshCurvature,
    output_dir: str | Path,
) -> dict[str, Path]:
    """Save a self-contained NPZ archive, clipped preview VTK, and summary JSON.

    ``result`` must be calculated from ``mesh``. The NPZ retains curvature,
    weights, and masks in the native mesh's vertex and face order. For display,
    VTK replicates each native face's mean of its three vertex values onto
    all clipped triangles from that face. The value is NaN when the native
    face or a contributing vertex is invalid. Use native arrays and area
    weights for statistics, since clipping may split one face into many.
    """
    mesh.validate()
    arrays = {}
    for name, count, dtype in (
        ("vertex_mean_curvature", len(mesh.vertices), np.float64),
        ("vertex_area", len(mesh.vertices), np.float64),
        ("vertex_valid", len(mesh.vertices), np.bool_),
        ("face_area", len(mesh.faces), np.float64),
        ("face_valid", len(mesh.faces), np.bool_),
    ):
        array = np.asarray(getattr(result, name), dtype=dtype)
        if array.shape != (count,):
            raise ValueError(f"{name} must have one value per corresponding native mesh element")
        arrays[name] = array
    vertex_curvature = arrays["vertex_mean_curvature"]
    vertex_valid = arrays["vertex_valid"]
    if not np.isfinite(vertex_curvature[vertex_valid]).all():
        raise ValueError("valid native vertices must have finite mean curvature")
    native_display_valid = arrays["face_valid"] & vertex_valid[mesh.faces].all(axis=1)
    source = mesh.clipped_face_sources
    clipped_valid = native_display_valid[source]
    clipped_curvature = np.full(len(source), np.nan, dtype=np.float64)
    clipped_curvature[clipped_valid] = vertex_curvature[mesh.faces[source[clipped_valid]]].mean(axis=1)

    destination = Path(output_dir)
    paths = {
        "mesh": destination / "defect_mesh_curvature.npz",
        "vtk": destination / "defect_mesh_curvature.vtk",
        "json": destination / "defect_mesh_curvature.json",
    }
    length_unit = str(mesh.parameters.get("length_unit", "input_length_unit"))
    metadata = _json_value({
        "format": _FORMAT,
        "version": 1,
        "summary": result.summary,
        "parameters": result.parameters,
        "mesh_parameters": mesh.parameters,
        "mesh_diagnostics": mesh.diagnostics,
        "length_unit": length_unit,
        "curvature_unit": f"1/{length_unit}",
        "area_unit": f"{length_unit}^2",
        "cell": mesh.cell,
        "origin": mesh.origin,
        "pbc": mesh.pbc,
        "array_conventions": {
            "cell": "Cartesian row vectors",
            "vertices_faces": "native mesh order; interpret positions with cell, origin, and pbc",
            "vertex_area": "one third of the area of each incident usable triangle, summed per native vertex",
            "vertex_valid": "native vertex curvature validity mask",
            "face_area": "native face area",
            "face_valid": "native face validity mask from the curvature analysis",
            "clipped_vertices_faces": "Cartesian display triangles cut at periodic boundaries, without caps",
            "clipped_face_sources": "zero-based native face index per clipped triangle",
            "vtk": "display only; each clipped triangle inherits mean curvature from its native face's three vertices, or NaN if invalid",
        },
        "files": {key: path.name for key, path in paths.items()},
    })
    serialized = json.dumps(metadata, ensure_ascii=False, indent=2, allow_nan=False)
    destination.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".defect_mesh_curvature-", dir=destination) as temporary:
        staged = {key: Path(temporary) / path.name for key, path in paths.items()}
        np.savez_compressed(
            staged["mesh"], vertices=mesh.vertices, faces=mesh.faces,
            clipped_vertices=mesh.clipped_vertices, clipped_faces=mesh.clipped_faces,
            clipped_face_sources=source, cell=mesh.cell, origin=mesh.origin,
            pbc=np.asarray(mesh.pbc, dtype=bool), metadata_json=np.asarray(serialized),
            **arrays,
        )
        _write_vtk(staged["vtk"], mesh, clipped_curvature, clipped_valid)
        staged["json"].write_text(serialized + "\n", encoding="utf-8")
        for key, path in paths.items():
            staged[key].replace(path)
    return paths
