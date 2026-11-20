"""Small output-contract cases independent of the curvature implementation."""

from __future__ import annotations

import csv
import json
from types import SimpleNamespace

import numpy as np

from crystalops.core import CrystalStructure
from crystalops.io.grain_boundary_results import write_grain_boundary_results


def _example_result():
    structure = CrystalStructure(
        positions=[[0.1, 0.2, 0.3], [0.7, 0.8, 0.9]],
        cell=np.eye(3), atom_types=[1, 1], type_elements={1: "Cu"},
        properties={"atom_ids": [10, 7], "grain_ids": [3, 8]},
        metadata={"length_unit": "angstrom", "source": "example.data"},
    )
    surface = SimpleNamespace(
        vertices=np.array([[0.0, 0.0, 0.5], [1.0, 0.0, 0.5], [0.0, 1.0, 0.5]]),
        faces=np.array([[0, 1, 2]]), grain_pair=(3, 8),
        face_mean_curvature=np.array([0.0]), face_area=np.array([0.5]),
        face_valid=np.array([True]),
    )
    result = SimpleNamespace(
        surfaces=[surface], atom_mean_curvature=np.array([0.0, np.nan]),
        atom_valid=np.array([True, False]), atom_grain_pairs=np.array([[3, 8], [-1, -1]]),
        atom_boundary_mask=np.array([True, False]), atom_junction_mask=np.array([False, False]),
        atom_surface_distance=np.array([0.2, np.nan]),
        summary={"area": np.float64(0.5), "undefined": np.nan},
        parameters={"smoothing_length": np.float64(0.2)},
    )
    return structure, result


def test_curvature_export_retains_ids_masks_and_strict_json(tmp_path) -> None:
    structure, result = _example_result()
    original_positions = structure.positions.copy()
    paths = write_grain_boundary_results(result, structure, tmp_path / "export")
    assert set(paths) == {"atoms", "surface_faces", "surface_mesh", "summary"}
    with paths["atoms"].open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert [int(row["atom_id"]) for row in rows] == [10, 7]
    assert [int(row["grain_id"]) for row in rows] == [3, 8]
    assert rows[0]["valid"] == "1" and rows[1]["valid"] == "0"
    assert rows[1]["mean_curvature"].lower() == "nan"
    payload = json.loads(paths["summary"].read_text(encoding="utf-8"))
    assert payload["summary"]["undefined"] is None
    assert payload["structure"]["coordinate_unit"] == "angstrom"
    assert payload["parameters"]["smoothing_length"] == 0.2
    assert np.array_equal(structure.positions, original_positions)
    vtk = paths["surface_mesh"].read_text(encoding="utf-8")
    assert "DATASET POLYDATA\nPOINTS 3 double" in vtk
    assert "POLYGONS 1 4\n3 0 1 2" in vtk
    assert "CELL_DATA 1" in vtk
    assert "SCALARS valid unsigned_char 1\nLOOKUP_TABLE default\n1" in vtk
    with paths["surface_faces"].open(encoding="utf-8", newline="") as handle:
        face = next(csv.DictReader(handle))
    assert int(face["grain_i"]) == 3 and int(face["grain_j"]) == 8
    assert float(face["area"]) == 0.5 and float(face["mean_curvature"]) == 0.0


def test_curvature_export_rejects_inconsistent_validity_before_writing(tmp_path) -> None:
    structure, result = _example_result()
    result.atom_valid[:] = True
    destination = tmp_path / "invalid_export"
    try:
        write_grain_boundary_results(result, structure, destination)
    except ValueError as error:
        assert "finite mean curvature" in str(error)
    else:
        raise AssertionError("Non-finite valid curvature must not be exported")
    assert not destination.exists()


def test_curvature_export_offsets_disconnected_surface_indices(tmp_path) -> None:
    structure, result = _example_result()
    first = result.surfaces[0]
    result.surfaces.append(SimpleNamespace(
        **{**vars(first), "vertices": first.vertices + 0.1, "grain_pair": (8, 9)}
    ))
    paths = write_grain_boundary_results(result, structure, tmp_path / "two_surfaces")
    vtk = paths["surface_mesh"].read_text(encoding="utf-8")
    assert "POINTS 6 double" in vtk
    assert "POLYGONS 2 8\n3 0 1 2\n3 3 4 5" in vtk


def test_grain_label_csv_aligns_by_atom_id_and_rejects_missing_rows(tmp_path) -> None:
    from example.crystalops.analysis.analyze_grain_boundary_curvature import load_grain_labels_csv

    structure, _ = _example_result()
    path = tmp_path / "labels.csv"
    path.write_text("atom_id,grain_id\n7,8\n10,3\n", encoding="utf-8")
    assert np.array_equal(load_grain_labels_csv(structure, path), [3, 8])
    path.write_text("atom_id,grain_id\n7,8\n", encoding="utf-8")
    try:
        load_grain_labels_csv(structure, path)
    except ValueError as error:
        assert "1 missing" in str(error)
    else:
        raise AssertionError("Missing grain-label rows must not be silently reordered or filled")
    path.write_text("atom_id,grain_id\n7,8\n7,8\n10,3\n", encoding="utf-8")
    try:
        load_grain_labels_csv(structure, path)
    except ValueError as error:
        assert "duplicate atom_id" in str(error)
    else:
        raise AssertionError("Duplicate grain-label rows must be rejected")
