"""DXA integration checks for empty and boundary-crossing defect meshes."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np

from crystalops.analysis import extract_dxa_defect_mesh
from crystalops.core import CrystalStructure
from crystalops.io import read_defect_mesh, write_defect_mesh
from crystalops.modeling import CrystalBuilder, FCC


def _require_ovito():
    if importlib.util.find_spec("ovito") is None:
        import pytest
        pytest.skip("optional OVITO backend not installed")


def _crystal_with_periodic_void():
    base = CrystalBuilder(FCC("Cu"), (6, 7, 8), 3.61).build()
    origin = np.array([2.0, -3.0, 7.0])
    # A void crossing the x seam: exporting native triangles naively would
    # create long spurious sheets across the cell.
    center = np.array([0.0, 3 * 3.61, 4 * 3.61])
    delta = base.positions - center
    lengths = np.diag(base.cell)
    delta -= np.rint(delta / lengths) * lengths
    keep = np.linalg.norm(delta, axis=1) > 3.2
    return CrystalStructure(
        positions=base.positions[keep] + origin, cell=base.cell, origin=origin,
        atom_types=base.atom_types[keep], type_elements=base.type_elements,
        metadata={"length_unit": "angstrom"},
    )


def test_perfect_crystal_empty_mesh_can_be_saved(tmp_path):
    _require_ovito()
    structure = CrystalBuilder(FCC("Cu"), (6, 7, 8), 3.61).build()
    result = extract_dxa_defect_mesh(structure)
    assert result.vertices.shape == result.faces.shape == (0, 3)
    assert result.clipped_vertices.shape == result.clipped_faces.shape == (0, 3)
    assert result.diagnostics["empty_mesh"]
    paths = write_defect_mesh(result, tmp_path)
    loaded = read_defect_mesh(paths["mesh"])
    assert loaded.faces.shape == (0, 3)
    assert "POINTS 0 double\nPOLYGONS 0 0" in paths["vtk"].read_text()


def test_periodic_void_clipping_and_portable_roundtrip(tmp_path):
    _require_ovito()
    structure = _crystal_with_periodic_void()
    before = structure.positions.copy()
    result = extract_dxa_defect_mesh(structure)
    assert len(result.faces) > 0
    assert result.diagnostics["excluded_cap_face_count"] > 0
    assert len(result.clipped_faces) > len(result.faces)
    fractional = np.linalg.solve(result.cell.T, (result.clipped_vertices - result.origin).T).T
    assert np.all(fractional >= -1e-12) and np.all(fractional <= 1 + 1e-12)
    triangles = result.clipped_vertices[result.clipped_faces]
    edges = np.linalg.norm(triangles - np.roll(triangles, 1, axis=1), axis=2)
    assert edges.max() < np.linalg.norm(result.cell, axis=1).min() / 2
    assert np.array_equal(structure.positions, before)
    assert np.array_equal(result.cell, structure.cell)
    assert np.array_equal(result.origin, structure.origin)
    assert result.pbc == structure.pbc

    paths = write_defect_mesh(result, tmp_path)
    loaded = read_defect_mesh(paths["mesh"])
    for name in ("vertices", "faces", "clipped_vertices", "clipped_faces", "clipped_face_sources", "cell", "origin"):
        assert np.array_equal(getattr(result, name), getattr(loaded, name)), name
    assert loaded.parameters == result.parameters
    assert loaded.diagnostics == result.diagnostics
    assert loaded.pbc == result.pbc
    payload = json.loads(paths["json"].read_text(encoding="utf-8"))
    assert payload["parameters"]["length_unit"] == "angstrom"
    assert payload["parameters"]["surface_smoothing"] == 0
    vtk = paths["vtk"].read_text(encoding="utf-8")
    assert f"POLYGONS {len(result.clipped_faces)} {4 * len(result.clipped_faces)}\n" in vtk
    assert "SCALARS source_face long 1" in vtk
    # VTK must reload with exactly the cap-free triangle count.
    from ovito.io import import_file
    imported = import_file(str(paths["vtk"])).compute()
    mesh = next(obj for obj in imported.objects if hasattr(obj, "face_count"))
    assert mesh.face_count == len(result.clipped_faces)


def test_extraction_matches_native_dxa_and_keeps_smoothing_explicit():
    _require_ovito()
    from ovito.data import DataCollection
    from ovito.modifiers import DislocationAnalysisModifier
    from ovito.pipeline import Pipeline, StaticSource

    structure = _crystal_with_periodic_void()
    result = extract_dxa_defect_mesh(structure, surface_smoothing=2)
    data = DataCollection()
    particles = data.create_particles(count=structure.number_of_atoms)
    particles.create_property("Position", data=structure.positions)
    data.create_cell(np.column_stack((structure.cell.T, structure.origin)), pbc=structure.pbc)
    pipeline = Pipeline(source=StaticSource(data=data))
    pipeline.modifiers.append(DislocationAnalysisModifier(defect_mesh_smoothing_level=2))
    native = pipeline.compute().surfaces["dxa-defect-mesh"]
    assert np.array_equal(result.vertices, native.get_vertices())
    assert np.array_equal(result.faces, native.get_faces())
    assert result.parameters["surface_smoothing"] == 2


def test_invalid_mesh_settings_fail_before_dxa():
    structure = _crystal_with_periodic_void()
    for options in ({"crystal_structure": "diamond_surface"}, {"surface_smoothing": -1},
                    {"surface_smoothing": 1.5}, {"surface_smoothing": True}):
        try:
            extract_dxa_defect_mesh(structure, **options)
        except ValueError:
            pass
        else:
            raise AssertionError(f"invalid settings were accepted: {options}")


def test_defect_mesh_imports_do_not_load_optional_backends():
    code = """
import sys
from crystalops.analysis import DefectMeshResult, extract_dxa_defect_mesh
from crystalops.io import read_defect_mesh, write_defect_mesh
from crystalops.visualization import plot_defect_mesh
for name in ('ovito', 'matplotlib', 'plotly', 'scipy', 'skimage'):
    assert name not in sys.modules, name
"""
    subprocess.run([sys.executable, "-c", code], check=True, capture_output=True,
                   text=True, cwd=Path(__file__).resolve().parents[3])
