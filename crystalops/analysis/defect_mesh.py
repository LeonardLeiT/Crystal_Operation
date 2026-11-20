"""Extract OVITO DXA defect meshes without fitting an interface or curvature.

The native mesh is periodic. A separate clipped triangle mesh is provided for
viewers that do not support periodic topology; artificial cell cap polygons
are deliberately omitted from that representation.
"""

from __future__ import annotations

from dataclasses import dataclass
from operator import index
from typing import Any

import numpy as np
from numpy.typing import NDArray

from crystalops.core import CrystalStructure

__all__ = ["DefectMeshResult", "extract_dxa_defect_mesh"]


def _check_mesh(vertices, faces, name):
    if vertices.ndim != 2 or vertices.shape[1:] != (3,) or not np.isfinite(vertices).all():
        raise ValueError(f"{name} vertices must be finite Cartesian coordinates with shape (N, 3)")
    if faces.ndim != 2 or faces.shape[1:] != (3,) or not np.issubdtype(faces.dtype, np.integer):
        raise ValueError(f"{name} faces must contain integer triangle indices with shape (M, 3)")
    if faces.size and (faces.min() < 0 or faces.max() >= len(vertices)):
        raise ValueError(f"{name} face index is out of bounds")


@dataclass(slots=True)
class DefectMeshResult:
    """Owned NumPy arrays for a DXA defect mesh, requiring no OVITO to reload.

    ``vertices``/``faces`` are the original periodic representation and must
    be interpreted together with the row-vector ``cell``, ``origin`` and
    ``pbc``. Do not render these triangles as non-periodic Cartesian faces.
    ``clipped_vertices``/``clipped_faces`` are cut at the periodic domain
    boundaries, without artificial caps. ``clipped_face_sources`` maps each
    clipped triangle back to its native face index. Clipping can split a face.

    This object represents DXA's non-dislocation defect envelope, not a
    grain-boundary middle surface. An empty mesh is a valid result.
    """

    vertices: NDArray[np.float64]
    faces: NDArray[np.int64]
    clipped_vertices: NDArray[np.float64]
    clipped_faces: NDArray[np.int64]
    clipped_face_sources: NDArray[np.int64]
    cell: NDArray[np.float64]
    origin: NDArray[np.float64]
    pbc: tuple[bool, bool, bool]
    parameters: dict[str, Any]
    diagnostics: dict[str, Any]

    def __post_init__(self):
        for name in ("vertices", "clipped_vertices", "cell", "origin"):
            setattr(self, name, np.asarray(getattr(self, name), dtype=np.float64).copy())
        # Check indices before conversion, so fractional values cannot be truncated.
        self.faces = np.asarray(self.faces).copy()
        self.clipped_faces = np.asarray(self.clipped_faces).copy()
        self.clipped_face_sources = np.asarray(self.clipped_face_sources).copy()
        self.pbc = tuple(bool(value) for value in self.pbc)
        self.parameters = dict(self.parameters)
        self.diagnostics = dict(self.diagnostics)
        self.validate()
        self.faces = self.faces.astype(np.int64)
        self.clipped_faces = self.clipped_faces.astype(np.int64)
        self.clipped_face_sources = self.clipped_face_sources.astype(np.int64)

    def validate(self):
        """Check geometry and source-face correspondence, including empty meshes."""
        _check_mesh(self.vertices, self.faces, "native")
        _check_mesh(self.clipped_vertices, self.clipped_faces, "clipped")
        sources = self.clipped_face_sources
        if sources.shape != (len(self.clipped_faces),) or not np.issubdtype(sources.dtype, np.integer):
            raise ValueError("clipped_face_sources must contain one integer per clipped triangle")
        if sources.size and (sources.min() < 0 or sources.max() >= len(self.faces)):
            raise ValueError("clipped source-face index is out of bounds")
        if self.cell.shape != (3, 3) or not np.isfinite(self.cell).all() or abs(np.linalg.det(self.cell)) <= np.finfo(float).eps:
            raise ValueError("cell must contain three finite, independent row vectors")
        if self.origin.shape != (3,) or not np.isfinite(self.origin).all():
            raise ValueError("origin must have shape (3,) and be finite")
        if len(self.pbc) != 3:
            raise ValueError("pbc must contain three flags")


def _integer(value, name, minimum):
    try:
        if isinstance(value, (bool, np.bool_)):
            raise TypeError
        value = index(value)
    except TypeError as exc:
        raise ValueError(f"{name} must be an integer >= {minimum}") from exc
    if value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def extract_dxa_defect_mesh(
    structure: CrystalStructure,
    *,
    crystal_structure: str = "fcc",
    surface_smoothing: int = 0,
    trial_circuit_length: int = 14,
    circuit_stretchability: int = 9,
) -> DefectMeshResult:
    """Run DXA on a complete structure and copy its defect mesh into arrays.

    ``crystal_structure`` describes the atomic lattice (``fcc``, ``hcp``,
    ``bcc``, ``cubic_diamond`` or ``hexagonal_diamond``), not the shape of the
    interface. No grain IDs are required. Coordinates, cell and PBC are passed
    to OVITO without box normalization or changes to the input structure.

    Surface smoothing defaults to zero to preserve the unfaired DXA mesh;
    set it explicitly to apply OVITO's surface fairing iterations (its native
    default is eight). Dislocation lines are used internally but not returned.
    No curvature, middle-surface extraction or feature fitting is performed.
    """
    lattices = {
        "fcc": "FCC", "hcp": "HCP", "bcc": "BCC",
        "cubic_diamond": "CubicDiamond", "hexagonal_diamond": "HexagonalDiamond",
    }
    if not isinstance(crystal_structure, str) or crystal_structure.lower() not in lattices:
        raise ValueError(f"crystal_structure must be one of {tuple(lattices)}")
    lattice = crystal_structure.lower()
    smoothing = _integer(surface_smoothing, "surface_smoothing", 0)
    circuit = _integer(trial_circuit_length, "trial_circuit_length", 3)
    stretch = _integer(circuit_stretchability, "circuit_stretchability", 0)
    structure.validate()
    if structure.number_of_atoms == 0:
        raise ValueError("DXA requires a nonempty complete structure")
    if not all(np.isfinite(a).all() for a in (structure.positions, structure.cell, structure.origin)):
        raise ValueError("DXA positions, cell and origin must be finite")
    if len(structure.pbc) != 3:
        raise ValueError("pbc must contain three flags")
    try:
        import ovito
        from ovito.data import DataCollection
        from ovito.modifiers import DislocationAnalysisModifier
        from ovito.pipeline import Pipeline, StaticSource
    except ImportError as exc:
        raise ImportError("DXA mesh extraction requires the optional 'ovito' package; install crystalops[defect-mesh]") from exc

    data = DataCollection()
    particles = data.create_particles(count=structure.number_of_atoms)
    particles.create_property("Position", data=structure.positions)
    particles.create_property("Particle Type", data=structure.atom_types)
    data.create_cell(
        np.column_stack((structure.cell.T, structure.origin)), pbc=tuple(structure.pbc),
    )
    modifier = DislocationAnalysisModifier(
        input_crystal_structure=getattr(DislocationAnalysisModifier.Lattice, lattices[lattice]),
        defect_mesh_smoothing_level=smoothing,
        trial_circuit_length=circuit,
        circuit_stretchability=stretch,
        only_selected=False,
        only_perfect_dislocations=False,
    )
    pipeline = Pipeline(source=StaticSource(data=data))
    pipeline.modifiers.append(modifier)
    output = pipeline.compute()
    key = next((key for key in ("dxa-defect-mesh", "defect-mesh") if key in output.surfaces), None)
    if key is None:
        raise RuntimeError(f"DXA produced no defect mesh object; available meshes: {list(output.surfaces)}")
    mesh = output.surfaces[key]
    triangles, caps, face_map = mesh.to_triangle_mesh()
    # OVITO may return a flat, untyped empty array for a defect-free crystal.
    vertices = np.asarray(mesh.get_vertices(), dtype=float).reshape(-1, 3).copy()
    faces = np.asarray(mesh.get_faces(), dtype=np.int64).reshape(-1, 3).copy()
    clipped_vertices = np.asarray(triangles.get_vertices(), dtype=float).reshape(-1, 3).copy()
    clipped_faces = np.asarray(triangles.get_faces(), dtype=np.int64).reshape(-1, 3).copy()
    parameters = {
        "method": "OVITO DXA defect mesh",
        "ovito_version": ovito.version_string,
        "crystal_structure": lattice,
        "surface_smoothing": smoothing,
        "trial_circuit_length": circuit,
        "circuit_stretchability": stretch,
        "only_selected": False,
        "only_perfect_dislocations": False,
        "mesh_identifier": key,
        "length_unit": str(structure.metadata.get("length_unit", "input_length_unit")),
        "source": str(structure.metadata.get("source", structure.name)),
        "interpretation": "non-dislocation defect envelope, not a grain-boundary middle surface",
        "clipped_mesh": "periodic domain clipping; artificial cap polygons excluded",
    }
    diagnostics = {
        "input_atom_count": structure.number_of_atoms,
        "native_vertex_count": len(vertices), "native_face_count": len(faces),
        "clipped_vertex_count": len(clipped_vertices), "clipped_face_count": len(clipped_faces),
        "excluded_cap_face_count": 0 if caps is None else int(caps.face_count),
        "empty_mesh": len(faces) == 0,
        "structure_type_counts": {
            str(key).removeprefix("DislocationAnalysis.counts."): int(value)
            for key, value in output.attributes.items() if str(key).startswith("DislocationAnalysis.counts.")
        },
    }
    return DefectMeshResult(
        vertices=vertices, faces=faces,
        clipped_vertices=clipped_vertices, clipped_faces=clipped_faces,
        clipped_face_sources=np.asarray(face_map),
        cell=structure.cell, origin=structure.origin, pbc=structure.pbc,
        parameters=parameters, diagnostics=diagnostics,
    )
