"""Read and write atomistic structures through a format registry."""

from .registry import (
    available_formats,
    read_structure,
    register_format,
    write_structure,
)
from .lammps_data import LammpsDataIO
from .grain_boundary_results import write_grain_boundary_results
from .defect_mesh import read_defect_mesh, write_defect_mesh
from .defect_mesh_curvature import write_defect_mesh_curvature
from .selection import write_selection_dump

register_format("lammps-data", LammpsDataIO, suffixes=(".data",), replace=True)

__all__ = [
    "LammpsDataIO",
    "available_formats",
    "read_structure",
    "register_format",
    "write_structure",
    "write_grain_boundary_results",
    "read_defect_mesh",
    "write_defect_mesh",
    "write_defect_mesh_curvature",
    "write_selection_dump",
]
