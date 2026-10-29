"""Read and write atomistic structures through a format registry."""

from .registry import (
    available_formats,
    read_structure,
    register_format,
    write_structure,
)
from .lammps_data import LammpsDataIO

register_format("lammps-data", LammpsDataIO, suffixes=(".data",), replace=True)

__all__ = [
    "LammpsDataIO",
    "available_formats",
    "read_structure",
    "register_format",
    "write_structure",
]
