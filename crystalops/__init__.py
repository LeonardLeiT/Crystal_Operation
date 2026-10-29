"""Build, transform, read, and write particle-level crystal structures.

Examples:
    >>> from crystalops import CrystalBuilder, FCC
    >>> structure = CrystalBuilder(FCC("Cu"), (4, 4, 4), 3.61).build()
    >>> structure.write("copper.data")
"""

from crystalops.builders import (
    CrystalBuilder,
    KelvinCrystal,
    SchwarzCrystal,
    StructureBuilder,
    TpmsBuilder,
    VoronoiPolycrystalBuilder,
    WeairePhelanCrystal,
    available_builders,
    create_structure,
)
from crystalops.core import CrystalStructure, Element, UnitCell, get_element
from crystalops.io import available_formats, read_structure, write_structure
from crystalops.lattices import (
    ABO,
    BCC,
    BCM,
    BCO,
    BCT,
    FCC,
    FCC_111,
    FCO,
    HCP,
    RHL,
    SC,
    SM,
    SO,
    ST,
    TRI,
)

__version__ = "0.1.0"

__all__ = [
    "ABO", "BCC", "BCM", "BCO", "BCT", "FCC", "FCC_111", "FCO",
    "HCP", "RHL", "SC", "SM", "SO", "ST", "TRI", "CrystalBuilder", "CrystalStructure",
    "Element", "KelvinCrystal", "SchwarzCrystal", "StructureBuilder",
    "TpmsBuilder",
    "UnitCell", "VoronoiPolycrystalBuilder", "available_builders",
    "WeairePhelanCrystal",
    "available_formats", "create_structure", "get_element", "read_structure",
    "write_structure",
]
