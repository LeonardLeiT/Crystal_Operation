"""Shared tools for atomistic modeling, analysis, visualization, and I/O.

The root namespace retains its original convenience exports. Implementations
are loaded on demand, so importing the data model does not load builders,
analysis routines, file-format registries, or optional visualization libraries.
"""

from importlib import import_module as _import_module

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

_EXPORTS = {
    **dict.fromkeys(
        ("CrystalStructure", "Element", "UnitCell", "get_element"),
        "crystalops.core",
    ),
    **dict.fromkeys(
        (
            "CrystalBuilder", "KelvinCrystal", "SchwarzCrystal",
            "StructureBuilder", "TpmsBuilder", "VoronoiPolycrystalBuilder",
            "WeairePhelanCrystal", "available_builders", "create_structure",
        ),
        "crystalops.modeling.builders",
    ),
    **dict.fromkeys(
        (
            "ABO", "BCC", "BCM", "BCO", "BCT", "FCC", "FCC_111", "FCO",
            "HCP", "RHL", "SC", "SM", "SO", "ST", "TRI",
        ),
        "crystalops.modeling.lattices",
    ),
    **dict.fromkeys(
        ("available_formats", "read_structure", "write_structure"),
        "crystalops.io",
    ),
}

_SUBMODULES = {
    "core", "modeling", "analysis", "visualization", "io",
    "builders", "operations", "lattices",
}


def __getattr__(name: str):
    if name in _SUBMODULES:
        value = _import_module(f"{__name__}.{name}")
    elif name in _EXPORTS:
        value = getattr(_import_module(_EXPORTS[name]), name)
    else:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__) | _SUBMODULES)
