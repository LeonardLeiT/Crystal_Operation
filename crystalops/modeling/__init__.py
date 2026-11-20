"""Structure builders, lattice factories, and structure operations.

Implementations are imported on demand; this namespace has no plotting dependency.
"""

from importlib import import_module as _import_module

_GROUPS = {
    "crystalops.modeling.builders": (
        "CrystalBuilder", "KelvinCrystal", "SchwarzCrystal", "StructureBuilder",
        "TpmsBuilder", "VoronoiPolycrystalBuilder", "WeairePhelanCrystal",
        "available_builders", "create_structure", "register_builder",
    ),
    "crystalops.modeling.operations": (
        "assign_atom_type", "merge", "remove_duplicates", "repeat", "rotate",
        "rotation_matrix", "select", "select_atoms", "select_sphere", "slice_by_planes",
        "translate", "wrap",
    ),
    "crystalops.modeling.lattices": (
        "BCC", "BCM", "BCO", "BCT", "FCC", "FCC_111", "FCO", "HCP", "RHL",
        "SC", "SM", "SO", "ST", "TRI", "ABO",
        "base_centered_monoclinic", "base_centered_orthorhombic",
        "body_centered_cubic", "body_centered_orthorhombic", "body_centered_tetragonal",
        "face_centered_cubic", "face_centered_cubic_111", "face_centered_orthorhombic",
        "hexagonal_close_packed", "rhombohedral", "simple_cubic", "simple_monoclinic",
        "simple_orthorhombic", "simple_tetragonal", "triclinic",
    ),
}
_EXPORTS = {name: module for module, names in _GROUPS.items() for name in names}
__all__ = list(_EXPORTS)


def __getattr__(name: str):
    """Load a public modeling object without importing unrelated modules."""
    try:
        module_name = _EXPORTS[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
    value = getattr(_import_module(module_name), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    """Expose lazy names to interactive completion."""
    return sorted(set(globals()) | set(__all__))
