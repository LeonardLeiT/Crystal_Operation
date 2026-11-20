"""Numerical structure and simulation-result analysis.

Analysis shares the data models in :mod:`crystalops.core` and has no plotting
dependency. Scientific backends are loaded only by the requested routines.
"""

from importlib import import_module as _import_module

_EXPORTS = {
    **dict.fromkeys(("aniso_alpha", "aniso_lattice"), ".anisotropy"),
    **dict.fromkeys((
        "mean_curvature_from_derivatives", "gaussian_curvature_from_derivatives",
        "principal_curvatures_from_derivatives",
    ), ".curvature"),
    **dict.fromkeys((
        "GrainBoundarySelection", "GrainBoundarySurface", "GrainBoundaryReconstruction",
        "GrainBoundaryCurvatureResult", "identify_grain_boundary_atoms",
        "reconstruct_grain_boundary_surfaces", "analyze_grain_boundary_curvature",
    ), ".grain_boundary"),
    **dict.fromkeys(("GrainSegmentationResult", "segment_grains_with_ptm"), ".grain_segmentation"),
    **dict.fromkeys(("DefectMeshResult", "extract_dxa_defect_mesh"), ".defect_mesh"),
    "select_dsc_control_band": ".dsc_control_band",
    **dict.fromkeys((
        "DefectMeshCurvature", "compute_defect_mesh_curvature",
    ), ".defect_mesh_curvature"),
}

__all__ = list(_EXPORTS)


def __getattr__(name: str):
    """Load numerical routines only when requested."""
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(_import_module(_EXPORTS[name], __name__), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    """Expose lazy public names to interactive completion."""
    return sorted(set(globals()) | set(__all__))
