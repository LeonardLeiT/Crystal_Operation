"""Visualize structures and numerical results with optional plotting tools.

Plotters are imported only when requested. Importing this namespace does not
load Matplotlib or scikit-image, or import project-specific plotting recipes.
"""

from __future__ import annotations

from importlib import import_module


__all__ = [
    "StructurePlotter", "TpmsPlotter", "VoronoiPlotter",
    "plot_grain_boundary_comparison", "plot_defect_mesh",
]

_PLOTTERS = {
    "StructurePlotter": ".structure",
    "TpmsPlotter": ".tpms",
    "VoronoiPlotter": ".voronoi",
    "plot_grain_boundary_comparison": ".grain_boundary",
    "plot_defect_mesh": ".defect_mesh",
}


def __getattr__(name: str):
    """Load one plotter without importing unrelated optional backends."""
    try:
        module_name = _PLOTTERS[name]
    except KeyError:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from None
    value = getattr(import_module(module_name, __name__), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    """Include lazy public names in interactive completion."""
    return sorted(set(globals()) | set(__all__))
