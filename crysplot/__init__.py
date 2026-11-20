"""Compatibility imports for :mod:`crystalops.visualization`.

New code should import plotters from ``crystalops.visualization``. Legacy
imports return the same plotter classes and retain lazy optional dependencies.
"""

from __future__ import annotations

from importlib import import_module


__all__ = ["StructurePlotter", "TpmsPlotter", "VoronoiPlotter"]


def __getattr__(name: str):
    """Forward public plotters to the canonical visualization package."""
    if name not in __all__:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module("crystalops.visualization"), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    """Include lazy compatibility exports in interactive completion."""
    return sorted(set(globals()) | set(__all__))
