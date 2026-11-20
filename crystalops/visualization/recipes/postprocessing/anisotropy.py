"""Compatibility entry point for :mod:`crystalops.analysis.anisotropy`."""

from __future__ import annotations

if __package__ in {None, ""}:
    from pathlib import Path as _Path
    import sys as _sys
    _sys.path.insert(0, str(_Path(__file__).resolve().parents[4]))

if __name__ == "__main__":
    from runpy import run_module as _run_module
    _run_module("crystalops.analysis.anisotropy", run_name="__main__")
else:
    from crystalops.analysis.anisotropy import aniso_alpha, aniso_lattice

    __all__ = ["aniso_alpha", "aniso_lattice"]
