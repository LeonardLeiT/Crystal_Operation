"""Compatibility entry point for :mod:`crystalops.visualization.recipes.figures.average_bar_column_single`."""

from __future__ import annotations

if __package__ in {None, ""}:
    from pathlib import Path as _Path
    import sys as _sys
    _sys.path.insert(0, str(_Path(__file__).resolve().parents[3]))

if __name__ == "__main__":
    from runpy import run_module as _run_module
    _run_module("crystalops.visualization.recipes.figures.average_bar_column_single", run_name="__main__")
else:
    from crystalops.visualization.recipes.figures.average_bar_column_single import *  # noqa: F401,F403
