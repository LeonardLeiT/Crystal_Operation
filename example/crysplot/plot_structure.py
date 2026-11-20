"""Compatibility entry point for :mod:`example.crystalops.visualization.plot_structure`."""

from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from example.crystalops.visualization.plot_structure import main

__all__ = ["main"]


if __name__ == "__main__":
    main()
