"""Compatibility entry point for :mod:`example.crystalops.modeling.builders.build_fcc`."""

from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from example.crystalops.modeling.builders.build_fcc import build_fcc, main

__all__ = ["build_fcc", "main"]


if __name__ == "__main__":
    main()
