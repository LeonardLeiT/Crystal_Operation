"""Compatibility entry point for :mod:`example.crystalops.modeling.builders.build_schwarz`."""

from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from example.crystalops.modeling.builders.build_schwarz import build_schwarz, main

__all__ = ["build_schwarz", "main"]


if __name__ == "__main__":
    main()
