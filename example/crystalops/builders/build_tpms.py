"""Compatibility entry point for :mod:`example.crystalops.modeling.builders.build_tpms`."""

from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from example.crystalops.modeling.builders.build_tpms import build_tpms, main

__all__ = ["build_tpms", "main"]


if __name__ == "__main__":
    main()
