"""Render structure and TPMS figures with the optional crysplot package."""

from __future__ import annotations

from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from crysplot import StructurePlotter, TpmsPlotter
from crystalops import CrystalBuilder, FCC


def main() -> None:
    """Generate two SVG files demonstrating the plotting API."""
    structure = CrystalBuilder(FCC("Cu"), (3, 3, 3), 3.61).build()
    StructurePlotter(structure).plot(output="structure_case.svg", show=False)
    TpmsPlotter("gyroid", resolution=40).plot(
        output="tpms_case.svg",
        show=False,
    )
    print("Wrote structure_case.svg and tpms_case.svg")


if __name__ == "__main__":
    main()
