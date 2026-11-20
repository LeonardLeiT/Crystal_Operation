"""Render structure and TPMS figures with the optional visualization module."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops.visualization import StructurePlotter, TpmsPlotter
from crystalops.modeling import CrystalBuilder, FCC


def main() -> None:
    """Generate two SVG files demonstrating the plotting API."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    options = parser.parse_args()
    options.output_dir.mkdir(parents=True, exist_ok=True)

    structure_output = options.output_dir / "structure_case.svg"
    tpms_output = options.output_dir / "tpms_case.svg"
    structure = CrystalBuilder(FCC("Cu"), (3, 3, 3), 3.61).build()
    StructurePlotter(structure).plot(output=structure_output, show=False)
    TpmsPlotter("gyroid", resolution=40).plot(
        output=tpms_output,
        show=False,
    )
    print(f"Wrote {structure_output} and {tpms_output}")


if __name__ == "__main__":
    main()
