"""Render a generated crystal with the shared visualization module."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops.visualization import StructurePlotter
from crystalops.modeling import CrystalBuilder, FCC


def main() -> None:
    """Build an FCC crystal and render it as an image."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("fcc.svg"))
    parser.add_argument("--show", action="store_true")
    options = parser.parse_args()

    structure = CrystalBuilder(FCC("Cu"), (3, 3, 3), 3.61).build()
    StructurePlotter(structure).plot(
        output=options.output,
        show=options.show,
    )
    print(f"Wrote {options.output}")


if __name__ == "__main__":
    main()
