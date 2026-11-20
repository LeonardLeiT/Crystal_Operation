"""Build an FCC supercell and write a LAMMPS data file."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

from crystalops.core import CrystalStructure
from crystalops.modeling import CrystalBuilder, FCC


def build_fcc(
    *,
    element: str = "Cu",
    repeats: tuple[int, int, int] = (10, 10, 10),
    lattice_constant: float = 3.61,
) -> CrystalStructure:
    """Build an FCC supercell.

    Args:
        element: Element occupying all FCC lattice sites.
        repeats: Unit-cell repetitions along the three cell vectors.
        lattice_constant: Cubic lattice constant in angstroms.

    Returns:
        The generated periodic crystal structure.
    """
    return CrystalBuilder(
        unit_cell=FCC(element),
        repeats=repeats,
        lattice_constant=lattice_constant,
    ).build()


def main() -> None:
    """Parse command-line options, build the crystal, and write it."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--element", default="Cu")
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--repeats", type=int, nargs=3, default=(10, 10, 10))
    parser.add_argument("--output", type=Path, default=Path("fcc.data"))
    options = parser.parse_args()

    structure = build_fcc(
        element=options.element,
        repeats=tuple(options.repeats),
        lattice_constant=options.lattice_constant,
    )
    structure.write(options.output)
    print(f"Wrote {structure.number_of_atoms} atoms to {options.output}")


if __name__ == "__main__":
    main()
