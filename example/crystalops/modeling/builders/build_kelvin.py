"""Build a reproducible Kelvin polycrystal."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

from crystalops.core import CrystalStructure
from crystalops.modeling import FCC, KelvinCrystal


def build_kelvin(
    *,
    element: str = "Cu",
    size: float = 15,
    lattice_constant: float = 3.61,
    preset: str = "full",
    random_seed: int = 42,
    label_grains: bool = False,
) -> CrystalStructure:
    """Build a seeded Kelvin polycrystal.

    Args:
        element: Element occupying the underlying FCC crystal.
        size: Half-length of the two-cell BCC grain arrangement in lattice units.
        lattice_constant: FCC lattice constant in angstroms.
        preset: Named fixed-orientation arrangement.
        random_seed: Seed controlling all remaining grain orientations.
        label_grains: Encode grain identifiers as LAMMPS atom types for
            direct coloring in visualization software.

    Returns:
        The generated 16-grain periodic polycrystal.
    """
    return KelvinCrystal(
        unit_cell=FCC(element),
        size=size,
        lattice_constant=lattice_constant,
        preset=preset,
        random_seed=random_seed,
        label_grains=label_grains,
    ).build()


def main() -> None:
    """Parse command-line options, build the polycrystal, and write it."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--element", default="Cu")
    parser.add_argument("--size", type=float, default=10)
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--preset", default="full")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--label-grains",
        action="store_true",
        help="encode the 16 grain IDs as LAMMPS atom types",
    )
    parser.add_argument("--output", type=Path, default=Path("kelvin.data"))
    options = parser.parse_args()

    structure = build_kelvin(
        element=options.element,
        size=options.size,
        lattice_constant=options.lattice_constant,
        preset=options.preset,
        random_seed=options.seed,
        label_grains=options.label_grains,
    )
    structure.write(options.output)
    fixed = structure.metadata["fixed_positions"]
    random = structure.metadata["random_positions"]
    outside = structure.metadata["outside_box_atom_count"]
    print(
        f"Wrote {structure.number_of_atoms} atoms in 16 grains to {options.output} "
        f"({len(fixed)} fixed orientations, {len(random)} random orientations, "
        f"{outside} atoms retained outside the periodic box)"
    )


if __name__ == "__main__":
    main()
