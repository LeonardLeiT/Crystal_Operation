"""Build a mathematical DSC or PSC Schwarz crystal."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys
from numpy.typing import ArrayLike

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops import CrystalStructure, SchwarzCrystal


def build_schwarz(
    *,
    size: int = 8,
    lattice_constant: float = 3.61,
    schwarz_type: str = "DSC",
    random_seed: int | None = None,
    orientation: ArrayLike | None = None,
    max_box_expansion: float = 0.25,
    max_orientation_index: int = 5,
) -> CrystalStructure:
    """Construct a DSC or PSC structure through ``SchwarzCrystal``."""
    return SchwarzCrystal(
        size=size,
        lattice_constant=lattice_constant,
        schwarz_type=schwarz_type,
        random_seed=random_seed,
        orientation=orientation,
        max_box_expansion=max_box_expansion,
        max_orientation_index=max_orientation_index,
    ).build()


def main() -> None:
    """Parse options, construct D-SCI, and write a LAMMPS data file."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=8)
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--type", choices=("DSC", "PSC"), default="PSC")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--max-box-expansion", type=float, default=0.25)
    parser.add_argument("--max-orientation-index", type=int, default=5)
    parser.add_argument("--output", type=Path, default=Path("schwarz.data"))
    options = parser.parse_args()

    structure = build_schwarz(
        size=options.size,
        lattice_constant=options.lattice_constant,
        schwarz_type=options.type,
        random_seed=options.seed,
        max_box_expansion=options.max_box_expansion,
        max_orientation_index=options.max_orientation_index,
    )
    structure.write(options.output)
    detail = structure.metadata.get(
        "domain_counts",
        structure.metadata.get("grain_counts"),
    )
    print(
        f"Wrote {structure.number_of_atoms} {options.type} atoms "
        f"({detail}) to {options.output}"
    )


if __name__ == "__main__":
    main()
