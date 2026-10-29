"""Build the Weaire-Phelan precursor from Yang et al. (2025)."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops import FCC, CrystalStructure, WeairePhelanCrystal


def build_weaire_phelan(
    *,
    size: float = 14.0,
    lattice_constant: float = 3.61,
    twin_type: str = "sigma3",
    dod_orientation_seed: int | None = None,
    label_grains: bool = False,
) -> CrystalStructure:
    """Build the article's eight-grain periodic W-P precursor."""
    return WeairePhelanCrystal(
        unit_cell=FCC("Cu"),
        size=size,
        lattice_constant=lattice_constant,
        twin_type=twin_type,
        dod_orientation_seed=dod_orientation_seed,
        label_grains=label_grains,
    ).build()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=float, default=14.0)
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument(
        "--twin-type",
        choices=("sigma3", "sigma5"),
        default="sigma3",
    )
    parser.add_argument("--label-grains", action="store_true")
    parser.add_argument(
        "--dod-orientation-seed",
        type=int,
        default=None,
        help="optional seed for experimental random Dod0/Dod1 orientations",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="output path; defaults to weaire_phelan_<twin-type>.data",
    )
    options = parser.parse_args()

    structure = build_weaire_phelan(
        size=options.size,
        lattice_constant=options.lattice_constant,
        twin_type=options.twin_type,
        dod_orientation_seed=options.dod_orientation_seed,
        label_grains=options.label_grains,
    )
    output = options.output or Path(
        f"weaire_phelan_{options.twin_type}.data"
    )
    structure.write(output)
    print(
        f"Wrote {structure.number_of_atoms} atoms in 8 grains to "
        f"{output} ({options.twin_type})."
    )


if __name__ == "__main__":
    main()
