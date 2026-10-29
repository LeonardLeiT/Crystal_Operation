"""Build a TPMS-selected FCC crystal."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops import CrystalStructure, FCC, TpmsBuilder


def build_tpms(
    *,
    surface: str = "gyroid",
    mode: str = "negative",
    repeats: tuple[int, int, int] = (20, 20, 20),
    lattice_constant: float = 3.61,
) -> CrystalStructure:
    """Build an FCC crystal selected by a TPMS implicit field."""
    return TpmsBuilder(
        unit_cell=FCC("Cu"),
        repeats=repeats,
        lattice_constant=lattice_constant,
        surface=surface,
        mode=mode,
    ).build()


def main() -> None:
    """Parse command-line options, build the TPMS crystal, and write it."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--surface",
        choices=("d-schwarz", "p-schwarz", "gyroid", "iwp"),
        default="gyroid",
    )
    parser.add_argument(
        "--mode",
        choices=("negative", "positive", "shell"),
        default="negative",
    )
    parser.add_argument("--repeats", type=int, nargs=3, default=(20, 20, 20))
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--output", type=Path, default=Path("gyroid.data"))
    options = parser.parse_args()

    structure = build_tpms(
        surface=options.surface,
        mode=options.mode,
        repeats=tuple(options.repeats),
        lattice_constant=options.lattice_constant,
    )
    structure.write(options.output)
    print(f"Wrote {structure.number_of_atoms} atoms to {options.output}")


if __name__ == "__main__":
    main()
