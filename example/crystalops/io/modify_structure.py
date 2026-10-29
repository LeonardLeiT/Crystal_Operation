"""Read, transform, and write one structure through the common API."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops import CrystalStructure
from crystalops.operations import rotation_matrix


def modify_structure(
    structure: CrystalStructure,
    *,
    repeats: tuple[int, int, int] = (2, 1, 1),
    z_angle: float = 30.0,
    translation: tuple[float, float, float] = (2.0, 0.0, 0.0),
) -> CrystalStructure:
    """Repeat, rotate, and translate a crystal structure."""
    result = structure.repeated(repeats)
    result = result.rotated(rotation_matrix(z_angle=z_angle))
    return result.translated(translation)


def main() -> None:
    """Read command-line paths, transform the input, and write the result."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, nargs="?", default=Path("fcc.data"))
    parser.add_argument("--output", type=Path, default=Path("fcc_modified.data"))
    parser.add_argument("--repeats", type=int, nargs=3, default=(2, 1, 1))
    parser.add_argument("--z-angle", type=float, default=30.0)
    parser.add_argument("--translation", type=float, nargs=3, default=(2.0, 0.0, 0.0))
    options = parser.parse_args()

    source = CrystalStructure.read(options.input)
    result = modify_structure(
        source,
        repeats=tuple(options.repeats),
        z_angle=options.z_angle,
        translation=tuple(options.translation),
    )
    result.write(options.output)
    print(f"Wrote {result.number_of_atoms} atoms to {options.output}")


if __name__ == "__main__":
    main()
