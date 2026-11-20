"""Mark atoms within a physical distance of the DSC builder's curved cuts."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

from crystalops.analysis import select_dsc_control_band
from crystalops.core import CrystalStructure
from crystalops.io import write_selection_dump


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "input", type=Path, nargs="?", default=Path("Schwarz"),
        help="a LAMMPS .data file or a directory of .data files",
    )
    parser.add_argument("--behind", type=float, default=2.0, help="depth on the negative side (Å)")
    parser.add_argument("--ahead", type=float, default=2.0, help="depth on the positive side (Å)")
    parser.add_argument("--lattice-constant", type=float, default=3.61, help="DSC lattice constant (Å)")
    parser.add_argument("--mesh-spacing", type=float, default=1.0, help="maximum mesh grid spacing (Å)")
    parser.add_argument(
        "--output-dir", type=Path, default=Path("analysis_results/schwarz_control_band")
    )
    options = parser.parse_args()

    if options.input.is_dir():
        files = sorted(options.input.glob("*.data"))
    elif options.input.is_file() and options.input.suffix.lower() == ".data":
        files = [options.input]
    else:
        files = []
    if not files:
        parser.error(f"no LAMMPS .data files found at {options.input}")

    for file in files:
        structure = CrystalStructure.read(file)
        marked = select_dsc_control_band(
            structure,
            behind=options.behind,
            ahead=options.ahead,
            lattice_constant=options.lattice_constant,
            mesh_spacing=options.mesh_spacing,
        )
        output = options.output_dir / f"{file.stem}_control_band.dump"
        write_selection_dump(marked, output, ovito_selection=True)
        selected = int(marked.properties["select"].sum())
        print(f"{file.name}: {selected}/{marked.number_of_atoms} selected -> {output}")


if __name__ == "__main__":
    main()
