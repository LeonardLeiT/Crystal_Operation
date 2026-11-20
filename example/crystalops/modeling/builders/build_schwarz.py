"""Build a six-domain Diamond Schwarz Crystal (DSC) or a Schwarz P crystal."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from numpy.typing import ArrayLike

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[4]))

from crystalops.core import CrystalStructure
from crystalops.modeling import SchwarzCrystal


DSC_SERIES = (1.0, 0.7, 1.2, 1.3, 1.4, 1.6)


def build_schwarz(
    *,
    size: int | None = None,
    lattice_constant: float = 3.61,
    schwarz_type: str = "DSC",
    random_seed: int | None = None,
    orientation: ArrayLike | None = None,
    max_box_expansion: float = 0.25,
    max_orientation_index: int = 5,
    z_ratio: float | None = None,
    refine_boundaries: bool = False,
) -> CrystalStructure:
    """Construct the selected structure through ``SchwarzCrystal``."""
    if size is None:
        size = 30 if z_ratio is not None else 10
    return SchwarzCrystal(
        size=size,
        lattice_constant=lattice_constant,
        schwarz_type=schwarz_type,
        random_seed=random_seed,
        orientation=orientation,
        max_box_expansion=max_box_expansion,
        max_orientation_index=max_orientation_index,
        z_ratio=z_ratio,
        refine_boundaries=refine_boundaries,
    ).build()


def write_dsc_labels(structure: CrystalStructure, path: Path) -> None:
    """Write both DSC domain and grain IDs alongside the LAMMPS atom IDs."""
    atom_ids = structure.properties.get("atom_ids")
    if atom_ids is None:
        atom_ids = range(1, structure.number_of_atoms + 1)
    labels_path = path.with_suffix(".grain_ids.csv")
    with labels_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("atom_id", "grain_id", "domain_id"))
        writer.writerows(zip(
            atom_ids, structure.properties["grain_ids"], structure.properties["domain_ids"]
        ))


def main() -> None:
    """Build one configuration or the six DSC box ratios."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, help="x/y length in conventional FCC cells")
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--type", choices=("DSC", "PSC"), default="DSC")
    ratio_options = parser.add_mutually_exclusive_group()
    ratio_options.add_argument("--z-ratio", type=float, help="build periodic DSC with this Lz/Lx")
    ratio_options.add_argument("--series", action="store_true", help="build all six DSC ratios")
    parser.add_argument(
        "--refine-boundaries", action="store_true",
        help="resolve competing DSC domains and repair detached islands in a noncubic box",
    )
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--max-box-expansion", type=float, default=0.25)
    parser.add_argument("--max-orientation-index", type=int, default=5)
    output_options = parser.add_mutually_exclusive_group()
    output_options.add_argument("--output", type=Path, help="path for one .data file")
    output_options.add_argument("--output-dir", type=Path, help="directory for --series")
    options = parser.parse_args()

    if (options.series or options.z_ratio is not None) and options.type != "DSC":
        parser.error("--series and --z-ratio require DSC")
    if options.output_dir is not None and not options.series:
        parser.error("--output-dir requires --series")
    if options.series and options.output is not None:
        parser.error("--series writes multiple files; use --output-dir")
    if options.series and options.size is not None and (options.size < 30 or options.size % 30):
        parser.error("--series requires --size to be a positive multiple of 30")
    if options.refine_boundaries and (options.type != "DSC" or options.z_ratio in (None, 1.0)):
        parser.error("--refine-boundaries requires a noncubic DSC specified with --z-ratio")

    output_dir = options.output_dir or Path("analysis_results/dsc_twin_facets")
    ratios = DSC_SERIES if options.series else (options.z_ratio,)
    for ratio in ratios:
        output = (
            output_dir / f"dsc_1-1-{ratio:g}.data"
            if options.series else options.output or Path("schwarz.data")
        )
        structure = build_schwarz(
            size=options.size,
            lattice_constant=options.lattice_constant,
            schwarz_type=options.type,
            random_seed=options.seed,
            max_box_expansion=options.max_box_expansion,
            max_orientation_index=options.max_orientation_index,
            z_ratio=ratio,
            refine_boundaries=options.refine_boundaries,
        )
        structure.write(output)
        if ratio is not None:
            write_dsc_labels(structure, output)
            print(
                f"Wrote {output}: {structure.number_of_atoms} DSC atoms at 1:1:{ratio:g}, "
                f"{structure.metadata['domain_counts']}; "
                f"removed {structure.metadata['removed_overlap_atoms']} overlapping atoms",
                flush=True,
            )
        else:
            detail = structure.metadata.get("domain_counts", structure.metadata.get("grain_counts"))
            print(f"Wrote {structure.number_of_atoms} {options.type} atoms ({detail}) to {output}")


if __name__ == "__main__":
    main()
