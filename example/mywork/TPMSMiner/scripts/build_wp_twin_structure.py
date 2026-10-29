"""Build one experimental W-P twin combination as a LAMMPS data file.

This is the TPMSMiner counterpart of
``example/crystalops/builders/build_weaire_phelan.py``.  It writes the actual
Cu atomistic configuration (``atom_style atomic``), not just a candidate
description or a geometry report.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


# Allow direct execution from a source checkout without installation.
REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from crystalops import CrystalStructure

from wp_twin_candidates import (
    DEFAULT_CONFIG,
    build_candidate,
    get_candidate,
    load_candidates,
)


DEFAULT_CANDIDATE = "wp_s5_s3_s3_open_x"


def build_wp_twin_structure(
    *,
    candidate_id: str = DEFAULT_CANDIDATE,
    size: float = 14.0,
    lattice_constant: float = 3.61,
    overlap_distance_factor: float = 0.5,
    config: Path = DEFAULT_CONFIG,
    label_grains: bool = False,
    dod_orientation_seed: int | None = None,
) -> CrystalStructure:
    """Build one complete atomistic W-P candidate for LAMMPS."""
    candidate = get_candidate(candidate_id, config)
    return build_candidate(
        candidate,
        size=size,
        lattice_constant=lattice_constant,
        overlap_distance_factor=overlap_distance_factor,
        label_grains=label_grains,
        dod_orientation_seed=dod_orientation_seed,
    )


def main() -> None:
    candidate_ids = [item["id"] for item in load_candidates(DEFAULT_CONFIG)]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidate",
        choices=candidate_ids,
        default=DEFAULT_CANDIDATE,
        help="constrained W-P twin combination to build",
    )
    parser.add_argument("--size", type=float, default=14.0)
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--overlap-distance-factor", type=float, default=0.5)
    parser.add_argument("--label-grains", action="store_true")
    parser.add_argument("--dod-orientation-seed", type=int, default=None)
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="LAMMPS data path; defaults to results/lammps/<candidate>.data",
    )
    options = parser.parse_args()

    structure = build_wp_twin_structure(
        candidate_id=options.candidate,
        size=options.size,
        lattice_constant=options.lattice_constant,
        overlap_distance_factor=options.overlap_distance_factor,
        label_grains=options.label_grains,
        dod_orientation_seed=options.dod_orientation_seed,
    )
    output = options.output
    if output is None:
        output = (
            Path(__file__).resolve().parents[1]
            / "results"
            / "lammps"
            / f"{options.candidate}.data"
        )
    structure.write(output, format="lammps-data", atom_style="atomic")

    pair_types = "/".join(structure.metadata["pair_types_xyz"])
    print(
        f"Wrote {structure.number_of_atoms} Cu atoms in 8 complete W-P grains "
        f"to {output.resolve()} (X/Y/Z={pair_types})."
    )


if __name__ == "__main__":
    main()
