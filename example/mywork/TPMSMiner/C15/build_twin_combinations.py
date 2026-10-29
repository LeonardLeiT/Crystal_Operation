"""Generate the constrained C15 twin-combination screening set."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np

from build_c15 import FIXED_SIZE, build_c15_laves, _jsonable


ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = ROOT / "twin_combinations.json"


def _cubic_symmetries() -> tuple[np.ndarray, ...]:
    operations: list[np.ndarray] = []
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product((-1.0, 1.0), repeat=3):
            operation = np.zeros((3, 3))
            operation[np.arange(3), permutation] = signs
            if np.linalg.det(operation) > 0.0:
                operations.append(operation)
    return tuple(operations)


CUBIC_SYMMETRIES = _cubic_symmetries()


def _disorientation_degrees(first: np.ndarray, second: np.ndarray) -> float:
    angles = []
    for symmetry in CUBIC_SYMMETRIES:
        relative = second @ symmetry @ first.T
        cosine = np.clip((np.trace(relative) - 1.0) / 2.0, -1.0, 1.0)
        angles.append(float(np.degrees(np.arccos(cosine))))
    return min(angles)


def _cna_counts(path: Path) -> dict[str, int]:
    try:
        import warnings

        warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
        from ovito.io import import_file
        from ovito.modifiers import CommonNeighborAnalysisModifier
    except ImportError as exc:
        raise RuntimeError("--cna requires the optional OVITO package") from exc
    pipeline = import_file(path)
    pipeline.modifiers.append(
        CommonNeighborAnalysisModifier(
            mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff
        )
    )
    data = pipeline.compute()
    values = np.asarray(data.particles["Structure Type"])
    counts = np.bincount(values, minlength=5)
    return {
        "other": int(counts[0]),
        "fcc": int(counts[1]),
        "hcp": int(counts[2]),
        "bcc": int(counts[3]),
        "icosahedral": int(counts[4]),
    }


def _load_config(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported twin-combination schema")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--overlap-distance-factor", type=float, default=0.5)
    parser.add_argument("--cna", action="store_true")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "output" / "twin_combinations",
    )
    options = parser.parse_args()

    payload = _load_config(options.config)
    options.output_dir.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "config": str(options.config.resolve()),
        "rules": payload["rules"],
        "fixed_size": FIXED_SIZE,
        "size_series_enabled": False,
        "lattice_constant": options.lattice_constant,
        "overlap_distance_factor": options.overlap_distance_factor,
        "combinations": [],
    }
    expected_angles = {"sigma3": 60.0, "sigma5": 36.86989764584402}
    for candidate in payload["combinations"]:
        structure = build_c15_laves(
            size=FIXED_SIZE,
            lattice_constant=options.lattice_constant,
            overlap_distance_factor=options.overlap_distance_factor,
            z16_twin_type=candidate["z16_twin_type"],
            z12_twin_type=candidate["z12_twin_type"],
            swap_classes=tuple(candidate["swap_classes"]),
            name=candidate["id"],
        )
        orientations = np.asarray(structure.metadata["grain_orientations"])
        pair_records = []
        for pair in structure.metadata["twin_pairs"]:
            first, second = (value - 1 for value in pair["pair"])
            measured = _disorientation_degrees(
                orientations[first], orientations[second]
            )
            expected = expected_angles[pair["twin_type"]]
            if not np.isclose(measured, expected, atol=1e-8):
                raise RuntimeError(
                    f"{candidate['id']} pair {pair['pair']} has angle "
                    f"{measured}, expected {expected}"
                )
            pair_records.append({**pair, "measured_disorientation": measured})

        output = options.output_dir / f"{candidate['id']}.data"
        structure.write(output, format="lammps-data", atom_style="atomic")
        record = {
            **candidate,
            "output": output.name,
            "number_of_atoms": structure.number_of_atoms,
            "sigma3_pair_count": structure.metadata["sigma3_pair_count"],
            "sigma5_pair_count": structure.metadata["sigma5_pair_count"],
            "overlap_atoms_removed": structure.metadata["overlap_atoms_removed"],
            "overlap_cutoff": structure.metadata["overlap_cutoff"],
            "grain_counts": structure.metadata["grain_counts"],
            "ordinary_grain_ids": structure.metadata["ordinary_grain_ids"],
            "ordinary_orientation_seed": structure.metadata[
                "background_orientation_seed"
            ],
            "ordinary_csl_exclusion_degrees": structure.metadata[
                "ordinary_csl_exclusion_degrees"
            ],
            "ordinary_minimum_csl_separation_degrees": structure.metadata[
                "ordinary_minimum_csl_separation_degrees"
            ],
            "ordinary_orientations_crystal_to_lab": orientations[
                [grain_id - 1 for grain_id in structure.metadata["ordinary_grain_ids"]]
            ].tolist(),
            "ordinary_lattice_anchors": np.asarray(
                structure.metadata["lattice_anchors"]
            )[[grain_id - 1 for grain_id in structure.metadata["ordinary_grain_ids"]]].tolist(),
            "pairs": pair_records,
        }
        if options.cna:
            record["cna"] = _cna_counts(output)
        manifest["combinations"].append(record)
        print(
            f"{candidate['id']}: {structure.number_of_atoms} atoms, "
            f"Sigma3={record['sigma3_pair_count']}, "
            f"Sigma5={record['sigma5_pair_count']}"
        )

    manifest_path = options.output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(_jsonable(manifest), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"Wrote {len(manifest['combinations'])} combinations to {options.output_dir.resolve()}")


if __name__ == "__main__":
    main()
