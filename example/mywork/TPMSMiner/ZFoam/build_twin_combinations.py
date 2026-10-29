"""Generate six fixed-size Z-foam twin/control combinations."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np

from build_z_foam import FIXED_SIZE, _jsonable, build_z_foam


ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = ROOT / "twin_combinations.json"
CLASSES = ("Z12", "Z14", "Z15")


def _cubic_symmetries() -> tuple[np.ndarray, ...]:
    operations = []
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product((-1.0, 1.0), repeat=3):
            operation = np.zeros((3, 3))
            operation[np.arange(3), permutation] = signs
            if np.linalg.det(operation) > 0.0:
                operations.append(operation)
    return tuple(operations)


CUBIC_SYMMETRIES = _cubic_symmetries()


def _disorientation(first: np.ndarray, second: np.ndarray) -> float:
    angles = []
    for symmetry in CUBIC_SYMMETRIES:
        relative = second @ symmetry @ first.T
        cosine = np.clip((np.trace(relative) - 1.0) / 2.0, -1.0, 1.0)
        angles.append(float(np.degrees(np.arccos(cosine))))
    return min(angles)


def _cna(path: Path) -> dict[str, int]:
    try:
        import warnings
        warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
        from ovito.io import import_file
        from ovito.modifiers import CommonNeighborAnalysisModifier
    except ImportError as exc:
        raise RuntimeError("--cna requires OVITO") from exc
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
        "other": int(counts[0]), "fcc": int(counts[1]), "hcp": int(counts[2]),
        "bcc": int(counts[3]), "icosahedral": int(counts[4]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--overlap-distance-factor", type=float, default=0.5)
    parser.add_argument("--cna", action="store_true")
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "output" / "twin_combinations"
    )
    options = parser.parse_args()
    payload = json.loads(options.config.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported Z-foam combination schema")
    options.output_dir.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "fixed_size": FIXED_SIZE,
        "config": str(options.config.resolve()),
        "rules": payload["rules"],
        "combinations": [],
    }
    expected = {"sigma3": 60.0, "sigma5": 36.86989764584402}
    for candidate in payload["combinations"]:
        sigma3 = set(candidate["sigma3_classes"])
        if sigma3 - set(CLASSES):
            raise ValueError(f"unknown classes: {sorted(sigma3-set(CLASSES))}")
        types = {key: "sigma3" if key in sigma3 else "sigma5" for key in CLASSES}
        structure = build_z_foam(
            lattice_constant=options.lattice_constant,
            overlap_distance_factor=options.overlap_distance_factor,
            z12_twin_type=types["Z12"],
            z14_twin_type=types["Z14"],
            z15_twin_type=types["Z15"],
            name=candidate["id"],
        )
        orientations = np.asarray(structure.metadata["grain_orientations"])
        pairs = []
        for pair in structure.metadata["twin_pairs"]:
            first, second = (value - 1 for value in pair["pair"])
            measured = _disorientation(orientations[first], orientations[second])
            if not np.isclose(measured, expected[pair["twin_type"]], atol=1e-8):
                raise RuntimeError(f"invalid angle in {candidate['id']} pair {pair['pair']}")
            pairs.append({**pair, "measured_disorientation": measured})
        output = options.output_dir / f"{candidate['id']}.data"
        structure.write(output, format="lammps-data", atom_style="atomic")
        record = {
            **candidate,
            "pair_types": types,
            "output": output.name,
            "number_of_atoms": structure.number_of_atoms,
            "box_vectors": structure.cell.tolist(),
            "sigma3_pair_count": structure.metadata["sigma3_pair_count"],
            "sigma5_pair_count": structure.metadata["sigma5_pair_count"],
            "overlap_atoms_removed": structure.metadata["overlap_atoms_removed"],
            "grain_counts": structure.metadata["grain_counts"],
            "pairs": pairs,
        }
        if options.cna:
            record["cna"] = _cna(output)
        manifest["combinations"].append(record)
        print(
            f"{candidate['id']}: {structure.number_of_atoms} atoms, "
            f"Sigma3={record['sigma3_pair_count']}, Sigma5={record['sigma5_pair_count']}"
        )
    (options.output_dir / "manifest.json").write_text(
        json.dumps(_jsonable(manifest), indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Wrote {len(manifest['combinations'])} combinations to {options.output_dir.resolve()}")


if __name__ == "__main__":
    main()
