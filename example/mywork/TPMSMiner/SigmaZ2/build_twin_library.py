"""Generate the independent SigmaZ2 twin-design library."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
from typing import Any

import numpy as np

from build_sigma_z2 import FIXED_SIZE, build_sigma_z2, _jsonable
from twin_designs import DEFAULT_CONFIG, load_design_library


ROOT = Path(__file__).resolve().parent


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
    values = np.asarray(pipeline.compute().particles["Structure Type"])
    counts = np.bincount(values, minlength=5)
    return {
        "other": int(counts[0]), "fcc": int(counts[1]), "hcp": int(counts[2]),
        "bcc": int(counts[3]), "icosahedral": int(counts[4]),
    }


def _build_record(
    candidate: dict[str, Any], twin_type: str, options: argparse.Namespace
) -> dict[str, Any]:
    pairs = tuple(candidate["pairs"])
    flat_grains = [grain for pair in pairs for grain in pair[:2]]
    if len(flat_grains) != len(set(flat_grains)):
        raise RuntimeError(f"non-disjoint design: {candidate['id']}")
    if candidate["inherits_sigma_pair_count"]:
        raise RuntimeError(f"design inherits an original Sigma face: {candidate['id']}")
    structure = build_sigma_z2(
        lattice_constant=options.lattice_constant,
        overlap_distance_factor=options.overlap_distance_factor,
        design_id=candidate["id"],
        twin_face_pairs=pairs,
        twin_type=twin_type,
    )
    orientations = np.asarray(structure.metadata["grain_orientations"])
    expected = 60.0 if twin_type == "sigma3" else 36.86989764584402
    pair_records = []
    for pair in structure.metadata["twin_pairs"]:
        first, second = (value - 1 for value in pair["pair"])
        measured = _disorientation(orientations[first], orientations[second])
        if not np.isclose(measured, expected, atol=1e-8):
            raise RuntimeError(f"invalid angle in {candidate['id']} pair {pair['pair']}")
        pair_records.append({**pair, "measured_disorientation": measured})
    output = options.output_dir / f"{candidate['id']}.data"
    structure.write(output, format="lammps-data", atom_style="atomic")
    record = {
        "id": candidate["id"],
        "hypothesis": candidate["hypothesis"],
        "twin_type": twin_type,
        "target_pairs": candidate["target_pairs"],
        "placement_counts": candidate["placement_counts"],
        "pair_class_counts": candidate["pair_class_counts"],
        "inherits_sigma_pair_count": candidate["inherits_sigma_pair_count"],
        "minimum_laguerre_face_margin": candidate["minimum_laguerre_face_margin"],
        "ordinary_grain_count": candidate["ordinary_grain_count"],
        "output": output.name,
        "number_of_atoms": structure.number_of_atoms,
        "box_lengths": np.diag(structure.cell).tolist(),
        "overlap_atoms_removed": structure.metadata["overlap_atoms_removed"],
        "pairs": pair_records,
    }
    if options.cna:
        record["cna"] = _cna(output)
    print(
        f"{candidate['id']}: {structure.number_of_atoms} atoms, "
        f"{len(pairs)} {twin_type} pairs, placements={candidate['placement_counts']}"
    )
    return record


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--overlap-distance-factor", type=float, default=0.5)
    parser.add_argument("--cna", action="store_true")
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "output" / "twin_library"
    )
    options = parser.parse_args()
    library = load_design_library(options.config)
    options.output_dir.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "fixed_size": FIXED_SIZE,
        "size_series_enabled": False,
        "config": str(options.config.resolve()),
        "available_non_sigma_true_faces": library["available_edge_count"],
        "rules": library["rules"],
        "designs": [],
    }
    for candidate in library["selected_designs"]:
        manifest["designs"].append(_build_record(candidate, "sigma3", options))
    source = next(
        design for design in library["selected_designs"]
        if design["id"] == library["control"]["source_design"]
    )
    control = {
        **source,
        "id": library["control"]["id"],
        "hypothesis": library["control"]["hypothesis"],
    }
    manifest["designs"].append(_build_record(control, "sigma5", options))
    (options.output_dir / "manifest.json").write_text(
        json.dumps(_jsonable(manifest), indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(f"Wrote {len(manifest['designs'])} designs to {options.output_dir.resolve()}")


if __name__ == "__main__":
    main()
