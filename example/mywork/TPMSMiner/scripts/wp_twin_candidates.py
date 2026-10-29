"""Generate constrained experimental W-P twin-network candidates.

The A15 geometry, equal-volume power cells, CTB anchors, and overlap handling
come from the validated :class:`crystalops.WeairePhelanCrystal` baseline.
Only the paper-defined Sigma3/Sigma5 orientation pair assigned to X, Y, and Z
and an optional T1/T2 exchange are varied here.  Experimental combinations are
kept out of the public CrystalOps builder registry.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

# Allow direct execution from a source checkout without requiring installation.
REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from crystalops import FCC, CrystalStructure, WeairePhelanCrystal
from crystalops.builders._lattice import orientation_from_rows


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "configs" / "wp_twin_candidates.json"
AXES = ("x", "y", "z")
PAIR_SLICES = {"x": slice(2, 4), "y": slice(4, 6), "z": slice(6, 8)}
EXPECTED_ANGLES = {"sigma3": 60.0, "sigma5": 36.86989764584402}


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


def _disorientation_angle(first: np.ndarray, second: np.ndarray) -> float:
    angles = []
    for symmetry in CUBIC_SYMMETRIES:
        relative = second @ symmetry @ first.T
        cosine = np.clip((np.trace(relative) - 1.0) / 2.0, -1.0, 1.0)
        angles.append(float(np.degrees(np.arccos(cosine))))
    return min(angles)


def compose_orientation_rows(candidate: dict[str, Any]) -> tuple[Any, ...]:
    """Compose eight Table-1 orientations from three constrained pair choices."""
    pair_types = candidate["pair_types"]
    swaps = set(candidate.get("swap_pairs", ()))
    unknown_axes = (set(pair_types) | swaps) - set(AXES)
    if unknown_axes:
        raise ValueError(f"unknown axes in {candidate['id']}: {sorted(unknown_axes)}")

    identity_rows = WeairePhelanCrystal.ORIENTATION_ROWS["sigma3"][:2]
    rows: list[Any] = list(identity_rows)
    for axis in AXES:
        family = pair_types[axis]
        if family not in EXPECTED_ANGLES:
            raise ValueError(
                f"{candidate['id']} axis {axis} must use sigma3 or sigma5"
            )
        pair = list(
            WeairePhelanCrystal.ORIENTATION_ROWS[family][PAIR_SLICES[axis]]
        )
        if axis in swaps:
            pair.reverse()
        rows.extend(pair)
    return tuple(rows)


def validate_orientation_rows(
    candidate: dict[str, Any],
    rows: tuple[Any, ...],
) -> dict[str, float]:
    """Verify every experimental Tet pair retains its requested CSL angle."""
    orientations = [orientation_from_rows(item) for item in rows]
    angles: dict[str, float] = {}
    for axis in AXES:
        pair_slice = PAIR_SLICES[axis]
        first, second = orientations[pair_slice]
        angle = _disorientation_angle(first, second)
        expected = EXPECTED_ANGLES[candidate["pair_types"][axis]]
        if not np.isclose(angle, expected, atol=1e-8):
            raise ValueError(
                f"{candidate['id']} has invalid {axis} pair angle "
                f"{angle:.8f}; expected {expected:.8f}"
            )
        angles[axis] = angle
    return angles


def build_candidate(
    candidate: dict[str, Any],
    *,
    size: float,
    lattice_constant: float,
    overlap_distance_factor: float = 0.5,
    label_grains: bool = False,
    dod_orientation_seed: int | None = None,
) -> CrystalStructure:
    """Build one locally valid but globally experimental twin combination."""
    rows = compose_orientation_rows(candidate)
    angles = validate_orientation_rows(candidate, rows)
    configured_class = type(
        f"_{candidate['id']}",
        (WeairePhelanCrystal,),
        {"ORIENTATION_ROWS": {"sigma3": rows, "sigma5": rows}},
    )
    structure = configured_class(
        unit_cell=FCC("Cu"),
        size=size,
        lattice_constant=lattice_constant,
        twin_type="sigma3",
        overlap_distance_factor=overlap_distance_factor,
        dod_orientation_seed=dod_orientation_seed,
        label_grains=label_grains,
        name=candidate["id"],
    ).build()
    orientations = np.asarray(structure.metadata["grain_orientations"])
    dod_minimum_csl_separation = min(
        abs(_disorientation_angle(orientations[dod], orientations[other]) - target)
        for dod in (0, 1)
        for other in range(len(orientations))
        if other != dod
        for target in EXPECTED_ANGLES.values()
    )
    pair_types = tuple(candidate["pair_types"][axis] for axis in AXES)
    structure.metadata.update(
        {
            "builder": "TPMSMiner.WPTwinCandidate",
            "candidate_id": candidate["id"],
            "experimental": True,
            "pair_types_xyz": pair_types,
            "swap_pairs": tuple(candidate.get("swap_pairs", ())),
            "sigma3_axis_count": pair_types.count("sigma3"),
            "anchor_rank": pair_types.count("sigma3"),
            "pair_disorientation_degrees": angles,
            "screening_tier": candidate["tier"],
            "hypothesis": candidate["hypothesis"],
            "dod_orientation_seed": dod_orientation_seed,
            "dod_orientation_sources": structure.metadata["orientation_sources"][:2],
            "dod_minimum_csl_separation_degrees": dod_minimum_csl_separation,
        }
    )
    return structure


def load_candidates(path: Path = DEFAULT_CONFIG) -> list[dict[str, Any]]:
    """Load the experimental W-P candidate definitions."""
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise ValueError("unsupported candidate schema")
    return list(payload["candidates"])


def get_candidate(
    candidate_id: str,
    path: Path = DEFAULT_CONFIG,
) -> dict[str, Any]:
    """Return one candidate definition by its stable identifier."""
    for candidate in load_candidates(path):
        if candidate["id"] == candidate_id:
            return candidate
    available = ", ".join(item["id"] for item in load_candidates(path))
    raise ValueError(
        f"unknown candidate ID {candidate_id!r}; available candidates: {available}"
    )


def _cna_counts(path: Path) -> dict[str, int]:
    try:
        import warnings

        warnings.filterwarnings("ignore", message=".*OVITO.*PyPI")
        from ovito.io import import_file
        from ovito.modifiers import CommonNeighborAnalysisModifier
    except ImportError as exc:
        raise RuntimeError("--cna requires the optional OVITO Python package") from exc

    pipeline = import_file(path)
    pipeline.modifiers.append(
        CommonNeighborAnalysisModifier(
            mode=CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff
        )
    )
    data = pipeline.compute()
    structure_types = np.asarray(data.particles["Structure Type"])
    counts = np.bincount(structure_types, minlength=5)
    return {
        "other": int(counts[0]),
        "fcc": int(counts[1]),
        "hcp": int(counts[2]),
        "bcc": int(counts[3]),
        "icosahedral": int(counts[4]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--size", type=float, default=14.0)
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--overlap-distance-factor", type=float, default=0.5)
    parser.add_argument(
        "--select",
        default="recommended",
        help="recommended, all, or a comma-separated list of candidate IDs",
    )
    parser.add_argument("--label-grains", action="store_true")
    parser.add_argument("--dod-orientation-seed", type=int, default=None)
    parser.add_argument("--cna", action="store_true")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "results" / "wp_twin_candidates",
    )
    options = parser.parse_args()

    candidates = load_candidates(options.config)
    if options.select == "recommended":
        candidates = [item for item in candidates if item["recommended"]]
    elif options.select != "all":
        requested = set(options.select.split(","))
        candidates = [item for item in candidates if item["id"] in requested]
        missing = requested - {item["id"] for item in candidates}
        if missing:
            raise ValueError(f"unknown candidate IDs: {sorted(missing)}")

    options.output_dir.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, Any] = {
        "config": str(options.config.resolve()),
        "size": options.size,
        "lattice_constant": options.lattice_constant,
        "grain_size": options.size * options.lattice_constant,
        "box_length": 2.0 * options.size * options.lattice_constant,
        "overlap_distance_factor": options.overlap_distance_factor,
        "candidates": [],
    }
    for candidate in candidates:
        structure = build_candidate(
            candidate,
            size=options.size,
            lattice_constant=options.lattice_constant,
            overlap_distance_factor=options.overlap_distance_factor,
            label_grains=options.label_grains,
            dod_orientation_seed=options.dod_orientation_seed,
        )
        output = options.output_dir / f"{candidate['id']}.data"
        structure.write(output)
        record = {
            "id": candidate["id"],
            "output": str(output.resolve()),
            "tier": candidate["tier"],
            "pair_types": candidate["pair_types"],
            "swap_pairs": candidate.get("swap_pairs", []),
            "anchor_rank": structure.metadata["anchor_rank"],
            "pair_disorientation_degrees": structure.metadata[
                "pair_disorientation_degrees"
            ],
            "number_of_atoms": structure.number_of_atoms,
            "grain_counts": structure.metadata["grain_counts"],
            "overlap_atoms_removed": structure.metadata[
                "overlap_atoms_removed"
            ],
            "overlap_cutoff": structure.metadata["overlap_cutoff"],
            "hypothesis": candidate["hypothesis"],
            "dod_orientation_seed": structure.metadata["dod_orientation_seed"],
            "dod_orientation_sources": structure.metadata["orientation_sources"][:2],
            "dod_orientations_crystal_to_lab": np.asarray(
                structure.metadata["grain_orientations"]
            )[:2].tolist(),
            "dod_lattice_origins_reduced": np.asarray(
                structure.metadata["lattice_origins_reduced"]
            )[:2].tolist(),
            "dod_minimum_csl_separation_degrees": structure.metadata[
                "dod_minimum_csl_separation_degrees"
            ],
        }
        if options.cna:
            record["cna"] = _cna_counts(output)
        manifest["candidates"].append(record)
        print(
            f"{candidate['id']}: {structure.number_of_atoms} atoms, "
            f"anchor rank {record['anchor_rank']}, wrote {output}"
        )

    manifest_path = options.output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(f"Wrote manifest to {manifest_path}")


if __name__ == "__main__":
    main()
