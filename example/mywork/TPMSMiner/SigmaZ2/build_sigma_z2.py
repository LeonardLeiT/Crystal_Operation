"""Build a Z-doubled Frank-Kasper sigma foam from an independent twin design."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
from numpy.typing import NDArray


REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
SIGMA_SOURCE = Path(__file__).resolve().parents[1] / "Sigma"
for search_path in (REPOSITORY_ROOT, SIGMA_SOURCE):
    if str(search_path) not in sys.path:
        sys.path.insert(0, str(search_path))

from crystalops import CrystalStructure, FCC
from crystalops.builders._lattice import candidate_lattice
from crystalops.builders.polycrystal import (
    _nearest_periodic_site_distance,
    _remove_periodic_overlaps,
)
from build_sigma import (
    COORDINATION_NUMBERS,
    C_OVER_A,
    EQUAL_VOLUME_WEIGHTS,
    FIXED_SIZE,
    GROUP_IDS,
    PERIODIC_SHIFTS,
    SIGMA_CENTERS,
    SITE_CLASSES,
    TWIN_SPECS,
    WYCKOFF_LABELS,
    _align_vector,
    _axis_angle,
    _jsonable,
)
from twin_designs import DEFAULT_CONFIG, load_design


FloatArray = NDArray[np.float64]
ROOT = Path(__file__).resolve().parent
RANDOM_SEED = 20260930

# Two copies of the original 30 sigma sites in a box twice as long in Z.
SIGMA_Z2_CENTERS = np.concatenate(
    (
        SIGMA_CENTERS * np.asarray([1.0, 1.0, 0.5]),
        SIGMA_CENTERS * np.asarray([1.0, 1.0, 0.5])
        + np.asarray([0.0, 0.0, 0.5]),
    )
)
Z2_SITE_CLASSES = SITE_CLASSES + SITE_CLASSES
Z2_WYCKOFF_LABELS = WYCKOFF_LABELS + WYCKOFF_LABELS
Z2_COORDINATION_NUMBERS = COORDINATION_NUMBERS + COORDINATION_NUMBERS
Z2_GROUP_IDS = np.tile(GROUP_IDS, 2)

def _random_rotation(rng: np.random.Generator) -> FloatArray:
    """Return a uniform proper rotation using a normalized quaternion."""
    quaternion = rng.normal(size=4)
    quaternion /= np.linalg.norm(quaternion)
    w, x, y, z = quaternion
    return np.asarray(
        [
            [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
            [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
            [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
        ],
        dtype=np.float64,
    )


def _orientations_and_anchors(
    cell: FloatArray,
    twin_face_pairs: tuple[tuple[int, int, tuple[float, float, float], str], ...],
    twin_type: str,
) -> tuple[FloatArray, FloatArray, FloatArray, list[dict[str, Any]]]:
    rng = np.random.default_rng(RANDOM_SEED)
    orientations = np.asarray([_random_rotation(rng) for _ in range(60)])
    anchors = SIGMA_Z2_CENTERS @ cell
    construction_centers = SIGMA_Z2_CENTERS.copy()
    records: list[dict[str, Any]] = []

    used: set[int] = set()
    crystal_axis, angle, sigma = TWIN_SPECS[twin_type]
    for pair_index, (first, second, vector_values, placement) in enumerate(twin_face_pairs):
        if first in used or second in used:
            raise ValueError("each SigmaZ2 grain may occur in at most one twin pair")
        used.update((first, second))
        vector = np.asarray(vector_values, dtype=np.float64)
        physical_vector = vector @ cell
        normal = physical_vector / np.linalg.norm(physical_vector)
        base = _align_vector(crystal_axis, normal)
        sign = 1.0 if pair_index % 2 == 0 else -1.0
        twin = _axis_angle(normal, sign * angle) @ base
        orientations[first] = base
        orientations[second] = twin

        anchor_reduced = SIGMA_Z2_CENTERS[first] + 0.5 * vector
        shared_anchor = anchor_reduced @ cell
        anchors[first] = shared_anchor
        anchors[second] = shared_anchor
        construction_centers[second] = SIGMA_Z2_CENTERS[first] + vector
        records.append(
            {
                "pair": [first + 1, second + 1],
                "placement": placement,
                "site_classes": [Z2_SITE_CLASSES[first], Z2_SITE_CLASSES[second]],
                "wyckoff_labels": [Z2_WYCKOFF_LABELS[first], Z2_WYCKOFF_LABELS[second]],
                "periodic_face_vector_reduced": vector.tolist(),
                "face_normal": normal.tolist(),
                "anchor_reduced_unwrapped": anchor_reduced.tolist(),
                "twin_type": twin_type,
                "sigma": sigma,
                "misorientation_degrees": angle,
                "coherent_phase_anchor": True,
            }
        )
    return orientations, anchors, construction_centers, records


def build_sigma_z2(
    *,
    lattice_constant: float = 3.61,
    overlap_distance_factor: float = 0.5,
    equal_volume: bool = True,
    label_grains: bool = False,
    design_id: str = "03_three_zone_balanced_15",
    design_config: Path = DEFAULT_CONFIG,
    twin_face_pairs: tuple[tuple[int, int, tuple[float, float, float], str], ...] | None = None,
    twin_type: str = "sigma3",
) -> CrystalStructure:
    """Build one fixed-size, 60-cell SigmaZ2 twin precursor."""
    if lattice_constant <= 0.0 or overlap_distance_factor <= 0.0:
        raise ValueError("lattice constant and overlap factor must be positive")
    if twin_type not in TWIN_SPECS:
        raise ValueError(f"unknown twin type {twin_type!r}")
    design = load_design(design_id, design_config) if twin_face_pairs is None else None
    if twin_face_pairs is None:
        twin_face_pairs = tuple(design["pairs"])

    mean_grain_length = FIXED_SIZE * lattice_constant
    a_length = np.cbrt(30.0 / C_OVER_A) * mean_grain_length
    cell = np.diag([a_length, a_length, 2.0 * C_OVER_A * a_length])
    orientations, anchors, construction_centers, twin_pairs = _orientations_and_anchors(
        cell, twin_face_pairs, twin_type
    )
    center_positions = construction_centers @ cell
    reduced_weights = EQUAL_VOLUME_WEIGHTS[Z2_GROUP_IDS] if equal_volume else np.zeros(60)
    weights = reduced_weights * a_length**2
    unit_cell = FCC("Cu")
    tolerance = max(a_length, 1.0) * 1e-10

    positions_by_grain: list[FloatArray] = []
    types_by_grain: list[NDArray[np.int64]] = []
    ids_by_grain: list[NDArray[np.int64]] = []
    for grain_index, (center, orientation, anchor) in enumerate(
        zip(center_positions, orientations, anchors, strict=True)
    ):
        normals: list[FloatArray] = []
        limits: list[float] = []
        for neighbor_index, neighbor in enumerate(SIGMA_Z2_CENTERS):
            for shift in PERIODIC_SHIFTS:
                if neighbor_index == grain_index and np.all(shift == 0.0):
                    continue
                relative = (neighbor + shift - construction_centers[grain_index]) @ cell
                normals.append(2.0 * relative)
                limits.append(
                    float(relative @ relative + weights[grain_index] - weights[neighbor_index])
                )
        candidates = candidate_lattice(unit_cell, lattice_constant, cell, orientation) + anchor
        relative_candidates = candidates - center
        # Do not materialize candidates x all-halfspaces at once: for the Z2
        # box that matrix is several GB.  Progressive blocks are equivalent
        # and quickly discard candidates outside the current Laguerre cell.
        normal_array = np.asarray(normals)
        limit_array = np.asarray(limits)
        inside = np.ones(len(relative_candidates), dtype=bool)
        for start in range(0, len(normal_array), 64):
            active = np.flatnonzero(inside)
            if not len(active):
                break
            stop = min(start + 64, len(normal_array))
            inside[active] = np.all(
                relative_candidates[active] @ normal_array[start:stop].T
                <= limit_array[None, start:stop] + tolerance,
                axis=1,
            )
        positions = candidates[inside]
        types = np.tile(unit_cell.atom_types, len(candidates) // unit_cell.number_of_sites)[inside]
        positions_by_grain.append(positions)
        types_by_grain.append(types)
        ids_by_grain.append(np.full(len(positions), grain_index + 1, dtype=np.int64))

    structure = CrystalStructure(
        positions=np.concatenate(positions_by_grain),
        cell=cell,
        atom_types=np.concatenate(types_by_grain),
        type_elements=unit_cell.type_elements,
        origin=np.zeros(3),
        pbc=(True, True, True),
        properties={"grain_ids": np.concatenate(ids_by_grain)},
        metadata={
            "builder": "TPMSMiner.SigmaZ2TwinLibrary",
            "experimental": True,
            "precursor": "Frank-Kasper-sigma-foam-Z2",
            "space_group": "P42/mnm",
            "grain_topology": "2 x (10xZ12 + 16xZ14 + 4xZ15)",
            "z_repetitions": 2,
            "c_over_a_original": C_OVER_A,
            "equal_volume": bool(equal_volume),
            "grain_weights_reduced": reduced_weights.copy(),
            "grain_centers_reduced": SIGMA_Z2_CENTERS.copy(),
            "construction_centers_reduced": construction_centers.copy(),
            "grain_centers": center_positions.copy(),
            "grain_orientations": orientations.copy(),
            "lattice_anchors": anchors.copy(),
            "coordination_numbers": Z2_COORDINATION_NUMBERS,
            "site_classes": Z2_SITE_CLASSES,
            "wyckoff_labels": Z2_WYCKOFF_LABELS,
            "twin_pairs": twin_pairs,
            "design_id": design_id,
            "twin_type": twin_type,
            "sigma3_pair_count": len(twin_pairs) if twin_type == "sigma3" else 0,
            "sigma5_pair_count": len(twin_pairs) if twin_type == "sigma5" else 0,
            "twin_grain_count": 2 * len(twin_pairs),
            "ordinary_grain_count": 60 - 2 * len(twin_pairs),
            "background_orientation_seed": RANDOM_SEED,
            "fixed_size": FIXED_SIZE,
            "mean_grain_length": float(mean_grain_length),
            "box_lengths": np.diag(cell).copy(),
            "lattice_constant": float(lattice_constant),
            "cropped_to_box": False,
        },
        name=f"sigma-z2-{design_id}-{twin_type}",
    )
    initial_count = structure.number_of_atoms
    nearest_distance = _nearest_periodic_site_distance(unit_cell, lattice_constant)
    overlap_cutoff = overlap_distance_factor * nearest_distance
    structure, removed = _remove_periodic_overlaps(structure, overlap_cutoff)
    ids, counts = np.unique(structure.properties["grain_ids"], return_counts=True)
    if not np.array_equal(ids, np.arange(1, 61)):
        raise RuntimeError("SigmaZ2 construction produced empty grains")
    lengths = np.diag(cell)
    outside = np.any(
        (structure.positions < -tolerance) | (structure.positions > lengths[None, :] + tolerance),
        axis=1,
    )
    structure.metadata.update(
        {
            "initial_atom_count": int(initial_count),
            "overlap_atoms_removed": int(removed),
            "overlap_cutoff": float(overlap_cutoff),
            "nearest_neighbor_distance": float(nearest_distance),
            "outside_box_atom_count": int(np.count_nonzero(outside)),
            "grain_counts": dict(zip(ids.tolist(), counts.tolist(), strict=True)),
        }
    )
    if label_grains:
        structure.atom_types = structure.properties["grain_ids"].copy()
        structure.type_elements = {grain_id: "Cu" for grain_id in range(1, 61)}
        structure.validate()
    return structure


def write_outputs(
    structure: CrystalStructure,
    output_directory: Path,
    *,
    filename: str = "01_sigma_z2_balanced15.data",
) -> dict[str, Any]:
    output_directory.mkdir(parents=True, exist_ok=True)
    data_path = output_directory / filename
    structure.write(data_path, format="lammps-data", atom_style="atomic")
    grain_ids = np.asarray(structure.properties["grain_ids"])
    twin_by_grain: dict[int, dict[str, Any]] = {}
    for twin in structure.metadata["twin_pairs"]:
        first, second = twin["pair"]
        twin_by_grain[first] = {**twin, "partner": second}
        twin_by_grain[second] = {**twin, "partner": first}

    cells: list[dict[str, Any]] = []
    for index in range(60):
        grain_id = index + 1
        twin = twin_by_grain.get(grain_id)
        cells.append(
            {
                "grain_id": grain_id,
                "z_copy": 1 if index < 30 else 2,
                "source_sigma_grain_id": index % 30 + 1,
                "site_class": Z2_SITE_CLASSES[index],
                "wyckoff_label": Z2_WYCKOFF_LABELS[index],
                "coordination_number": Z2_COORDINATION_NUMBERS[index],
                "atom_count": int(np.count_nonzero(grain_ids == grain_id)),
                "center_reduced": SIGMA_Z2_CENTERS[index].tolist(),
                "center_cartesian_unwrapped": np.asarray(structure.metadata["grain_centers"])[index].tolist(),
                "orientation_matrix_crystal_to_lab": np.asarray(structure.metadata["grain_orientations"])[index].tolist(),
                "boundary_role": (
                    f"designed_{twin['twin_type']}_twin" if twin else "ordinary_grain"
                ),
                "twin_partner_grain_id": twin["partner"] if twin else None,
                "twin_placement": twin["placement"] if twin else None,
                "shared_phase_anchor": twin["anchor_reduced_unwrapped"] if twin else None,
            }
        )

    manifest = {
        "schema_version": 1,
        "configuration": structure.metadata["design_id"],
        "output": data_path.name,
        "fixed_size": FIXED_SIZE,
        "random_seed": RANDOM_SEED,
        "z_repetitions": 2,
        "number_of_cells": 60,
        "number_of_atoms": structure.number_of_atoms,
        "topology": "20xZ12 + 32xZ14 + 8xZ15",
        "box_lengths": np.diag(structure.cell).tolist(),
        "twin_type": structure.metadata["twin_type"],
        "designed_twin_pair_count": len(structure.metadata["twin_pairs"]),
        "designed_twin_grain_count": 2 * len(structure.metadata["twin_pairs"]),
        "ordinary_grain_count": structure.metadata["ordinary_grain_count"],
        "design_rule": (
            "disjoint true-face pairs selected from the periodic 60-cell graph after "
            "excluding every modulo-30 pair used by the original Sigma library"
        ),
        "twin_pairs": structure.metadata["twin_pairs"],
        "cells": cells,
    }
    (output_directory / "cell_manifest.json").write_text(
        json.dumps(_jsonable(manifest), indent=2, ensure_ascii=False), encoding="utf-8"
    )

    fields = [
        "grain_id", "z_copy", "source_sigma_grain_id", "site_class",
        "wyckoff_label", "coordination_number", "atom_count", "center_fx",
        "center_fy", "center_fz", "boundary_role", "twin_partner_grain_id",
        "twin_placement",
    ]
    with (output_directory / "cell_manifest.csv").open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for cell_record in cells:
            center = cell_record["center_reduced"]
            writer.writerow(
                {
                    **{key: cell_record[key] for key in fields if not key.startswith("center_f")},
                    "center_fx": center[0],
                    "center_fy": center[1],
                    "center_fz": center[2],
                }
            )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--overlap-distance-factor", type=float, default=0.5)
    parser.add_argument("--raw-voronoi", action="store_true")
    parser.add_argument("--label-grains", action="store_true")
    parser.add_argument("--design-id", default="03_three_zone_balanced_15")
    parser.add_argument("--design-config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--twin-type", choices=sorted(TWIN_SPECS), default="sigma3")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output")
    options = parser.parse_args()
    structure = build_sigma_z2(
        lattice_constant=options.lattice_constant,
        overlap_distance_factor=options.overlap_distance_factor,
        equal_volume=not options.raw_voronoi,
        label_grains=options.label_grains,
        design_id=options.design_id,
        design_config=options.design_config,
        twin_type=options.twin_type,
    )
    manifest = write_outputs(structure, options.output_dir)
    print(
        f"Wrote {structure.number_of_atoms} atoms in 60 cells to "
        f"{(options.output_dir / manifest['output']).resolve()}"
    )
    print(
        f"box={manifest['box_lengths']} A; {manifest['twin_type']} pairs="
        f"{manifest['designed_twin_pair_count']}; ordinary grains="
        f"{manifest['ordinary_grain_count']}"
    )


if __name__ == "__main__":
    main()
