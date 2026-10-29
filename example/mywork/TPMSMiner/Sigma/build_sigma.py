"""Build the tetragonal 30-cell Frank-Kasper sigma foam precursor."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np
from numpy.typing import NDArray


REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from crystalops import CrystalStructure, FCC
from crystalops.builders._lattice import candidate_lattice
from crystalops.builders.polycrystal import (
    _nearest_periodic_site_distance,
    _remove_periodic_overlaps,
)
from crystalops.operations import select


FloatArray = NDArray[np.float64]
ROOT = Path(__file__).resolve().parent
FIXED_SIZE = 8.6
C_OVER_A = 131.0 / 253.0

Z12_8I = np.asarray(
    [
        [0.071, 0.257, 0.5], [0.243, 0.571, 0.0],
        [0.257, 0.071, 0.5], [0.429, 0.757, 0.0],
        [0.571, 0.243, 0.0], [0.743, 0.929, 0.5],
        [0.757, 0.429, 0.0], [0.929, 0.743, 0.5],
    ],
    dtype=np.float64,
)
Z12_2B = np.asarray([[0.0, 0.0, 0.5], [0.5, 0.5, 0.0]], dtype=np.float64)
Z14_8I = np.asarray(
    [
        [0.031, 0.371, 0.0], [0.129, 0.531, 0.5],
        [0.371, 0.031, 0.0], [0.469, 0.871, 0.5],
        [0.531, 0.129, 0.5], [0.629, 0.969, 0.0],
        [0.871, 0.469, 0.5], [0.969, 0.629, 0.0],
    ],
    dtype=np.float64,
)
Z14_8J = np.asarray(
    [
        [0.181, 0.819, 0.25], [0.181, 0.819, 0.75],
        [0.319, 0.319, 0.25], [0.319, 0.319, 0.75],
        [0.681, 0.681, 0.25], [0.681, 0.681, 0.75],
        [0.819, 0.181, 0.25], [0.819, 0.181, 0.75],
    ],
    dtype=np.float64,
)
Z15_4F = np.asarray(
    [[0.104, 0.104, 0.0], [0.396, 0.604, 0.5],
     [0.604, 0.396, 0.5], [0.896, 0.896, 0.0]],
    dtype=np.float64,
)
SIGMA_CENTERS = np.concatenate((Z12_8I, Z12_2B, Z14_8I, Z14_8J, Z15_4F))
COORDINATION_NUMBERS = (12,) * 10 + (14,) * 16 + (15,) * 4
SITE_CLASSES = ("Z12",) * 10 + ("Z14",) * 16 + ("Z15",) * 4
WYCKOFF_LABELS = (
    ("8i-Z12",) * 8 + ("2b-Z12",) * 2
    + ("8i-Z14",) * 8 + ("8j-Z14",) * 8 + ("4f-Z15",) * 4
)
GROUP_IDS = np.repeat(np.arange(5), (8, 2, 8, 8, 4))

# Gauge: Z12-8i weight is zero.  These dimensionless weights use a=1 and
# c=131/253 and make all five Wyckoff-cell volumes exactly det(cell)/30.
EQUAL_VOLUME_WEIGHTS = np.asarray(
    [
        0.0,
        0.0009810127647068775,
        -0.002252513628455159,
        -0.0014524364976177743,
        -0.003137917345335853,
    ],
    dtype=np.float64,
)

PERIODIC_SHIFTS = np.asarray(
    list(itertools.product((-1.0, 0.0, 1.0), repeat=3)), dtype=np.float64
)

# Perfect matching over true periodic Voronoi faces. Indices are zero-based;
# vectors are fractional displacements to the chosen periodic image.
TWIN_FACE_PAIRS: tuple[tuple[int, int, tuple[float, float, float]], ...] = (
    (0, 2, (0.186, -0.186, 0.0)),
    (1, 3, (0.186, 0.186, 0.0)),
    (4, 6, (0.186, 0.186, 0.0)),
    (5, 8, (0.257, 0.071, 0.0)),
    (7, 22, (-0.248, -0.062, -0.25)),
    (9, 23, (0.181, 0.181, -0.25)),
    (10, 17, (-0.062, 0.258, 0.0)),
    (11, 16, (-0.258, -0.062, 0.0)),
    (12, 13, (0.098, -0.160, 0.5)),
    (14, 15, (0.098, -0.160, -0.5)),
    (18, 19, (0.0, 0.0, -0.5)),
    (20, 21, (0.0, 0.0, -0.5)),
    (24, 25, (0.0, 0.0, -0.5)),
    (26, 29, (-0.208, -0.208, 0.0)),
    (27, 28, (0.208, -0.208, 0.0)),
)

TWIN_SPECS: dict[str, tuple[FloatArray, float, int]] = {
    "sigma3": (np.asarray([1.0, 1.0, 1.0]), 60.0, 3),
    "sigma5": (np.asarray([1.0, 0.0, 0.0]), 36.86989764584402, 5),
}


def _axis_angle(axis: FloatArray, degrees: float) -> FloatArray:
    axis = np.asarray(axis, dtype=np.float64)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    skew = np.asarray([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    angle = np.radians(degrees)
    return np.eye(3) + np.sin(angle) * skew + (1.0 - np.cos(angle)) * (skew @ skew)


def _align_vector(source: FloatArray, target: FloatArray) -> FloatArray:
    source = np.asarray(source, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    source /= np.linalg.norm(source)
    target /= np.linalg.norm(target)
    cross = np.cross(source, target)
    sine = np.linalg.norm(cross)
    cosine = float(np.clip(source @ target, -1.0, 1.0))
    if sine < 1e-12:
        if cosine > 0.0:
            return np.eye(3)
        trial = np.asarray([1.0, 0.0, 0.0])
        if abs(source @ trial) > 0.9:
            trial = np.asarray([0.0, 1.0, 0.0])
        return _axis_angle(np.cross(source, trial), 180.0)
    x, y, z = cross / sine
    skew = np.asarray([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) + sine * skew + (1.0 - cosine) * (skew @ skew)


def _pair_group(first: int, second: int) -> str:
    classes = sorted((SITE_CLASSES[first], SITE_CLASSES[second]))
    return f"{classes[0]}-{classes[1]}"


def _orientations_and_anchors(
    cell: FloatArray,
    pair_types: dict[str, str],
    swap_groups: frozenset[str],
) -> tuple[FloatArray, FloatArray, FloatArray, list[dict[str, Any]]]:
    orientations = np.empty((30, 3, 3), dtype=np.float64)
    anchors = np.empty((30, 3), dtype=np.float64)
    construction_centers = SIGMA_CENTERS.copy()
    records: list[dict[str, Any]] = []
    for pair_index, (first, second, vector_values) in enumerate(TWIN_FACE_PAIRS):
        vector = np.asarray(vector_values, dtype=np.float64)
        physical_vector = vector @ cell
        normal = physical_vector / np.linalg.norm(physical_vector)
        pair_group = _pair_group(first, second)
        twin_type = pair_types.get(pair_group, "sigma3")
        crystal_axis, angle, sigma = TWIN_SPECS[twin_type]
        base = _align_vector(crystal_axis, normal)
        sign = 1.0 if pair_index % 2 == 0 else -1.0
        twin = _axis_angle(normal, sign * angle) @ base
        swapped = pair_group in swap_groups
        orientations[first] = twin if swapped else base
        orientations[second] = base if swapped else twin
        anchor_reduced = SIGMA_CENTERS[first] + 0.5 * vector
        anchor = anchor_reduced @ cell
        anchors[first] = anchor
        anchors[second] = anchor
        construction_centers[first] = SIGMA_CENTERS[first]
        construction_centers[second] = SIGMA_CENTERS[first] + vector
        records.append(
            {
                "pair": [first + 1, second + 1],
                "pair_group": pair_group,
                "site_classes": [SITE_CLASSES[first], SITE_CLASSES[second]],
                "wyckoff_labels": [WYCKOFF_LABELS[first], WYCKOFF_LABELS[second]],
                "twin_type": twin_type,
                "sigma": sigma,
                "angle_degrees": angle,
                "face_normal": normal.tolist(),
                "anchor_reduced": anchor_reduced.tolist(),
                "swapped": swapped,
            }
        )
    return orientations, anchors, construction_centers, records


def build_sigma_foam(
    *,
    size: float = FIXED_SIZE,
    lattice_constant: float = 3.61,
    overlap_distance_factor: float = 0.5,
    equal_volume: bool = True,
    label_grains: bool = False,
    pair_types: dict[str, str] | None = None,
    swap_groups: tuple[str, ...] = (),
    name: str = "sigma-foam-sigma3-reference",
) -> CrystalStructure:
    """Build a fixed-topology tetragonal sigma-foam atomistic precursor."""
    if size <= 0.0 or lattice_constant <= 0.0:
        raise ValueError("size and lattice_constant must be positive")
    if overlap_distance_factor <= 0.0:
        raise ValueError("overlap_distance_factor must be positive")
    pair_type_map = {} if pair_types is None else dict(pair_types)
    unknown_types = set(pair_type_map.values()) - set(TWIN_SPECS)
    if unknown_types:
        raise ValueError(f"unsupported twin types: {sorted(unknown_types)}")
    valid_groups = {_pair_group(first, second) for first, second, _ in TWIN_FACE_PAIRS}
    if set(pair_type_map) - valid_groups:
        raise ValueError(f"unsupported pair groups: {sorted(set(pair_type_map)-valid_groups)}")
    swap_set = frozenset(swap_groups)
    if swap_set - valid_groups:
        raise ValueError(f"unsupported swap groups: {sorted(swap_set-valid_groups)}")

    mean_grain_length = size * lattice_constant
    a_length = np.cbrt(30.0 / C_OVER_A) * mean_grain_length
    c_length = C_OVER_A * a_length
    cell = np.diag([a_length, a_length, c_length])
    orientations, anchors, construction_centers, twin_pairs = (
        _orientations_and_anchors(cell, pair_type_map, swap_set)
    )
    center_positions = construction_centers @ cell
    reduced_weights = EQUAL_VOLUME_WEIGHTS[GROUP_IDS] if equal_volume else np.zeros(30)
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
        for neighbor_index, neighbor in enumerate(SIGMA_CENTERS):
            for shift in PERIODIC_SHIFTS:
                if neighbor_index == grain_index and np.all(shift == 0.0):
                    continue
                relative = (
                    neighbor + shift - construction_centers[grain_index]
                ) @ cell
                normals.append(2.0 * relative)
                limits.append(
                    float(
                        relative @ relative
                        + weights[grain_index]
                        - weights[neighbor_index]
                    )
                )
        candidates = candidate_lattice(
            unit_cell, lattice_constant, cell, orientation
        ) + anchor
        relative_candidates = candidates - center
        inside = np.all(
            relative_candidates @ np.asarray(normals).T
            <= np.asarray(limits)[None, :] + tolerance,
            axis=1,
        )
        positions = candidates[inside]
        types = np.tile(
            unit_cell.atom_types,
            len(candidates) // unit_cell.number_of_sites,
        )[inside]
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
            "builder": "TPMSMiner.SigmaFoamTwinPrecursor",
            "experimental": True,
            "precursor": "Frank-Kasper-sigma-foam",
            "space_group": "P42/mnm",
            "grain_topology": "10xZ12 + 16xZ14 + 4xZ15",
            "c_over_a": C_OVER_A,
            "equal_volume": bool(equal_volume),
            "grain_weights_reduced": reduced_weights.copy(),
            "grain_weights": weights.copy(),
            "grain_centers_reduced": SIGMA_CENTERS.copy(),
            "construction_centers_reduced": construction_centers.copy(),
            "grain_centers": center_positions.copy(),
            "grain_orientations": orientations.copy(),
            "lattice_anchors": anchors.copy(),
            "coordination_numbers": COORDINATION_NUMBERS,
            "site_classes": SITE_CLASSES,
            "wyckoff_labels": WYCKOFF_LABELS,
            "twin_pairs": twin_pairs,
            "sigma3_pair_count": sum(item["sigma"] == 3 for item in twin_pairs),
            "sigma5_pair_count": sum(item["sigma"] == 5 for item in twin_pairs),
            "pair_types": pair_type_map,
            "swap_groups": tuple(sorted(swap_set)),
            "fixed_size": float(size),
            "mean_grain_length": float(mean_grain_length),
            "box_lengths": np.diag(cell).copy(),
            "lattice_constant": float(lattice_constant),
            "cropped_to_box": False,
        },
        name=name,
    )
    initial_count = structure.number_of_atoms
    nearest_distance = _nearest_periodic_site_distance(unit_cell, lattice_constant)
    overlap_cutoff = overlap_distance_factor * nearest_distance
    structure, removed = _remove_periodic_overlaps(structure, overlap_cutoff)
    ids, counts = np.unique(structure.properties["grain_ids"], return_counts=True)
    if not np.array_equal(ids, np.arange(1, 31)):
        raise RuntimeError("sigma-foam construction produced empty grains")
    lengths = np.diag(cell)
    outside = np.any(
        (structure.positions < -tolerance)
        | (structure.positions > lengths[None, :] + tolerance),
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
        structure.type_elements = {grain_id: "Cu" for grain_id in range(1, 31)}
        structure.validate()
    return structure


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_jsonable(item) for item in value]
    return value


def write_sigma_outputs(structure: CrystalStructure, output_directory: Path) -> dict[str, Any]:
    """Write full, representative, class-only structures and manifests."""
    output_directory.mkdir(parents=True, exist_ok=True)
    grain_ids = np.asarray(structure.properties["grain_ids"])

    def one_type(candidate: CrystalStructure, name: str) -> CrystalStructure:
        candidate.name = name
        candidate.atom_types = np.ones(candidate.number_of_atoms, dtype=np.int64)
        candidate.type_elements = {1: "Cu"}
        candidate.validate()
        return candidate

    configuration_files: dict[str, str] = {}
    full_path = output_directory / "01_sigma_all.data"
    structure.write(full_path, format="lammps-data", atom_style="atomic")
    configuration_files["sigma_all"] = full_path.name

    representatives = {"Z12": 0, "Z14": 10, "Z15": 26}
    for output_index, (site_class, grain_index) in enumerate(representatives.items(), start=2):
        representative = one_type(
            select(structure, grain_ids == grain_index + 1),
            f"representative-{site_class.lower()}-cell-at-origin",
        )
        representative.positions -= np.asarray(structure.metadata["grain_centers"])[grain_index]
        representative.origin = -0.5 * np.diag(representative.cell)
        representative.validate()
        path = output_directory / f"{output_index:02d}_{site_class.lower()}_cell_origin.data"
        representative.write(path, format="lammps-data", atom_style="atomic")
        configuration_files[f"{site_class.lower()}_representative_at_origin"] = path.name

    for output_index, site_class in enumerate(("Z12", "Z14", "Z15"), start=5):
        mask = np.asarray(SITE_CLASSES) == site_class
        allowed_ids = np.flatnonzero(mask) + 1
        candidate = one_type(
            select(structure, np.isin(grain_ids, allowed_ids)),
            f"sigma-{site_class.lower()}-only",
        )
        path = output_directory / f"{output_index:02d}_sigma_{site_class.lower()}_only.data"
        candidate.write(path, format="lammps-data", atom_style="atomic")
        configuration_files[f"{site_class.lower()}_only"] = path.name

    partner_records: dict[int, dict[str, Any]] = {}
    for pair in structure.metadata["twin_pairs"]:
        first, second = pair["pair"]
        partner_records[first] = {**pair, "partner": second}
        partner_records[second] = {
            **pair,
            "partner": first,
            "face_normal": (-np.asarray(pair["face_normal"])).tolist(),
        }

    class_indices = {"Z12": 0, "Z14": 0, "Z15": 0}
    cells: list[dict[str, Any]] = []
    for index in range(30):
        grain_id = index + 1
        site_class = SITE_CLASSES[index]
        class_indices[site_class] += 1
        twin = partner_records[grain_id]
        cells.append(
            {
                "grain_id": grain_id,
                "cell_name": f"{site_class}_{class_indices[site_class]:02d}_grain_{grain_id:02d}",
                "site_class": site_class,
                "wyckoff_label": WYCKOFF_LABELS[index],
                "coordination_number": COORDINATION_NUMBERS[index],
                "atom_count": int(np.count_nonzero(grain_ids == grain_id)),
                "center_reduced": SIGMA_CENTERS[index].tolist(),
                "construction_center_reduced_unwrapped": np.asarray(
                    structure.metadata["construction_centers_reduced"]
                )[index].tolist(),
                "center_cartesian_unwrapped": np.asarray(
                    structure.metadata["grain_centers"]
                )[index].tolist(),
                "orientation_matrix_crystal_to_lab": np.asarray(
                    structure.metadata["grain_orientations"]
                )[index].tolist(),
                "twin_partner_grain_id": int(twin["partner"]),
                "twin_pair_group": twin["pair_group"],
                "twin_face_normal_outward": twin["face_normal"],
                "shared_anchor_reduced": twin["anchor_reduced"],
                "sigma": int(twin["sigma"]),
                "misorientation_degrees": float(twin["angle_degrees"]),
            }
        )

    manifest = {
        "schema_version": 1,
        "fixed_size": FIXED_SIZE,
        "size_series_enabled": False,
        "topology": "10xZ12 + 16xZ14 + 4xZ15",
        "space_group": "P42/mnm",
        "c_over_a": C_OVER_A,
        "box_lengths": np.diag(structure.cell).tolist(),
        "configuration_files": configuration_files,
        "coordinate_convention": {
            "center_reduced": "fractional coordinate in the primary tetragonal box",
            "orientation_matrix_crystal_to_lab": "v_lab = O @ v_crystal",
            "representative_cells": "geometric cell center translated to (0,0,0)",
            "class_only_files": "retain global sigma-cell relative positions",
        },
        "cells": cells,
    }
    (output_directory / "cell_manifest.json").write_text(
        json.dumps(_jsonable(manifest), indent=2, ensure_ascii=False), encoding="utf-8"
    )

    csv_fields = [
        "grain_id", "cell_name", "site_class", "wyckoff_label",
        "coordination_number", "atom_count", "center_fx", "center_fy", "center_fz",
        "center_x_A", "center_y_A", "center_z_A", "twin_partner_grain_id",
        "twin_pair_group", "sigma", "misorientation_degrees",
    ] + [f"orientation_{row}{column}" for row in range(3) for column in range(3)]
    with (output_directory / "cell_manifest.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=csv_fields)
        writer.writeheader()
        for record in cells:
            center = record["center_reduced"]
            cartesian = record["center_cartesian_unwrapped"]
            row = {
                key: record[key]
                for key in (
                    "grain_id", "cell_name", "site_class", "wyckoff_label",
                    "coordination_number", "atom_count", "twin_partner_grain_id",
                    "twin_pair_group", "sigma", "misorientation_degrees",
                )
            }
            row.update(
                {
                    "center_fx": center[0], "center_fy": center[1], "center_fz": center[2],
                    "center_x_A": cartesian[0], "center_y_A": cartesian[1],
                    "center_z_A": cartesian[2],
                }
            )
            orientation = record["orientation_matrix_crystal_to_lab"]
            row.update(
                {
                    f"orientation_{r}{c}": orientation[r][c]
                    for r in range(3) for c in range(3)
                }
            )
            writer.writerow(row)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--overlap-distance-factor", type=float, default=0.5)
    parser.add_argument("--raw-voronoi", action="store_true")
    parser.add_argument("--label-grains", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output")
    options = parser.parse_args()
    structure = build_sigma_foam(
        size=FIXED_SIZE,
        lattice_constant=options.lattice_constant,
        overlap_distance_factor=options.overlap_distance_factor,
        equal_volume=not options.raw_voronoi,
        label_grains=options.label_grains,
    )
    manifest = write_sigma_outputs(structure, options.output_dir)
    print(
        f"Wrote {structure.number_of_atoms} atoms in 30 sigma-foam cells to "
        f"{(options.output_dir / manifest['configuration_files']['sigma_all']).resolve()}."
    )
    print(
        f"Topology: 10xZ12 + 16xZ14 + 4xZ15; box="
        f"{np.diag(structure.cell).tolist()} A; fixed size={FIXED_SIZE}."
    )
    print(f"Wrote seven configurations and manifests to {options.output_dir.resolve()}.")


if __name__ == "__main__":
    main()
