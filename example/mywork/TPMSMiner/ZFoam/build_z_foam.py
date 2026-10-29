"""Build the hexagonal hP7 Frank-Kasper Z-phase / HS-I-dual foam."""

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
from scipy.spatial import Voronoi, cKDTree


REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from crystalops import CrystalStructure, FCC
from crystalops.builders._lattice import candidate_lattice
from crystalops.builders.polycrystal import (
    _nearest_periodic_site_distance,
)
from crystalops.operations import select


FloatArray = NDArray[np.float64]
ROOT = Path(__file__).resolve().parent
FIXED_SIZE = 13.9
C_OVER_A = 5.390 / 5.433
RANDOM_SEED = 20260930

# hP7 Zr4Al3 prototype, P6/mmm (#191).  The ideal/experimental 2e internal
# parameter is z=1/4.  Ordering by FK coordination makes the topology explicit.
Z12_3F = np.asarray(
    [[0.5, 0.0, 0.0], [0.0, 0.5, 0.0], [0.5, 0.5, 0.0]], dtype=np.float64
)
Z15_2D = np.asarray(
    [[1.0 / 3.0, 2.0 / 3.0, 0.5], [2.0 / 3.0, 1.0 / 3.0, 0.5]],
    dtype=np.float64,
)
Z14_2E = np.asarray([[0.0, 0.0, 0.25], [0.0, 0.0, 0.75]], dtype=np.float64)
Z_CENTERS = np.concatenate((Z12_3F, Z15_2D, Z14_2E))
SITE_CLASSES = ("Z12",) * 3 + ("Z15",) * 2 + ("Z14",) * 2
WYCKOFF_LABELS = ("3f-Z12",) * 3 + ("2d-Z15",) * 2 + ("2e-Z14",) * 2
COORDINATION_NUMBERS = (12, 12, 12, 15, 15, 14, 14)
GROUP_IDS = np.asarray([0, 0, 0, 1, 1, 2, 2], dtype=np.int64)

# Gauge Z12=0. These dimensionless power weights use a=1 and c/a above and
# make all seven periodic Laguerre cells exactly det(cell)/7.
EQUAL_VOLUME_WEIGHTS = np.asarray(
    [0.0, -0.022579650470292777, -0.011872417238690782], dtype=np.float64
)
PERIODIC_SHIFTS = np.asarray(
    list(itertools.product((-1.0, 0.0, 1.0), repeat=3)), dtype=np.float64
)

# One true same-class face per class. One of the three Z12 cells remains an
# ordinary independently oriented grain because the primitive cell is odd.
TWIN_FACE_PAIRS: tuple[tuple[int, int, tuple[float, float, float]], ...] = (
    (0, 1, (0.5, 0.5, 0.0)),
    (3, 4, (-2.0 / 3.0, -1.0 / 3.0, 0.0)),
    (5, 6, (0.0, 0.0, -0.5)),
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


def _random_rotation(seed: int) -> FloatArray:
    quaternion = np.random.default_rng(seed).normal(size=4)
    quaternion /= np.linalg.norm(quaternion)
    w, x, y, z = quaternion
    return np.asarray(
        [
            [1 - 2 * (y*y + z*z), 2 * (x*y - z*w), 2 * (x*z + y*w)],
            [2 * (x*y + z*w), 1 - 2 * (x*x + z*z), 2 * (y*z - x*w)],
            [2 * (x*z - y*w), 2 * (y*z + x*w), 1 - 2 * (x*x + y*y)],
        ]
    )


def _hexagonal_cell(a_length: float) -> FloatArray:
    return np.asarray(
        [
            [a_length, 0.0, 0.0],
            [-0.5 * a_length, np.sqrt(3.0) * 0.5 * a_length, 0.0],
            [0.0, 0.0, C_OVER_A * a_length],
        ],
        dtype=np.float64,
    )


def z_coordination_numbers() -> tuple[int, ...]:
    """Recover the seven FK coordination numbers from periodic Voronoi faces."""
    cell = _hexagonal_cell(1.0)
    shifts = np.asarray(
        list(itertools.product(range(-2, 3), repeat=3)), dtype=np.float64
    )
    points: list[FloatArray] = []
    labels: list[tuple[int, tuple[int, int, int]]] = []
    for shift in shifts:
        for site_index, center in enumerate(Z_CENTERS):
            points.append((center + shift) @ cell)
            labels.append((site_index, tuple(shift.astype(int).tolist())))
    tessellation = Voronoi(np.asarray(points))
    neighbors: dict[int, set[tuple[int, tuple[int, int, int]]]] = {
        index: set() for index in range(7)
    }
    for (first, second), vertices in zip(
        tessellation.ridge_points, tessellation.ridge_vertices, strict=True
    ):
        if -1 in vertices:
            continue
        first_site, first_shift = labels[first]
        second_site, second_shift = labels[second]
        if first_shift == (0, 0, 0):
            neighbors[first_site].add((second_site, second_shift))
        if second_shift == (0, 0, 0):
            neighbors[second_site].add((first_site, first_shift))
    return tuple(len(neighbors[index]) for index in range(7))


def _remove_hexagonal_periodic_overlaps(
    structure: CrystalStructure, cutoff: float
) -> tuple[CrystalStructure, int]:
    """Remove overlaps under a general triclinic periodic cell."""
    inverse_cell = np.linalg.inv(structure.cell)
    fractional = (structure.positions - structure.origin) @ inverse_cell
    wrapped = np.mod(fractional, 1.0) @ structure.cell
    shifts = np.asarray(
        list(itertools.product((-1.0, 0.0, 1.0), repeat=3)), dtype=np.float64
    )
    replicas = np.concatenate([wrapped + shift @ structure.cell for shift in shifts])
    tree = cKDTree(replicas)
    neighbor_lists = tree.query_ball_point(wrapped, np.nextafter(cutoff, 0.0))
    atom_count = structure.number_of_atoms
    retained = np.ones(atom_count, dtype=bool)
    for atom_index, neighbors in enumerate(neighbor_lists):
        periodic_indices = {int(index % atom_count) for index in neighbors}
        if any(
            other < atom_index and retained[other]
            for other in periodic_indices
            if other != atom_index
        ):
            retained[atom_index] = False
    selected = np.flatnonzero(retained)
    return select(structure, selected), atom_count - len(selected)


def _orientations_and_anchors(
    cell: FloatArray,
    pair_types: dict[str, str],
) -> tuple[FloatArray, FloatArray, FloatArray, list[dict[str, Any]]]:
    orientations = np.asarray([_random_rotation(RANDOM_SEED + i) for i in range(7)])
    anchors = Z_CENTERS @ cell
    construction_centers = Z_CENTERS.copy()
    records: list[dict[str, Any]] = []
    for pair_index, (first, second, vector_values) in enumerate(TWIN_FACE_PAIRS):
        site_class = SITE_CLASSES[first]
        vector = np.asarray(vector_values, dtype=np.float64)
        physical_vector = vector @ cell
        normal = physical_vector / np.linalg.norm(physical_vector)
        twin_type = pair_types[site_class]
        crystal_axis, angle, sigma = TWIN_SPECS[twin_type]
        base = _align_vector(crystal_axis, normal)
        sign = 1.0 if pair_index % 2 == 0 else -1.0
        twin = _axis_angle(normal, sign * angle) @ base
        orientations[first] = base
        orientations[second] = twin
        anchor_reduced = Z_CENTERS[first] + 0.5 * vector
        anchors[first] = anchor_reduced @ cell
        anchors[second] = anchor_reduced @ cell
        construction_centers[second] = Z_CENTERS[first] + vector
        records.append(
            {
                "pair": [first + 1, second + 1],
                "site_class": site_class,
                "twin_type": twin_type,
                "sigma": sigma,
                "angle_degrees": angle,
                "periodic_face_vector_reduced": vector.tolist(),
                "face_normal": normal.tolist(),
                "anchor_reduced_unwrapped": anchor_reduced.tolist(),
                "coherent_phase_anchor": True,
            }
        )
    return orientations, anchors, construction_centers, records


def build_z_foam(
    *,
    lattice_constant: float = 3.61,
    overlap_distance_factor: float = 0.5,
    equal_volume: bool = True,
    label_grains: bool = False,
    z12_twin_type: str = "sigma3",
    z14_twin_type: str = "sigma3",
    z15_twin_type: str = "sigma3",
    name: str = "z-foam-all-sigma3",
) -> CrystalStructure:
    """Build one fixed-size seven-cell FK Z-phase atomistic precursor."""
    if lattice_constant <= 0.0 or overlap_distance_factor <= 0.0:
        raise ValueError("lattice constant and overlap factor must be positive")
    pair_types = {
        "Z12": z12_twin_type.lower(), "Z14": z14_twin_type.lower(),
        "Z15": z15_twin_type.lower(),
    }
    unknown = set(pair_types.values()) - set(TWIN_SPECS)
    if unknown:
        raise ValueError(f"unsupported twin types: {sorted(unknown)}")
    coordination = z_coordination_numbers()
    if coordination != COORDINATION_NUMBERS:
        raise RuntimeError(f"invalid Z-phase topology: {coordination}")

    mean_grain_length = FIXED_SIZE * lattice_constant
    dimensionless_volume = np.sqrt(3.0) * 0.5 * C_OVER_A
    a_length = np.cbrt(7.0 * mean_grain_length**3 / dimensionless_volume)
    cell = _hexagonal_cell(a_length)
    orientations, anchors, construction_centers, twin_pairs = (
        _orientations_and_anchors(cell, pair_types)
    )
    center_positions = construction_centers @ cell
    reduced_weights = (
        EQUAL_VOLUME_WEIGHTS[GROUP_IDS] if equal_volume else np.zeros(7)
    )
    weights = reduced_weights * a_length**2
    unit_cell = FCC("Cu")
    tolerance = max(a_length, 1.0) * 1e-10

    position_groups: list[FloatArray] = []
    type_groups: list[NDArray[np.int64]] = []
    id_groups: list[NDArray[np.int64]] = []
    for grain_index, (center, orientation, anchor) in enumerate(
        zip(center_positions, orientations, anchors, strict=True)
    ):
        normals: list[FloatArray] = []
        limits: list[float] = []
        for neighbor_index, neighbor in enumerate(Z_CENTERS):
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
        inside = np.all(
            relative_candidates @ np.asarray(normals).T
            <= np.asarray(limits)[None, :] + tolerance,
            axis=1,
        )
        positions = candidates[inside]
        types = np.tile(unit_cell.atom_types, len(candidates) // unit_cell.number_of_sites)[inside]
        position_groups.append(positions)
        type_groups.append(types)
        id_groups.append(np.full(len(positions), grain_index + 1, dtype=np.int64))

    structure = CrystalStructure(
        positions=np.concatenate(position_groups),
        cell=cell,
        atom_types=np.concatenate(type_groups),
        type_elements=unit_cell.type_elements,
        origin=np.zeros(3),
        pbc=(True, True, True),
        properties={"grain_ids": np.concatenate(id_groups)},
        metadata={
            "builder": "TPMSMiner.ZFoamTwinPrecursor",
            "experimental": True,
            "precursor": "Frank-Kasper-Z-phase-HS-I-dual",
            "prototype": "hP7-Zr4Al3",
            "space_group": "P6/mmm",
            "grain_topology": "3xZ12 + 2xZ14 + 2xZ15",
            "c_over_a": C_OVER_A,
            "equal_volume": bool(equal_volume),
            "grain_weights_reduced": reduced_weights.copy(),
            "grain_centers_reduced": Z_CENTERS.copy(),
            "construction_centers_reduced": construction_centers.copy(),
            "grain_centers": center_positions.copy(),
            "grain_orientations": orientations.copy(),
            "lattice_anchors": anchors.copy(),
            "coordination_numbers": coordination,
            "site_classes": SITE_CLASSES,
            "wyckoff_labels": WYCKOFF_LABELS,
            "twin_pairs": twin_pairs,
            "pair_types_by_class": pair_types,
            "sigma3_pair_count": sum(item["sigma"] == 3 for item in twin_pairs),
            "sigma5_pair_count": sum(item["sigma"] == 5 for item in twin_pairs),
            "ordinary_grain_ids": (3,),
            "fixed_size": FIXED_SIZE,
            "mean_grain_length": float(mean_grain_length),
            "box_vectors": cell.copy(),
            "lattice_constant": float(lattice_constant),
            "cropped_to_box": False,
        },
        name=name,
    )
    initial_count = structure.number_of_atoms
    nearest_distance = _nearest_periodic_site_distance(unit_cell, lattice_constant)
    overlap_cutoff = overlap_distance_factor * nearest_distance
    structure, removed = _remove_hexagonal_periodic_overlaps(structure, overlap_cutoff)
    ids, counts = np.unique(structure.properties["grain_ids"], return_counts=True)
    if not np.array_equal(ids, np.arange(1, 8)):
        raise RuntimeError("Z-foam construction produced empty grains")
    structure.metadata.update(
        {
            "initial_atom_count": int(initial_count),
            "overlap_atoms_removed": int(removed),
            "overlap_cutoff": float(overlap_cutoff),
            "nearest_neighbor_distance": float(nearest_distance),
            "grain_counts": dict(zip(ids.tolist(), counts.tolist(), strict=True)),
        }
    )
    if label_grains:
        structure.atom_types = structure.properties["grain_ids"].copy()
        structure.type_elements = {grain_id: "Cu" for grain_id in range(1, 8)}
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


def write_z_outputs(structure: CrystalStructure, output_directory: Path) -> dict[str, Any]:
    """Write full, representative and class-only configurations."""
    output_directory.mkdir(parents=True, exist_ok=True)
    grain_ids = np.asarray(structure.properties["grain_ids"])

    def single_type(candidate: CrystalStructure, name: str) -> CrystalStructure:
        candidate.name = name
        candidate.atom_types = np.ones(candidate.number_of_atoms, dtype=np.int64)
        candidate.type_elements = {1: "Cu"}
        candidate.validate()
        return candidate

    files: dict[str, str] = {}
    full_path = output_directory / "01_z_foam_all.data"
    structure.write(full_path, format="lammps-data", atom_style="atomic")
    files["z_foam_all"] = full_path.name

    representatives = {"Z12": 0, "Z14": 5, "Z15": 3}
    for output_index, (site_class, grain_index) in enumerate(representatives.items(), start=2):
        candidate = single_type(
            select(structure, grain_ids == grain_index + 1),
            f"representative-{site_class.lower()}-cell-at-origin",
        )
        candidate.positions -= np.asarray(structure.metadata["grain_centers"])[grain_index]
        candidate.origin = -0.5 * np.diag(candidate.cell)
        path = output_directory / f"{output_index:02d}_{site_class.lower()}_cell_origin.data"
        candidate.write(path, format="lammps-data", atom_style="atomic")
        files[f"{site_class.lower()}_representative_at_origin"] = path.name

    for output_index, site_class in enumerate(("Z12", "Z14", "Z15"), start=5):
        allowed = np.flatnonzero(np.asarray(SITE_CLASSES) == site_class) + 1
        candidate = single_type(
            select(structure, np.isin(grain_ids, allowed)), f"z-foam-{site_class.lower()}-only"
        )
        path = output_directory / f"{output_index:02d}_z_foam_{site_class.lower()}_only.data"
        candidate.write(path, format="lammps-data", atom_style="atomic")
        files[f"{site_class.lower()}_only"] = path.name

    twin_by_grain: dict[int, dict[str, Any]] = {}
    for twin in structure.metadata["twin_pairs"]:
        first, second = twin["pair"]
        twin_by_grain[first] = {**twin, "partner": second}
        twin_by_grain[second] = {**twin, "partner": first}
    cells = []
    for index in range(7):
        twin = twin_by_grain.get(index + 1)
        cells.append(
            {
                "grain_id": index + 1,
                "site_class": SITE_CLASSES[index],
                "wyckoff_label": WYCKOFF_LABELS[index],
                "coordination_number": COORDINATION_NUMBERS[index],
                "atom_count": int(np.count_nonzero(grain_ids == index + 1)),
                "center_reduced": Z_CENTERS[index].tolist(),
                "construction_center_reduced_unwrapped": np.asarray(
                    structure.metadata["construction_centers_reduced"]
                )[index].tolist(),
                "center_cartesian_unwrapped": np.asarray(
                    structure.metadata["grain_centers"]
                )[index].tolist(),
                "orientation_matrix_crystal_to_lab": np.asarray(
                    structure.metadata["grain_orientations"]
                )[index].tolist(),
                "boundary_role": "designed_csl_pair" if twin else "ordinary_grain",
                "twin_partner_grain_id": twin["partner"] if twin else None,
                "sigma": twin["sigma"] if twin else None,
                "misorientation_degrees": twin["angle_degrees"] if twin else None,
                "shared_anchor_reduced": twin["anchor_reduced_unwrapped"] if twin else None,
            }
        )
    manifest = {
        "schema_version": 1,
        "topology": "3xZ12 + 2xZ14 + 2xZ15",
        "prototype": "hP7-Zr4Al3",
        "space_group": "P6/mmm",
        "fixed_size": FIXED_SIZE,
        "number_of_atoms": structure.number_of_atoms,
        "box_vectors": structure.cell.tolist(),
        "configuration_files": files,
        "cells": cells,
    }
    (output_directory / "cell_manifest.json").write_text(
        json.dumps(_jsonable(manifest), indent=2, ensure_ascii=False), encoding="utf-8"
    )
    fields = [
        "grain_id", "site_class", "wyckoff_label", "coordination_number",
        "atom_count", "center_fx", "center_fy", "center_fz", "boundary_role",
        "twin_partner_grain_id", "sigma", "misorientation_degrees",
    ]
    with (output_directory / "cell_manifest.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for record in cells:
            center = record["center_reduced"]
            row = {key: record[key] for key in fields if not key.startswith("center_f")}
            row.update(center_fx=center[0], center_fy=center[1], center_fz=center[2])
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
    structure = build_z_foam(
        lattice_constant=options.lattice_constant,
        overlap_distance_factor=options.overlap_distance_factor,
        equal_volume=not options.raw_voronoi,
        label_grains=options.label_grains,
    )
    manifest = write_z_outputs(structure, options.output_dir)
    print(
        f"Wrote {structure.number_of_atoms} atoms in 7 Z-phase cells to "
        f"{(options.output_dir / manifest['configuration_files']['z_foam_all']).resolve()}"
    )
    print(
        f"Topology: 3xZ12 + 2xZ14 + 2xZ15; box vectors="
        f"{np.asarray(structure.cell).tolist()}; Sigma3 pairs=3; ordinary grains=1"
    )


if __name__ == "__main__":
    main()
