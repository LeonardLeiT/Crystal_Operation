"""Build an experimental C15-Laves Sigma3 twin precursor for LAMMPS.

The geometry is the cubic C15 (MgCu2 / clathrate-II-dual) point topology:
eight Z16 sites and sixteen Z12 sites per conventional periodic cell.  This
is a different periodic tessellation from both BCC/Kelvin and A15/W-P.

Each site is converted to its complete periodic Voronoi cell.  Sites of the
same C15 class are paired across a real Voronoi face.  The two FCC crystals
in every pair share a lattice-site anchor at the face midpoint and differ by
60 degrees about the face normal, producing a geometrically registered
Sigma3 {111} candidate boundary.  The cells are not cropped or wrapped into
the primary box before writing.
"""

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
FIXED_SIZE = 9.25
RANDOM_SEED = 20260930
ORDINARY_CSL_EXCLUSION_DEGREES = 1.0

# C15 8a sites: Z16 Friauf polyhedra (large-site sublattice).
Z16_CENTERS = np.asarray(
    [
        [1 / 8, 1 / 8, 1 / 8],
        [1 / 8, 5 / 8, 5 / 8],
        [5 / 8, 1 / 8, 5 / 8],
        [5 / 8, 5 / 8, 1 / 8],
        [7 / 8, 7 / 8, 7 / 8],
        [7 / 8, 3 / 8, 3 / 8],
        [3 / 8, 7 / 8, 3 / 8],
        [3 / 8, 3 / 8, 7 / 8],
    ],
    dtype=np.float64,
)

# C15 16d sites: pyrochlore sublattice, Z12 polyhedra.
_FCC_TRANSLATIONS = np.asarray(
    [[0, 0, 0], [0, 1 / 2, 1 / 2], [1 / 2, 0, 1 / 2], [1 / 2, 1 / 2, 0]],
    dtype=np.float64,
)
_PYROCHLORE_BASIS = np.asarray(
    [[1 / 2, 1 / 2, 1 / 2], [1 / 2, 1 / 4, 1 / 4],
     [1 / 4, 1 / 2, 1 / 4], [1 / 4, 1 / 4, 1 / 2]],
    dtype=np.float64,
)
Z12_CENTERS = np.unique(
    np.round(
        np.mod(
            (
                _FCC_TRANSLATIONS[:, None, :]
                + _PYROCHLORE_BASIS[None, :, :]
            ).reshape(-1, 3),
            1.0,
        ),
        12,
    ),
    axis=0,
)
C15_CENTERS = np.concatenate((Z16_CENTERS, Z12_CENTERS), axis=0)
SITE_CLASSES = ("Z16",) * 8 + ("Z12",) * 16

# Dimensionless Z16-minus-Z12 power weight.  With Z12 fixed to the zero
# gauge, this value makes every periodic C15 Laguerre cell exactly 1/24 of
# the cubic box volume.  It was solved from the half-space polyhedron volume.
Z16_EQUAL_VOLUME_WEIGHT = -0.009748098845136293

TWIN_SPECS: dict[str, tuple[FloatArray, float, int]] = {
    "sigma3": (np.asarray([1.0, 1.0, 1.0]), 60.0, 3),
    "sigma5": (np.asarray([1.0, 0.0, 0.0]), 36.86989764584402, 5),
}
CUBIC_SYMMETRIES = np.asarray(
    [
        operation
        for permutation in itertools.permutations(range(3))
        for signs in itertools.product((-1.0, 1.0), repeat=3)
        for operation in [
            np.eye(3, dtype=np.float64)[list(permutation)]
            * np.asarray(signs, dtype=np.float64)[:, None]
        ]
        if np.linalg.det(operation) > 0.0
    ]
)

PERIODIC_SHIFTS = np.asarray(
    list(itertools.product((-1.0, 0.0, 1.0), repeat=3)),
    dtype=np.float64,
)


def _axis_angle(axis: FloatArray, degrees: float) -> FloatArray:
    axis = np.asarray(axis, dtype=np.float64)
    axis /= np.linalg.norm(axis)
    x, y, z = axis
    skew = np.asarray([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    angle = np.radians(degrees)
    return np.eye(3) + np.sin(angle) * skew + (1.0 - np.cos(angle)) * (skew @ skew)


def _align_vector(source: FloatArray, target: FloatArray) -> FloatArray:
    """Return a proper rotation mapping ``source`` onto ``target``."""
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


def _cubic_disorientation(first: FloatArray, second: FloatArray) -> float:
    """Return the minimum misorientation angle under cubic symmetry."""
    relative = np.einsum(
        "ij,sjk,kl->sil", second, CUBIC_SYMMETRIES, first.T
    )
    cosines = np.clip(
        (np.trace(relative, axis1=1, axis2=2) - 1.0) / 2.0, -1.0, 1.0
    )
    return float(np.min(np.degrees(np.arccos(cosines))))


def _random_non_csl_rotation(
    rng: np.random.Generator,
    existing: list[FloatArray],
) -> FloatArray:
    """Sample a random ordinary orientation away from near-Sigma3/5 angles."""
    targets = tuple(specification[1] for specification in TWIN_SPECS.values())
    for _ in range(20000):
        candidate = _random_rotation(rng)
        angles = [_cubic_disorientation(candidate, item) for item in existing]
        if all(
            abs(angle - target) >= ORDINARY_CSL_EXCLUSION_DEGREES
            for angle in angles
            for target in targets
        ):
            return candidate
    raise RuntimeError("could not sample a non-CSL C15 ordinary orientation")


def _periodic_vector(first: FloatArray, second: FloatArray) -> FloatArray:
    vector = np.asarray(second) - np.asarray(first)
    vector -= np.round(vector)
    return vector


def _random_rotation(rng: np.random.Generator) -> FloatArray:
    quaternion = rng.normal(size=4)
    quaternion /= np.linalg.norm(quaternion)
    w, x, y, z = quaternion
    return np.asarray(
        [
            [1 - 2 * (y*y + z*z), 2 * (x*y - z*w), 2 * (x*z + y*w)],
            [2 * (x*y + z*w), 1 - 2 * (x*x + z*z), 2 * (y*z - x*w)],
            [2 * (x*z - y*w), 2 * (y*z + x*w), 1 - 2 * (x*x + y*y)],
        ],
        dtype=np.float64,
    )


def c15_coordination_numbers() -> tuple[int, ...]:
    """Compute C15 Z16/Z12 counts from the exact periodic distance shells."""
    counts: list[int] = []
    z16_cutoff = np.sqrt(3.0) / 4.0 + 1e-10
    z12_cutoff = np.sqrt(11.0) / 8.0 + 1e-10
    for index, center in enumerate(C15_CENTERS):
        cutoff = z16_cutoff if index < 8 else z12_cutoff
        count = 0
        for neighbor_index, neighbor in enumerate(C15_CENTERS):
            for shift in PERIODIC_SHIFTS:
                if neighbor_index == index and np.all(shift == 0.0):
                    continue
                distance = np.linalg.norm(neighbor + shift - center)
                if distance <= cutoff:
                    count += 1
        counts.append(count)
    return tuple(counts)


def _same_class_perfect_matching() -> list[tuple[int, int, FloatArray]]:
    """Return the most directionally balanced same-class perfect matching."""
    pairs: list[tuple[int, int, FloatArray]] = []
    for group, neighbor_distance in (
        (tuple(range(0, 8)), np.sqrt(3.0) / 4.0),
        (tuple(range(8, 24)), np.sqrt(2.0) / 4.0),
    ):
        edges: dict[int, list[tuple[int, FloatArray]]] = {
            index: [] for index in group
        }
        for first_position, first in enumerate(group):
            for second in group[first_position + 1 :]:
                vector = _periodic_vector(
                    C15_CENTERS[first], C15_CENTERS[second]
                )
                if np.isclose(np.linalg.norm(vector), neighbor_distance, atol=1e-12):
                    edges[first].append((second, vector))
                    edges[second].append((first, -vector))

        best_score: tuple[int, int] | None = None
        best_pairs: list[tuple[int, int, FloatArray]] | None = None

        def direction_key(vector: FloatArray) -> tuple[int, int, int]:
            scaled = vector / np.min(np.abs(vector[np.abs(vector) > 1e-12]))
            rounded = np.rint(scaled).astype(int)
            first_nonzero = rounded[np.flatnonzero(rounded)[0]]
            if first_nonzero < 0:
                rounded = -rounded
            return tuple(rounded.tolist())

        def search(
            unmatched: frozenset[int],
            selected: list[tuple[int, int, FloatArray]],
        ) -> None:
            nonlocal best_score, best_pairs
            if not unmatched:
                histogram: dict[tuple[int, int, int], int] = {}
                for _, _, vector in selected:
                    key = direction_key(vector)
                    histogram[key] = histogram.get(key, 0) + 1
                score = (len(histogram), -sum(value * value for value in histogram.values()))
                if best_score is None or score > best_score:
                    best_score = score
                    best_pairs = list(selected)
                return
            first = min(unmatched)
            for second, vector in edges[first]:
                if second not in unmatched:
                    continue
                search(
                    unmatched - {first, second},
                    selected + [(first, second, vector)],
                )

        search(frozenset(group), [])
        if best_pairs is None:
            raise RuntimeError(f"C15 sites {group} have no perfect face matching")
        pairs.extend(best_pairs)
    return pairs


def _single_layer_twin_pairs() -> list[tuple[int, int, FloatArray]]:
    """Keep only one spatial twin layer for each unoriented face normal."""
    selected: list[tuple[int, int, FloatArray]] = []
    directions: set[tuple[int, int, int, int]] = set()
    for first, second, vector in _same_class_perfect_matching():
        scaled = vector / np.min(np.abs(vector[np.abs(vector) > 1e-12]))
        rounded = np.rint(scaled).astype(int)
        first_nonzero = rounded[np.flatnonzero(rounded)[0]]
        if first_nonzero < 0:
            rounded = -rounded
        direction = tuple(rounded.tolist())
        class_direction = (0 if SITE_CLASSES[first] == "Z16" else 1, *direction)
        if class_direction in directions:
            continue
        directions.add(class_direction)
        selected.append((first, second, vector))
    return selected


def _orientations_and_anchors(
    box_length: float,
    pair_types: dict[str, str],
    swap_classes: frozenset[str],
) -> tuple[FloatArray, FloatArray, FloatArray, list[dict[str, Any]]]:
    orientations = np.empty((24, 3, 3), dtype=np.float64)
    anchors = C15_CENTERS * box_length
    construction_centers = C15_CENTERS.copy()
    records: list[dict[str, Any]] = []
    selected_pairs = _single_layer_twin_pairs()
    for pair_index, (first, second, vector) in enumerate(selected_pairs):
        site_class = SITE_CLASSES[first]
        twin_type = pair_types[site_class]
        crystal_axis, angle, sigma = TWIN_SPECS[twin_type]
        normal = vector / np.linalg.norm(vector)
        base = _align_vector(crystal_axis, normal)
        sign = 1.0 if pair_index % 2 == 0 else -1.0
        twin = _axis_angle(normal, sign * angle) @ base
        anchor_reduced = C15_CENTERS[first] + 0.5 * vector
        anchor = anchor_reduced * box_length
        # Keep a boundary-crossing pair in one unwrapped coordinate image.
        # Otherwise its anchor can be near zero while the second center remains
        # near one, and a centered candidate lattice will not cover that cell.
        construction_centers[first] = C15_CENTERS[first]
        construction_centers[second] = C15_CENTERS[first] + vector
        swapped = site_class in swap_classes
        orientations[first] = twin if swapped else base
        orientations[second] = base if swapped else twin
        anchors[first] = anchor
        anchors[second] = anchor
        records.append(
            {
                "pair": [first + 1, second + 1],
                "site_class": site_class,
                "twin_type": twin_type,
                "face_normal": normal.tolist(),
                "anchor_reduced": anchor_reduced.tolist(),
                "sigma": sigma,
                "angle_degrees": angle,
                "swapped": swapped,
            }
        )
    # Build a common exclusion catalog containing both supported CSL choices
    # on every retained face.  Ordinary orientations are therefore identical
    # across all Sigma3/Sigma5 and swap controls, rather than being a hidden
    # extra variable in the screening set.
    controlled_catalog: list[FloatArray] = []
    for pair_index, (_, _, vector) in enumerate(selected_pairs):
        normal = vector / np.linalg.norm(vector)
        sign = 1.0 if pair_index % 2 == 0 else -1.0
        for crystal_axis, angle, _ in TWIN_SPECS.values():
            base = _align_vector(crystal_axis, normal)
            controlled_catalog.extend((base, _axis_angle(normal, sign * angle) @ base))
    paired = {grain for pair in selected_pairs for grain in pair[:2]}
    ordinary_indices = [index for index in range(24) if index not in paired]
    rng = np.random.default_rng(RANDOM_SEED)
    existing = list(controlled_catalog)
    for index in ordinary_indices:
        orientations[index] = _random_non_csl_rotation(rng, existing)
        existing.append(orientations[index])
    return orientations, anchors, construction_centers, records


def build_c15_laves(
    *,
    size: float = FIXED_SIZE,
    lattice_constant: float = 3.61,
    overlap_distance_factor: float = 0.5,
    label_grains: bool = False,
    equal_volume: bool = True,
    z16_twin_type: str = "sigma3",
    z12_twin_type: str = "sigma3",
    swap_classes: tuple[str, ...] = (),
    name: str | None = None,
) -> CrystalStructure:
    """Build one constrained 24-grain C15-Laves twin precursor.

    ``size * lattice_constant`` is the cube root of the mean grain volume.
    Consequently the conventional C15 box edge is ``24**(1/3)`` times this
    length, allowing direct size comparison with the A15/W-P builder.
    """
    if size <= 0.0 or lattice_constant <= 0.0:
        raise ValueError("size and lattice_constant must be positive")
    if overlap_distance_factor <= 0.0:
        raise ValueError("overlap_distance_factor must be positive")
    pair_types = {"Z16": z16_twin_type.lower(), "Z12": z12_twin_type.lower()}
    unknown_twin_types = set(pair_types.values()) - set(TWIN_SPECS)
    if unknown_twin_types:
        raise ValueError(
            f"unsupported twin types: {sorted(unknown_twin_types)}; "
            f"choose from {sorted(TWIN_SPECS)}"
        )
    swap_set = frozenset(swap_classes)
    unknown_swap_classes = swap_set - {"Z16", "Z12"}
    if unknown_swap_classes:
        raise ValueError(
            f"unsupported swap classes: {sorted(unknown_swap_classes)}"
        )
    coordination = c15_coordination_numbers()
    if coordination != (16,) * 8 + (12,) * 16:
        raise RuntimeError(f"invalid C15 topology: {coordination}")

    mean_grain_length = size * lattice_constant
    box_length = np.cbrt(24.0) * mean_grain_length
    cell = np.eye(3, dtype=np.float64) * box_length
    orientations, lattice_anchors, construction_centers, twin_pairs = (
        _orientations_and_anchors(box_length, pair_types, swap_set)
    )
    twin_grain_ids = {item for pair in twin_pairs for item in pair["pair"]}
    ordinary_grain_ids = tuple(sorted(set(range(1, 25)) - twin_grain_ids))
    ordinary_indices = [grain_id - 1 for grain_id in ordinary_grain_ids]
    ordinary_minimum_csl_separation = min(
        abs(_cubic_disorientation(orientations[index], orientations[other]) - angle)
        for index in ordinary_indices
        for other in range(24)
        if other != index
        for _, angle, _ in TWIN_SPECS.values()
    )
    center_positions = construction_centers * box_length
    reduced_weights = np.zeros(24, dtype=np.float64)
    if equal_volume:
        reduced_weights[:8] = Z16_EQUAL_VOLUME_WEIGHT
    weights = reduced_weights * box_length**2
    unit_cell = FCC("Cu")
    tolerance = max(box_length, 1.0) * 1e-10

    grain_positions: list[FloatArray] = []
    grain_types: list[NDArray[np.int64]] = []
    grain_ids: list[NDArray[np.int64]] = []
    for grain_index, (center, orientation, anchor) in enumerate(
        zip(center_positions, orientations, lattice_anchors, strict=True)
    ):
        normals: list[FloatArray] = []
        limits: list[float] = []
        for neighbor_index, neighbor in enumerate(C15_CENTERS):
            for shift in PERIODIC_SHIFTS:
                if neighbor_index == grain_index and np.all(shift == 0.0):
                    continue
                relative = (
                    neighbor + shift - construction_centers[grain_index]
                ) * box_length
                normals.append(2.0 * relative)
                limits.append(float(relative @ relative))
                limits[-1] += float(
                    weights[grain_index] - weights[neighbor_index]
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
        grain_positions.append(positions)
        grain_types.append(types)
        grain_ids.append(
            np.full(len(positions), grain_index + 1, dtype=np.int64)
        )

    structure = CrystalStructure(
        positions=np.concatenate(grain_positions),
        cell=cell,
        atom_types=np.concatenate(grain_types),
        type_elements=unit_cell.type_elements,
        origin=np.zeros(3),
        pbc=(True, True, True),
        properties={"grain_ids": np.concatenate(grain_ids)},
        metadata={
            "builder": "TPMSMiner.C15LavesTwinPrecursor",
            "experimental": True,
            "precursor": "C15-Laves",
            "space_group": "Fd-3m",
            "grain_topology": "8xZ16 + 16xZ12",
            "pair_types_by_class": pair_types.copy(),
            "swap_classes": tuple(sorted(swap_set)),
            "sigma3_pair_count": sum(
                record["sigma"] == 3 for record in twin_pairs
            ),
            "sigma5_pair_count": sum(
                record["sigma"] == 5 for record in twin_pairs
            ),
            "equal_volume": bool(equal_volume),
            "grain_weights_reduced": reduced_weights.copy(),
            "grain_weights": weights.copy(),
            "grain_centers_reduced": C15_CENTERS.copy(),
            "construction_centers_reduced": construction_centers.copy(),
            "grain_centers": center_positions.copy(),
            "grain_orientations": orientations.copy(),
            "lattice_anchors": lattice_anchors.copy(),
            "twin_pairs": twin_pairs,
            "twin_layer_rule": "one spatial layer per site-class face-normal direction",
            "twin_grain_count": 2 * len(twin_pairs),
            "ordinary_grain_ids": ordinary_grain_ids,
            "background_orientation_seed": RANDOM_SEED,
            "orientation_sources": tuple(
                "random-ordinary" if grain_id in ordinary_grain_ids else "controlled-CSL"
                for grain_id in range(1, 25)
            ),
            "ordinary_csl_exclusion_degrees": ORDINARY_CSL_EXCLUSION_DEGREES,
            "ordinary_minimum_csl_separation_degrees": ordinary_minimum_csl_separation,
            "coordination_numbers": coordination,
            "mean_grain_length": float(mean_grain_length),
            "box_length": float(box_length),
            "lattice_constant": float(lattice_constant),
            "cropped_to_box": False,
        },
        name=name or f"c15-laves-{pair_types['Z16']}-{pair_types['Z12']}",
    )
    initial_atom_count = structure.number_of_atoms
    nearest_distance = _nearest_periodic_site_distance(
        unit_cell, lattice_constant
    )
    overlap_cutoff = overlap_distance_factor * nearest_distance
    structure, removed = _remove_periodic_overlaps(structure, overlap_cutoff)
    ids, counts = np.unique(structure.properties["grain_ids"], return_counts=True)
    if not np.array_equal(ids, np.arange(1, 25)):
        raise RuntimeError("C15 construction produced one or more empty grains")
    outside = np.any(
        (structure.positions < -tolerance)
        | (structure.positions > box_length + tolerance),
        axis=1,
    )
    structure.metadata.update(
        {
            "initial_atom_count": int(initial_atom_count),
            "overlap_atoms_removed": int(removed),
            "overlap_cutoff": float(overlap_cutoff),
            "nearest_neighbor_distance": float(nearest_distance),
            "outside_box_atom_count": int(np.count_nonzero(outside)),
            "grain_counts": dict(zip(ids.tolist(), counts.tolist(), strict=True)),
        }
    )
    if label_grains:
        structure.atom_types = structure.properties["grain_ids"].copy()
        structure.type_elements = {grain_id: "Cu" for grain_id in range(1, 25)}
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


def write_c15_outputs(
    structure: CrystalStructure,
    output_directory: Path,
) -> dict[str, Any]:
    """Write five C15 configurations and the 24-cell position manifest."""
    output_directory.mkdir(parents=True, exist_ok=True)
    grain_ids = np.asarray(structure.properties["grain_ids"])

    def as_single_type(candidate: CrystalStructure, name: str) -> CrystalStructure:
        candidate.name = name
        candidate.atom_types = np.ones(candidate.number_of_atoms, dtype=np.int64)
        candidate.type_elements = {1: "Cu"}
        candidate.validate()
        return candidate

    combined_path = output_directory / "01_c15_all.data"
    structure.write(combined_path, format="lammps-data", atom_style="atomic")

    z16_only = as_single_type(
        select(structure, grain_ids <= 8), "c15-z16-only"
    )
    z12_only = as_single_type(
        select(structure, grain_ids >= 9), "c15-z12-only"
    )
    z16_only_path = output_directory / "04_c15_z16_only.data"
    z12_only_path = output_directory / "05_c15_z12_only.data"
    z16_only.write(z16_only_path, format="lammps-data", atom_style="atomic")
    z12_only.write(z12_only_path, format="lammps-data", atom_style="atomic")

    representative_paths: dict[str, Path] = {}
    for site_class, grain_index, filename in (
        ("Z16", 0, "02_z16_cell_origin.data"),
        ("Z12", 8, "03_z12_cell_origin.data"),
    ):
        grain_id = grain_index + 1
        representative = as_single_type(
            select(structure, grain_ids == grain_id),
            f"representative-{site_class.lower()}-cell-at-origin",
        )
        representative.positions -= np.asarray(
            structure.metadata["grain_centers"][grain_index]
        )
        representative.origin = -0.5 * np.diag(representative.cell)
        representative.metadata.update(
            {
                "representative_site_class": site_class,
                "source_grain_id": grain_id,
                "cell_center": np.zeros(3),
            }
        )
        representative.validate()
        path = output_directory / filename
        representative.write(path, format="lammps-data", atom_style="atomic")
        representative_paths[site_class] = path

    pair_by_grain: dict[int, dict[str, Any]] = {}
    for record in structure.metadata["twin_pairs"]:
        first, second = record["pair"]
        pair_by_grain[first] = {**record, "partner": second}
        pair_by_grain[second] = {
            **record,
            "partner": first,
            "face_normal": (-np.asarray(record["face_normal"])).tolist(),
        }

    class_indices = {"Z16": 0, "Z12": 0}
    cell_records: list[dict[str, Any]] = []
    for grain_index in range(24):
        grain_id = grain_index + 1
        site_class = SITE_CLASSES[grain_index]
        class_indices[site_class] += 1
        class_index = class_indices[site_class]
        cell_name = f"{site_class}_{class_index:02d}_grain_{grain_id:02d}"

        twin = pair_by_grain.get(grain_id)
        cell_records.append(
            {
                "grain_id": grain_id,
                "cell_name": cell_name,
                "site_class": site_class,
                "coordination_number": int(
                    structure.metadata["coordination_numbers"][grain_index]
                ),
                "atom_count": int(np.count_nonzero(grain_ids == grain_id)),
                "class_configuration_file": (
                    z16_only_path.name if site_class == "Z16" else z12_only_path.name
                ),
                "center_reduced": C15_CENTERS[grain_index].tolist(),
                "construction_center_reduced_unwrapped": structure.metadata[
                    "construction_centers_reduced"
                ][grain_index].tolist(),
                "center_cartesian_unwrapped": structure.metadata[
                    "grain_centers"
                ][grain_index].tolist(),
                "orientation_matrix_crystal_to_lab": structure.metadata[
                    "grain_orientations"
                ][grain_index].tolist(),
                "orientation_source": structure.metadata["orientation_sources"][
                    grain_index
                ],
                "lattice_anchor_cartesian": structure.metadata["lattice_anchors"][
                    grain_index
                ].tolist(),
                "boundary_role": "designed_csl_pair" if twin else "ordinary_grain",
                "twin_partner_grain_id": int(twin["partner"]) if twin else None,
                "twin_face_normal_outward": twin["face_normal"] if twin else None,
                "shared_anchor_reduced": twin["anchor_reduced"] if twin else None,
                "sigma": int(twin["sigma"]) if twin else None,
                "misorientation_degrees": float(twin["angle_degrees"]) if twin else None,
            }
        )

    manifest = {
        "schema_version": 1,
        "coordinate_convention": {
            "center_reduced": "fractional position in the primary C15 box",
            "construction_center_reduced_unwrapped": (
                "periodic image used to construct the complete unwrapped cell"
            ),
            "center_cartesian_unwrapped": "angstrom",
            "orientation_matrix_crystal_to_lab": (
                "column-vector convention: v_lab = O @ v_crystal"
            ),
            "representative_cells": "cell center translated to Cartesian (0,0,0)",
            "class_only_files": "retain global C15 box and relative positions",
        },
        "configuration_files": {
            "c15_all": combined_path.name,
            "z16_representative_at_origin": representative_paths["Z16"].name,
            "z12_representative_at_origin": representative_paths["Z12"].name,
            "z16_only": z16_only_path.name,
            "z12_only": z12_only_path.name,
        },
        "number_of_cells": 24,
        "site_class_counts": {"Z16": 8, "Z12": 16},
        "twin_layer_rule": structure.metadata["twin_layer_rule"],
        "designed_csl_pair_count": len(structure.metadata["twin_pairs"]),
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
        "box_length": float(structure.metadata["box_length"]),
        "lattice_constant": float(structure.metadata["lattice_constant"]),
        "equal_volume": bool(structure.metadata["equal_volume"]),
        "cells": cell_records,
    }
    manifest_path = output_directory / "cell_manifest.json"
    manifest_path.write_text(
        json.dumps(_jsonable(manifest), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    csv_path = output_directory / "cell_manifest.csv"
    csv_fields = [
        "grain_id", "cell_name", "site_class", "coordination_number",
        "atom_count", "class_configuration_file",
        "center_fx", "center_fy", "center_fz",
        "construction_fx", "construction_fy", "construction_fz",
        "center_x_A", "center_y_A", "center_z_A", "twin_partner_grain_id",
        "normal_x", "normal_y", "normal_z", "anchor_fx", "anchor_fy",
        "anchor_fz", "boundary_role", "orientation_source", "sigma",
        "misorientation_degrees",
    ] + [f"orientation_{row}{column}" for row in range(3) for column in range(3)]
    with csv_path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=csv_fields)
        writer.writeheader()
        for record in cell_records:
            center = record["center_reduced"]
            construction = record["construction_center_reduced_unwrapped"]
            cartesian = record["center_cartesian_unwrapped"]
            normal = record["twin_face_normal_outward"] or (None, None, None)
            anchor = record["shared_anchor_reduced"] or (None, None, None)
            row = {
                key: record[key]
                for key in (
                    "grain_id", "cell_name", "site_class",
                    "coordination_number", "atom_count",
                    "class_configuration_file",
                    "boundary_role", "orientation_source",
                    "twin_partner_grain_id", "sigma",
                    "misorientation_degrees",
                )
            }
            row.update(
                {
                    "center_fx": center[0], "center_fy": center[1],
                    "center_fz": center[2],
                    "construction_fx": construction[0],
                    "construction_fy": construction[1],
                    "construction_fz": construction[2],
                    "center_x_A": cartesian[0], "center_y_A": cartesian[1],
                    "center_z_A": cartesian[2],
                    "normal_x": normal[0], "normal_y": normal[1],
                    "normal_z": normal[2],
                    "anchor_fx": anchor[0], "anchor_fy": anchor[1],
                    "anchor_fz": anchor[2],
                }
            )
            orientation = record["orientation_matrix_crystal_to_lab"]
            row.update(
                {
                    f"orientation_{r}{c}": orientation[r][c]
                    for r in range(3)
                    for c in range(3)
                }
            )
            writer.writerow(row)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lattice-constant", type=float, default=3.61)
    parser.add_argument("--overlap-distance-factor", type=float, default=0.5)
    parser.add_argument("--label-grains", action="store_true")
    parser.add_argument(
        "--raw-voronoi",
        action="store_true",
        help="disable the default equal-volume C15 Laguerre weighting",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "output",
        help="directory for five configurations and the position manifests",
    )
    options = parser.parse_args()

    structure = build_c15_laves(
        size=FIXED_SIZE,
        lattice_constant=options.lattice_constant,
        overlap_distance_factor=options.overlap_distance_factor,
        label_grains=options.label_grains,
        equal_volume=not options.raw_voronoi,
    )
    manifest = write_c15_outputs(structure, options.output_dir)
    print(
        f"Wrote {structure.number_of_atoms} Cu atoms in 24 C15 grains to "
        f"{(options.output_dir / manifest['configuration_files']['c15_all']).resolve()}."
    )
    print(
        f"Topology: 8xZ16 + 16xZ12; "
        f"Sigma3 pairs: {structure.metadata['sigma3_pair_count']}; "
        f"Sigma5 pairs: {structure.metadata['sigma5_pair_count']}; "
        f"box edge: {structure.metadata['box_length']:.6g} A."
    )
    print(f"Fixed size: {FIXED_SIZE}; size-series screening is disabled.")
    print(
        f"Wrote five configurations and JSON/CSV manifests to "
        f"{options.output_dir.resolve()}."
    )


if __name__ == "__main__":
    main()
