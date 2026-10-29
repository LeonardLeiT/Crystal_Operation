"""Voronoi and Kelvin polycrystal builders."""

from __future__ import annotations

from dataclasses import dataclass
import itertools
from typing import ClassVar

import numpy as np
from numpy.typing import ArrayLike, NDArray

from crystalops.builders._lattice import candidate_lattice, orientation_from_rows
from crystalops.builders.base import StructureBuilder
from crystalops.core import CrystalStructure, UnitCell
from crystalops.lattices import body_centered_cubic
from crystalops.operations import remove_duplicates, rotation_matrix, select


FloatArray = NDArray[np.float64]


def _nearest_periodic_site_distance(
    unit_cell: UnitCell,
    lattice_constant: float,
) -> float:
    """Return the shortest non-zero site distance of a periodic unit cell."""
    translations = np.stack(
        np.meshgrid((-1, 0, 1), (-1, 0, 1), (-1, 0, 1), indexing="ij"),
        axis=-1,
    ).reshape(-1, 3)
    primitive = unit_cell.cell * lattice_constant
    shortest = np.inf
    for first_index, first in enumerate(unit_cell.fractional_positions):
        for second_index, second in enumerate(unit_cell.fractional_positions):
            differences = second + translations - first
            distances = np.linalg.norm(differences @ primitive, axis=1)
            if first_index == second_index:
                distances[np.all(translations == 0, axis=1)] = np.inf
            shortest = min(shortest, float(np.min(distances)))
    if not np.isfinite(shortest) or shortest <= 0.0:
        raise ValueError("unit cell has no positive periodic site distance")
    return shortest


def _remove_periodic_overlaps(
    structure: CrystalStructure,
    cutoff: float,
) -> tuple[CrystalStructure, int]:
    """Remove periodic overlaps while retaining original unwrapped positions."""
    if cutoff <= 0.0:
        raise ValueError("overlap cutoff must be positive")
    if not np.allclose(structure.cell, np.diag(np.diag(structure.cell))):
        raise ValueError("periodic overlap removal currently requires an orthogonal cell")
    lengths = np.diag(structure.cell)
    if np.any(lengths <= 0.0):
        raise ValueError("cell lengths must be positive")

    wrapped = np.mod(structure.positions - structure.origin, lengths)
    wrapped[np.isclose(wrapped, lengths[None, :], atol=1e-12)] = 0.0
    bins_per_axis = np.maximum(1, np.floor(lengths / cutoff).astype(np.int64))
    bin_indices = np.floor(
        wrapped / lengths[None, :] * bins_per_axis[None, :]
    ).astype(np.int64)
    bin_indices %= bins_per_axis[None, :]

    neighbor_offsets = np.stack(
        np.meshgrid((-1, 0, 1), (-1, 0, 1), (-1, 0, 1), indexing="ij"),
        axis=-1,
    ).reshape(-1, 3)
    buckets: dict[tuple[int, int, int], list[int]] = {}
    retained: list[int] = []
    cutoff_squared = cutoff**2
    for atom_index, bin_index in enumerate(bin_indices):
        overlaps = False
        for offset in neighbor_offsets:
            neighbor_bin = tuple(
                ((bin_index + offset) % bins_per_axis).tolist()
            )
            for retained_index in buckets.get(neighbor_bin, ()):
                difference = wrapped[atom_index] - wrapped[retained_index]
                difference -= np.round(difference / lengths) * lengths
                if difference @ difference < cutoff_squared:
                    overlaps = True
                    break
            if overlaps:
                break
        if not overlaps:
            retained.append(atom_index)
            key = tuple(bin_index.tolist())
            buckets.setdefault(key, []).append(atom_index)

    removed = structure.number_of_atoms - len(retained)
    return select(structure, np.asarray(retained, dtype=np.int64)), removed


@dataclass(slots=True)
class VoronoiPolycrystalBuilder(StructureBuilder):
    """Fill periodic Voronoi grains with independently oriented crystals.

    Grain centers are given in fractional coordinates of the final simulation
    box. Candidate atoms are generated around every center, wrapped into the
    box, and retained only when that center is their nearest periodic seed.

    Examples:
        >>> import numpy as np
        >>> from crystalops.lattices import face_centered_cubic
        >>> builder = VoronoiPolycrystalBuilder(
        ...     unit_cell=face_centered_cubic("Cu"),
        ...     cell=np.eye(3) * 10.0,
        ...     grain_centers=[[0.25, 0.5, 0.5], [0.75, 0.5, 0.5]],
        ...     orientations=[np.eye(3), np.eye(3)],
        ...     lattice_constant=3.61,
        ... )
        >>> builder.build().number_of_atoms > 0
        True
    """

    unit_cell: UnitCell
    cell: ArrayLike
    grain_centers: ArrayLike
    orientations: ArrayLike
    lattice_constant: float = 1.0
    origin: ArrayLike = (0.0, 0.0, 0.0)
    tolerance: float = 1e-7
    name: str = "voronoi-polycrystal"
    grain_weights: ArrayLike | None = None

    def build(self) -> CrystalStructure:
        """Build the periodic polycrystal.

        Returns:
            A structure containing all generated Voronoi grains.
        """
        cell = np.asarray(self.cell, dtype=np.float64)
        centers = np.asarray(self.grain_centers, dtype=np.float64)
        orientations = np.asarray(self.orientations, dtype=np.float64)
        weights = (
            np.zeros(len(centers), dtype=np.float64)
            if self.grain_weights is None
            else np.asarray(self.grain_weights, dtype=np.float64)
        )
        origin = np.asarray(self.origin, dtype=np.float64)
        if cell.shape != (3, 3):
            raise ValueError("cell must have shape (3, 3)")
        if centers.ndim != 2 or centers.shape[1] != 3 or len(centers) == 0:
            raise ValueError("grain_centers must have shape (number_of_grains, 3)")
        if orientations.shape != (len(centers), 3, 3):
            raise ValueError("orientations must have shape (number_of_grains, 3, 3)")
        if weights.shape != (len(centers),):
            raise ValueError("grain_weights must have shape (number_of_grains,)")
        if not np.all(np.isfinite(weights)):
            raise ValueError("grain_weights must contain only finite values")
        if origin.shape != (3,):
            raise ValueError("origin must have shape (3,)")
        if self.lattice_constant <= 0:
            raise ValueError("lattice_constant must be positive")
        for orientation in orientations:
            if not np.allclose(orientation @ orientation.T, np.eye(3), atol=1e-10):
                raise ValueError("each orientation must be orthogonal")

        grain_positions: list[FloatArray] = []
        grain_types: list[NDArray[np.int64]] = []
        grain_ids: list[NDArray[np.int64]] = []
        center_cartesian = centers @ cell + origin
        inverse_cell = np.linalg.inv(cell)
        for grain_index, (center, orientation) in enumerate(
            zip(center_cartesian, orientations, strict=True)
        ):
            candidates = candidate_lattice(
                self.unit_cell,
                self.lattice_constant,
                cell,
                orientation,
            )
            candidates += center
            fractional = (candidates - origin) @ inverse_cell
            inside_box = np.all(
                (fractional >= -self.tolerance)
                & (fractional < 1.0 - self.tolerance),
                axis=1,
            )
            fractional = fractional[inside_box]

            differences = fractional[:, None, :] - centers[None, :, :]
            differences -= np.round(differences)
            cartesian_differences = differences @ cell
            distances_squared = np.einsum(
                "ngk,ngk->ng",
                cartesian_differences,
                cartesian_differences,
            )
            power_distances = distances_squared - weights[None, :]
            selected = np.argmin(power_distances, axis=1) == grain_index
            positions = fractional[selected] @ cell + origin
            site_types = np.tile(
                self.unit_cell.atom_types,
                len(candidates) // self.unit_cell.number_of_sites,
            )[inside_box][selected]
            grain_positions.append(positions)
            grain_types.append(site_types)
            grain_ids.append(
                np.full(len(positions), grain_index + 1, dtype=np.int64)
            )

        structure = CrystalStructure(
            positions=np.concatenate(grain_positions, axis=0),
            cell=cell,
            atom_types=np.concatenate(grain_types),
            type_elements=self.unit_cell.type_elements,
            origin=origin,
            properties={"grain_ids": np.concatenate(grain_ids)},
            metadata={
                "builder": type(self).__name__,
                "grain_centers": centers.copy(),
                "grain_weights": weights.copy(),
                "grain_orientations": orientations.copy(),
                "lattice_constant": float(self.lattice_constant),
            },
            name=self.name,
        )
        return remove_duplicates(structure, tolerance=self.tolerance)


@dataclass(slots=True)
class WeairePhelanCrystal(StructureBuilder):
    """Build the A15 Weaire-Phelan precursor used by Yang et al. (2025).

    The cubic periodic cell contains the two dodecahedral grains ``Dod0``
    and ``Dod1`` and three perpendicular pairs of tetrakaidecahedral grains.
    Their fractional centers and crystallographic orientations reproduce
    Section 2 and Table 1 of the paper.  A small Laguerre weight on the two
    Dod sites makes all eight planar precursor cells exactly equal in volume.

    ``size`` is the nominal grain size ``d`` in lattice units.  Following the
    paper's convention ``a = 2d``, the cubic simulation-cell edge is
    ``2 * size * lattice_constant``.
    """

    unit_cell: UnitCell
    size: float
    lattice_constant: float = 1.0
    twin_type: str = "sigma3"
    label_grains: bool = False
    overlap_distance_factor: float = 0.5
    dod_orientation_seed: int | None = None
    tolerance: float = 1e-7
    name: str = "weaire-phelan-polycrystal"

    # Fractional A15 sites in the order used in Table 1: the two Dods,
    # followed by Tet X, Y, and Z pairs (T1, T2).
    CENTERS: ClassVar[FloatArray] = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.5],
            [0.75, 0.0, 0.5],
            [0.25, 0.0, 0.5],
            [0.5, 0.75, 0.0],
            [0.5, 0.25, 0.0],
            [0.0, 0.5, 0.75],
            [0.0, 0.5, 0.25],
        ],
        dtype=np.float64,
    )
    CELL_LABELS: ClassVar[tuple[str, ...]] = (
        "Dod0", "Dod1", "TetX-T1", "TetX-T2",
        "TetY-T1", "TetY-T2", "TetZ-T1", "TetZ-T2",
    )

    # Dimensionless power weight for the two Dod sites when the cubic box
    # edge is one.  It is the root for which every periodic A15 power cell
    # has volume 1/8; Tet weights use the zero gauge.
    DOD_EQUAL_VOLUME_WEIGHT: ClassVar[float] = 0.002480262473718248
    MAX_FACE_NEIGHBOR_DISTANCE: ClassVar[float] = float(np.sqrt(3.0 / 8.0))
    DOD_CSL_EXCLUSION_DEGREES: ClassVar[float] = 3.0

    ORIENTATION_ROWS: ClassVar[
        dict[str, tuple[tuple[tuple[int, int, int], ...], ...]]
    ] = {
        "sigma3": (
            ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            ((-1, -1, -1), (1, 1, -2), (1, -1, 0)),
            ((1, 1, -1), (1, 1, 2), (1, -1, 0)),
            ((0, 1, -1), (-1, -1, -1), (-2, 1, 1)),
            ((0, 1, -1), (-1, 1, 1), (2, 1, 1)),
            ((1, -2, 1), (-1, 0, 1), (-1, -1, -1)),
            ((1, 2, 1), (-1, 0, 1), (1, -1, 1)),
        ),
        "sigma5": (
            ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            ((1, 0, 0), (0, 1, 0), (0, 0, 1)),
            ((1, -2, 0), (2, 1, 0), (0, 0, 1)),
            ((2, -1, 0), (1, 2, 0), (0, 0, 1)),
            ((1, 0, 0), (0, 1, -2), (0, 2, 1)),
            ((1, 0, 0), (0, 2, -1), (0, 1, 2)),
            ((1, 0, -2), (0, 1, 0), (2, 0, 1)),
            ((2, 0, -1), (0, 1, 0), (1, 0, 2)),
        ),
    }

    @staticmethod
    def _random_rotation(rng: np.random.Generator) -> FloatArray:
        """Return a uniformly sampled proper rotation from a unit quaternion."""
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

    @staticmethod
    def _cubic_disorientation(first: FloatArray, second: FloatArray) -> float:
        """Return the minimum misorientation angle under cubic symmetry."""
        angles: list[float] = []
        for permutation in itertools.permutations(range(3)):
            for signs in itertools.product((-1.0, 1.0), repeat=3):
                symmetry = np.zeros((3, 3), dtype=np.float64)
                symmetry[np.arange(3), permutation] = signs
                if np.linalg.det(symmetry) <= 0.0:
                    continue
                relative = second @ symmetry @ first.T
                cosine = np.clip((np.trace(relative) - 1.0) / 2.0, -1.0, 1.0)
                angles.append(float(np.degrees(np.arccos(cosine))))
        return min(angles)

    def _random_non_csl_rotation(
        self,
        rng: np.random.Generator,
        existing: list[FloatArray],
    ) -> FloatArray:
        """Sample a random orientation separated from near-Sigma3/5 relations."""
        csl_angles = (60.0, 36.86989764584402)
        for _ in range(10000):
            candidate = self._random_rotation(rng)
            angles = [self._cubic_disorientation(candidate, item) for item in existing]
            if all(
                abs(angle - target) >= self.DOD_CSL_EXCLUSION_DEGREES
                for angle in angles
                for target in csl_angles
            ):
                return candidate
        raise RuntimeError("could not sample non-CSL Dod orientation")

    def build(self) -> CrystalStructure:
        """Build the eight-grain periodic W-P precursor."""
        if self.size <= 0:
            raise ValueError("size must be positive")
        if self.lattice_constant <= 0:
            raise ValueError("lattice_constant must be positive")
        if self.overlap_distance_factor <= 0:
            raise ValueError("overlap_distance_factor must be positive")

        twin_key = self.twin_type.lower().replace("-", "")
        aliases = {"sigma3": "sigma3", "s3": "sigma3", "3": "sigma3",
                   "sigma5": "sigma5", "s5": "sigma5", "5": "sigma5"}
        try:
            canonical_twin = aliases[twin_key]
        except KeyError as exc:
            raise ValueError("twin_type must be 'sigma3' or 'sigma5'") from exc

        orientations = np.asarray(
            [
                orientation_from_rows(rows)
                for rows in self.ORIENTATION_ROWS[canonical_twin]
            ]
        )
        orientation_sources = [f"paper-{canonical_twin}"] * len(orientations)
        if self.dod_orientation_seed is not None:
            rng = np.random.default_rng(self.dod_orientation_seed)
            # Use the complete paper orientation catalog, rather than only the
            # currently selected family.  Dod0/1 are therefore identical in
            # every Sigma3/Sigma5 screening candidate and cannot confound the
            # comparison between Tet networks.
            existing = [
                orientation_from_rows(rows)
                for family in WeairePhelanCrystal.ORIENTATION_ROWS.values()
                for rows in family[2:]
            ]
            first_dod = self._random_non_csl_rotation(rng, existing)
            second_dod = self._random_non_csl_rotation(rng, [*existing, first_dod])
            orientations[:2] = np.asarray([first_dod, second_dod])
            orientation_sources[:2] = ["random-Dod0", "random-Dod1"]
        box_length = 2.0 * self.size * self.lattice_constant
        cell = np.eye(3, dtype=np.float64) * box_length
        weights = np.zeros(len(self.CENTERS), dtype=np.float64)
        weights[:2] = self.DOD_EQUAL_VOLUME_WEIGHT * box_length**2

        # Construct every complete power cell around its own center, then
        # translate it to the corresponding A15 site.  In particular, Dod0
        # and the boundary-crossing Tets are not clipped to or wrapped into
        # the primary simulation box.  This mirrors KelvinCrystal and keeps
        # every precursor grain visually and atomistically intact.
        periodic_shifts = np.stack(
            np.meshgrid((-1, 0, 1), (-1, 0, 1), (-1, 0, 1), indexing="ij"),
            axis=-1,
        ).reshape(-1, 3)
        center_positions = self.CENTERS * box_length
        lattice_origins_reduced = np.asarray(
            [
                self.CENTERS[0],
                self.CENTERS[1],
                [0.5, 0.0, 0.5],
                [0.5, 0.0, 0.5],
                [0.5, 0.5, 0.0],
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.5],
                [0.0, 0.5, 0.5],
            ],
            dtype=np.float64,
        )
        lattice_origins = lattice_origins_reduced * box_length
        grain_positions: list[FloatArray] = []
        grain_types: list[NDArray[np.int64]] = []
        grain_ids_list: list[NDArray[np.int64]] = []

        for grain_index, (center, orientation, lattice_origin) in enumerate(
            zip(center_positions, orientations, lattice_origins, strict=True)
        ):
            neighbor_vectors: list[FloatArray] = []
            limits: list[float] = []
            for neighbor_index, neighbor_center in enumerate(self.CENTERS):
                for shift in periodic_shifts:
                    if neighbor_index == grain_index and np.all(shift == 0):
                        continue
                    relative = (
                        neighbor_center + shift - self.CENTERS[grain_index]
                    ) * box_length
                    if (
                        np.linalg.norm(relative) / box_length
                        > self.MAX_FACE_NEIGHBOR_DISTANCE + 1e-12
                    ):
                        continue
                    neighbor_vectors.append(2.0 * relative)
                    limits.append(
                        float(
                            relative @ relative
                            + weights[grain_index]
                            - weights[neighbor_index]
                        )
                    )

            halfspace_normals = np.asarray(neighbor_vectors)
            halfspace_limits = np.asarray(limits)
            candidates = candidate_lattice(
                self.unit_cell,
                self.lattice_constant,
                cell,
                orientation,
            )
            absolute_candidates = candidates + lattice_origin
            relative_to_cell = absolute_candidates - center
            inside = np.all(
                relative_to_cell @ halfspace_normals.T
                <= halfspace_limits[None, :] + self.tolerance,
                axis=1,
            )
            positions = absolute_candidates[inside]
            site_types = np.tile(
                self.unit_cell.atom_types,
                len(candidates) // self.unit_cell.number_of_sites,
            )[inside]
            grain_positions.append(positions)
            grain_types.append(site_types)
            grain_ids_list.append(
                np.full(len(positions), grain_index + 1, dtype=np.int64)
            )

        origin = np.zeros(3, dtype=np.float64)
        structure = CrystalStructure(
            positions=np.concatenate(grain_positions),
            cell=cell,
            atom_types=np.concatenate(grain_types),
            type_elements=self.unit_cell.type_elements,
            origin=origin,
            pbc=(True, True, True),
            properties={"grain_ids": np.concatenate(grain_ids_list)},
            metadata={
                "grain_centers": self.CENTERS.copy(),
                "grain_centers_cartesian": center_positions.copy(),
                "lattice_origins_reduced": lattice_origins_reduced.copy(),
                "lattice_origins": lattice_origins.copy(),
                "grain_weights": weights.copy(),
                "grain_orientations": orientations.copy(),
                "orientation_sources": tuple(orientation_sources),
                "dod_orientation_seed": self.dod_orientation_seed,
                "dod_csl_exclusion_degrees": self.DOD_CSL_EXCLUSION_DEGREES,
                "lattice_constant": float(self.lattice_constant),
                "cropped_to_box": False,
            },
            name=self.name,
        )
        initial_atom_count = structure.number_of_atoms
        nearest_neighbor_distance = _nearest_periodic_site_distance(
            self.unit_cell,
            self.lattice_constant,
        )
        overlap_cutoff = (
            self.overlap_distance_factor * nearest_neighbor_distance
        )
        structure, overlap_atoms_removed = _remove_periodic_overlaps(
            structure,
            overlap_cutoff,
        )
        outside_box = np.any(
            (structure.positions < origin - self.tolerance)
            | (structure.positions > box_length + self.tolerance),
            axis=1,
        )
        structure.metadata["outside_box_atom_count"] = int(
            np.count_nonzero(outside_box)
        )
        grain_ids = structure.properties["grain_ids"]
        unique_grains, grain_counts = np.unique(grain_ids, return_counts=True)
        expected_grains = np.arange(1, len(self.CENTERS) + 1, dtype=np.int64)
        if not np.array_equal(unique_grains, expected_grains):
            missing = sorted(set(expected_grains) - set(unique_grains))
            raise RuntimeError(
                f"Weaire-Phelan construction produced empty grains {missing}; "
                "increase size"
            )

        structure.metadata.update(
            {
                "builder": type(self).__name__,
                "precursor": "A15-Weaire-Phelan",
                "paper_model": "Yang et al., Acta Materialia 296 (2025) 121289",
                "twin_type": canonical_twin,
                "paper_reference_orientation": self.dod_orientation_seed is None,
                "available_twin_types": ("sigma3", "sigma5"),
                "cell_labels": self.CELL_LABELS,
                "cell_types": ("Dod", "Dod", "Tet", "Tet", "Tet", "Tet", "Tet", "Tet"),
                "nominal_grain_size": float(self.size * self.lattice_constant),
                "supercell_size": float(box_length),
                "initial_atom_count": int(initial_atom_count),
                "overlap_atoms_removed": int(overlap_atoms_removed),
                "overlap_cutoff": float(overlap_cutoff),
                "nearest_neighbor_distance": float(nearest_neighbor_distance),
                "grain_counts": dict(
                    zip(unique_grains.tolist(), grain_counts.tolist(), strict=True)
                ),
            }
        )
        if self.label_grains:
            elements = set(self.unit_cell.type_elements.values())
            if len(elements) != 1:
                raise ValueError(
                    "label_grains requires a single-element unit cell because "
                    "LAMMPS atom types cannot represent element and grain separately"
                )
            element = next(iter(elements))
            structure.atom_types = grain_ids.copy()
            structure.type_elements = {
                int(grain_id): element for grain_id in unique_grains
            }
            structure.validate()
        return structure


@dataclass(slots=True)
class KelvinCrystal(StructureBuilder):
    """Arrange 16 complete truncated-octahedral grains on BCC centers.

    The named orientation presets retain the intent of the original CrystalOps
    DSC structures, while ``random_seed`` makes all remaining grain
    orientations reproducible. Each grain is clipped only by the six square
    and eight hexagonal faces of its own truncated octahedron. Grains are then
    translated to a 2 x 2 x 2 BCC arrangement without clipping them against a
    global periodic box.

    Examples:
        >>> from crystalops.builders import KelvinCrystal
        >>> from crystalops.lattices import face_centered_cubic
        >>> structure = KelvinCrystal(
        ...     unit_cell=face_centered_cubic("Cu"),
        ...     size=3.0,
        ...     lattice_constant=3.61,
        ...     preset="DSC1",
        ...     random_seed=7,
        ... ).build()
        >>> structure.number_of_atoms > 0
        True
    """

    unit_cell: UnitCell
    size: float
    lattice_constant: float = 1.0
    preset: str = "DSC1"
    random_seed: int | None = None
    grain_cell: UnitCell | None = None
    label_grains: bool = False
    name: str = "kelvin-polycrystal"

    PRESETS: ClassVar[dict[str, set[int]]] = {
        "full": {0, 1, 2, 4, 6, 8, 10, 12, 14, 15},
        "DSC1": {1, 2, 4, 8, 14, 15},
        "DSC2": {0, 1, 6, 10, 12, 15},
        "FCC": {1, 6, 10, 12, 15},
        "anti-FCC": {1, 2, 4, 8, 15},
        "NB": {1, 2, 4, 14, 15},
        "anti-NB": {0, 1, 10, 12, 15},
        "NC": {1, 2, 8, 14, 15},
        "anti-NC": {0, 1, 6, 12, 15},
        "ND": {1, 4, 8, 14, 15},
        "anti-ND": {0, 1, 6, 10, 15},
    }

    def build(self) -> CrystalStructure:
        """Build 16 complete, oriented truncated-octahedral grains."""
        if self.size <= 0:
            raise ValueError("size must be positive")
        preset_key = self.preset.replace("_", "-")
        aliases = {key.lower(): key for key in self.PRESETS}
        try:
            canonical_preset = aliases[preset_key.lower()]
        except KeyError as exc:
            raise ValueError(
                f"unknown preset {self.preset!r}; available presets: "
                f"{tuple(self.PRESETS)}"
            ) from exc

        grain_cell = self.grain_cell or body_centered_cubic("Cu")
        translations = np.stack(
            np.meshgrid(np.arange(2), np.arange(2), np.arange(2), indexing="ij"),
            axis=-1,
        ).reshape(-1, 3)
        center_indices = (
            translations[:, None, :] + grain_cell.fractional_positions[None, :, :]
        ).reshape(-1, 3)
        centers = center_indices / 2.0

        fixed_orientations = [
            np.eye(3),
            orientation_from_rows([[-1, 2, 2], [2, -1, 2], [2, 2, -1]]),
            orientation_from_rows([[-1, -2, -2], [-2, -1, 2], [-2, 2, -1]]),
            orientation_from_rows([[-1, -2, 2], [-2, -1, -2], [2, -2, -1]]),
            orientation_from_rows([[-1, 2, -2], [2, -1, -2], [-2, -2, -1]]),
        ]
        position_orientation = {
            1: 0, 15: 0, 0: 1, 14: 1, 8: 2,
            6: 2, 10: 3, 4: 3, 2: 4, 12: 4,
        }
        rng = np.random.default_rng(self.random_seed)
        fixed_positions = self.PRESETS[canonical_preset]
        orientations: list[FloatArray] = []
        orientation_sources: list[str] = []
        for index in range(len(centers)):
            if index in fixed_positions and index in position_orientation:
                orientation_index = position_orientation[index]
                orientations.append(fixed_orientations[orientation_index])
                orientation_sources.append(f"fixed-{orientation_index}")
            else:
                orientations.append(
                    rotation_matrix(*rng.uniform(0.0, 360.0, size=3))
                )
                orientation_sources.append("random")

        if not np.allclose(grain_cell.cell, np.eye(3), atol=1e-12):
            raise ValueError(
                "KelvinCrystal requires a cubic BCC grain cell"
            )

        grain_spacing = self.size * self.lattice_constant
        half_square_face = 0.5 * grain_spacing
        hexagonal_face = 0.75 * grain_spacing
        tolerance = max(grain_spacing, 1.0) * 1e-10
        center_positions = center_indices * grain_spacing
        grain_positions: list[FloatArray] = []
        grain_types: list[NDArray[np.int64]] = []
        grain_ids_list: list[NDArray[np.int64]] = []
        shape_box = np.eye(3) * grain_spacing

        for grain_index, (center, orientation) in enumerate(
            zip(center_positions, orientations, strict=True)
        ):
            candidates = candidate_lattice(
                self.unit_cell,
                self.lattice_constant,
                shape_box,
                orientation,
            )
            absolute_coordinates = np.abs(candidates)
            inside = (
                np.all(
                    absolute_coordinates <= half_square_face + tolerance,
                    axis=1,
                )
                & (
                    np.sum(absolute_coordinates, axis=1)
                    <= hexagonal_face + tolerance
                )
            )
            positions = candidates[inside] + center
            site_types = np.tile(
                self.unit_cell.atom_types,
                len(candidates) // self.unit_cell.number_of_sites,
            )[inside]
            grain_positions.append(positions)
            grain_types.append(site_types)
            grain_ids_list.append(
                np.full(len(positions), grain_index + 1, dtype=np.int64)
            )

        origin = np.zeros(3, dtype=np.float64)
        box = grain_cell.cell * (2.0 * grain_spacing)
        structure = CrystalStructure(
            positions=np.concatenate(grain_positions, axis=0),
            cell=box,
            atom_types=np.concatenate(grain_types),
            type_elements=self.unit_cell.type_elements,
            origin=origin,
            pbc=(True, True, True),
            properties={"grain_ids": np.concatenate(grain_ids_list)},
            metadata={
                "grain_centers": center_positions.copy(),
                "grain_centers_reduced": centers.copy(),
                "grain_orientations": np.asarray(orientations).copy(),
                "grain_shape": "truncated-octahedron",
                "grain_spacing": float(grain_spacing),
                "cropped_to_box": False,
                "lattice_constant": float(self.lattice_constant),
            },
            name=self.name,
        )
        outside_box = np.any(
            (structure.positions < structure.origin - tolerance)
            | (
                structure.positions
                > structure.origin + np.sum(structure.cell, axis=0) + tolerance
            ),
            axis=1,
        )
        structure.metadata["outside_box_atom_count"] = int(
            np.count_nonzero(outside_box)
        )
        grain_ids = structure.properties["grain_ids"]
        unique_grains, grain_counts = np.unique(grain_ids, return_counts=True)
        expected_grains = np.arange(1, len(centers) + 1, dtype=np.int64)
        if not np.array_equal(unique_grains, expected_grains):
            missing = sorted(set(expected_grains) - set(unique_grains))
            raise RuntimeError(
                "Kelvin construction produced empty grains "
                f"{missing}; increase size or lattice resolution"
            )

        structure.metadata.update(
            {
                "builder": type(self).__name__,
                "preset": canonical_preset,
                "random_seed": self.random_seed,
                "orientation_sources": tuple(orientation_sources),
                "fixed_positions": tuple(sorted(fixed_positions)),
                "random_positions": tuple(
                    index
                    for index, source in enumerate(orientation_sources)
                    if source == "random"
                ),
                "grain_counts": dict(
                    zip(unique_grains.tolist(), grain_counts.tolist(), strict=True)
                ),
            }
        )
        if self.label_grains:
            elements = set(self.unit_cell.type_elements.values())
            if len(elements) != 1:
                raise ValueError(
                    "label_grains requires a single-element unit cell because "
                    "LAMMPS atom types cannot represent element and grain separately"
                )
            element = next(iter(elements))
            structure.atom_types = grain_ids.copy()
            structure.type_elements = {
                int(grain_id): element for grain_id in unique_grains
            }
            structure.validate()
        return structure
