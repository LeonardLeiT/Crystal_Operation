"""Mathematical D- and P-Schwarz crystal builder."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from crystalops.builders._lattice import candidate_lattice, orientation_from_rows
from crystalops.builders.base import StructureBuilder
from crystalops.core import CrystalStructure, UnitCell
from crystalops.lattices import face_centered_cubic
from crystalops.operations import rotation_matrix


FloatArray = NDArray[np.float64]


_ORIENTATION_ROWS = (
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[-1, 2, 2], [2, -1, 2], [2, 2, -1]],
    [[-1, 2, -2], [2, -1, -2], [-2, -2, -1]],
    [[-1, -2, -2], [-2, -1, 2], [-2, 2, -1]],
    [[-1, -2, 2], [-2, -1, -2], [2, -2, -1]],
)


def _half_space(
    positions: FloatArray,
    point: Sequence[float],
    normal: Sequence[float],
    sense: str,
) -> NDArray[np.bool_]:
    values = (positions - np.asarray(point, dtype=np.float64)) @ np.asarray(
        normal, dtype=np.float64
    )
    if sense == "ge":
        return values >= 0.0
    if sense == "le":
        return values <= 0.0
    if sense == "gt":
        return values > 0.0
    if sense == "lt":
        return values < 0.0
    raise ValueError(f"unknown half-space sense {sense!r}")


def _middle_region(
    positions: FloatArray,
    size: FloatArray,
    upper: bool,
) -> NDArray[np.bool_]:
    """Select either of the two central polyhedra from the legacy model."""
    x, y, z = size
    if not upper:
        planes = (
            ((0, 0, z / 8), (1, 1, -1), "ge"),
            ((z, 0, 0), (1, 1, -1), "le"),
            ((x / 8, 0, 0), (-1, 1, 1), "ge"),
            ((0, z, 0), (-1, 1, 1), "le"),
            ((0, y / 8, 0), (1, -1, 1), "ge"),
            ((z, 0, 0), (1, -1, 1), "le"),
            ((0, z, z / 8), (1, 1, 1), "le"),
            ((0, 0, 0), (1, 1, 1), "ge"),
        )
    else:
        planes = (
            ((x, y, 7 * z / 8), (1, 1, -1), "le"),
            ((x, 0, z), (1, 1, -1), "ge"),
            ((7 * x / 8, y, z), (-1, 1, 1), "le"),
            ((x, y, 0), (-1, 1, 1), "ge"),
            ((x, 7 * y / 8, z), (1, -1, 1), "le"),
            ((x, y, 0), (1, -1, 1), "ge"),
            ((x, 7 * y / 8, 0), (1, 1, 1), "ge"),
            ((x, y, z), (1, 1, 1), "le"),
        )
    mask = np.ones(len(positions), dtype=bool)
    for point, normal, sense in planes:
        mask &= _half_space(positions, point, normal, sense)
    return mask


def _inside_box(positions: FloatArray, size: FloatArray) -> NDArray[np.bool_]:
    return np.all((positions >= 0.0) & (positions < size), axis=1)


def _d_schwarz_region(
    positions: FloatArray,
    size: FloatArray,
    *,
    positive: bool,
    angle: float = 0.0,
) -> NDArray[np.bool_]:
    """Evaluate the second D-Schwarz equation used by the source model."""
    radians = np.radians(angle)
    cosine, sine = np.cos(radians), np.sin(radians)
    rotation = np.array(
        [[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]]
    )
    rotated = positions @ rotation.T
    a, b, c = (rotated + 0.5).T
    x, y, z = size
    field = (
        np.sin(a * np.pi / x) * np.sin(b * np.pi / y) * np.sin(c * np.pi / z)
        + np.sin(a * np.pi / x) * np.cos(b * np.pi / y) * np.cos(c * np.pi / z)
        + np.cos(a * np.pi / x) * np.sin(b * np.pi / y) * np.cos(c * np.pi / z)
        + np.cos(a * np.pi / x) * np.cos(b * np.pi / y) * np.sin(c * np.pi / z)
    )
    return field > 0.0 if positive else field < 0.0


def _transition_region(
    positions: FloatArray,
    size: FloatArray,
    domain: int,
) -> NDArray[np.bool_]:
    """Select one of the four oriented transition regions."""
    x, y, z = size
    definitions = {
        3: (((x, 7 * y / 8, 0), (1, 1, 1), "lt"),
            ((0, z, z / 8), (1, 1, 1), "gt"), True, 0.0),
        4: (((x, 7 * y / 8, z), (1, 1, -1), "lt"),
            ((z / 8, 0, 0), (1, 1, -1), "gt"), True, 180.0),
        5: (((7 * x / 8, 0, z), (-1, 1, 1), "gt"),
            ((y / 8, 0, z), (-1, 1, 1), "lt"), False, -90.0),
        6: (((x, y / 8, 0), (1, -1, 1), "lt"),
            ((x, 7 * y / 8, 0), (1, -1, 1), "gt"), False, 90.0),
    }
    first, second, positive, angle = definitions[domain]
    mask = _half_space(positions, *first)
    mask &= _half_space(positions, *second)
    mask &= _inside_box(positions, size)
    mask &= _d_schwarz_region(positions, size, positive=positive, angle=angle)
    return mask


@dataclass(slots=True)
class SchwarzCrystal(StructureBuilder):
    """Build a mathematical DSC or PSC Schwarz crystal.

    ``schwarz_type="DSC"`` builds the existing six-domain construction with
    two central FCC regions and four oriented D-Schwarz transition regions.
    ``schwarz_type="PSC"`` creates two FCC grains separated by the P-Schwarz
    field. Grain 1 occupies ``P <= 0`` and grain 2 occupies ``P > 0``. The
    caller may provide only grain 1's orientation; grain 2 is chosen
    automatically to share the same periodic box. If omitted, grain 1 is also
    random. Coordinates are generated in lattice units and finally scaled by
    ``lattice_constant``.

    Examples:
        >>> from crystalops.builders import SchwarzCrystal
        >>> structure = SchwarzCrystal(size=8).build()
        >>> set(structure.properties["domain_ids"])
        {1, 2, 3, 4, 5, 6}
        >>> psc = SchwarzCrystal(
        ...     size=8, schwarz_type="PSC", random_seed=7
        ... ).build()
        >>> set(psc.properties["grain_ids"])
        {1, 2}
    """

    size: int = 8
    lattice_constant: float = 3.61
    element: str = "Cu"
    atomic_mass: float = 63.546
    schwarz_type: str = "DSC"
    orientation: ArrayLike | None = None
    random_seed: int | None = None
    enforce_pbc: bool = True
    max_box_expansion: float = 0.25
    max_orientation_index: int = 5
    center: tuple[float, float, float] = (0.5, 0.5, 0.5)
    name: str = "schwarz-crystal"

    def build(self) -> CrystalStructure:
        """Construct and return the selected Schwarz crystal type."""
        if isinstance(self.size, bool) or int(self.size) != self.size or self.size <= 0:
            raise ValueError("size must be a positive integer")
        if self.lattice_constant <= 0:
            raise ValueError("lattice_constant must be positive")
        if self.max_box_expansion < 0:
            raise ValueError("max_box_expansion must be non-negative")
        if self.max_orientation_index < 2:
            raise ValueError("max_orientation_index must be at least 2")

        unit_cell = face_centered_cubic(self.element)
        schwarz_type = self.schwarz_type.upper()
        if schwarz_type == "DSC":
            size = np.full(3, float(self.size))
            box = np.diag(size)
            return self._build_dsc(size, box, unit_cell)
        if schwarz_type == "PSC":
            return self._build_psc(int(self.size), unit_cell)
        raise ValueError("schwarz_type must be 'DSC' or 'PSC'")

    def _build_dsc(
        self,
        size: FloatArray,
        box: FloatArray,
        unit_cell: UnitCell,
    ) -> CrystalStructure:
        """Build the six oriented domains of the D-Schwarz crystal."""
        domain_positions: list[FloatArray] = []
        domain_ids: list[NDArray[np.int64]] = []
        domain_counts: dict[int, int] = {}
        translations = {
            4: np.array([-0.5, -0.5, 0.0]) * size,
            5: np.array([0.0, -0.5, -0.5]) * size,
            6: np.array([-0.5, 0.0, -0.5]) * size,
        }

        for domain, rows in enumerate(_ORIENTATION_ROWS, start=1):
            orientation = orientation_from_rows(rows)
            candidates = candidate_lattice(unit_cell, 1.0, box, orientation)
            candidates += size / 2.0
            if domain == 1:
                mask = _middle_region(candidates, size, upper=False)
            elif domain == 2:
                mask = _middle_region(candidates, size, upper=True)
            else:
                mask = _transition_region(candidates, size, domain)
            selected = candidates[mask]
            if domain in translations:
                selected = selected + translations[domain]
            if len(selected) == 0:
                raise RuntimeError(
                    f"D-SCI domain {domain} is empty; increase size"
                )
            domain_positions.append(selected)
            domain_ids.append(np.full(len(selected), domain, dtype=np.int64))
            domain_counts[domain] = len(selected)

        positions_lattice = np.concatenate(domain_positions)
        identifiers = np.concatenate(domain_ids)
        return CrystalStructure(
            positions=positions_lattice * self.lattice_constant,
            cell=np.diag(size * self.lattice_constant),
            atom_types=np.ones(len(positions_lattice), dtype=np.int64),
            type_elements={1: self.element},
            properties={"domain_ids": identifiers},
            metadata={
                "builder": type(self).__name__,
                "size": int(size[0]),
                "lattice_constant": float(self.lattice_constant),
                "atomic_mass": float(self.atomic_mass),
                "type_masses": {1: float(self.atomic_mass)},
                "domain_counts": domain_counts,
                "cropped_to_box": False,
                "schwarz_type": "DSC",
            },
            name=self.name,
        )

    def _build_psc(
        self,
        requested_size: int,
        unit_cell: UnitCell,
    ) -> CrystalStructure:
        """Build an FCC crystal selected by the P-Schwarz implicit field."""
        center = np.asarray(self.center, dtype=np.float64)
        if center.shape != (3,):
            raise ValueError("center must contain three fractional values")
        effective_size, orientations, orientation_sources = (
            self._psc_periodic_geometry(requested_size)
        )
        size = np.full(3, float(effective_size))
        box = np.diag(size)
        grain_positions: list[FloatArray] = []
        grain_fields: list[FloatArray] = []
        grain_ids_list: list[NDArray[np.int64]] = []
        grain_counts: dict[int, int] = {}
        for grain_id, orientation in enumerate(orientations, start=1):
            candidates = candidate_lattice(unit_cell, 1.0, box, orientation)
            candidates += size / 2.0
            candidates = candidates[_inside_box(candidates, size)]
            field = _p_schwarz_field(candidates, size, center)
            mask = field <= 0.0 if grain_id == 1 else field > 0.0
            selected = candidates[mask]
            if len(selected) == 0:
                raise RuntimeError(
                    f"PSC grain {grain_id} is empty; increase size"
                )
            grain_positions.append(selected)
            grain_fields.append(field[mask])
            grain_ids_list.append(
                np.full(len(selected), grain_id, dtype=np.int64)
            )
            grain_counts[grain_id] = len(selected)

        positions_lattice = np.concatenate(grain_positions)
        grain_ids = np.concatenate(grain_ids_list)
        fields = np.concatenate(grain_fields)

        return CrystalStructure(
            positions=positions_lattice * self.lattice_constant,
            cell=np.diag(size * self.lattice_constant),
            atom_types=np.ones(len(positions_lattice), dtype=np.int64),
            type_elements={1: self.element},
            properties={
                "grain_ids": grain_ids,
                "schwarz_field": fields,
            },
            metadata={
                "builder": type(self).__name__,
                "requested_size": requested_size,
                "effective_size": effective_size,
                "size": effective_size,
                "box_expansion": effective_size / requested_size - 1.0,
                "lattice_constant": float(self.lattice_constant),
                "atomic_mass": float(self.atomic_mass),
                "type_masses": {1: float(self.atomic_mass)},
                "schwarz_type": "PSC",
                "grain_counts": grain_counts,
                "grain_regions": {1: "P <= 0", 2: "P > 0"},
                "grain_orientations": orientations.copy(),
                "orientation_sources": orientation_sources,
                "random_seed": self.random_seed,
                "enforce_pbc": self.enforce_pbc,
                "periodicity_residuals": tuple(
                    _periodicity_residual(orientation, effective_size)
                    for orientation in orientations
                ),
                "center": center.copy(),
                "cropped_to_box": True,
            },
            name=self.name,
        )

    def _psc_periodic_geometry(
        self,
        requested_size: int,
    ) -> tuple[int, FloatArray, tuple[str, str]]:
        """Resolve two orientations and a shared commensurate cubic box."""
        rng = np.random.default_rng(self.random_seed)
        supplied: FloatArray | None = None
        if self.orientation is not None:
            matrix = np.asarray(self.orientation, dtype=np.float64)
            if matrix.shape != (3, 3):
                raise ValueError("orientation must have shape (3, 3)")
            supplied = orientation_from_rows(matrix)

        if not self.enforce_pbc:
            first = (
                rotation_matrix(*rng.uniform(0.0, 360.0, size=3))
                if supplied is None
                else supplied
            )
            second = rotation_matrix(*rng.uniform(0.0, 360.0, size=3))
            first_source = "random" if supplied is None else "provided"
            return requested_size, np.asarray([first, second]), (
                first_source,
                "random",
            )

        maximum_size = int(
            np.floor(requested_size * (1.0 + self.max_box_expansion) + 1e-12)
        )
        catalog = _commensurate_orientation_catalog(
            self.max_orientation_index
        )
        if not catalog:
            raise RuntimeError("no commensurate orientation candidates found")

        # When grain 1 is unspecified, choose it first and adjust the box to
        # its smallest allowed commensurate length. Grain 2 is then selected
        # only from orientations compatible with that same box.
        if supplied is None:
            order = rng.permutation(len(catalog))
            candidate_pairs = []
            for index in order:
                first = catalog[index]
                for box_size in range(requested_size, maximum_size + 1):
                    if _is_periodic_orientation(first, box_size):
                        candidate_pairs.append((box_size, first))
                        break
            for box_size, first in candidate_pairs:
                compatible = [
                    matrix
                    for matrix in catalog
                    if _is_periodic_orientation(matrix, box_size)
                    and not np.allclose(matrix, first, atol=1e-12)
                ]
                if compatible:
                    second = compatible[int(rng.integers(len(compatible)))]
                    return (
                        box_size,
                        np.asarray([first, second]),
                        ("random-commensurate", "random-commensurate"),
                    )
        else:
            for box_size in range(requested_size, maximum_size + 1):
                if not _is_periodic_orientation(supplied, box_size):
                    continue
                compatible = [
                    matrix
                    for matrix in catalog
                    if _is_periodic_orientation(matrix, box_size)
                    and not np.allclose(matrix, supplied, atol=1e-12)
                ]
                if compatible:
                    second = compatible[int(rng.integers(len(compatible)))]
                    return (
                        box_size,
                        np.asarray([supplied, second]),
                        ("provided", "random-commensurate"),
                    )

        raise ValueError(
            "the requested PSC orientation cannot share a periodic box "
            f"between size {requested_size} and {maximum_size}; increase "
            "max_box_expansion or provide commensurate orientations"
        )


def _commensurate_orientation_catalog(maximum_index: int) -> list[FloatArray]:
    """Return non-cubic rational rotations from primitive quaternions."""
    component_limit = int(np.sqrt(maximum_index)) + 1
    rotations: dict[tuple[float, ...], FloatArray] = {}
    values = range(-component_limit, component_limit + 1)
    for quaternion in product(values, repeat=4):
        vector = np.asarray(quaternion, dtype=np.int64)
        norm = int(vector @ vector)
        if norm < 2 or norm > maximum_index:
            continue
        nonzero = vector[vector != 0]
        if nonzero.size == 0 or nonzero[0] < 0:
            continue
        if int(np.gcd.reduce(np.abs(nonzero))) != 1:
            continue
        w, x, y, z = (int(value) for value in vector)
        numerator = np.array(
            [
                [
                    w * w + x * x - y * y - z * z,
                    2 * (x * y - w * z),
                    2 * (x * z + w * y),
                ],
                [
                    2 * (x * y + w * z),
                    w * w - x * x + y * y - z * z,
                    2 * (y * z - w * x),
                ],
                [
                    2 * (x * z - w * y),
                    2 * (y * z + w * x),
                    w * w - x * x - y * y + z * z,
                ],
            ],
            dtype=np.float64,
        )
        rotation = numerator / norm
        if np.allclose(np.abs(rotation), np.round(np.abs(rotation)), atol=1e-12):
            continue
        key = tuple(np.round(rotation, 12).ravel())
        rotations[key] = rotation
    return list(rotations.values())


def _periodicity_residual(orientation: FloatArray, box_size: int) -> float:
    """Return the largest deviation of a box translation from the lattice."""
    crystal_translations = (np.eye(3) * box_size) @ orientation
    return float(
        np.max(np.abs(crystal_translations - np.round(crystal_translations)))
    )


def _is_periodic_orientation(
    orientation: FloatArray,
    box_size: int,
    *,
    tolerance: float = 1e-10,
) -> bool:
    """Return whether a cubic box is an exact supercell of an orientation."""
    return _periodicity_residual(orientation, box_size) <= tolerance


def _p_schwarz_field(
    positions: FloatArray,
    size: FloatArray,
    center: FloatArray,
) -> FloatArray:
    """Evaluate the P-Schwarz field from the original PSC definition."""
    shifted = positions - (center - 0.5)
    x, y, z = size
    a, b, c = shifted.T
    return (
        np.cos(2.0 * a * np.pi / x)
        + np.cos(2.0 * b * np.pi / y)
        + np.cos(2.0 * c * np.pi / z)
    )


__all__ = ["SchwarzCrystal"]
