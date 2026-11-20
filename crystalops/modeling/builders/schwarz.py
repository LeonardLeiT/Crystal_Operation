"""Mathematical D- and P-Schwarz crystal builder."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from crystalops.modeling.builders._lattice import candidate_lattice, orientation_from_rows
from crystalops.modeling.builders.base import StructureBuilder
from crystalops.modeling.builders.polycrystal import _remove_periodic_overlaps
from crystalops.core import CrystalStructure, UnitCell
from crystalops.modeling.lattices import face_centered_cubic
from crystalops.modeling.operations import rotation_matrix


FloatArray = NDArray[np.float64]


_ORIENTATION_ROWS = (
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[-1, 2, 2], [2, -1, 2], [2, 2, -1]],
    [[-1, 2, -2], [2, -1, -2], [-2, -2, -1]],
    [[-1, -2, -2], [-2, -1, 2], [-2, 2, -1]],
    [[-1, -2, 2], [-2, -1, -2], [2, -2, -1]],
)

# All four DSC twin orientations repeat after three conventional lattice
# spacings along each Cartesian box axis.
_DSC_REPEAT = 3
_DSC_TWIN_NORMALS = np.array(
    [[1, 1, 1], [-1, -1, 1], [1, -1, -1], [-1, 1, -1]], dtype=np.float64
)
_DSC_TWIN_ANGLES = (0.0, 180.0, -90.0, 90.0)
_DSC_TWIN_SIGNS = (1, 1, -1, -1)


def _dsc_translations(size: FloatArray) -> dict[int, FloatArray]:
    return {
        3: np.zeros(3),
        4: np.array([-0.5, -0.5, 0.0]) * size,
        5: np.array([0.0, -0.5, -0.5]) * size,
        6: np.array([-0.5, 0.0, -0.5]) * size,
    }


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


def _d_schwarz_field(
    positions: FloatArray,
    size: FloatArray,
    *,
    angle: float = 0.0,
) -> FloatArray:
    """Evaluate the second D-Schwarz equation used by the source model."""
    radians = np.radians(angle)
    cosine, sine = np.cos(radians), np.sin(radians)
    rotation = np.array(
        [[cosine, -sine, 0.0], [sine, cosine, 0.0], [0.0, 0.0, 1.0]]
    )
    rotated = positions @ rotation.T
    a, b, c = (rotated + 0.5).T
    x, y, z = size
    return (
        np.sin(a * np.pi / x) * np.sin(b * np.pi / y) * np.sin(c * np.pi / z)
        + np.sin(a * np.pi / x) * np.cos(b * np.pi / y) * np.cos(c * np.pi / z)
        + np.cos(a * np.pi / x) * np.sin(b * np.pi / y) * np.cos(c * np.pi / z)
        + np.cos(a * np.pi / x) * np.cos(b * np.pi / y) * np.sin(c * np.pi / z)
    )


def _dsc_transition_planes(size: FloatArray, domain: int):
    """Return the two legacy planar cuts in the local chart."""
    x, y, z = size
    definitions = {
        3: (((x, 7 * y / 8, 0), (1, 1, 1), "lt"),
            ((0, z, z / 8), (1, 1, 1), "gt")),
        4: (((x, 7 * y / 8, z), (1, 1, -1), "lt"),
            ((z / 8, 0, 0), (1, 1, -1), "gt")),
        5: (((7 * x / 8, 0, z), (-1, 1, 1), "gt"),
            ((y / 8, 0, z), (-1, 1, 1), "lt")),
        6: (((x, y / 8, 0), (1, -1, 1), "lt"),
            ((x, 7 * y / 8, 0), (1, -1, 1), "gt")),
    }
    return definitions[domain]


def _transition_region(
    positions: FloatArray, size: FloatArray, domain: int
) -> NDArray[np.bool_]:
    """Select one of the four oriented transition regions."""
    first, second = _dsc_transition_planes(size, domain)
    mask = _half_space(positions, *first)
    mask &= _half_space(positions, *second)
    mask &= _inside_box(positions, size)
    index = domain - 3
    field = _d_schwarz_field(positions, size, angle=_DSC_TWIN_ANGLES[index])
    mask &= _DSC_TWIN_SIGNS[index] * field > 0.0
    return mask


def _dsc_twin_layers(desired: FloatArray) -> NDArray[np.int64]:
    """Snap two faces to parent layers with the same FCC stacking registry."""
    near = np.floor(desired[:, None]).astype(np.int64) + np.arange(-2, 4)
    low, high = np.meshgrid(*near, indexing="ij")
    allowed = (high > low) & ((high - low) % 3 == 0)
    low, high = low[allowed], high[allowed]
    error_low, error_high = low - desired[0], high - desired[1]
    closest = np.lexsort((
        low,
        np.abs(error_low + error_high),
        np.maximum(np.abs(error_low), np.abs(error_high)),
    ))[0]
    return np.array([low[closest], high[closest]])


def _dsc_projection_limits(projection: FloatArray) -> FloatArray:
    """Bound x, y, z by adding pairs of twin-normal projections."""
    return (projection[0] + projection[[2, 3, 1]]) / 2.0


def _dsc_in_parent(
    projections: FloatArray, low: FloatArray, high: FloatArray
) -> NDArray[np.bool_]:
    return np.all(
        (projections >= low - 1e-10) & (projections <= high + 1e-10), axis=1
    )


@dataclass(slots=True)
class _DscPartition:
    """Assign a single DSC domain at each point of the periodic cell."""

    size: FloatArray
    reference_size: FloatArray
    offsets: FloatArray
    low: FloatArray
    high: FloatArray
    lower_back: FloatArray
    upper_low: FloatArray
    upper_back: FloatArray

    def assign(
        self, points: FloatArray
    ) -> tuple[NDArray[np.int8], NDArray[np.bool_], NDArray[np.bool_]]:
        normals = _DSC_TWIN_NORMALS
        points = points % self.size
        owners = np.zeros(len(points), dtype=np.int8)
        projected = points @ normals.T
        for label, low, high in (
            (1, self.lower_back, self.low),
            (2, self.upper_low, self.upper_back),
        ):
            low_image = np.floor(_dsc_projection_limits(low) / self.size).astype(int)
            high_image = np.floor(_dsc_projection_limits(high) / self.size).astype(int)
            for image in product(*(range(a, b + 1) for a, b in zip(low_image, high_image))):
                shift = (self.size * image) @ normals.T
                in_parent = _dsc_in_parent(projected, low - shift, high - shift)
                owners[in_parent & (owners == 0)] = label

        outside = np.flatnonzero(owners == 0)
        eligible = np.zeros((len(points), 4), dtype=bool)
        best = np.full(len(outside), -np.inf)
        collar_owned = np.zeros(len(outside), dtype=bool)
        aspect = self.size / self.reference_size
        for index, (normal, offset) in enumerate(zip(normals, self.offsets)):
            raw = (points[outside] - offset) % self.size
            chart = raw + offset
            projected = chart @ normal
            from_low = projected - self.low[index]
            from_high = self.high[index] - projected
            in_slab = (from_low > 1e-10) & (from_high > 1e-10)
            eligible[outside, index] = in_slab
            score = _d_schwarz_field(
                raw / aspect, self.reference_size, angle=_DSC_TWIN_ANGLES[index]
            )
            score *= _DSC_TWIN_SIGNS[index]

            # The first two layers next to either parent facet should use
            # that facet's own mirror lattice even when slivers overlap.
            collar = np.zeros(len(outside), dtype=bool)
            near_low = in_slab & (from_low <= 2.0 + 1e-10)
            on_low = chart[near_low] - from_low[near_low, None] * normal / 3.0
            collar[near_low] = _dsc_in_parent(on_low @ normals.T, self.lower_back, self.low)
            near_high = in_slab & (from_high <= 2.0 + 1e-10)
            on_high = chart[near_high] + from_high[near_high, None] * normal / 3.0
            collar[near_high] |= _dsc_in_parent(
                (on_high - 2.0 * offset) @ normals.T, self.upper_low, self.upper_back
            )

            choose = in_slab & (
                (collar & ~collar_owned) | ((collar == collar_owned) & (score > best))
            )
            owners[outside[choose]] = index + 3
            best[choose] = score[choose]
            collar_owned[choose] = collar[choose]

        protected = (owners == 1) | (owners == 2)
        protected[outside] = collar_owned
        return owners, eligible, protected

    def refine(self, domains: list[FloatArray]) -> tuple[list[FloatArray], int]:
        """Keep exact atom assignments, with periodic corrections for islands."""
        from ._dsc_connectivity import reassign_dsc_islands

        grid_size = 96
        # Offset the samples from voxel centers so narrow junctions on
        # symmetry planes do not inherit an artificial one-voxel bridge.
        phase = 1.0 / 3.0
        grid_shape = (grid_size,) * 3
        grid = (np.indices(grid_shape).reshape(3, -1).T + phase) * (self.size / grid_size)
        owners, eligible, protected = self.assign(grid)
        sampled_owners = owners.reshape(grid_shape)
        repaired, moved = reassign_dsc_islands(
            sampled_owners,
            eligible.reshape((*grid_shape, 4)),
            protected.reshape(grid_shape),
        )
        refined = []
        for domain, points in enumerate(domains, start=1):
            exact, eligible, protected = self.assign(points)
            reduced = (points % self.size) * (grid_size / self.size)
            index = np.floor(reduced + (0.5 - phase)).astype(int)
            index = tuple((index % grid_size).T)
            sampled = repaired[index]
            changed = (sampled_owners[index] == exact) & (sampled != exact)
            allowed = eligible[np.arange(len(points)), (sampled - 3).clip(0, 3)]
            replace = changed & allowed & ~protected & (sampled >= 3)
            exact[replace] = sampled[replace]
            refined.append(points[exact == domain])
        return refined, moved


def _dsc_coherent_geometry(
    size: FloatArray,
    reference_size: FloatArray,
    translations: dict[int, FloatArray],
) -> tuple[_DscPartition, NDArray[np.int64]]:
    """Locate the shared {111} cuts for a rectangular DSC."""
    normals = _DSC_TWIN_NORMALS
    parent_origin = size / 2.0
    center_lower, center_upper = size / 4.0, 3.0 * size / 4.0
    half_width = size.sum() / 8.0
    ideal_low = normals @ (center_lower - parent_origin) + half_width
    layers = np.stack([
        _dsc_twin_layers(np.array([level, level + 2.0 * half_width]))
        for level in ideal_low
    ])
    low, high = (layers + (normals @ parent_origin)[:, None]).T
    # Each upper central facet is the second face of a translated slab.
    offsets = np.stack([translations[domain] for domain in range(3, 7)])
    upper_low = high - 2.0 * np.einsum("ij,ij->i", normals, offsets)
    lower_back = normals @ center_lower - 2.0 * half_width
    upper_back = normals @ center_upper + 2.0 * half_width
    partition = _DscPartition(size, reference_size, offsets, low, high, lower_back, upper_low, upper_back)
    return partition, layers


def _dsc_coherent_domains(
    unit_cell: UnitCell,
    size: FloatArray,
    reference_size: FloatArray,
    translations: dict[int, FloatArray],
    *,
    refine_boundaries: bool = False,
) -> tuple[list[FloatArray], NDArray[np.int64], int]:
    """Place the rectangular DSC facets on shared parent/twin {111} layers."""
    partition, layers = _dsc_coherent_geometry(size, reference_size, translations)
    normals, parent_origin = _DSC_TWIN_NORMALS, size / 2.0
    low, high = partition.low, partition.high
    lower_back, upper_low, upper_back = (
        partition.lower_back, partition.upper_low, partition.upper_back
    )
    extent = np.maximum(
        parent_origin - _dsc_projection_limits(np.minimum(lower_back, upper_low)),
        _dsc_projection_limits(np.maximum(low, upper_back)) - parent_origin,
    )
    parent = candidate_lattice(unit_cell, 1.0, np.diag(2.0 * extent), np.eye(3))
    parent += parent_origin
    projection = parent @ normals.T
    domains = [
        parent[np.all((projection >= lower_back) & (projection <= low), axis=1)],
        parent[np.all((projection >= upper_low) & (projection <= upper_back), axis=1)],
    ]
    aspect = size / reference_size
    for index, (normal, offset) in enumerate(zip(normals, partition.offsets)):
        rotation = orientation_from_rows(_ORIENTATION_ROWS[index + 2])
        # The twin lattice is the mirror image about the lower {111} face.
        # The upper face shares its registry because the layers differ by 3n.
        mirror_shift = (2.0 / 3.0) * layers[index, 0] * normal
        displacement = mirror_shift - offset
        displacement -= _DSC_REPEAT * np.rint(displacement / _DSC_REPEAT)
        box = np.diag(size + 2.0 * np.abs(displacement))
        candidates = candidate_lattice(unit_cell, 1.0, box, rotation)
        candidates += parent_origin + displacement
        candidates[np.abs(candidates) < 1e-10] = 0.0
        candidates = candidates[_inside_box(candidates, size - 1e-10)]
        projected = (candidates + offset) @ normal
        between = (projected > low[index] + 1e-10) & (projected < high[index] - 1e-10)
        candidates = candidates[between]
        if not refine_boundaries:
            field = _d_schwarz_field(
                candidates / aspect, reference_size, angle=_DSC_TWIN_ANGLES[index]
            )
            candidates = candidates[_DSC_TWIN_SIGNS[index] * field > 0.0]
        domains.append(candidates + offset)
    moved = 0
    if refine_boundaries:
        domains, moved = partition.refine(domains)
    return domains, layers, moved


@dataclass(slots=True)
class SchwarzCrystal(StructureBuilder):
    """Build a mathematical DSC or PSC Schwarz crystal.

    ``schwarz_type="DSC"`` builds the existing six-domain construction with
    two central FCC regions and four Σ3 twin variants relative to the
    center. Set ``z_ratio`` to build the DSC in a periodic rectangular box.
    The noncubic geometry adjusts the physical {111} twin facets and the
    FCC stacking registry for each box. For noncubic DSC, ``refine_boundaries``
    resolves twin-region overlap and reconnects detached regions while
    preserving the first two atomic layers next to the mother grain.
    This option requires SciPy. In the periodic DSC mode,
    ``size`` and ``size * z_ratio`` must be multiples of 3.
    ``schwarz_type="PSC"`` creates two FCC grains separated by the P-Schwarz
    field. Grain 1 occupies ``P <= 0`` and grain 2 occupies ``P > 0``. The
    caller may provide only grain 1's orientation; grain 2 is chosen
    automatically to share the same periodic box. If omitted, grain 1 is also
    random. Coordinates are generated in lattice units and finally scaled by
    ``lattice_constant``.

    Examples:
        >>> from crystalops.modeling.builders import SchwarzCrystal
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
    z_ratio: float | None = None
    refine_boundaries: bool = False

    def build(self) -> CrystalStructure:
        """Construct and return the selected Schwarz crystal type."""
        if isinstance(self.size, bool) or int(self.size) != self.size or self.size <= 0:
            raise ValueError("size must be a positive integer")
        if not np.isfinite(self.lattice_constant) or self.lattice_constant <= 0:
            raise ValueError("lattice_constant must be finite and positive")
        if self.max_box_expansion < 0:
            raise ValueError("max_box_expansion must be non-negative")
        if self.max_orientation_index < 2:
            raise ValueError("max_orientation_index must be at least 2")

        unit_cell = face_centered_cubic(self.element)
        schwarz_type = self.schwarz_type.upper()
        if schwarz_type == "DSC":
            return self._build_dsc(unit_cell)
        if self.z_ratio is not None or self.refine_boundaries:
            raise ValueError("z_ratio and refine_boundaries require schwarz_type='DSC'")
        if schwarz_type == "PSC":
            return self._build_psc(int(self.size), unit_cell)
        raise ValueError("schwarz_type must be 'DSC' or 'PSC'")

    def _build_dsc(self, unit_cell: UnitCell) -> CrystalStructure:
        """Build the six legacy domains, optionally in a rectangular period cell."""
        reference_size = np.full(3, float(self.size))
        size = reference_size.copy()
        periodic = self.z_ratio is not None
        if self.z_ratio is not None:
            z_cells = self.size * self.z_ratio
            if not np.isfinite(z_cells) or z_cells <= 0:
                raise ValueError("z_ratio must be finite and positive")
            repeat = _DSC_REPEAT
            z_periods = round(z_cells / repeat)
            if (
                self.size % repeat
                or z_periods < 1
                or not np.isclose(z_cells, repeat * z_periods, rtol=0.0, atol=1e-10)
            ):
                raise ValueError(
                    "Periodic DSC requires size and size * z_ratio to be positive "
                    "multiples of 3; use size=30 for the six example ratios"
                )
            size[2] = repeat * z_periods
        translations = _dsc_translations(size)

        twin_layers = None
        relocated_cells = 0
        if self.refine_boundaries and size[2] == size[0]:
            raise ValueError("refine_boundaries requires a noncubic DSC box")
        if size[2] != size[0]:
            domain_positions, twin_layers, relocated_cells = _dsc_coherent_domains(
                unit_cell, size, reference_size, translations,
                refine_boundaries=self.refine_boundaries,
            )
        else:
            box = np.diag(size)
            domain_positions = []
            for domain, rows in enumerate(_ORIENTATION_ROWS, start=1):
                orientation = orientation_from_rows(rows)
                candidates = candidate_lattice(unit_cell, 1.0, box, orientation)
                candidates += size / 2.0
                if domain in (1, 2):
                    mask = _middle_region(candidates, size, upper=domain == 2)
                    selected = candidates[mask]
                else:
                    mask = _transition_region(candidates, size, domain)
                    selected = candidates[mask] + translations[domain]
                domain_positions.append(selected)

        domain_counts = {
            domain: len(positions) for domain, positions in enumerate(domain_positions, start=1)
        }
        if not all(domain_counts.values()):
            raise RuntimeError("a DSC domain is empty; increase size")

        positions_lattice = np.concatenate(domain_positions)
        identifiers = np.repeat(np.arange(1, 7), list(domain_counts.values()))
        properties = {"domain_ids": identifiers}
        if periodic:
            positions_lattice %= size
            near_face = np.isclose(positions_lattice, size, rtol=0.0, atol=1e-10)
            positions_lattice[near_face] = 0.0
            properties["grain_ids"] = np.maximum(1, identifiers - 1)
        structure = CrystalStructure(
            positions=positions_lattice * self.lattice_constant,
            cell=np.diag(size * self.lattice_constant),
            atom_types=np.ones(len(positions_lattice), dtype=np.int64),
            type_elements={1: self.element},
            properties=properties,
            metadata={
                "builder": type(self).__name__,
                "size": int(size[0]),
                "lattice_constant": float(self.lattice_constant),
                "atomic_mass": float(self.atomic_mass),
                "type_masses": {1: float(self.atomic_mass)},
                "domain_counts": domain_counts,
                "cropped_to_box": periodic,
                "schwarz_type": "DSC",
            },
            name=self.name,
        )
        if periodic:
            cutoff = self.lattice_constant / (2.0 * np.sqrt(2.0))
            structure, removed = _remove_periodic_overlaps(structure, cutoff)
            remaining = np.bincount(structure.properties["domain_ids"], minlength=7)
            if not np.all(remaining[1:]):
                raise RuntimeError("DSC overlap removal emptied a domain; increase size")
            structure.metadata.update(
                z_ratio=float(size[2] / size[0]),
                domain_counts={domain: int(remaining[domain]) for domain in range(1, 7)},
                overlap_cutoff=float(cutoff),
                removed_overlap_atoms=removed,
            )
        if twin_layers is not None:
            structure.metadata["twin_plane_layers"] = twin_layers
        if self.refine_boundaries:
            structure.metadata["reassigned_grid_cells"] = relocated_cells
        return structure

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
