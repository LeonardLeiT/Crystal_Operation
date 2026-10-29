"""Geometry and atom-set operations."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from numpy.typing import ArrayLike, NDArray

from crystalops.core import CrystalStructure


BoolArray = NDArray[np.bool_]
FloatArray = NDArray[np.float64]


def _clone(
    structure: CrystalStructure,
    *,
    positions: ArrayLike | None = None,
    cell: ArrayLike | None = None,
    atom_types: ArrayLike | None = None,
    type_elements: dict[int, str] | None = None,
    origin: ArrayLike | None = None,
    properties: dict[str, ArrayLike] | None = None,
    name: str | None = None,
) -> CrystalStructure:
    """Create a changed copy while retaining unspecified structure fields."""
    return CrystalStructure(
        positions=structure.positions if positions is None else positions,
        cell=structure.cell if cell is None else cell,
        atom_types=structure.atom_types if atom_types is None else atom_types,
        type_elements=(
            structure.type_elements if type_elements is None else type_elements
        ),
        origin=structure.origin if origin is None else origin,
        pbc=structure.pbc,
        properties=structure.properties if properties is None else properties,
        metadata=structure.metadata,
        name=structure.name if name is None else name,
    )


def rotation_matrix(
    x_angle: float = 0.0,
    y_angle: float = 0.0,
    z_angle: float = 0.0,
    *,
    degrees: bool = True,
) -> FloatArray:
    """Return an active XYZ Euler rotation matrix.

    Args:
        x_angle: Rotation around the Cartesian x axis.
        y_angle: Rotation around the Cartesian y axis.
        z_angle: Rotation around the Cartesian z axis.
        degrees: Interpret angles as degrees when true.

    Returns:
        A ``(3, 3)`` matrix applying x, then y, then z rotation.
    """
    angles = np.array([x_angle, y_angle, z_angle], dtype=np.float64)
    if degrees:
        angles = np.radians(angles)
    x, y, z = angles
    rx = np.array(
        [[1.0, 0.0, 0.0], [0.0, np.cos(x), -np.sin(x)], [0.0, np.sin(x), np.cos(x)]]
    )
    ry = np.array(
        [[np.cos(y), 0.0, np.sin(y)], [0.0, 1.0, 0.0], [-np.sin(y), 0.0, np.cos(y)]]
    )
    rz = np.array(
        [[np.cos(z), -np.sin(z), 0.0], [np.sin(z), np.cos(z), 0.0], [0.0, 0.0, 1.0]]
    )
    return rz @ ry @ rx


def translate(structure: CrystalStructure, vector: ArrayLike) -> CrystalStructure:
    """Translate all atoms without moving the simulation cell.

    Args:
        structure: Source structure.
        vector: Cartesian translation vector.

    Returns:
        A translated copy.
    """
    displacement = np.asarray(vector, dtype=np.float64)
    if displacement.shape != (3,):
        raise ValueError("vector must have shape (3,)")
    return _clone(structure, positions=structure.positions + displacement)


def rotate(
    structure: CrystalStructure,
    rotation: ArrayLike,
    *,
    center: ArrayLike | None = None,
    rotate_cell: bool = False,
) -> CrystalStructure:
    """Rotate atom positions around a Cartesian point.

    Args:
        structure: Source structure.
        rotation: Orthogonal ``(3, 3)`` active rotation matrix.
        center: Cartesian rotation center. Defaults to the cell center.
        rotate_cell: Rotate the cell vectors and origin as a rigid body.

    Returns:
        A rotated copy.
    """
    matrix = np.asarray(rotation, dtype=np.float64)
    if matrix.shape != (3, 3):
        raise ValueError("rotation must have shape (3, 3)")
    if not np.allclose(matrix @ matrix.T, np.eye(3), atol=1e-10):
        raise ValueError("rotation must be orthogonal")

    pivot = (
        structure.origin + 0.5 * np.sum(structure.cell, axis=0)
        if center is None
        else np.asarray(center, dtype=np.float64)
    )
    if pivot.shape != (3,):
        raise ValueError("center must have shape (3,)")

    positions = (structure.positions - pivot) @ matrix.T + pivot
    if not rotate_cell:
        return _clone(structure, positions=positions)

    cell = structure.cell @ matrix.T
    origin = (structure.origin - pivot) @ matrix.T + pivot
    return _clone(structure, positions=positions, cell=cell, origin=origin)


def repeat(
    structure: CrystalStructure,
    repeats: tuple[int, int, int],
) -> CrystalStructure:
    """Construct a supercell by repeating a structure.

    Args:
        structure: Source structure.
        repeats: Positive repeat counts along the three cell vectors.

    Returns:
        The repeated supercell.
    """
    counts = np.asarray(repeats, dtype=np.int64)
    if counts.shape != (3,) or np.any(counts < 1):
        raise ValueError("repeats must contain three positive integers")

    indices = np.stack(
        np.meshgrid(
            np.arange(counts[0]),
            np.arange(counts[1]),
            np.arange(counts[2]),
            indexing="ij",
        ),
        axis=-1,
    ).reshape(-1, 3)
    translations = indices @ structure.cell
    positions = (
        structure.positions[None, :, :] + translations[:, None, :]
    ).reshape(-1, 3)
    copies = len(translations)
    atom_types = np.tile(structure.atom_types, copies)
    properties = {
        name: np.concatenate([values] * copies, axis=0)
        for name, values in structure.properties.items()
    }
    if "atom_ids" in properties:
        properties["atom_ids"] = np.arange(1, len(positions) + 1, dtype=np.int64)
    cell = structure.cell * counts[:, None]
    return _clone(
        structure,
        positions=positions,
        cell=cell,
        atom_types=atom_types,
        properties=properties,
    )


def select(structure: CrystalStructure, mask: ArrayLike) -> CrystalStructure:
    """Return only atoms selected by a Boolean mask or integer indices."""
    selection = np.asarray(mask)
    if selection.dtype == np.bool_ and selection.shape != (structure.number_of_atoms,):
        raise ValueError("Boolean mask must contain one value per atom")
    properties = {
        name: values[selection]
        for name, values in structure.properties.items()
    }
    return _clone(
        structure,
        positions=structure.positions[selection],
        atom_types=structure.atom_types[selection],
        properties=properties,
    )


def select_sphere(
    structure: CrystalStructure,
    center: ArrayLike,
    radius: float,
    *,
    inside: bool = True,
) -> CrystalStructure:
    """Select atoms inside or outside a Cartesian sphere."""
    if radius < 0:
        raise ValueError("radius must be non-negative")
    sphere_center = np.asarray(center, dtype=np.float64)
    if sphere_center.shape != (3,):
        raise ValueError("center must have shape (3,)")
    distances = np.linalg.norm(structure.positions - sphere_center, axis=1)
    mask = distances <= radius
    return select(structure, mask if inside else ~mask)


def slice_by_planes(
    structure: CrystalStructure,
    points: ArrayLike,
    normals: ArrayLike,
    *,
    tolerance: float = 1e-10,
) -> CrystalStructure:
    """Keep atoms on the non-positive side of every oriented plane.

    Args:
        structure: Source structure.
        points: One Cartesian point on each plane, shape ``(M, 3)``.
        normals: Outward normal of each plane, shape ``(M, 3)``.
        tolerance: Positive numerical allowance at plane boundaries.

    Returns:
        Structure containing the atoms inside all half-spaces.
    """
    plane_points = np.atleast_2d(np.asarray(points, dtype=np.float64))
    plane_normals = np.atleast_2d(np.asarray(normals, dtype=np.float64))
    if plane_points.shape != plane_normals.shape or plane_points.shape[1:] != (3,):
        raise ValueError("points and normals must have matching shape (M, 3)")
    signed = np.einsum(
        "nmk,mk->nm",
        structure.positions[:, None, :] - plane_points[None, :, :],
        plane_normals,
    )
    return select(structure, np.all(signed <= tolerance, axis=1))


def wrap(structure: CrystalStructure) -> CrystalStructure:
    """Wrap periodic fractional coordinates into the primary simulation cell."""
    fractional = structure.fractional_positions
    for axis, periodic in enumerate(structure.pbc):
        if periodic:
            fractional[:, axis] %= 1.0
    positions = fractional @ structure.cell + structure.origin
    return _clone(structure, positions=positions)


def assign_atom_type(
    structure: CrystalStructure,
    mask: ArrayLike,
    type_id: int,
    element: str,
) -> CrystalStructure:
    """Assign an atom type and element mapping to selected atoms."""
    if type_id < 1:
        raise ValueError("type_id must be positive")
    selection = np.asarray(mask)
    if selection.dtype == np.bool_ and selection.shape != (structure.number_of_atoms,):
        raise ValueError("Boolean mask must contain one value per atom")
    atom_types = structure.atom_types.copy()
    atom_types[selection] = type_id
    type_elements = dict(structure.type_elements)
    type_elements[type_id] = element
    return _clone(
        structure,
        atom_types=atom_types,
        type_elements=type_elements,
    )


def merge(structures: Iterable[CrystalStructure], *, name: str = "merged") -> CrystalStructure:
    """Merge structures that share the same simulation cell and properties."""
    items = list(structures)
    if not items:
        raise ValueError("at least one structure is required")
    first = items[0]
    for item in items[1:]:
        if not np.allclose(item.cell, first.cell) or not np.allclose(
            item.origin,
            first.origin,
        ):
            raise ValueError("all structures must use the same cell and origin")
        if item.properties.keys() != first.properties.keys():
            raise ValueError("all structures must define the same properties")

    type_elements: dict[int, str] = {}
    for item in items:
        for type_id, element in item.type_elements.items():
            if type_id in type_elements and type_elements[type_id] != element:
                raise ValueError(f"conflicting element mapping for type {type_id}")
            type_elements[type_id] = element

    properties = {
        property_name: np.concatenate(
            [item.properties[property_name] for item in items],
            axis=0,
        )
        for property_name in first.properties
    }
    if "atom_ids" in properties:
        properties["atom_ids"] = np.arange(
            1,
            sum(item.number_of_atoms for item in items) + 1,
            dtype=np.int64,
        )
    return CrystalStructure(
        positions=np.concatenate([item.positions for item in items], axis=0),
        cell=first.cell,
        atom_types=np.concatenate([item.atom_types for item in items]),
        type_elements=type_elements,
        origin=first.origin,
        pbc=first.pbc,
        properties=properties,
        metadata=first.metadata,
        name=name,
    )


def remove_duplicates(
    structure: CrystalStructure,
    *,
    tolerance: float = 1e-8,
) -> CrystalStructure:
    """Remove later atoms whose Cartesian positions are equal within tolerance."""
    if tolerance <= 0:
        raise ValueError("tolerance must be positive")
    keys = np.round(structure.positions / tolerance).astype(np.int64)
    _, indices = np.unique(keys, axis=0, return_index=True)
    return select(structure, np.sort(indices))
