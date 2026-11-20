"""Cases for immutable structure operations."""

import numpy as np

from crystalops.modeling import CrystalBuilder, FCC
from crystalops.modeling.operations import (
    remove_duplicates,
    rotation_matrix,
    select_sphere,
    wrap,
)


def test_transform_and_wrap_structure() -> None:
    """Translate, rotate, and wrap atoms without changing the source object."""
    source = CrystalBuilder(FCC("Cu"), (2, 2, 2), 3.61).build()
    transformed = source.translated([8.0, -1.0, 0.5])
    transformed = transformed.rotated(
        rotation_matrix(z_angle=30.0),
        center=transformed.origin,
    )
    transformed = wrap(transformed)

    assert transformed is not source
    assert not np.shares_memory(transformed.positions, source.positions)
    periodic_axes = np.asarray(transformed.pbc)
    fractional = transformed.fractional_positions[:, periodic_axes]
    assert np.all(fractional >= 0.0)
    assert np.all(fractional < 1.0)


def test_select_and_remove_duplicate_atoms() -> None:
    """Select a spherical region and remove an intentionally duplicated atom."""
    structure = CrystalBuilder(FCC("Cu"), (3, 3, 3), 3.61).build()
    center = structure.origin + 0.5 * structure.cell.sum(axis=0)
    selected = select_sphere(structure, center, radius=4.0)

    duplicated = selected.copy()
    duplicated.positions = np.vstack([selected.positions, selected.positions[:1]])
    duplicated.atom_types = np.concatenate([selected.atom_types, selected.atom_types[:1]])
    duplicated.validate()
    cleaned = remove_duplicates(duplicated)

    assert 0 < selected.number_of_atoms < structure.number_of_atoms
    assert cleaned.number_of_atoms == selected.number_of_atoms
