"""Cases for unit cells and the canonical crystal structure model."""

import numpy as np

from crystalops.modeling import BCC, CrystalBuilder, FCC


def test_build_fcc_supercell() -> None:
    """Build a 3 x 2 x 1 FCC copper supercell."""
    structure = CrystalBuilder(
        unit_cell=FCC("Cu"),
        repeats=(3, 2, 1),
        lattice_constant=3.61,
    ).build()

    assert structure.number_of_atoms == 24
    assert np.allclose(np.diag(structure.cell), [10.83, 7.22, 3.61])
    assert structure.type_elements == {1: "Cu"}
    assert np.isclose(structure.volume, 10.83 * 7.22 * 3.61)


def test_build_bcc_supercell() -> None:
    """Build a 2 x 2 x 2 BCC iron supercell."""
    structure = CrystalBuilder(
        unit_cell=BCC("Fe"),
        repeats=(2, 2, 2),
        lattice_constant=2.86,
    ).build()

    assert structure.number_of_atoms == 16
    assert structure.type_elements == {1: "Fe"}
    assert np.all(structure.fractional_positions >= 0.0)
    assert np.all(structure.fractional_positions < 1.0)
