"""Cases for reading and writing LAMMPS data files."""

import numpy as np

from crystalops.core import CrystalStructure
from crystalops.modeling import CrystalBuilder, FCC
from crystalops.modeling.operations import assign_atom_type


def test_lammps_data_round_trip(tmp_path) -> None:
    """Write a binary alloy and read it back into the common structure class."""
    structure = CrystalBuilder(FCC("Cu"), (2, 2, 2), 3.61).build()
    nickel_sites = np.arange(structure.number_of_atoms) % 4 == 0
    structure = assign_atom_type(structure, nickel_sites, 2, "Ni")

    output = tmp_path / "cu_ni.data"
    structure.write(output)
    loaded = CrystalStructure.read(output)

    assert loaded.number_of_atoms == structure.number_of_atoms
    assert loaded.type_elements == {1: "Cu", 2: "Ni"}
    assert np.array_equal(loaded.atom_types, structure.atom_types)
    assert np.allclose(loaded.positions, structure.positions)
    assert np.allclose(loaded.cell, structure.cell)


def test_triclinic_lammps_data_round_trip(tmp_path) -> None:
    """Preserve a restricted triclinic box and non-zero origin."""
    structure = CrystalStructure.from_fractional(
        fractional_positions=[[0.1, 0.2, 0.3], [0.7, 0.5, 0.9]],
        cell=[[5.0, 0.0, 0.0], [1.2, 6.0, 0.0], [0.4, 0.8, 7.0]],
        atom_types=[1, 1],
        type_elements={1: "Cu"},
        origin=[-2.0, 1.0, 3.0],
        name="triclinic-case",
    )

    output = tmp_path / "triclinic.data"
    structure.write(output)
    loaded = CrystalStructure.read(output)

    assert np.allclose(loaded.origin, structure.origin)
    assert np.allclose(loaded.cell, structure.cell)
    assert np.allclose(loaded.positions, structure.positions)
