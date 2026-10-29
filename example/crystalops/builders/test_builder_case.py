"""Cases for specialized structure builders."""

import numpy as np

from crystalops import (
    FCC,
    KelvinCrystal,
    SchwarzCrystal,
    TpmsBuilder,
    WeairePhelanCrystal,
)


def test_kelvin_polycrystal_is_reproducible() -> None:
    """Build the same Kelvin structure twice with one random seed."""
    options = {
        "unit_cell": FCC("Cu"),
        "size": 2.0,
        "lattice_constant": 1.0,
        "preset": "DSC1",
        "random_seed": 42,
    }
    first = KelvinCrystal(**options).build()
    second = KelvinCrystal(**options).build()

    assert first.number_of_atoms > 0
    assert np.array_equal(
        np.unique(first.properties["grain_ids"]),
        np.arange(1, 17),
    )
    assert len(first.metadata["fixed_positions"]) == 6
    assert len(first.metadata["random_positions"]) == 10
    assert len(first.metadata["grain_counts"]) == 16
    assert first.metadata["grain_shape"] == "truncated-octahedron"
    assert first.metadata["cropped_to_box"] is False
    assert first.pbc == (True, True, True)
    assert np.allclose(
        np.diag(first.cell),
        np.full(3, 2.0 * first.metadata["grain_spacing"]),
    )

    grain_ids = first.properties["grain_ids"]
    centers = first.metadata["grain_centers"]
    relative = first.positions - centers[grain_ids - 1]
    spacing = first.metadata["grain_spacing"]
    assert np.all(np.max(np.abs(relative), axis=1) <= 0.5 * spacing + 1e-10)
    assert np.all(np.sum(np.abs(relative), axis=1) <= 0.75 * spacing + 1e-10)
    assert np.array_equal(first.atom_types, second.atom_types)
    assert np.allclose(first.positions, second.positions)


def test_build_weaire_phelan_sigma3_precursor() -> None:
    """Reproduce the eight A15 sites and Sigma3 orientations from Yang et al."""
    structure = WeairePhelanCrystal(
        unit_cell=FCC("Cu"),
        size=4.0,
        lattice_constant=1.0,
        twin_type="sigma3",
    ).build()

    assert structure.number_of_atoms > 0
    assert set(structure.properties["grain_ids"]) == set(range(1, 9))
    assert structure.metadata["precursor"] == "A15-Weaire-Phelan"
    assert structure.metadata["twin_type"] == "sigma3"
    assert structure.metadata["cropped_to_box"] is False
    assert structure.metadata["outside_box_atom_count"] > 0
    assert structure.metadata["overlap_atoms_removed"] > 0
    assert structure.metadata["initial_atom_count"] > structure.number_of_atoms
    assert structure.metadata["cell_labels"] == (
        "Dod0", "Dod1", "TetX-T1", "TetX-T2",
        "TetY-T1", "TetY-T2", "TetZ-T1", "TetZ-T2",
    )
    assert np.allclose(np.diag(structure.cell), 8.0)
    assert np.allclose(
        structure.metadata["grain_centers"],
        WeairePhelanCrystal.CENTERS,
    )
    lattice_origins = structure.metadata["lattice_origins_reduced"]
    assert np.allclose(lattice_origins[2], lattice_origins[3])
    assert np.allclose(lattice_origins[4], lattice_origins[5])
    assert np.allclose(lattice_origins[6], lattice_origins[7])
    orientations = structure.metadata["grain_orientations"]
    assert structure.metadata["orientation_sources"][:2] == (
        "paper-sigma3", "paper-sigma3"
    )
    assert structure.metadata["dod_orientation_seed"] is None
    assert structure.metadata["paper_reference_orientation"] is True
    assert np.allclose(orientations[0], np.eye(3))
    assert np.allclose(orientations[1], np.eye(3))

    # Dod0 is centered at the cell corner but remains one complete cell;
    # it therefore contains coordinates on both sides of the primary box.
    dod0 = structure.positions[structure.properties["grain_ids"] == 1]
    assert np.any(dod0 < 0.0)
    assert np.any(dod0 > 0.0)


def test_build_weaire_phelan_sigma5_precursor() -> None:
    """Build the paper's Sigma5 control with cubic Dod orientations."""
    structure = WeairePhelanCrystal(
        unit_cell=FCC("Cu"),
        size=4.0,
        lattice_constant=1.0,
        twin_type="sigma5",
    ).build()

    assert set(structure.properties["grain_ids"]) == set(range(1, 9))
    assert structure.metadata["twin_type"] == "sigma5"
    assert structure.metadata["available_twin_types"] == ("sigma3", "sigma5")
    assert structure.metadata["orientation_sources"] == (
        "paper-sigma5",
        "paper-sigma5",
        "paper-sigma5",
        "paper-sigma5",
        "paper-sigma5",
        "paper-sigma5",
        "paper-sigma5",
        "paper-sigma5",
    )
    assert np.allclose(structure.metadata["grain_orientations"][:2], np.eye(3))


def test_weaire_phelan_paper_twin_pairs_have_expected_angles() -> None:
    """Check the crystallographic disorientation of the three Tet pairs."""
    expected = {"sigma3": 60.0, "sigma5": 36.86989764584402}

    # The 24 proper signed permutation matrices are the cubic point group.
    cubic_symmetries = []
    from itertools import permutations, product

    for permutation in permutations(range(3)):
        for signs in product((-1.0, 1.0), repeat=3):
            symmetry = np.zeros((3, 3))
            symmetry[np.arange(3), permutation] = signs
            if np.linalg.det(symmetry) > 0.0:
                cubic_symmetries.append(symmetry)

    for twin_type, expected_angle in expected.items():
        orientations = [
            np.asarray(rows, dtype=float)
            / np.linalg.norm(rows, axis=1, keepdims=True)
            for rows in WeairePhelanCrystal.ORIENTATION_ROWS[twin_type]
        ]
        for first, second in ((2, 3), (4, 5), (6, 7)):
            angles = []
            for symmetry in cubic_symmetries:
                relative = orientations[second] @ symmetry @ orientations[first].T
                cosine = np.clip((np.trace(relative) - 1.0) / 2.0, -1.0, 1.0)
                angles.append(np.degrees(np.arccos(cosine)))
            assert np.isclose(min(angles), expected_angle)


def test_build_gyroid_crystal() -> None:
    """Use the gyroid field to retain one labyrinth of an FCC crystal."""
    structure = TpmsBuilder(
        unit_cell=FCC("Cu"),
        repeats=(8, 8, 8),
        lattice_constant=3.61,
        surface="gyroid",
        mode="negative",
    ).build()

    assert 0 < structure.number_of_atoms < 4 * 8**3
    assert structure.metadata["surface"] == "gyroid"
    assert structure.metadata["mode"] == "negative"


def test_build_schwarz_crystal() -> None:
    """Build all six regions of the mathematical D-SCI structure."""
    structure = SchwarzCrystal(
        size=8,
        lattice_constant=3.61,
    ).build()

    assert structure.number_of_atoms > 0
    assert set(structure.properties["domain_ids"]) == {1, 2, 3, 4, 5, 6}
    assert set(structure.metadata["domain_counts"]) == {1, 2, 3, 4, 5, 6}
    assert np.allclose(np.diag(structure.cell), np.full(3, 8 * 3.61))
    assert structure.metadata["schwarz_type"] == "DSC"
    assert structure.metadata["type_masses"] == {1: 63.546}
    assert isinstance(structure.metadata["builder"], str)
    assert structure.metadata["builder"] == "SchwarzCrystal"


def test_build_p_schwarz_crystal() -> None:
    """Build two reproducible grains separated by the P-Schwarz field."""
    options = dict(
        size=8,
        schwarz_type="PSC",
        random_seed=42,
    )
    first = SchwarzCrystal(**options).build()
    second = SchwarzCrystal(**options).build()

    grain_ids = first.properties["grain_ids"]
    field = first.properties["schwarz_field"]
    assert first.metadata["schwarz_type"] == "PSC"
    assert set(grain_ids) == {1, 2}
    assert set(first.atom_types) == {1}
    assert first.type_elements == {1: "Cu"}
    assert np.all(field[grain_ids == 1] <= 0.0)
    assert np.all(field[grain_ids == 2] > 0.0)
    assert first.metadata["orientation_sources"] == (
        "random-commensurate",
        "random-commensurate",
    )
    assert max(first.metadata["periodicity_residuals"]) <= 1e-10
    assert 8 <= first.metadata["effective_size"] <= 10
    assert np.allclose(first.positions, second.positions)


def test_build_p_schwarz_with_given_orientation() -> None:
    """Use one provided orientation and derive a periodic second grain."""
    orientation = np.array(
        [[-1, 2, 2], [2, -1, 2], [2, 2, -1]],
        dtype=float,
    )
    structure = SchwarzCrystal(
        size=8,
        schwarz_type="PSC",
        orientation=orientation,
    ).build()

    assert structure.metadata["orientation_sources"] == (
        "provided",
        "random-commensurate",
    )
    assert structure.metadata["effective_size"] == 9
    assert max(structure.metadata["periodicity_residuals"]) <= 1e-10

