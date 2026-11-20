"""Factories for common crystallographic unit cells."""

from __future__ import annotations

import numpy as np

from crystalops.core import UnitCell


def _single_element_cell(
    name: str,
    fractional_positions: list[list[float]],
    cell: list[list[float]] | np.ndarray,
    element: str,
) -> UnitCell:
    """Build a single-element unit cell from fractional basis positions."""
    return UnitCell(
        name=name,
        fractional_positions=fractional_positions,
        cell=cell,
        atom_types=np.ones(len(fractional_positions), dtype=np.int64),
        type_elements={1: element},
    )


def simple_cubic(element: str = "Cu") -> UnitCell:
    """Return a primitive simple-cubic unit cell."""
    return _single_element_cell("SC", [[0.0, 0.0, 0.0]], np.eye(3), element)


def body_centered_cubic(element: str = "Cu") -> UnitCell:
    """Return a conventional body-centered-cubic unit cell."""
    return _single_element_cell(
        "BCC",
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        np.eye(3),
        element,
    )


def face_centered_cubic(element: str = "Cu") -> UnitCell:
    """Return a conventional face-centered-cubic unit cell."""
    return _single_element_cell(
        "FCC",
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ],
        np.eye(3),
        element,
    )


def face_centered_cubic_111(element: str = "Cu") -> UnitCell:
    """Return the six-site orthogonal FCC cell aligned with the [111] axis."""
    return _single_element_cell(
        "FCC_111",
        [
            [0.0, 0.0, 0.0],
            [2.0 / 3.0, 2.0 / 3.0, 0.0],
            [1.0 / 3.0, 1.0 / 3.0, 0.0],
            [0.5, 0.0, 0.5],
            [5.0 / 6.0, 1.0 / 3.0, 0.5],
            [1.0 / 6.0, 2.0 / 3.0, 0.5],
        ],
        np.diag([np.sqrt(1.5), np.sqrt(3.0), np.sqrt(0.5)]),
        element,
    )


def hexagonal_close_packed(
    element: str = "Cu",
    *,
    c_over_a: float = np.sqrt(8.0 / 3.0),
) -> UnitCell:
    """Return a primitive hexagonal-close-packed unit cell.

    Args:
        element: Element occupying both basis sites.
        c_over_a: Ratio between the ``c`` and ``a`` lattice constants.
    """
    cell = np.array(
        [
            [1.0, 0.0, 0.0],
            [-0.5, np.sqrt(3.0) / 2.0, 0.0],
            [0.0, 0.0, c_over_a],
        ]
    )
    return _single_element_cell(
        "HCP",
        [[0.0, 0.0, 0.0], [2.0 / 3.0, 1.0 / 3.0, 0.5]],
        cell,
        element,
    )


def simple_tetragonal(element: str = "Cu", *, c_over_a: float = 1.5) -> UnitCell:
    """Return a primitive simple-tetragonal unit cell."""
    return _single_element_cell(
        "ST", [[0.0, 0.0, 0.0]], np.diag([1.0, 1.0, c_over_a]), element
    )


def body_centered_tetragonal(
    element: str = "Cu",
    *,
    c_over_a: float = 1.5,
) -> UnitCell:
    """Return a conventional body-centered-tetragonal unit cell."""
    return _single_element_cell(
        "BCT",
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        np.diag([1.0, 1.0, c_over_a]),
        element,
    )


def simple_orthorhombic(
    element: str = "Cu",
    *,
    b_over_a: float = 2.0,
    c_over_a: float = 3.0,
) -> UnitCell:
    """Return a primitive simple-orthorhombic unit cell."""
    return _single_element_cell(
        "SO", [[0.0, 0.0, 0.0]], np.diag([1.0, b_over_a, c_over_a]), element
    )


def base_centered_orthorhombic(
    element: str = "Cu",
    *,
    b_over_a: float = 2.0,
    c_over_a: float = 3.0,
) -> UnitCell:
    """Return a conventional base-centered-orthorhombic unit cell."""
    return _single_element_cell(
        "ABO",
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]],
        np.diag([1.0, b_over_a, c_over_a]),
        element,
    )


def face_centered_orthorhombic(
    element: str = "Cu",
    *,
    b_over_a: float = 2.0,
    c_over_a: float = 3.0,
) -> UnitCell:
    """Return a conventional face-centered-orthorhombic unit cell."""
    return _single_element_cell(
        "FCO",
        [
            [0.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
            [0.0, 0.5, 0.5],
        ],
        np.diag([1.0, b_over_a, c_over_a]),
        element,
    )


def body_centered_orthorhombic(
    element: str = "Cu",
    *,
    b_over_a: float = 2.0,
    c_over_a: float = 3.0,
) -> UnitCell:
    """Return a conventional body-centered-orthorhombic unit cell."""
    return _single_element_cell(
        "BCO",
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
        np.diag([1.0, b_over_a, c_over_a]),
        element,
    )


def simple_monoclinic(
    element: str = "Cu",
    *,
    b_over_a: float = 2.0,
    c_over_a: float = 3.0,
    beta: float = 100.0,
) -> UnitCell:
    """Return a primitive monoclinic unit cell with beta given in degrees."""
    angle = np.radians(beta)
    cell = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, b_over_a, 0.0],
            [c_over_a * np.cos(angle), 0.0, c_over_a * np.sin(angle)],
        ]
    )
    return _single_element_cell("SM", [[0.0, 0.0, 0.0]], cell, element)


def base_centered_monoclinic(
    element: str = "Cu",
    *,
    b_over_a: float = 2.0,
    c_over_a: float = 3.0,
    beta: float = 100.0,
) -> UnitCell:
    """Return a conventional base-centered-monoclinic unit cell."""
    base = simple_monoclinic(
        element,
        b_over_a=b_over_a,
        c_over_a=c_over_a,
        beta=beta,
    )
    return _single_element_cell(
        "BCM",
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.0]],
        base.cell,
        element,
    )


def triclinic(
    element: str = "Cu",
    *,
    lengths: tuple[float, float, float] = (1.0, 1.1, 1.2),
    angles: tuple[float, float, float] = (80.0, 85.0, 75.0),
) -> UnitCell:
    """Return a primitive triclinic unit cell.

    ``angles`` contains alpha, beta, and gamma in degrees.
    """
    a, b, c = lengths
    alpha, beta, gamma = np.radians(angles)
    cell = np.array(
        [
            [a, 0.0, 0.0],
            [b * np.cos(gamma), b * np.sin(gamma), 0.0],
            [
                c * np.cos(beta),
                c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma),
                0.0,
            ],
        ]
    )
    cell[2, 2] = np.sqrt(c**2 - cell[2, 0] ** 2 - cell[2, 1] ** 2)
    return _single_element_cell("TRI", [[0.0, 0.0, 0.0]], cell, element)


def rhombohedral(element: str = "Cu", *, alpha: float = 75.0) -> UnitCell:
    """Return a primitive rhombohedral unit cell with equal edge lengths."""
    cosine = np.cos(np.radians(alpha))
    sine = np.sin(np.radians(alpha))
    third_x = cosine
    third_y = (cosine - cosine**2) / sine
    third_z = np.sqrt(1.0 - third_x**2 - third_y**2)
    cell = np.array(
        [[1.0, 0.0, 0.0], [cosine, sine, 0.0], [third_x, third_y, third_z]]
    )
    return _single_element_cell("RHL", [[0.0, 0.0, 0.0]], cell, element)


SC = simple_cubic
BCC = body_centered_cubic
FCC = face_centered_cubic
FCC_111 = face_centered_cubic_111
HCP = hexagonal_close_packed
ST = simple_tetragonal
BCT = body_centered_tetragonal
SO = simple_orthorhombic
ABO = base_centered_orthorhombic
FCO = face_centered_orthorhombic
BCO = body_centered_orthorhombic
SM = simple_monoclinic
BCM = base_centered_monoclinic
TRI = triclinic
RHL = rhombohedral

__all__ = [
    "BCC",
    "BCM",
    "BCO",
    "BCT",
    "FCC",
    "FCC_111",
    "FCO",
    "HCP",
    "RHL",
    "SC",
    "SM",
    "SO",
    "ST",
    "TRI",
    "ABO",
    "base_centered_monoclinic",
    "base_centered_orthorhombic",
    "body_centered_cubic",
    "body_centered_orthorhombic",
    "body_centered_tetragonal",
    "face_centered_cubic",
    "face_centered_cubic_111",
    "face_centered_orthorhombic",
    "hexagonal_close_packed",
    "rhombohedral",
    "simple_cubic",
    "simple_monoclinic",
    "simple_orthorhombic",
    "simple_tetragonal",
    "triclinic",
]
