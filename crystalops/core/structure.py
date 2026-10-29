"""Canonical atomistic structure containers."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def _copy_properties(
    properties: Mapping[str, ArrayLike] | None,
) -> dict[str, NDArray[Any]]:
    """Return owned NumPy arrays for a per-atom property mapping."""
    if properties is None:
        return {}
    return {name: np.asarray(values).copy() for name, values in properties.items()}


@dataclass(slots=True)
class CrystalStructure:
    """Store a complete particle-level crystal structure.

    Positions always use Cartesian coordinates. The three rows of ``cell``
    are the simulation-cell vectors, and atom type identifiers start at one
    to match LAMMPS conventions. Every entry in ``properties`` must have one
    value per atom along its first dimension.

    Examples:
        >>> import numpy as np
        >>> structure = CrystalStructure(
        ...     positions=np.array([[0.0, 0.0, 0.0]]),
        ...     cell=np.eye(3) * 3.61,
        ...     atom_types=np.array([1]),
        ...     type_elements={1: "Cu"},
        ... )
        >>> structure.number_of_atoms
        1
        >>> structure.fractional_positions.tolist()
        [[0.0, 0.0, 0.0]]
    """

    positions: ArrayLike
    cell: ArrayLike
    atom_types: ArrayLike
    type_elements: Mapping[int, str]
    origin: ArrayLike = field(default_factory=lambda: np.zeros(3))
    pbc: tuple[bool, bool, bool] = (True, True, True)
    properties: Mapping[str, ArrayLike] = field(default_factory=dict)
    metadata: Mapping[str, Any] = field(default_factory=dict)
    name: str = "structure"

    def __post_init__(self) -> None:
        """Normalize input arrays and validate all structural invariants."""
        self.positions = np.asarray(self.positions, dtype=np.float64).copy()
        self.cell = np.asarray(self.cell, dtype=np.float64).copy()
        self.atom_types = np.asarray(self.atom_types, dtype=np.int64).copy()
        self.type_elements = {
            int(type_id): str(symbol)
            for type_id, symbol in self.type_elements.items()
        }
        self.origin = np.asarray(self.origin, dtype=np.float64).copy()
        self.properties = _copy_properties(self.properties)
        self.metadata = dict(self.metadata)
        self.validate()

    @property
    def number_of_atoms(self) -> int:
        """Return the number of atoms in the structure."""
        return int(self.positions.shape[0])

    @property
    def volume(self) -> float:
        """Return the absolute volume of the simulation cell."""
        return float(abs(np.linalg.det(self.cell)))

    @property
    def fractional_positions(self) -> FloatArray:
        """Return positions expressed in fractional cell coordinates."""
        return np.linalg.solve(
            self.cell.T,
            (self.positions - self.origin).T,
        ).T

    def validate(self) -> None:
        """Validate shapes, type mappings, and cell invertibility.

        Raises:
            ValueError: If an array shape is invalid, the cell is singular,
                or an atom type has no element mapping.
        """
        if self.positions.ndim != 2 or self.positions.shape[1] != 3:
            raise ValueError("positions must have shape (number_of_atoms, 3)")
        if self.cell.shape != (3, 3):
            raise ValueError("cell must have shape (3, 3)")
        if self.origin.shape != (3,):
            raise ValueError("origin must have shape (3,)")
        if self.atom_types.shape != (self.number_of_atoms,):
            raise ValueError("atom_types must contain one value per atom")
        if self.volume <= np.finfo(float).eps:
            raise ValueError("cell vectors must define a non-zero volume")
        if self.number_of_atoms and np.any(self.atom_types < 1):
            raise ValueError("atom type identifiers must start at one")

        missing_types = set(np.unique(self.atom_types)) - set(self.type_elements)
        if missing_types:
            raise ValueError(
                f"missing element mappings for atom types {sorted(missing_types)}"
            )

        for name, values in self.properties.items():
            if values.ndim == 0 or len(values) != self.number_of_atoms:
                raise ValueError(
                    f"property {name!r} must contain one value per atom"
                )

    def copy(self, *, name: str | None = None) -> CrystalStructure:
        """Return a deep copy of the structure.

        Args:
            name: Optional name for the copied structure.

        Returns:
            An independent :class:`CrystalStructure` instance.
        """
        return CrystalStructure(
            positions=self.positions,
            cell=self.cell,
            atom_types=self.atom_types,
            type_elements=self.type_elements,
            origin=self.origin,
            pbc=self.pbc,
            properties=self.properties,
            metadata=self.metadata,
            name=self.name if name is None else name,
        )

    @classmethod
    def from_fractional(
        cls,
        fractional_positions: ArrayLike,
        cell: ArrayLike,
        atom_types: ArrayLike,
        type_elements: Mapping[int, str],
        *,
        origin: ArrayLike | None = None,
        **kwargs: Any,
    ) -> CrystalStructure:
        """Create a structure from fractional coordinates.

        Args:
            fractional_positions: Fractional coordinates with shape ``(N, 3)``.
            cell: Three cell vectors stored as rows.
            atom_types: Positive integer type identifier for each atom.
            type_elements: Mapping from type identifiers to element symbols.
            origin: Cartesian cell origin. Defaults to zero.
            **kwargs: Additional :class:`CrystalStructure` fields.

        Returns:
            A validated structure with Cartesian positions.
        """
        cell_array = np.asarray(cell, dtype=np.float64)
        origin_array = (
            np.zeros(3, dtype=np.float64)
            if origin is None
            else np.asarray(origin, dtype=np.float64)
        )
        positions = np.asarray(fractional_positions, dtype=np.float64) @ cell_array
        positions += origin_array
        return cls(
            positions=positions,
            cell=cell_array,
            atom_types=atom_types,
            type_elements=type_elements,
            origin=origin_array,
            **kwargs,
        )

    @classmethod
    def read(
        cls,
        filename: str | Path,
        *,
        format: str | None = None,
        **kwargs: Any,
    ) -> CrystalStructure:
        """Read an external structure file.

        Args:
            filename: Input file path.
            format: Explicit format name. The suffix is used when omitted.
            **kwargs: Format-specific reader options.

        Returns:
            The structure decoded by the selected reader.
        """
        from crystalops.io import read_structure

        return read_structure(filename, format=format, **kwargs)

    def write(
        self,
        filename: str | Path,
        *,
        format: str | None = None,
        **kwargs: Any,
    ) -> None:
        """Write the structure to an external file.

        Args:
            filename: Output file path.
            format: Explicit format name. The suffix is used when omitted.
            **kwargs: Format-specific writer options.
        """
        from crystalops.io import write_structure

        write_structure(self, filename, format=format, **kwargs)

    def translated(self, vector: ArrayLike) -> CrystalStructure:
        """Return a copy whose atom positions are translated by ``vector``."""
        from crystalops.operations import translate

        return translate(self, vector)

    def rotated(
        self,
        rotation: ArrayLike,
        *,
        center: ArrayLike | None = None,
        rotate_cell: bool = False,
    ) -> CrystalStructure:
        """Return a rotated copy of the structure."""
        from crystalops.operations import rotate

        return rotate(self, rotation, center=center, rotate_cell=rotate_cell)

    def repeated(self, repeats: tuple[int, int, int]) -> CrystalStructure:
        """Return a supercell repeated along all three cell vectors."""
        from crystalops.operations import repeat

        return repeat(self, repeats)


@dataclass(slots=True)
class UnitCell:
    """Describe a reusable crystallographic unit-cell template.

    ``fractional_positions`` are dimensionless, while ``cell`` describes the
    cell shape relative to a lattice constant supplied by a builder.

    Examples:
        >>> from crystalops.lattices import face_centered_cubic
        >>> cell = face_centered_cubic("Cu")
        >>> cell.number_of_sites
        4
    """

    name: str
    fractional_positions: ArrayLike
    cell: ArrayLike
    atom_types: ArrayLike
    type_elements: Mapping[int, str]

    def __post_init__(self) -> None:
        """Normalize arrays and validate the unit-cell definition."""
        self.fractional_positions = np.asarray(
            self.fractional_positions,
            dtype=np.float64,
        ).copy()
        self.cell = np.asarray(self.cell, dtype=np.float64).copy()
        self.atom_types = np.asarray(self.atom_types, dtype=np.int64).copy()
        self.type_elements = {
            int(type_id): str(symbol)
            for type_id, symbol in self.type_elements.items()
        }
        self.validate()

    @property
    def number_of_sites(self) -> int:
        """Return the number of basis sites in the unit cell."""
        return int(self.fractional_positions.shape[0])

    def validate(self) -> None:
        """Validate coordinates, cell vectors, and atom type mappings."""
        if (
            self.fractional_positions.ndim != 2
            or self.fractional_positions.shape[1] != 3
        ):
            raise ValueError("fractional_positions must have shape (N, 3)")
        if self.cell.shape != (3, 3):
            raise ValueError("cell must have shape (3, 3)")
        if abs(np.linalg.det(self.cell)) <= np.finfo(float).eps:
            raise ValueError("cell vectors must define a non-zero volume")
        if self.atom_types.shape != (self.number_of_sites,):
            raise ValueError("atom_types must contain one value per site")
        if self.number_of_sites and np.any(self.atom_types < 1):
            raise ValueError("atom type identifiers must start at one")
        missing_types = set(np.unique(self.atom_types)) - set(self.type_elements)
        if missing_types:
            raise ValueError(
                f"missing element mappings for atom types {sorted(missing_types)}"
            )
