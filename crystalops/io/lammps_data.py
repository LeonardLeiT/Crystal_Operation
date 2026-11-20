"""LAMMPS data reader and writer for the atomic atom style."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray

from crystalops.core import CrystalStructure, get_element
from crystalops.core.elements import infer_element
from crystalops.io.base import StructureIO


FloatArray = NDArray[np.float64]

_FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
_SECTIONS = {
    "Masses",
    "Atoms",
    "Velocities",
    "Bonds",
    "Angles",
    "Dihedrals",
    "Impropers",
    "Pair Coeffs",
    "Bond Coeffs",
    "Angle Coeffs",
    "Dihedral Coeffs",
    "Improper Coeffs",
}


class LammpsDataIO(StructureIO):
    """Read and write LAMMPS data files using ``atom_style atomic``.

    Orthogonal and restricted triclinic boxes are supported. General cells
    are rotated to LAMMPS restricted-triclinic form during writing while
    preserving all fractional coordinates.

    Examples:
        >>> from crystalops.modeling.builders import CrystalBuilder
        >>> from crystalops.modeling.lattices import face_centered_cubic
        >>> structure = CrystalBuilder(face_centered_cubic(), (2, 2, 2), 3.61).build()
        >>> LammpsDataIO.write(structure, Path("copper.data"))
        >>> loaded = LammpsDataIO.read(Path("copper.data"))
        >>> loaded.number_of_atoms
        32
    """

    @classmethod
    def read(
        cls,
        filename: Path,
        *,
        atom_style: str = "atomic",
        encoding: str = "utf-8",
        **_: Any,
    ) -> CrystalStructure:
        """Read a LAMMPS data file.

        Args:
            filename: Path to the data file.
            atom_style: Layout of the ``Atoms`` section. Only ``atomic`` is
                currently supported.
            encoding: Text encoding used to open the file.

        Returns:
            A validated crystal structure sorted by atom identifier.
        """
        if atom_style != "atomic":
            raise NotImplementedError("only LAMMPS atom_style atomic is supported")
        lines = filename.read_text(encoding=encoding).splitlines()
        atom_count = cls._header_integer(lines, "atoms")
        atom_type_count = cls._header_integer(lines, "atom types")
        origin, cell = cls._read_box(lines)
        masses = cls._read_masses(lines)
        atom_ids, atom_types, positions = cls._read_atoms(lines, atom_style)

        if atom_count != len(atom_ids):
            raise ValueError(
                f"header declares {atom_count} atoms but Atoms contains {len(atom_ids)}"
            )
        if atom_type_count < int(atom_types.max(initial=0)):
            raise ValueError("an atom type exceeds the declared atom type count")

        order = np.argsort(atom_ids)
        atom_ids = atom_ids[order]
        atom_types = atom_types[order]
        positions = positions[order]
        type_elements = {
            type_id: symbol
            for type_id in range(1, atom_type_count + 1)
            for symbol in [cls._element_for_type(type_id, masses)]
        }
        return CrystalStructure(
            positions=positions,
            cell=cell,
            atom_types=atom_types,
            type_elements=type_elements,
            origin=origin,
            properties={"atom_ids": atom_ids},
            metadata={
                "source": str(filename),
                "format": "lammps-data",
                "atom_style": atom_style,
                "type_masses": {
                    type_id: mass for type_id, (mass, _) in masses.items()
                },
            },
            name=filename.stem,
        )

    @classmethod
    def write(
        cls,
        structure: CrystalStructure,
        filename: Path,
        *,
        atom_style: str = "atomic",
        encoding: str = "utf-8",
        **_: Any,
    ) -> None:
        """Write a structure as a LAMMPS data file.

        Args:
            structure: Structure to encode.
            filename: Destination path.
            atom_style: LAMMPS atom style. Only ``atomic`` is supported.
            encoding: Text encoding used to create the file.
        """
        if atom_style != "atomic":
            raise NotImplementedError("only LAMMPS atom_style atomic is supported")
        if structure.number_of_atoms == 0:
            raise ValueError("cannot write an empty structure")

        maximum_type = int(structure.atom_types.max())
        missing = set(range(1, maximum_type + 1)) - set(structure.type_elements)
        if missing:
            raise ValueError(
                "LAMMPS requires contiguous type metadata; missing mappings for "
                f"{sorted(missing)}"
            )

        restricted_cell, positions = cls._to_restricted_cell(structure)
        origin = structure.origin
        lx, ly, lz = np.diag(restricted_cell)
        xy = restricted_cell[1, 0]
        xz = restricted_cell[2, 0]
        yz = restricted_cell[2, 1]

        atom_ids = structure.properties.get(
            "atom_ids",
            np.arange(1, structure.number_of_atoms + 1, dtype=np.int64),
        )
        atom_ids = np.asarray(atom_ids, dtype=np.int64)
        if len(np.unique(atom_ids)) != structure.number_of_atoms or np.any(atom_ids < 1):
            raise ValueError("atom_ids must be unique positive integers")

        output: list[str] = [
            f"# CrystalOps structure: {structure.name}",
            "",
            f"{structure.number_of_atoms} atoms",
            f"{maximum_type} atom types",
            "",
            f"{origin[0]:.16g} {origin[0] + lx:.16g} xlo xhi",
            f"{origin[1]:.16g} {origin[1] + ly:.16g} ylo yhi",
            f"{origin[2]:.16g} {origin[2] + lz:.16g} zlo zhi",
            f"{xy:.16g} {xz:.16g} {yz:.16g} xy xz yz",
            "",
            "Masses",
            "",
        ]
        type_masses = structure.metadata.get("type_masses", {})
        for type_id in range(1, maximum_type + 1):
            symbol = structure.type_elements[type_id]
            if type_id in type_masses:
                mass = float(type_masses[type_id])
            else:
                try:
                    mass = get_element(symbol).atomic_mass
                except KeyError:
                    raise ValueError(
                        f"no atomic mass is available for type {type_id} ({symbol})"
                    ) from None
            output.append(f"{type_id} {mass:.10g} # {symbol}")

        output.extend(["", "Atoms # atomic", ""])
        for atom_id, type_id, position in zip(
            atom_ids,
            structure.atom_types,
            positions,
            strict=True,
        ):
            output.append(
                f"{atom_id} {type_id} "
                f"{position[0]:.16g} {position[1]:.16g} {position[2]:.16g}"
            )
        output.append("")
        filename.parent.mkdir(parents=True, exist_ok=True)
        filename.write_text("\n".join(output), encoding=encoding)

    @staticmethod
    def _header_integer(lines: list[str], label: str) -> int:
        """Extract an integer count from the LAMMPS header."""
        pattern = re.compile(rf"^\s*(\d+)\s+{re.escape(label)}\s*$")
        for line in lines:
            match = pattern.match(line.split("#", 1)[0])
            if match:
                return int(match.group(1))
        raise ValueError(f"LAMMPS header is missing the {label!r} count")

    @staticmethod
    def _read_box(lines: list[str]) -> tuple[FloatArray, FloatArray]:
        """Read orthogonal or restricted-triclinic box parameters."""
        bounds: dict[str, tuple[float, float]] = {}
        tilt = (0.0, 0.0, 0.0)
        bound_pattern = re.compile(
            rf"^\s*({_FLOAT})\s+({_FLOAT})\s+([xyz])lo\s+\3hi\s*$"
        )
        tilt_pattern = re.compile(
            rf"^\s*({_FLOAT})\s+({_FLOAT})\s+({_FLOAT})\s+xy\s+xz\s+yz\s*$"
        )
        for line in lines:
            content = line.split("#", 1)[0]
            bound_match = bound_pattern.match(content)
            if bound_match:
                bounds[bound_match.group(3)] = (
                    float(bound_match.group(1)),
                    float(bound_match.group(2)),
                )
            tilt_match = tilt_pattern.match(content)
            if tilt_match:
                tilt = tuple(float(value) for value in tilt_match.groups())
        if bounds.keys() != {"x", "y", "z"}:
            raise ValueError("LAMMPS header must define x, y, and z box bounds")
        origin = np.array([bounds[axis][0] for axis in "xyz"], dtype=np.float64)
        lengths = [bounds[axis][1] - bounds[axis][0] for axis in "xyz"]
        xy, xz, yz = tilt
        cell = np.array(
            [
                [lengths[0], 0.0, 0.0],
                [xy, lengths[1], 0.0],
                [xz, yz, lengths[2]],
            ],
            dtype=np.float64,
        )
        return origin, cell

    @classmethod
    def _section(cls, lines: list[str], name: str) -> list[str]:
        """Return non-empty data lines belonging to a named section."""
        start: int | None = None
        for index, line in enumerate(lines):
            heading = line.split("#", 1)[0].strip()
            if heading == name:
                start = index + 1
                break
        if start is None:
            raise ValueError(f"LAMMPS data file is missing the {name!r} section")

        result: list[str] = []
        for line in lines[start:]:
            heading = line.split("#", 1)[0].strip()
            if heading in _SECTIONS:
                break
            if line.strip():
                result.append(line)
        return result

    @classmethod
    def _read_masses(cls, lines: list[str]) -> dict[int, tuple[float, str | None]]:
        """Read atom masses and optional element-symbol comments."""
        masses: dict[int, tuple[float, str | None]] = {}
        for line in cls._section(lines, "Masses"):
            values, _, comment = line.partition("#")
            fields = values.split()
            if len(fields) < 2:
                continue
            symbol = comment.strip().split()[0] if comment.strip() else None
            masses[int(fields[0])] = (float(fields[1]), symbol)
        return masses

    @classmethod
    def _read_atoms(
        cls,
        lines: list[str],
        atom_style: str,
    ) -> tuple[NDArray[np.int64], NDArray[np.int64], FloatArray]:
        """Read IDs, types, and Cartesian positions from ``Atoms``."""
        atom_lines = cls._section(lines, "Atoms")
        records: list[tuple[int, int, float, float, float]] = []
        for line in atom_lines:
            fields = line.split("#", 1)[0].split()
            if not fields:
                continue
            if atom_style == "atomic" and len(fields) < 5:
                raise ValueError("atomic-style Atoms lines require id type x y z")
            records.append(
                (
                    int(fields[0]),
                    int(fields[1]),
                    float(fields[2]),
                    float(fields[3]),
                    float(fields[4]),
                )
            )
        array = np.asarray(records, dtype=np.float64)
        if array.size == 0:
            return (
                np.empty(0, dtype=np.int64),
                np.empty(0, dtype=np.int64),
                np.empty((0, 3), dtype=np.float64),
            )
        return (
            array[:, 0].astype(np.int64),
            array[:, 1].astype(np.int64),
            array[:, 2:5],
        )

    @staticmethod
    def _element_for_type(
        type_id: int,
        masses: dict[int, tuple[float, str | None]],
    ) -> str:
        """Resolve an element symbol from a Masses record."""
        if type_id not in masses:
            return f"X{type_id}"
        mass, symbol = masses[type_id]
        return symbol or infer_element(mass) or f"X{type_id}"

    @staticmethod
    def _to_restricted_cell(
        structure: CrystalStructure,
    ) -> tuple[FloatArray, FloatArray]:
        """Rotate a general cell and its atoms into LAMMPS restricted form."""
        a, b, c = structure.cell
        lx = np.linalg.norm(a)
        a_hat = a / lx
        xy = float(np.dot(b, a_hat))
        b_perpendicular = b - xy * a_hat
        ly = np.linalg.norm(b_perpendicular)
        if ly <= np.finfo(float).eps:
            raise ValueError("first and second cell vectors are collinear")
        b_hat = b_perpendicular / ly
        c_hat = np.cross(a_hat, b_hat)
        xz = float(np.dot(c, a_hat))
        yz = float(np.dot(c, b_hat))
        lz = float(np.dot(c, c_hat))
        if lz < 0:
            c_hat = -c_hat
            lz = -lz
        restricted = np.array(
            [[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]],
            dtype=np.float64,
        )
        positions = structure.fractional_positions @ restricted + structure.origin
        return restricted, positions
