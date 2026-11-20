"""Export a complete structure with its per-atom selection flag."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from crystalops.core import CrystalStructure


def write_selection_dump(
    structure: CrystalStructure,
    path: str | Path,
    *,
    ovito_selection: bool = False,
) -> None:
    """Write one LAMMPS dump frame with ``id type x y z select`` columns.

    The ``select`` property must contain one integer (0 or 1) per atom.
    All atoms are written in structure order. Their original IDs are used
    when present. Orthogonal and restricted triclinic cells preserve the
    Cartesian coordinates, cell origin, and periodic boundary flags.
    Set ``ovito_selection=True`` to also write the native OVITO Selection column.
    """
    structure.validate()
    count = structure.number_of_atoms
    selection = structure.properties.get("select")
    if (
        selection is None
        or selection.shape != (count,)
        or not np.issubdtype(selection.dtype, np.integer)
        or np.any((selection != 0) & (selection != 1))
    ):
        raise ValueError("the 'select' property must contain one integer 0 or 1 per atom")

    atom_ids = structure.properties.get("atom_ids")
    if atom_ids is None:
        atom_ids = np.arange(1, count + 1, dtype=np.int64)
    elif (
        atom_ids.shape != (count,)
        or not np.issubdtype(atom_ids.dtype, np.integer)
        or np.any(atom_ids <= 0)
        or int(atom_ids.max(initial=0)) > np.iinfo(np.int64).max
        or len(np.unique(atom_ids)) != count
    ):
        raise ValueError("the 'atom_ids' property must contain unique positive 64-bit integers")

    cell = structure.cell
    origin = structure.origin
    if not all(np.isfinite(array).all() for array in (structure.positions, cell, origin)):
        raise ValueError("positions, cell, and origin must be finite")
    if len(structure.pbc) != 3:
        raise ValueError("pbc must contain three boundary flags")
    lengths = np.diag(cell)
    if np.any(np.triu(cell, k=1) != 0) or np.any(lengths <= 0):
        raise ValueError("selection dump requires an axis-aligned orthogonal or restricted triclinic cell")

    flags = " ".join("pp" if flag else "ff" for flag in structure.pbc)
    box_header = f"ITEM: BOX BOUNDS {flags}"
    xy, xz, yz = cell[1, 0], cell[2, 0], cell[2, 1]
    lower = origin.copy()
    upper = origin + lengths
    tilts = None
    if xy or xz or yz:
        # Dump files contain enclosing Cartesian bounds, unlike data files.
        box_header = f"ITEM: BOX BOUNDS xy xz yz {flags}"
        lower += (min(0, xy, xz, xy + xz), min(0, yz), 0)
        upper += (max(0, xy, xz, xy + xz), max(0, yz), 0)
        tilts = (xy, xz, yz)
    if not np.isfinite(lower).all() or not np.isfinite(upper).all() or np.any(upper <= lower):
        raise ValueError("the simulation box has unrepresentable dump bounds")

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as file:
        file.write(f"ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS\n{count}\n{box_header}\n")
        for axis, (lo, hi) in enumerate(zip(lower, upper, strict=True)):
            tilt = "" if tilts is None else f" {tilts[axis]:.17g}"
            file.write(f"{lo:.17g} {hi:.17g}{tilt}\n")
        columns = "id type x y z select"
        if ovito_selection:
            columns += " selection"
        file.write(f"ITEM: ATOMS {columns}\n")
        for atom_id, atom_type, (x, y, z), flag in zip(
            atom_ids, structure.atom_types, structure.positions, selection, strict=True,
        ):
            native_flag = f" {flag}" if ovito_selection else ""
            file.write(f"{atom_id} {atom_type} {x:.17g} {y:.17g} {z:.17g} {flag}{native_flag}\n")
