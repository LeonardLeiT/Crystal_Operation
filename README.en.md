# CrystalOps

CrystalOps builds, transforms, reads, and writes particle-level crystal
structures. Structures generated in Python and structures loaded from external
files share the same `CrystalStructure` representation.

## Installation

```bash
pip install -e .
```

Install optional plotting dependencies with:

```bash
pip install -e ".[plot]"
```

## Example

```python
from crystalops import CrystalBuilder, FCC

structure = CrystalBuilder(
    unit_cell=FCC("Cu"),
    repeats=(10, 10, 10),
    lattice_constant=3.61,
).build()
structure.write("fcc.data")
```

The numerical package is organized into `core`, `operations`, `builders`, and
`io` modules. Visualization lives in the separate `crysplot` package. LAMMPS
data files using `atom_style atomic` are currently supported for both reading
and writing.
