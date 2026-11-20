# CrystalOps

CrystalOps provides particle-level structure modeling, analysis, visualization,
and file exchange. Structures generated in Python and structures loaded from
external files share the same `CrystalStructure` representation.

## Organization

- `crystalops.core`: shared `CrystalStructure`, `UnitCell`, and `Element` data
  models, including coordinates, simulation cells, periodic boundaries, and
  per-atom properties.
- `crystalops.modeling`: builders, structure operations, and lattice factories.
- `crystalops.analysis`: numerical analysis, including anisotropy, periodic
  grain-boundary identification, Gaussian occupancy-interface reconstruction,
  analytic local curvature, and optional PTM grain-label recovery.
- `crystalops.visualization`: structure plotting and the existing plotting
  recipes, previously provided by `crysplot`.
- `crystalops.io`: structure-format readers and writers.

Modeling, analysis, and I/O share the core data models. Plotting dependencies
are optional and are not required to import the core package.
The retained recipes are standalone scripts that require their own input
data and additional dependencies when run.

Legacy imports from `crystalops.builders`, `crystalops.operations`,
`crystalops.lattices`, and `crysplot` remain compatibility entry points.
Top-level imports such as `from crystalops import CrystalBuilder` also remain
available. New code should use the modules listed above.

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
from crystalops.modeling import CrystalBuilder, FCC

structure = CrystalBuilder(
    unit_cell=FCC("Cu"),
    repeats=(10, 10, 10),
    lattice_constant=3.61,
).build()
structure.write("fcc.data")
```

LAMMPS data files using `atom_style atomic` are currently supported for both
reading and writing. This format handler does not preserve arbitrary per-atom
properties such as `grain_ids`. Keep the structure returned by a builder for
analysis requiring these labels, or load labels explicitly and align them by
atom ID.

Examples follow the same module layout; see [the example guide](example/README.md).
