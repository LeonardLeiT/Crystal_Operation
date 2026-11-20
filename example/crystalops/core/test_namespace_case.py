"""Compatibility and dependency-boundary cases for the unified namespace."""

from __future__ import annotations

import importlib
from pathlib import Path
import pickle
import subprocess
import sys

import numpy as np


def test_core_import_does_not_load_feature_modules() -> None:
    """Data-model users need neither builders nor optional scientific tools."""
    code = """
import sys
import crystalops
assert 'crystalops.modeling' not in sys.modules
assert 'crystalops.io' not in sys.modules
from crystalops.core import CrystalStructure, UnitCell, Element
assert CrystalStructure.__module__ == 'crystalops.core.structure'
assert UnitCell.__module__ == 'crystalops.core.unit_cell'
for name in ('crystalops.modeling', 'crystalops.analysis',
             'crystalops.visualization', 'crystalops.io',
             'matplotlib', 'scipy', 'skimage', 'pandas', 'seaborn'):
    assert name not in sys.modules, name
"""
    subprocess.run(
        [sys.executable, "-c", code],
        cwd=Path(__file__).resolve().parents[3],
        check=True,
        capture_output=True,
        text=True,
    )


def test_analysis_and_visualization_namespaces_are_lightweight() -> None:
    """Importing namespaces must not require or initialize a plotting backend."""
    code = """
import sys
import crystalops.analysis as analysis
import crystalops.visualization as visualization
assert callable(analysis.aniso_lattice)
assert 'StructurePlotter' in dir(visualization)
for name in ('matplotlib', 'scipy', 'skimage', 'pandas', 'seaborn'):
    assert name not in sys.modules, name
"""
    subprocess.run(
        [sys.executable, "-c", code],
        cwd=Path(__file__).resolve().parents[3],
        check=True,
        capture_output=True,
        text=True,
    )


def test_old_and_new_modeling_exports_share_implementations() -> None:
    """Compatibility paths preserve class identity and private helper access."""
    import crystalops
    from crystalops import modeling

    for package in ("builders", "operations", "lattices"):
        legacy = importlib.import_module(f"crystalops.{package}")
        canonical = importlib.import_module(f"crystalops.modeling.{package}")
        for name in canonical.__all__:
            assert getattr(legacy, name) is getattr(canonical, name), name
            assert getattr(modeling, name) is getattr(canonical, name), name
            if name in crystalops.__all__:
                assert getattr(crystalops, name) is getattr(canonical, name), name

    for module in ("base", "crystal", "polycrystal", "schwarz", "tpms", "registry", "_lattice"):
        assert importlib.import_module(f"crystalops.builders.{module}") is (
            importlib.import_module(f"crystalops.modeling.builders.{module}")
        )
    assert importlib.import_module("crystalops.operations.geometry") is (
        importlib.import_module("crystalops.modeling.operations.geometry")
    )


def test_builder_registry_is_shared_across_import_paths() -> None:
    """A builder registered through the old API is usable through the new one."""
    from crystalops.builders import register_builder
    from crystalops.modeling import CrystalBuilder, FCC, StructureBuilder, create_structure
    from crystalops.modeling.builders import registry

    class ExampleBuilder(StructureBuilder):
        def build(self):
            return CrystalBuilder(FCC("Cu"), (1, 1, 1), 3.61).build()

    name = "namespace-compatibility-case"
    register_builder(name, ExampleBuilder)
    try:
        assert name in registry.available_builders()
        assert create_structure(name).number_of_atoms == 4
    finally:
        registry._BUILDERS.pop(name, None)


def test_unit_cell_and_structure_keep_one_shared_identity() -> None:
    """Models survive serialization and are shared by old and new factories."""
    from crystalops.core import CrystalStructure, UnitCell
    from crystalops.core.structure import UnitCell as LegacyUnitCell
    from crystalops.core.unit_cell import UnitCell as CanonicalUnitCell
    from crystalops.modeling import CrystalBuilder, FCC

    assert UnitCell is LegacyUnitCell is CanonicalUnitCell
    unit_cell = pickle.loads(pickle.dumps(FCC("Cu")))
    assert isinstance(unit_cell, UnitCell)
    structure = CrystalBuilder(unit_cell, (1, 1, 1), 3.61).build()
    structure.properties["grain_ids"] = np.ones(structure.number_of_atoms, dtype=int)
    restored = pickle.loads(pickle.dumps(structure))
    assert isinstance(restored, CrystalStructure)
    assert np.array_equal(restored.positions, structure.positions)
    assert np.array_equal(restored.properties["grain_ids"], structure.properties["grain_ids"])
    assert isinstance(restored.translated((1, 0, 0)), CrystalStructure)
    assert isinstance(restored.repeated((2, 1, 1)), CrystalStructure)
