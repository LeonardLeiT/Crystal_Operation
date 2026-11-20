"""Abstract interfaces for structure construction."""

from __future__ import annotations

from abc import ABC, abstractmethod

from crystalops.core import CrystalStructure


class StructureBuilder(ABC):
    """Build a complete :class:`CrystalStructure` without writing files."""

    @abstractmethod
    def build(self) -> CrystalStructure:
        """Construct and return a new crystal structure."""
