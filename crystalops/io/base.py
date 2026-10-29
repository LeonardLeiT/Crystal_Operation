"""Interfaces implemented by structure file formats."""

from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any

from crystalops.core import CrystalStructure


class StructureIO(ABC):
    """Decode and encode one atomistic structure file format."""

    @classmethod
    @abstractmethod
    def read(cls, filename: Path, **kwargs: Any) -> CrystalStructure:
        """Read ``filename`` into a :class:`CrystalStructure`."""

    @classmethod
    @abstractmethod
    def write(
        cls,
        structure: CrystalStructure,
        filename: Path,
        **kwargs: Any,
    ) -> None:
        """Write ``structure`` to ``filename``."""
