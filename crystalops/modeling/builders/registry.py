"""Named registry for reusable structure builders."""

from __future__ import annotations

from typing import Any

from crystalops.modeling.builders.base import StructureBuilder
from crystalops.core import CrystalStructure


_BUILDERS: dict[str, type[StructureBuilder]] = {}


def register_builder(
    name: str,
    builder: type[StructureBuilder],
    *,
    replace: bool = False,
) -> None:
    """Register a builder class under a normalized name."""
    key = name.strip().lower().replace("_", "-")
    if not key:
        raise ValueError("builder name cannot be empty")
    if key in _BUILDERS and not replace:
        raise ValueError(f"builder {key!r} is already registered")
    _BUILDERS[key] = builder


def available_builders() -> tuple[str, ...]:
    """Return registered builder names in deterministic order."""
    return tuple(sorted(_BUILDERS))


def create_structure(name: str, **kwargs: Any) -> CrystalStructure:
    """Instantiate a registered builder and return its structure."""
    key = name.strip().lower().replace("_", "-")
    try:
        builder_class = _BUILDERS[key]
    except KeyError as exc:
        raise ValueError(
            f"unknown builder {name!r}; available builders: {available_builders()}"
        ) from exc
    return builder_class(**kwargs).build()
