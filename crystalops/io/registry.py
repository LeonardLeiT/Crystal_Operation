"""Registry and public dispatch functions for structure formats."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from crystalops.core import CrystalStructure
from crystalops.io.base import StructureIO


_FORMATS: dict[str, type[StructureIO]] = {}
_SUFFIXES: dict[str, str] = {}


def register_format(
    name: str,
    handler: type[StructureIO],
    *,
    suffixes: tuple[str, ...] = (),
    replace: bool = False,
) -> None:
    """Register a file-format handler and optional filename suffixes."""
    key = name.strip().lower().replace("_", "-")
    if not key:
        raise ValueError("format name cannot be empty")
    if key in _FORMATS and not replace:
        raise ValueError(f"format {key!r} is already registered")
    _FORMATS[key] = handler
    for suffix in suffixes:
        normalized = suffix.lower()
        if not normalized.startswith("."):
            normalized = f".{normalized}"
        _SUFFIXES[normalized] = key


def available_formats() -> tuple[str, ...]:
    """Return registered format names in deterministic order."""
    return tuple(sorted(_FORMATS))


def _resolve_format(filename: Path, format: str | None) -> type[StructureIO]:
    """Resolve an explicit format or infer one from a filename suffix."""
    if format is None:
        try:
            key = _SUFFIXES[filename.suffix.lower()]
        except KeyError as exc:
            raise ValueError(
                f"cannot infer format from suffix {filename.suffix!r}; "
                "pass format explicitly"
            ) from exc
    else:
        key = format.strip().lower().replace("_", "-")
        aliases = {"lammps": "lammps-data", "data": "lammps-data"}
        key = aliases.get(key, key)
    try:
        return _FORMATS[key]
    except KeyError as exc:
        raise ValueError(
            f"unknown format {key!r}; available formats: {available_formats()}"
        ) from exc


def read_structure(
    filename: str | Path,
    *,
    format: str | None = None,
    **kwargs: Any,
) -> CrystalStructure:
    """Read one structure from an external file."""
    path = Path(filename)
    return _resolve_format(path, format).read(path, **kwargs)


def write_structure(
    structure: CrystalStructure,
    filename: str | Path,
    *,
    format: str | None = None,
    **kwargs: Any,
) -> None:
    """Write a structure to an external file."""
    path = Path(filename)
    _resolve_format(path, format).write(structure, path, **kwargs)
