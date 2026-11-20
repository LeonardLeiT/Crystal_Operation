"""Compatibility alias for :mod:`crystalops.modeling.operations.geometry`."""

from importlib import import_module as _import_module
import sys as _sys

_sys.modules[__name__] = _import_module("crystalops.modeling.operations.geometry")
