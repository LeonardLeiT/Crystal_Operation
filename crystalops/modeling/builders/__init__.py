"""Structure builders and the public builder registry."""

from .base import StructureBuilder
from .crystal import CrystalBuilder
from .polycrystal import (
    KelvinCrystal,
    VoronoiPolycrystalBuilder,
    WeairePhelanCrystal,
)
from .registry import available_builders, create_structure, register_builder
from .schwarz import SchwarzCrystal
from .tpms import TpmsBuilder

register_builder("crystal", CrystalBuilder, replace=True)
register_builder("voronoi", VoronoiPolycrystalBuilder, replace=True)
register_builder("kelvin", KelvinCrystal, replace=True)
register_builder("weaire-phelan", WeairePhelanCrystal, replace=True)
register_builder("tpms", TpmsBuilder, replace=True)
register_builder("schwarz", SchwarzCrystal, replace=True)

__all__ = [
    "CrystalBuilder",
    "KelvinCrystal",
    "SchwarzCrystal",
    "StructureBuilder",
    "TpmsBuilder",
    "VoronoiPolycrystalBuilder",
    "WeairePhelanCrystal",
    "available_builders",
    "create_structure",
    "register_builder",
]
