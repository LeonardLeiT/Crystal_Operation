"""Recover current grain labels from local lattice orientations using OVITO.

OVITO is optional and is imported only when segmentation is requested.  This
module does not infer interfaces from a sample's name or an ideal surface.
Returned arrays follow the input structure's atom order.  Grain IDs are positive
integers and the reserved value -1 consistently denotes an unassigned atom,
independently of the installed OVITO version.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import NDArray

from crystalops.core import CrystalStructure


@dataclass(slots=True)
class GrainSegmentationResult:
    """PTM diagnostics, pre-adoption grain cores, and final grain assignment."""

    atom_ids: NDArray[np.int64]
    grain_ids: NDArray[np.int64]
    core_grain_ids: NDArray[np.int64]
    structure_types: NDArray[np.int64]
    rmsd: NDArray[np.float64]
    orientations: NDArray[np.float64]
    interatomic_distances: NDArray[np.float64]
    assignment_distances: NDArray[np.float64]
    assignment_distance_margins: NDArray[np.float64]
    diagnostics: dict[str, Any]


def _counts(values: NDArray) -> dict[str, int]:
    labels, counts = np.unique(values, return_counts=True)
    return {str(int(label)): int(count) for label, count in zip(labels, counts)}


def _nearest_grain_core(
    structure: CrystalStructure,
    core_labels: NDArray[np.int64],
) -> tuple[NDArray[np.int64], NDArray[np.float64], NDArray[np.float64], int]:
    """Assign orphans to the closest periodic core; resolve distance ties by ID."""
    from itertools import product
    from scipy.spatial import cKDTree

    cell = np.asarray(structure.cell, dtype=float)
    lengths = np.linalg.norm(cell, axis=1)
    basis = cell / lengths[:, None]
    if not np.allclose(basis @ basis.T, np.eye(3), atol=1e-10):
        raise ValueError("nearest_core adoption currently requires an orthogonal cell; use orphan_method='ovito' for a triclinic cell")
    coordinates = (np.asarray(structure.positions) - structure.origin) @ basis.T
    for axis, periodic in enumerate(structure.pbc):
        if periodic:
            coordinates[:, axis] %= lengths[axis]
    orphan_rows = np.flatnonzero(core_labels < 0)
    labels = core_labels.copy()
    distances = np.zeros(len(labels), dtype=float)
    margins = np.full(len(labels), np.nan)
    if not len(orphan_rows):
        return labels, distances, margins, 0
    grains = np.unique(core_labels[core_labels > 0])
    if not len(grains):
        raise ValueError("no grain cores are available for orphan adoption")
    candidates = np.empty((len(orphan_rows), len(grains)))
    for index, grain in enumerate(grains):
        seeds = coordinates[core_labels == grain]
        if all(structure.pbc):
            tree = cKDTree(seeds, boxsize=lengths)
        else:
            shifts = np.asarray(list(product(*[(-1, 0, 1) if p else (0,) for p in structure.pbc])))
            tree = cKDTree((seeds[None, :, :] + shifts[:, None, :] * lengths).reshape(-1, 3))
        candidates[:, index] = tree.query(coordinates[orphan_rows], k=1)[0]
    minimum = candidates.min(axis=1)
    tolerance = 1e-8
    tied = candidates <= minimum[:, None] + tolerance
    choice = tied.argmax(axis=1)
    labels[orphan_rows] = grains[choice]
    distances[orphan_rows] = candidates[np.arange(len(orphan_rows)), choice]
    if len(grains) > 1:
        ordered = np.sort(candidates, axis=1)
        margins[orphan_rows] = ordered[:, 1] - ordered[:, 0]
    return labels, distances, margins, int(np.count_nonzero(tied.sum(axis=1) > 1))


def segment_grains_with_ptm(
    structure: CrystalStructure,
    *,
    rmsd_cutoff: float = 0.1,
    min_grain_size: int = 50,
    handle_stacking_faults: bool = True,
    orphan_adoption: bool = True,
    orphan_method: str = "nearest_core",
    algorithm: str = "graph_auto",
    merging_threshold: float = 0.0,
) -> GrainSegmentationResult:
    """Infer grain labels using PTM orientations and graph segmentation.

    The complete atomistic structure is analyzed with its actual cell vectors,
    origin, and PBC flags.  Before orphan adoption, noncrystalline boundary atoms
    remain unassigned.  The optional adoption pass assigns these atoms to nearby
    grains.  Its default periodic nearest-core rule resolves distance ties by
    grain ID; ``orphan_method='ovito'`` selects OVITO's native adoption instead.
    Original core labels and nearest-distinct-grain distance margins are kept.

    ``algorithm`` is ``graph_auto``, ``graph_manual``, or
    ``minimum_spanning_tree``.  ``merging_threshold`` is ignored in auto mode;
    its interpretation depends on the selected OVITO algorithm.  No file is
    modified, and grain labels are not written back into ``structure``.
    """
    if not np.isfinite(rmsd_cutoff) or rmsd_cutoff <= 0:
        raise ValueError("rmsd_cutoff must be a finite positive number")
    if isinstance(min_grain_size, bool) or int(min_grain_size) != min_grain_size or min_grain_size < 1:
        raise ValueError("min_grain_size must be a positive integer")
    if not np.isfinite(merging_threshold):
        raise ValueError("merging_threshold must be finite")
    if orphan_method not in {"nearest_core", "ovito"}:
        raise ValueError("orphan_method must be 'nearest_core' or 'ovito'")
    if structure.number_of_atoms == 0:
        raise ValueError("grain segmentation requires a nonempty structure")
    try:
        import ovito
        from ovito.data import DataCollection
        from ovito.modifiers import GrainSegmentationModifier, PolyhedralTemplateMatchingModifier
        from ovito.pipeline import Pipeline, StaticSource
    except ImportError as exc:
        raise ImportError("Grain recovery requires the optional 'ovito' package") from exc

    algorithms = {
        "graph_auto": GrainSegmentationModifier.Algorithm.GraphClusteringAuto,
        "graph_manual": GrainSegmentationModifier.Algorithm.GraphClusteringManual,
        "minimum_spanning_tree": GrainSegmentationModifier.Algorithm.MinimumSpanningTree,
    }
    if algorithm not in algorithms:
        raise ValueError(f"algorithm must be one of {tuple(algorithms)}")

    atom_ids = np.asarray(
        structure.properties.get("atom_ids", np.arange(1, structure.number_of_atoms + 1)),
        dtype=np.int64,
    )
    if atom_ids.shape != (structure.number_of_atoms,) or len(np.unique(atom_ids)) != len(atom_ids):
        raise ValueError("atom_ids must contain one unique identifier per atom")

    data = DataCollection()
    particles = data.create_particles(count=structure.number_of_atoms)
    particles.create_property("Position", data=np.asarray(structure.positions))
    particles.create_property("Particle Identifier", data=atom_ids)
    particles.create_property("Particle Type", data=np.asarray(structure.atom_types))
    matrix = np.column_stack((np.asarray(structure.cell).T, np.asarray(structure.origin)))
    data.create_cell(matrix, pbc=tuple(structure.pbc))

    pipeline = Pipeline(source=StaticSource(data=data))
    pipeline.modifiers.append(
        PolyhedralTemplateMatchingModifier(
            output_orientation=True,
            output_rmsd=True,
            output_interatomic_distance=True,
            rmsd_cutoff=float(rmsd_cutoff),
        )
    )
    segmentation = GrainSegmentationModifier(
        algorithm=algorithms[algorithm],
        min_grain_size=int(min_grain_size),
        handle_stacking_faults=bool(handle_stacking_faults),
        orphan_adoption=False,
        merging_threshold=float(merging_threshold),
        color_particles=False,
    )
    pipeline.modifiers.append(segmentation)
    core = pipeline.compute()
    core_labels = np.asarray(core.particles["Grain"], dtype=np.int64).copy()
    out_ids = np.asarray(core.particles["Particle Identifier"], dtype=np.int64).copy()
    sorted_rows = np.argsort(out_ids)
    rows = sorted_rows[np.searchsorted(out_ids[sorted_rows], atom_ids)]
    if not np.array_equal(out_ids[rows], atom_ids):
        raise RuntimeError("OVITO output atom IDs differ from the input structure")

    if orphan_adoption and orphan_method == "ovito":
        segmentation.orphan_adoption = True
        final = pipeline.compute()
    else:
        final = core
    final_labels = np.asarray(final.particles["Grain"], dtype=np.int64).copy()
    version = tuple(int(v) for v in ovito.version)
    if version < (3, 17, 0):
        core_labels = np.where(core_labels == 0, -1, core_labels)
        final_labels = np.where(final_labels == 0, -1, final_labels)
        native_unknown = 0
    else:
        core_labels = np.where(core_labels < 0, -1, core_labels + 1)
        final_labels = np.where(final_labels < 0, -1, final_labels + 1)
        native_unknown = -1

    core_labels = core_labels[rows]
    final_labels = final_labels[rows]
    structure_types = np.asarray(final.particles["Structure Type"], dtype=np.int64)[rows].copy()
    # Naming anchors stabilize arbitrary size-ranked IDs without changing cores.
    grains = np.unique(core_labels[core_labels > 0])
    def anchor(grain):
        members = core_labels == grain
        fcc = members & (structure_types == 1)
        return int(atom_ids[fcc if np.any(fcc) else members].min())
    grains = sorted(grains.tolist(), key=anchor)
    mapping = {int(grain): index + 1 for index, grain in enumerate(grains)}
    old_core = core_labels.copy()
    old_final = final_labels.copy()
    for previous, current in mapping.items():
        core_labels[old_core == previous] = current
        final_labels[old_final == previous] = current
    assignment_distances = np.where(core_labels < 0, np.nan, 0.0)
    assignment_margins = np.full(len(core_labels), np.nan)
    ties = 0
    if orphan_adoption and orphan_method == "nearest_core":
        final_labels, assignment_distances, assignment_margins, ties = _nearest_grain_core(structure, core_labels)
    grain_count = int(final.attributes["GrainSegmentation.grain_count"])
    if grain_count < 1:
        raise ValueError("PTM grain segmentation found no grains; inspect the structure and segmentation settings")
    diagnostics = {
        "method": "OVITO PTM orientations + grain segmentation",
        "ovito_version": ovito.version_string,
        "parameters": {
            "rmsd_cutoff": float(rmsd_cutoff),
            "min_grain_size": int(min_grain_size),
            "algorithm": algorithm,
            "merging_threshold": float(merging_threshold),
            "handle_stacking_faults": bool(handle_stacking_faults),
            "orphan_adoption": bool(orphan_adoption),
            "orphan_method": orphan_method,
        },
        "atom_count": int(structure.number_of_atoms),
        "grain_count": grain_count,
        "core_grain_counts": _counts(core_labels),
        "final_grain_counts": _counts(final_labels),
        "ptm_structure_counts": _counts(structure_types),
        "unassigned_core_count": int(np.count_nonzero(core_labels < 0)),
        "unassigned_final_count": int(np.count_nonzero(final_labels < 0)),
        "native_unassigned_grain_id": native_unknown,
        "returned_unassigned_grain_id": -1,
        "returned_valid_grain_ids": "positive integers ordered by the smallest core atom ID (FCC naming anchors preferred)",
        "normalized_native_to_returned_grain_ids": mapping,
        "orphan_nearest_core_distance_quantiles": [float(v) for v in np.quantile(assignment_distances[core_labels < 0], [0, 0.5, 0.95, 1])] if orphan_adoption and orphan_method == "nearest_core" and np.any(core_labels < 0) else [],
        "orphan_nearest_distinct_grain_margin_quantiles": [float(v) for v in np.quantile(assignment_margins[np.isfinite(assignment_margins)], [0, 0.5, 0.95, 1])] if np.any(np.isfinite(assignment_margins)) else [],
        "orphan_equidistant_tie_count": ties,
        "orphan_tie_tolerance_A": 1e-8,
        "pbc": [bool(value) for value in structure.pbc],
        "source": structure.metadata.get("source"),
        "surface_geometry_used_for_segmentation": False,
    }
    table = core.tables["grains"]
    table_ids = np.asarray(table["Grain Identifier"], dtype=np.int64)
    if version >= (3, 17, 0):
        table_ids = table_ids + 1
    diagnostics["grain_core_orientations_xyzw"] = {
        str(mapping[int(grain)]): np.asarray(table["Orientation"])[index].tolist()
        for index, grain in enumerate(table_ids)
        if int(grain) in mapping
    }
    return GrainSegmentationResult(
        atom_ids=atom_ids.copy(),
        grain_ids=final_labels,
        core_grain_ids=core_labels,
        structure_types=structure_types,
        rmsd=np.asarray(final.particles["RMSD"], dtype=float)[rows].copy(),
        orientations=np.asarray(final.particles["Orientation"], dtype=float)[rows].copy(),
        interatomic_distances=np.asarray(final.particles["Interatomic Distance"], dtype=float)[rows].copy(),
        assignment_distances=assignment_distances,
        assignment_distance_margins=assignment_margins,
        diagnostics=diagnostics,
    )


__all__ = ["GrainSegmentationResult", "segment_grains_with_ptm"]
