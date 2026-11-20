"""Read a labeled periodic structure, reconstruct its interfaces, and export curvature."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops.core import CrystalStructure
from crystalops.io.grain_boundary_results import write_grain_boundary_results


def load_grain_labels_csv(structure: CrystalStructure, filename: str | Path) -> np.ndarray:
    """Align ``atom_id,grain_id`` CSV labels by IDs, never by file row order."""
    by_id = {}
    with Path(filename).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not {"atom_id", "grain_id"} <= set(reader.fieldnames):
            raise ValueError("grain-label CSV must contain atom_id and grain_id columns")
        for row in reader:
            try:
                atom_id = int(row["atom_id"])
                grain_id = int(row["grain_id"])
            except (ValueError, TypeError) as error:
                raise ValueError("grain-label atom_id and grain_id entries must be integers") from error
            if atom_id in by_id:
                raise ValueError(f"duplicate atom_id in grain-label CSV: {atom_id}")
            by_id[atom_id] = grain_id
    atom_ids = np.asarray(
        structure.properties.get("atom_ids", np.arange(1, structure.number_of_atoms + 1)),
        dtype=np.int64,
    )
    if atom_ids.shape != (structure.number_of_atoms,) or len(np.unique(atom_ids)) != len(atom_ids):
        raise ValueError("input structure must contain one unique atom ID per atom")
    input_ids = set(int(atom_id) for atom_id in atom_ids)
    missing = input_ids - set(by_id)
    extra = set(by_id) - input_ids
    if missing or extra:
        raise ValueError(
            f"grain-label CSV IDs do not match structure: {len(missing)} missing, {len(extra)} extra"
        )
    return np.array([by_id[int(atom_id)] for atom_id in atom_ids], dtype=np.int64)


def analyze_grain_boundary_file(
    filename: str | Path,
    *,
    output_dir: str | Path,
    smoothing_length: float,
    neighbor_cutoff: float,
    sampling_spacing: float | None = None,
    grain_labels: str | Path | None = None,
    segment_grains: bool = False,
    rmsd_cutoff: float = 0.1,
    min_grain_size: int = 50,
    grain_algorithm: str = "graph_auto",
    merging_threshold: float = 0.0,
    handle_stacking_faults: bool = True,
    junction_threshold: float = 0.1,
    kernel_cutoff: float = 5.0,
    length_unit: str = "angstrom",
):
    """Execute one reproducible analysis without rewriting input structure/label files.

    Length parameters use the input coordinates' unit; ``length_unit`` labels
    outputs and performs no coordinate conversion. The supplied Cu examples use
    angstroms. Provide labels explicitly or opt into PTM grain segmentation.
    """
    from crystalops.analysis import analyze_grain_boundary_curvature

    if grain_labels is not None and segment_grains:
        raise ValueError("choose either grain_labels or segment_grains")
    structure = CrystalStructure.read(filename)
    structure.metadata["length_unit"] = length_unit
    segmentation = None
    if grain_labels is not None:
        labels = load_grain_labels_csv(structure, grain_labels)
        source = f"atom-ID aligned CSV: {Path(grain_labels).resolve()}"
    elif segment_grains:
        from crystalops.analysis.grain_segmentation import segment_grains_with_ptm

        segmentation = segment_grains_with_ptm(
            structure, rmsd_cutoff=rmsd_cutoff, min_grain_size=min_grain_size,
            algorithm=grain_algorithm, merging_threshold=merging_threshold,
            handle_stacking_faults=handle_stacking_faults, orphan_adoption=True,
        )
        labels = segmentation.grain_ids
        structure.metadata["grain_segmentation"] = segmentation.diagnostics
        source = "OVITO PTM orientations + grain segmentation with orphan adoption"
    elif "grain_ids" in structure.properties:
        labels = np.asarray(structure.properties["grain_ids"])
        source = "input structure grain_ids property"
    else:
        raise ValueError(
            "input has no grain_ids; provide --grain-labels CSV or request --segment-grains"
        )
    structure.properties["grain_ids"] = np.asarray(labels).copy()
    structure.metadata["grain_id_source"] = source
    result = analyze_grain_boundary_curvature(
        structure, neighbor_cutoff=neighbor_cutoff, smoothing_length=smoothing_length,
        grain_ids=labels, sampling_spacing=sampling_spacing,
        junction_threshold=junction_threshold, kernel_cutoff=kernel_cutoff,
        grain_id_source=source,
    )
    paths = write_grain_boundary_results(result, structure, output_dir, grain_ids=labels)
    # Save the actual labels used, including pre-adoption cores when recovered.
    label_path = Path(output_dir) / "grain_labels.csv"
    atom_ids = structure.properties.get("atom_ids", np.arange(1, structure.number_of_atoms + 1))
    with label_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        if segmentation is None:
            writer.writerow(("atom_id", "grain_id"))
            writer.writerows(zip(atom_ids, labels))
        else:
            writer.writerow(("atom_id", "grain_id", "core_grain_id", "ptm_structure_type", "ptm_rmsd"))
            writer.writerows(zip(
                atom_ids, labels, segmentation.core_grain_ids,
                segmentation.structure_types, segmentation.rmsd,
            ))
    paths["grain_labels"] = label_path
    payload = json.loads(paths["summary"].read_text(encoding="utf-8"))
    payload["files"]["grain_labels"] = label_path.name
    paths["summary"].write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return result, paths


def main() -> None:
    """Parse one physical scale and write audited numeric/mesh results."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="input structure, currently LAMMPS atomic data")
    labels = parser.add_mutually_exclusive_group()
    labels.add_argument("--grain-labels", type=Path, help="CSV with atom_id,grain_id columns")
    labels.add_argument("--segment-grains", action="store_true", help="recover current grains with OVITO PTM")
    parser.add_argument("--sigma", type=float, required=True, help="Gaussian smoothing length in coordinate units")
    parser.add_argument("--neighbor-cutoff", type=float, required=True, help="unlike-grain contact distance")
    parser.add_argument("--spacing", type=float, default=None, help="auxiliary surface-sampling spacing; default sigma/3")
    parser.add_argument("--junction-threshold", type=float, default=0.1)
    parser.add_argument("--kernel-cutoff", type=float, default=5.0, help="Gaussian support radius in multiples of sigma")
    parser.add_argument("--rmsd-cutoff", type=float, default=0.1)
    parser.add_argument("--min-grain-size", type=int, default=50)
    parser.add_argument(
        "--grain-algorithm", choices=("graph_auto", "graph_manual", "minimum_spanning_tree"),
        default="graph_auto",
    )
    parser.add_argument("--merging-threshold", type=float, default=0.0)
    parser.add_argument("--no-handle-stacking-faults", dest="handle_stacking_faults", action="store_false")
    parser.add_argument(
        "--length-unit", default="angstrom",
        help="output unit label, no conversion; default angstrom for the supplied Cu samples",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    options = parser.parse_args()
    try:
        result, paths = analyze_grain_boundary_file(
            options.input, output_dir=options.output_dir,
            smoothing_length=options.sigma, neighbor_cutoff=options.neighbor_cutoff,
            sampling_spacing=options.spacing, grain_labels=options.grain_labels,
            segment_grains=options.segment_grains, rmsd_cutoff=options.rmsd_cutoff,
            min_grain_size=options.min_grain_size, grain_algorithm=options.grain_algorithm,
            merging_threshold=options.merging_threshold,
            handle_stacking_faults=options.handle_stacking_faults,
            junction_threshold=options.junction_threshold, kernel_cutoff=options.kernel_cutoff,
            length_unit=options.length_unit,
        )
    except (ValueError, ImportError, NotImplementedError) as error:
        parser.error(str(error))
    summary = result.summary
    print(
        f"Grains: {summary['grain_count']}; interface pairs: {summary['interface_pair_count']}; "
        f"valid area fraction: {summary['valid_area_fraction']:.3%}"
    )
    for key in ("mean_absolute_H", "rms_H"):
        value = summary[key]
        print(f"{key}: {value:.8g} {options.length_unit}^-1" if value is not None else f"{key}: undefined")
    signed = summary["mean_signed_H"]
    print(
        f"Network signed area average H: {signed:.8g} {options.length_unit}^-1 "
        "(normal-convention dependent; not a minimal-surface test)"
        if signed is not None else "Network signed area average H: undefined"
    )
    cancellation = summary.get("signed_cancellation_fraction")
    if cancellation is not None:
        print(f"Positive/negative curvature cancellation: {cancellation:.2%}")
    print("VTK curvature is a face CELL_DATA attribute; NaN is invalid, not zero. Select valid=1.")
    print(f"Results: {paths['summary'].resolve()}")


if __name__ == "__main__":
    main()
