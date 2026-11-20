"""Analyze labeled samples at several physical scales and sampling spacings."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops.core import CrystalStructure
from crystalops.analysis import analyze_grain_boundary_curvature
from crystalops.io import write_grain_boundary_results
from example.crystalops.analysis.analyze_grain_boundary_curvature import load_grain_labels_csv


def save_surface_arrays(result, filename):
    """Retain every face for independent plotting and numerical inspection."""
    arrays = {"grain_pairs": np.asarray([s.grain_pair for s in result.surfaces], dtype=np.int64)}
    for i, surface in enumerate(result.surfaces):
        for name in (
            "vertices", "faces", "face_centers", "face_mean_curvature",
            "face_gaussian_curvature", "face_valid", "face_area", "face_junction",
            "face_projection_distance", "face_projection_residual",
        ):
            arrays[f"surface_{i}_{name}"] = getattr(surface, name)
    np.savez_compressed(filename, **arrays)


def run_sample_comparison(
    input_files, *, labels_dir, output_dir, sigmas=(1.8, 2.2, 2.6),
    spacings=(1.0, 0.75, 0.5), reference_sigma=2.2, reference_spacing=0.75,
    neighbor_cutoff=3.2, dsc_six_grain_labels=None, make_plots=True,
):
    """Run one-factor sensitivity checks, retaining fixed label files per sample."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    rows, baselines = [], {}
    configurations = [(reference_sigma, reference_spacing)]
    configurations += [(float(s), reference_spacing) for s in sigmas]
    configurations += [(reference_sigma, float(h)) for h in spacings]
    configurations = list(dict.fromkeys(configurations))
    comparison_path = output / "comparison.csv"

    def save_rows():
        with comparison_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        (output / "comparison.json").write_text(
            json.dumps(rows, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8"
        )

    for filename in input_files:
        filename = Path(filename)
        structure = CrystalStructure.read(filename)
        structure.metadata["length_unit"] = "angstrom"
        label_file = Path(labels_dir) / f"{filename.stem}_grain_labels.csv"
        labels = load_grain_labels_csv(structure, label_file)
        diagnostics_file = Path(labels_dir) / f"{filename.stem}_grain_diagnostics.json"
        if diagnostics_file.exists():
            structure.metadata["grain_segmentation"] = json.loads(diagnostics_file.read_text(encoding="utf-8"))
        recovered_algorithm = structure.metadata.get("grain_segmentation", {}).get("parameters", {}).get("algorithm")
        primary_scheme = "primary_mst" if recovered_algorithm == "minimum_spanning_tree" else "primary_provided"
        scenarios = [(primary_scheme, labels, label_file, configurations)]
        if filename.stem.upper() == "DSC" and dsc_six_grain_labels:
            alternative = Path(dsc_six_grain_labels)
            scenarios.append((
                "graph_auto_6grain_sensitivity", load_grain_labels_csv(structure, alternative),
                alternative, [(reference_sigma, reference_spacing)],
            ))
        for scheme, current_labels, current_label_file, settings in scenarios:
            scenario_diagnostics = current_label_file.with_name(
                current_label_file.name.replace("_labels.csv", "_diagnostics.json")
            )
            structure.metadata["grain_segmentation"] = (
                json.loads(scenario_diagnostics.read_text(encoding="utf-8"))
                if scenario_diagnostics.exists() else {
                    "method": "provided atom-ID aligned label file",
                    "source": str(current_label_file.resolve()),
                    "grain_count": int(len(np.unique(current_labels))),
                }
            )
            for sigma, spacing in settings:
                run_name = f"sigma_{sigma:g}_spacing_{spacing:g}".replace(".", "p")
                sample_name = filename.stem if scheme.startswith("primary_") else filename.stem + "_graph_auto_6grain"
                directory = output / sample_name / run_name
                print(f"START {sample_name}: sigma={sigma:g} A, spacing={spacing:g} A", flush=True)
                started = time.perf_counter()
                source = f"atom-ID aligned {current_label_file.resolve()} ({scheme})"
                structure.metadata["grain_id_source"] = source
                result = analyze_grain_boundary_curvature(
                    structure, grain_ids=current_labels, neighbor_cutoff=neighbor_cutoff,
                    smoothing_length=sigma, sampling_spacing=spacing, grain_id_source=source,
                )
                result.parameters.update(
                    input_sha256=hashlib.sha256(filename.read_bytes()).hexdigest(),
                    grain_labels_sha256=hashlib.sha256(current_label_file.read_bytes()).hexdigest(),
                    grain_scheme=scheme,
                )
                paths = write_grain_boundary_results(result, structure, directory, grain_ids=current_labels)
                save_surface_arrays(result, directory / "surface_arrays.npz")
                row = {
                    "structure": filename.stem, "grain_scheme": scheme,
                    "grain_count": result.summary["grain_count"],
                    "smoothing_length": float(sigma), "sampling_spacing": float(spacing),
                    **{key: value for key, value in result.summary.items() if key != "per_pair"},
                    "seconds": round(time.perf_counter() - started, 3),
                    "summary_file": str(paths["summary"].relative_to(output)),
                }
                rows.append(row)
                save_rows()
                if scheme.startswith("primary_") and (sigma, spacing) == (reference_sigma, reference_spacing):
                    baselines[filename.stem] = result
                absolute = result.summary['mean_absolute_H']
                rms = result.summary['rms_H']
                absolute_text = f"{absolute:.7g}" if absolute is not None else "undefined"
                rms_text = f"{rms:.7g}" if rms is not None else "undefined"
                print(
                    f"DONE {sample_name}: <|H|>={absolute_text} A^-1, "
                    f"RMS={rms_text} A^-1, "
                    f"valid area={result.summary['valid_area_fraction']:.2%}, {row['seconds']:.1f}s",
                    flush=True,
                )
    if make_plots:
        from crystalops.visualization.grain_boundary import plot_grain_boundary_comparison

        figure, _ = plot_grain_boundary_comparison(
            baselines, rows, output_prefix=output / "curvature_comparison",
            smoothing_length=reference_sigma, sampling_spacing=reference_spacing,
        )
        import matplotlib.pyplot as plt
        plt.close(figure)
    return baselines, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--labels-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--sigmas", nargs="+", type=float, default=(1.8, 2.2, 2.6))
    parser.add_argument("--spacings", nargs="+", type=float, default=(1.0, 0.75, 0.5))
    parser.add_argument("--reference-sigma", type=float, default=2.2)
    parser.add_argument("--reference-spacing", type=float, default=0.75)
    parser.add_argument("--neighbor-cutoff", type=float, default=3.2)
    parser.add_argument("--dsc-six-grain-labels", type=Path)
    parser.add_argument("--skip-plots", action="store_true")
    options = parser.parse_args()
    run_sample_comparison(
        options.inputs, labels_dir=options.labels_dir, output_dir=options.output_dir,
        sigmas=options.sigmas, spacings=options.spacings,
        reference_sigma=options.reference_sigma, reference_spacing=options.reference_spacing,
        neighbor_cutoff=options.neighbor_cutoff, dsc_six_grain_labels=options.dsc_six_grain_labels,
        make_plots=not options.skip_plots,
    )


if __name__ == "__main__":
    main()
