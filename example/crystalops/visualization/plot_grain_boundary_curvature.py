"""Plot full grain-boundary meshes and scale sensitivity from saved analysis results."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))


def load_grain_boundary_plot_result(directory: str | Path):
    """Load saved measured arrays without recomputing labels or curvature."""
    directory = Path(directory)
    payload = json.loads((directory / "summary.json").read_text(encoding="utf-8"))
    parameters = dict(payload["parameters"])
    structure = payload.get("structure", {})
    parameters.setdefault("periodic_cell", structure.get("cell_vectors"))
    parameters.setdefault("periodic_origin", structure.get("origin", [0.0, 0.0, 0.0]))
    parameters.setdefault("length_unit", structure.get("coordinate_unit", "input_length_unit"))
    surfaces = []
    with np.load(directory / "surface_arrays.npz", allow_pickle=False) as arrays:
        pairs = np.asarray(arrays["grain_pairs"])
        if pairs.ndim != 2 or pairs.shape[1:] != (2,):
            raise ValueError("grain_pairs must have shape (K, 2)")
        for index, pair in enumerate(pairs):
            prefix = f"surface_{index}_"
            surfaces.append(SimpleNamespace(
                grain_pair=tuple(int(value) for value in pair),
                vertices=np.array(arrays[prefix + "vertices"], copy=True),
                faces=np.array(arrays[prefix + "faces"], copy=True),
                face_mean_curvature=np.array(arrays[prefix + "face_mean_curvature"], copy=True),
                face_valid=np.array(arrays[prefix + "face_valid"], dtype=bool, copy=True),
                face_area=np.array(arrays[prefix + "face_area"], copy=True),
            ))
    return SimpleNamespace(surfaces=surfaces, summary=payload["summary"], parameters=parameters)


def main() -> None:
    """Read two baseline runs and their actual sensitivity comparison table."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dsc-dir", type=Path, required=True, help="baseline DSC result directory")
    parser.add_argument("--psc-dir", type=Path, required=True, help="baseline PSC result directory")
    parser.add_argument("--comparison", type=Path, required=True, help="measured comparison.csv table")
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--sigma", type=float, default=2.2, help="baseline smoothing length")
    parser.add_argument("--spacing", type=float, default=0.75, help="baseline requested sampling spacing")
    options = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from crystalops.visualization.grain_boundary import plot_grain_boundary_comparison

    try:
        results = {
            "DSC": load_grain_boundary_plot_result(options.dsc_dir),
            "PSC": load_grain_boundary_plot_result(options.psc_dir),
        }
        with options.comparison.open("r", encoding="utf-8-sig", newline="") as handle:
            rows = list(csv.DictReader(handle))
        figure, paths = plot_grain_boundary_comparison(
            results, rows, output_prefix=options.output_prefix,
            smoothing_length=options.sigma, sampling_spacing=options.spacing,
        )
    except (ValueError, KeyError, OSError) as error:
        parser.error(str(error))
    plt.close(figure)
    for kind, path in paths.items():
        print(f"{kind}: {path.resolve()}")


if __name__ == "__main__":
    main()
