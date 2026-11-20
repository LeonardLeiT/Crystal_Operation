"""Calculate area-weighted mean curvature from a saved DXA defect mesh."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops.analysis import compute_defect_mesh_curvature
from crystalops.io import read_defect_mesh, write_defect_mesh_curvature


def analyze_defect_mesh_file(
    filename: str | Path, *, output_dir: str | Path, fit_radius: float | None = None,
):
    """Use the saved periodic geometry and export curvature and valid areas."""
    mesh = read_defect_mesh(filename)
    result = compute_defect_mesh_curvature(mesh, fit_radius=fit_radius)
    return result, write_defect_mesh_curvature(mesh, result, output_dir)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="defect_mesh.npz exported by the DXA example")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--fit-radius", type=float,
        help="graph-geodesic radius of the local quadratic fit, in the input mesh length unit",
    )
    options = parser.parse_args()
    try:
        result, paths = analyze_defect_mesh_file(
            options.input, output_dir=options.output_dir, fit_radius=options.fit_radius,
        )
    except (ImportError, OSError, ValueError) as exc:
        parser.error(str(exc))

    summary = result.summary
    if summary["rms_H"] is None:
        print("Curvature is undefined because the mesh has no valid integration area.")
    else:
        unit = result.parameters["curvature_unit"]
        print(f"Mean absolute H: {summary['mean_abs_H']:.8g} {unit}")
        print(f"RMS H: {summary['rms_H']:.8g} {unit}")
        print(f"Valid area: {summary['valid_area_fraction']:.2%}")
    for kind, path in paths.items():
        print(f"{kind}: {path.resolve()}")


if __name__ == "__main__":
    main()
