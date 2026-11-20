"""Calibrate reconstruction using one fixed labeled catenoid point cloud.

This is a numerical control, not experimental data. The central evaluated
patch is farther than the Gaussian support radius from every periodic seam.
The undeformed catenoid is exactly minimal; the same points and labels are
then stretched along y. Reconstruction error is reported, never set to zero.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from crystalops.analysis import analyze_grain_boundary_curvature, mean_curvature_from_derivatives
from crystalops.core import CrystalStructure


def audit_catenoid_reconstruction(output_dir, *, sigma=1.2, spacing=0.9):
    """Measure local reconstruction error before and after a 1.4x y stretch."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    neck_radius, length, step = 4.0, 24.0, 0.8
    axis = np.arange(-length / 2 + 0.17, length / 2, step)
    points = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    labels = np.where(
        points[:, 0]**2 + points[:, 1]**2
        < neck_radius**2 * np.cosh(points[:, 2] / neck_radius)**2, 1, 2,
    )
    cases = []
    for stretch in (1.0, 1.4):
        started = time.perf_counter()
        scale = np.array([1.0, stretch, 1.0])
        structure = CrystalStructure(
            positions=points * scale, cell=np.diag(length * scale), origin=-length / 2 * scale,
            atom_types=np.ones(len(points), dtype=int), type_elements={1: "Cu"},
            properties={"grain_ids": labels.copy()},
        )
        result = analyze_grain_boundary_curvature(
            structure, neighbor_cutoff=1.6, smoothing_length=sigma, sampling_spacing=spacing,
        )
        surface = result.surfaces[0]
        seam_distance = np.min(np.minimum(
            surface.face_centers - structure.origin,
            structure.origin + np.diag(structure.cell) - surface.face_centers,
        ), axis=1)
        mask = (surface.face_valid & (np.abs(surface.face_centers[:, 2]) < 1.2)
                & (seam_distance > 5 * sigma))
        measured_points = surface.face_centers[mask]
        measured = surface.face_mean_curvature[mask]
        area = surface.face_area[mask]
        if not len(measured):
            raise ValueError("no central control patch remains at this sigma/spacing")
        u = np.arctan2(measured_points[:, 1] / stretch, measured_points[:, 0])
        v = measured_points[:, 2] / neck_radius
        reference_points = np.column_stack((
            neck_radius * np.cosh(v) * np.cos(u),
            stretch * neck_radius * np.cosh(v) * np.sin(u), measured_points[:, 2],
        ))
        gradient = np.column_stack((
            2 * reference_points[:, 0], 2 * reference_points[:, 1] / stretch**2,
            -neck_radius * np.sinh(2 * v),
        ))
        hessian = np.zeros((len(area), 3, 3))
        hessian[:, 0, 0] = 2
        hessian[:, 1, 1] = 2 / stretch**2
        hessian[:, 2, 2] = -2 * np.cosh(2 * v)
        reference = mean_curvature_from_derivatives(gradient, hessian)
        rms = lambda values: float(np.sqrt(np.average(values**2, weights=area)))
        case = {
            "stretch_y": stretch, "atoms": len(points), "neck_radius_A": neck_radius,
            "sigma_A": sigma, "sampling_spacing_A": spacing, "face_count": len(area),
            "reference_rms_H": rms(reference), "reconstructed_rms_H": rms(measured),
            "reference_mean_absolute_H": float(np.average(np.abs(reference), weights=area)),
            "reconstructed_mean_absolute_H": float(np.average(np.abs(measured), weights=area)),
            "rms_local_H_error": rms(measured - reference),
            "mean_H_error": float(np.average(measured - reference, weights=area)),
            "local_H_correlation": float(np.corrcoef(measured, reference)[0, 1]) if stretch != 1 else None,
            "area_A2": float(area.sum()), "seconds": time.perf_counter() - started,
        }
        np.savez_compressed(
            output / f"catenoid_y{stretch:g}.npz", positions=measured_points,
            reference_positions=reference_points, area=area, H=measured, reference_H=reference,
        )
        cases.append(case)
        print(json.dumps(case), flush=True)
    payload = {
        "control": "same point identities and labels; y affine stretch only; no resegmentation",
        "reference": "exact catenoid and its stretched parametric image; implicit derivatives independently checked against fundamental forms",
        "scope": "central patch only, excluding periodic seams by more than 5 sigma",
        "conclusion_limit": "a correct curvature kernel does not certify accurate reconstruction or minimality of supplied atomic structures",
        "cases": cases,
    }
    (output / "catenoid_reconstruction.json").write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("analysis_results/self_audit/end_to_end"))
    parser.add_argument("--sigma", type=float, default=1.2)
    parser.add_argument("--spacing", type=float, default=0.9)
    options = parser.parse_args()
    audit_catenoid_reconstruction(options.output_dir, sigma=options.sigma, spacing=options.spacing)


if __name__ == "__main__":
    main()
