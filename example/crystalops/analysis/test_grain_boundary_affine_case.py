"""Whole-reconstruction metric and deformation checks on a fixed labeled cloud."""

import numpy as np

from crystalops.analysis import analyze_grain_boundary_curvature
from crystalops.core import CrystalStructure


def _sphere(scale):
    axis = np.arange(0.23, 16, 0.8)
    positions = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    labels = np.where(np.linalg.norm(positions - 8, axis=1) < 4, 1, 2)
    scale = np.asarray(scale, dtype=float)
    return CrystalStructure(
        positions=positions * scale, cell=np.diag(16 * scale),
        atom_types=np.ones(len(positions), dtype=int), type_elements={1: "Cu"},
        properties={"grain_ids": labels},
    )


def _caps(result, scale):
    surface = result.surfaces[0]
    scale = np.asarray(scale)
    direction = (surface.face_centers - 8 * scale) / (4 * scale)
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    values = []
    for axis in (0, 1):
        selected = surface.face_valid & (np.abs(direction[:, axis]) > 0.97)
        values.append(float(np.average(surface.face_mean_curvature[selected], weights=surface.face_area[selected])))
    return values


def test_uniform_scaling_preserves_the_full_reconstruction_metric():
    base = analyze_grain_boundary_curvature(
        _sphere([1, 1, 1]), neighbor_cutoff=1.4, smoothing_length=0.7, sampling_spacing=0.6,
    )
    factor = 1.4
    scaled = analyze_grain_boundary_curvature(
        _sphere([factor] * 3), neighbor_cutoff=1.4 * factor,
        smoothing_length=0.7 * factor, sampling_spacing=0.6 * factor,
    )
    assert np.isclose(scaled.summary["rms_H"] * factor, base.summary["rms_H"], rtol=1e-6)
    assert np.isclose(scaled.summary["total_interface_area"] / factor**2,
                      base.summary["total_interface_area"], rtol=1e-6)


def test_uniaxial_stretch_changes_local_curvature_in_opposite_directions():
    base = analyze_grain_boundary_curvature(
        _sphere([1, 1, 1]), neighbor_cutoff=1.4, smoothing_length=0.7, sampling_spacing=0.6,
    )
    stretched = analyze_grain_boundary_curvature(
        _sphere([1, 1.4, 1]), neighbor_cutoff=1.4, smoothing_length=0.7, sampling_spacing=0.6,
    )
    x0, y0 = _caps(base, [1, 1, 1])
    x1, y1 = _caps(stretched, [1, 1.4, 1])
    # This checks deformation sensitivity, not absolute cap accuracy of a
    # voxel-sampled dividing sphere. The audit separately records that bias.
    assert x1 < 0.9 * x0
    assert y1 > 1.25 * y0
    assert abs(x1 - y1) > 0.1
