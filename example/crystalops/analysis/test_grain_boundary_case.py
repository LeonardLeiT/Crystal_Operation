"""End-to-end checks for periodic labeled interfaces and Gaussian derivatives."""

import numpy as np

from crystalops.core import CrystalStructure
from crystalops.analysis import analyze_grain_boundary_curvature
from crystalops.analysis.grain_boundary import _PeriodicGaussianField


def _labeled_cube(length=12, step=1.0):
    axis = np.arange(0.23, length, step)
    positions = np.stack(np.meshgrid(axis, axis, axis, indexing="ij"), axis=-1).reshape(-1, 3)
    labels = np.where(positions[:, 0] < length / 2, 1, 2)
    return CrystalStructure(
        positions=positions, cell=np.eye(3) * length,
        atom_types=np.ones(len(positions), dtype=int), type_elements={1: "Cu"},
        properties={"grain_ids": labels}, metadata={"length_unit": "angstrom"},
    )


def test_periodic_planar_interfaces_have_correct_area_and_zero_curvature():
    structure = _labeled_cube()
    original = structure.positions.copy()
    result = analyze_grain_boundary_curvature(
        structure, neighbor_cutoff=1.1, smoothing_length=0.7, sampling_spacing=0.6,
    )
    assert len(result.surfaces) == 1
    assert abs(result.summary["total_interface_area"] - 2 * 12**2) < 1e-3
    assert result.summary["valid_area_fraction"] > 0.999
    assert result.summary["rms_H"] < 1e-7
    assert np.array_equal(original, structure.positions)
    assert np.all(np.isnan(result.atom_mean_curvature[~result.atom_boundary_mask]))
    assert np.any(result.atom_valid)


def test_periodic_field_derivatives_and_box_translations():
    structure = _labeled_cube()
    field = _PeriodicGaussianField(structure, structure.properties["grain_ids"], 0.7)
    points = np.array([[5.33, 2.17, 7.31], [-0.41, 8.13, 1.53]])
    f, g, hessian = field.pair_derivatives(points, (1, 2))
    shifted = points + np.array([12, -12, 24])
    for left, right in zip((f, g, hessian), field.pair_derivatives(shifted, (1, 2))):
        assert np.allclose(left, right, atol=1e-11)
    epsilon = 1e-4
    for axis in range(3):
        offset = np.eye(3)[axis] * epsilon
        plus = field.pair_derivatives(points + offset, (1, 2))
        minus = field.pair_derivatives(points - offset, (1, 2))
        assert np.allclose((plus[0] - minus[0]) / (2 * epsilon), g[:, axis], atol=2e-6)
        assert np.allclose((plus[1] - minus[1]) / (2 * epsilon), hessian[:, :, axis], atol=2e-6)


def test_spherical_labeled_interface_recovers_curvature_without_fitting_a_sphere():
    structure = _labeled_cube(length=16, step=0.8)
    radius = 4.0
    structure.properties["grain_ids"] = np.where(
        np.linalg.norm(structure.positions - 8, axis=1) < radius, 1, 2,
    )
    result = analyze_grain_boundary_curvature(
        structure, neighbor_cutoff=1.0, smoothing_length=0.7, sampling_spacing=0.6,
    )
    # The occupancy interface differs slightly from the prescribed atom-center
    # dividing sphere. This checks the full reconstruction error, not only H(F).
    assert abs(result.summary["mean_signed_H"] - 1 / radius) < 0.15 / radius
    assert abs(result.summary["total_interface_area"] / (4 * np.pi * radius**2) - 1) < 0.15
    assert result.summary["valid_area_fraction"] > 0.999
