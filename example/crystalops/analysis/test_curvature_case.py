"""Analytical and invariance checks for implicit-surface curvature."""

import numpy as np

from crystalops.analysis.curvature import (
    gaussian_curvature_from_derivatives,
    mean_curvature_from_derivatives,
    principal_curvatures_from_derivatives,
)


def test_plane_has_zero_curvature() -> None:
    gradient = np.array([[1.0, 2.0, -3.0], [0.0, 0.0, 1.0]])
    hessian = np.zeros((2, 3, 3))
    assert np.array_equal(mean_curvature_from_derivatives(gradient, hessian), [0.0, 0.0])
    assert np.array_equal(gaussian_curvature_from_derivatives(gradient, hessian), [0.0, 0.0])
    assert np.array_equal(principal_curvatures_from_derivatives(gradient, hessian), np.zeros((2, 2)))


def test_sphere_matches_radius_at_arbitrary_normals() -> None:
    radii = np.array([0.8, 2.0, 7.5, 15.0])
    direction = np.array([[1, 0, 0], [1, 2, 3], [-2, 3, 1], [0, 0, -1]], dtype=float)
    direction /= np.linalg.norm(direction, axis=-1)[:, None]
    position = radii[:, None] * direction
    # F = |x|^2 - R^2, evaluated on each sphere.
    gradient = 2.0 * position
    hessian = np.broadcast_to(2.0 * np.eye(3), (len(radii), 3, 3))
    assert np.allclose(mean_curvature_from_derivatives(gradient, hessian), 1.0 / radii)
    assert np.allclose(gaussian_curvature_from_derivatives(gradient, hessian), 1.0 / radii**2)
    assert np.allclose(principal_curvatures_from_derivatives(gradient, hessian), (1.0 / radii)[:, None])


def test_cylinder_has_one_zero_principal_curvature() -> None:
    radius = 3.7
    angle = np.linspace(0.0, 2.0 * np.pi, 9, endpoint=False)
    gradient = np.stack((2 * radius * np.cos(angle), 2 * radius * np.sin(angle), np.zeros_like(angle)), axis=-1)
    hessian = np.broadcast_to(np.diag([2.0, 2.0, 0.0]), (len(angle), 3, 3))
    assert np.allclose(mean_curvature_from_derivatives(gradient, hessian), 1.0 / (2.0 * radius))
    assert np.allclose(gaussian_curvature_from_derivatives(gradient, hessian), 0.0, atol=1e-14)
    assert np.allclose(principal_curvatures_from_derivatives(gradient, hessian), [0.0, 1.0 / radius])


def test_saddle_has_zero_mean_but_negative_gaussian_curvature() -> None:
    # F = z - (x^2 - y^2)/2, evaluated at the origin.
    gradient = np.array([0.0, 0.0, 1.0])
    hessian = np.diag([-1.0, 1.0, 0.0])
    assert np.isclose(mean_curvature_from_derivatives(gradient, hessian), 0.0)
    assert np.isclose(gaussian_curvature_from_derivatives(gradient, hessian), -1.0)
    assert np.allclose(principal_curvatures_from_derivatives(gradient, hessian), [-1.0, 1.0])


def test_field_scaling_and_normal_reversal() -> None:
    gradient = np.array([[1.0, -2.0, 0.5], [-0.2, 0.5, 3.0]])
    hessian = np.array([[[2, 0.1, 0.2], [0.1, -0.5, 0.4], [0.2, 0.4, 1.3]], np.diag([1.0, 2.0, 3.0])])
    mean = mean_curvature_from_derivatives(gradient, hessian)
    gaussian = gaussian_curvature_from_derivatives(gradient, hessian)
    principal = principal_curvatures_from_derivatives(gradient, hessian)
    for factor in (1e-7, 2.3, 1e7):
        assert np.allclose(mean_curvature_from_derivatives(factor * gradient, factor * hessian), mean)
        assert np.allclose(gaussian_curvature_from_derivatives(factor * gradient, factor * hessian), gaussian)
        assert np.allclose(principal_curvatures_from_derivatives(factor * gradient, factor * hessian), principal)
    assert np.allclose(mean_curvature_from_derivatives(-gradient, -hessian), -mean)
    assert np.allclose(gaussian_curvature_from_derivatives(-gradient, -hessian), gaussian)
    assert np.allclose(principal_curvatures_from_derivatives(-gradient, -hessian), -principal[:, ::-1])


def test_field_reparameterization_preserves_surface_curvature() -> None:
    # Replacing F by exp(F)-1 at F=0 adds g g^T to its Hessian.
    gradient = np.array([[2.0, 1.0, -0.5], [0.0, 0.0, 3.0]])
    hessian = np.broadcast_to(np.diag([1.0, 2.0, 3.0]), (2, 3, 3)).copy()
    changed = hessian + np.einsum("ni,nj->nij", gradient, gradient)
    assert np.allclose(mean_curvature_from_derivatives(gradient, changed), mean_curvature_from_derivatives(gradient, hessian))
    assert np.allclose(gaussian_curvature_from_derivatives(gradient, changed), gaussian_curvature_from_derivatives(gradient, hessian))


def test_batch_shape_and_invalid_points() -> None:
    gradient = np.zeros((2, 2, 3))
    gradient[0, 0] = [1.0, 0.0, 0.0]
    gradient[0, 1] = [1e-14, 0.0, 0.0]
    gradient[1, 0] = [np.nan, 0.0, 0.0]
    gradient[1, 1] = [1.0, 0.0, 0.0]
    hessian = np.zeros((2, 2, 3, 3))
    hessian[1, 1, 0, 0] = np.inf
    for function in (mean_curvature_from_derivatives, gaussian_curvature_from_derivatives):
        result = function(gradient, hessian)
        assert result.shape == (2, 2)
        assert result[0, 0] == 0.0
        assert np.isnan(result[0, 1])
        assert np.all(np.isnan(result[1]))
    principal = principal_curvatures_from_derivatives(gradient, hessian)
    assert principal.shape == (2, 2, 2)
    assert np.array_equal(principal[0, 0], [0.0, 0.0])
    assert np.all(np.isnan(principal[0, 1]))
    assert mean_curvature_from_derivatives([0, 0, 1], np.eye(3)).shape == ()
    assert mean_curvature_from_derivatives(np.empty((0, 3)), np.empty((0, 3, 3))).shape == (0,)


def test_invalid_shapes_and_thresholds_raise() -> None:
    invalid = [
        (1.0, np.eye(3), 1e-12),
        ([1.0, 2.0], np.eye(3), 1e-12),
        (np.ones((2, 3)), np.eye(3), 1e-12),
        (np.ones(3), np.ones((3, 2)), 1e-12),
        (np.ones(3), np.eye(3), -1.0),
        (np.ones(3), np.eye(3), np.nan),
        (np.ones(3), np.eye(3), [1e-12, 1e-12]),
    ]
    for function in (mean_curvature_from_derivatives, gaussian_curvature_from_derivatives, principal_curvatures_from_derivatives):
        for gradient, hessian, threshold in invalid:
            try:
                function(gradient, hessian, gradient_tolerance=threshold)
            except ValueError:
                pass
            else:
                raise AssertionError("invalid derivatives or tolerance were accepted")
