"""Exact minimal-surface and affine-stretch curvature regression cases.

These tests use a true catenoid, not a trigonometric TPMS approximation.
The independent reference uses parametric first/second fundamental forms.
No reference calculation calls the implicit gradient/Hessian curvature kernel.
"""

from __future__ import annotations

import numpy as np

from crystalops.analysis.curvature import (
    gaussian_curvature_from_derivatives,
    mean_curvature_from_derivatives,
    principal_curvatures_from_derivatives,
)


def _catenoid_samples(affine=None, *, a=10.0, nu=144, nv=64):
    """Sample an affine image of X=(a cosh(v) cos(u),a cosh(v) sin(u),a v)."""
    affine = np.eye(3) if affine is None else np.asarray(affine, dtype=float)
    u = (np.arange(nu) + 0.5) * 2 * np.pi / nu
    v = -1.2 + (np.arange(nv) + 0.5) * 2.4 / nv
    u, v = np.meshgrid(u, v, indexing="ij")
    c, s, co, si = np.cosh(v), np.sinh(v), np.cos(u), np.sin(u)
    zero = np.zeros_like(u)
    position = np.stack((a * c * co, a * c * si, a * v), axis=-1)

    # F0=x^2+y^2-a^2*cosh(z/a)^2; F_A(x)=F0(A^-1 x).
    base_gradient = np.stack((2 * position[..., 0], 2 * position[..., 1], -a * np.sinh(2 * v)), axis=-1)
    base_hessian = np.zeros(u.shape + (3, 3))
    base_hessian[..., 0, 0] = 2
    base_hessian[..., 1, 1] = 2
    base_hessian[..., 2, 2] = -2 * np.cosh(2 * v)
    inverse = np.linalg.inv(affine)
    gradient = base_gradient @ inverse
    hessian = np.einsum("ai,...ij,jb->...ab", inverse.T, base_hessian, inverse)

    # Independently differentiate the parametric surface, then transform tangents.
    xu = np.stack((-a * c * si, a * c * co, zero), axis=-1) @ affine.T
    xv = np.stack((a * s * co, a * s * si, np.full_like(u, a)), axis=-1) @ affine.T
    xuu = np.stack((-a * c * co, -a * c * si, zero), axis=-1) @ affine.T
    xuv = np.stack((-a * s * si, a * s * co, zero), axis=-1) @ affine.T
    xvv = np.stack((a * c * co, a * c * si, zero), axis=-1) @ affine.T
    cross = np.cross(xu, xv)
    jacobian = np.linalg.norm(cross, axis=-1)
    normal = cross / jacobian[..., None]

    def dot(left, right):
        return np.sum(left * right, axis=-1)

    first_e, first_f, first_g = dot(xu, xu), dot(xu, xv), dot(xv, xv)
    second_e, second_f, second_g = dot(normal, xuu), dot(normal, xuv), dot(normal, xvv)
    # With outward n, n.dot(X_uu) gives the opposite sign to H=div(n)/2.
    reference_h = -(second_e * first_g - 2 * second_f * first_f + second_g * first_e) / (
        2 * (first_e * first_g - first_f**2)
    )
    reference_k = (second_e * second_g - second_f**2) / (first_e * first_g - first_f**2)
    return {
        "position": position @ affine.T,
        "gradient": gradient,
        "hessian": hessian,
        "area_weight": jacobian,
        "parametric_mean_curvature": reference_h,
        "parametric_gaussian_curvature": reference_k,
    }


def _weighted_metrics(values, weights):
    return {
        "mean_signed_H": float(np.average(values, weights=weights)),
        "mean_absolute_H": float(np.average(np.abs(values), weights=weights)),
        "rms_H": float(np.sqrt(np.average(values**2, weights=weights))),
        "maximum_absolute_H": float(np.max(np.abs(values))),
    }


def test_exact_catenoid_is_locally_minimal():
    samples = _catenoid_samples()
    mean = mean_curvature_from_derivatives(samples["gradient"], samples["hessian"])
    assert np.max(np.abs(mean)) < 1e-14
    assert np.allclose(mean, samples["parametric_mean_curvature"], atol=1e-14)


def test_anisotropic_stretch_of_catenoid_is_not_forced_to_zero():
    samples = _catenoid_samples(np.diag([1, 1.4, 1]))
    mean = mean_curvature_from_derivatives(samples["gradient"], samples["hessian"])
    gaussian = gaussian_curvature_from_derivatives(samples["gradient"], samples["hessian"])
    metrics = _weighted_metrics(mean, samples["area_weight"])
    assert metrics["rms_H"] > 0.005
    assert metrics["mean_absolute_H"] > 0.005
    assert np.allclose(mean, samples["parametric_mean_curvature"], atol=1e-13, rtol=1e-11)
    assert np.allclose(gaussian, samples["parametric_gaussian_curvature"], atol=1e-13, rtol=1e-11)


def test_uniform_scaling_of_nonminimal_surface_rescales_curvatures():
    stretch = np.diag([1, 1.4, 1])
    scale = 1.7
    original = _catenoid_samples(stretch)
    scaled = _catenoid_samples(scale * stretch)
    for function, exponent in (
        (mean_curvature_from_derivatives, 1),
        (principal_curvatures_from_derivatives, 1),
        (gaussian_curvature_from_derivatives, 2),
    ):
        before = function(original["gradient"], original["hessian"])
        after = function(scaled["gradient"], scaled["hessian"])
        assert np.allclose(after, before / scale**exponent, atol=1e-13, rtol=1e-11)


def test_general_affine_gradient_transform_matches_parametric_geometry():
    affine = np.array([[1.1, 0.15, 0.0], [0.0, 1.4, 0.1], [0.0, 0.0, 0.8]])
    samples = _catenoid_samples(affine)
    mean = mean_curvature_from_derivatives(samples["gradient"], samples["hessian"])
    gaussian = gaussian_curvature_from_derivatives(samples["gradient"], samples["hessian"])
    assert np.allclose(mean, samples["parametric_mean_curvature"], atol=1e-13, rtol=1e-11)
    assert np.allclose(gaussian, samples["parametric_gaussian_curvature"], atol=1e-13, rtol=1e-11)


def test_signed_area_mean_can_vanish_for_stretched_nonminimal_periodic_surface():
    # z=A*cos(k*x) is periodic in x and extruded in y. Its H has positive and
    # negative lobes; the area integral cancels over a period at any x stretch.
    theta = (np.arange(1024) + 0.5) * 2 * np.pi / 1024
    amplitude = 2.0
    for stretch in (1.0, 1.4):
        k = 0.2 / stretch
        slope = amplitude * k * np.sin(theta)
        gradient = np.stack((slope, np.zeros_like(theta), np.ones_like(theta)), axis=-1)
        hessian = np.zeros((len(theta), 3, 3))
        hessian[:, 0, 0] = amplitude * k**2 * np.cos(theta)
        mean = mean_curvature_from_derivatives(gradient, hessian)
        reference_h = amplitude * k**2 * np.cos(theta) / (2 * (1 + slope**2)**1.5)
        weights = np.sqrt(1 + slope**2)
        metrics = _weighted_metrics(mean, weights)
        assert np.allclose(mean, reference_h, atol=1e-14)
        assert abs(metrics["mean_signed_H"]) < 1e-14
        assert metrics["rms_H"] > 0.01
        assert metrics["mean_absolute_H"] > 0.01
