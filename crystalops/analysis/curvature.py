"""Differential curvature of regular implicit surfaces using NumPy.

For an implicit surface ``F(x) = 0``, the normal is ``grad(F)/|grad(F)|``.
We use ``H = div(normal)/2``: a sphere with outward normal has ``H = +1/R``.
Changing ``F`` to ``-F`` reverses signed mean and principal curvatures while
leaving Gaussian curvature unchanged. These routines evaluate the chosen
implicit surface; they do not reconstruct or smooth an atomic point cloud.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]

__all__ = [
    "mean_curvature_from_derivatives",
    "gaussian_curvature_from_derivatives",
    "principal_curvatures_from_derivatives",
]


def _prepare_derivatives(
    gradients: ArrayLike,
    hessians: ArrayLike,
    gradient_tolerance: float,
) -> tuple[FloatArray, FloatArray, FloatArray, NDArray[np.bool_], tuple[int, ...]]:
    """Validate matching batches and isolate regular, finite sample points."""
    gradient = np.asarray(gradients, dtype=np.float64)
    hessian = np.asarray(hessians, dtype=np.float64)
    if gradient.ndim < 1 or gradient.shape[-1] != 3:
        raise ValueError("gradients must have shape (..., 3)")
    if hessian.shape != gradient.shape[:-1] + (3, 3):
        raise ValueError("hessians must have matching shape (..., 3, 3)")
    tolerance = np.asarray(gradient_tolerance, dtype=np.float64)
    if tolerance.ndim != 0 or not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError("gradient_tolerance must be a finite non-negative scalar")

    leading_shape = gradient.shape[:-1]
    gradient = gradient.reshape((-1, 3))
    hessian = hessian.reshape((-1, 3, 3))
    with np.errstate(over="ignore", invalid="ignore"):
        norms = np.linalg.norm(gradient, axis=-1)
    valid = (
        np.all(np.isfinite(gradient), axis=-1)
        & np.all(np.isfinite(hessian), axis=(-2, -1))
        & np.isfinite(norms)
        & (norms > tolerance)
    )
    # An exact Hessian is symmetric. This removes floating-point asymmetry
    # from approximate derivative evaluations without introducing a new scale.
    with np.errstate(over="ignore", invalid="ignore"):
        hessian = 0.5 * hessian + 0.5 * np.swapaxes(hessian, -1, -2)
    return gradient, hessian, norms, valid, leading_shape


def mean_curvature_from_derivatives(
    gradients: ArrayLike,
    hessians: ArrayLike,
    *,
    gradient_tolerance: float = 1e-12,
) -> FloatArray:
    """Compute signed mean curvature from implicit-field derivatives.

    Args:
        gradients: Cartesian gradients, with shape ``(..., 3)``.
        hessians: Cartesian Hessians with exactly matching leading dimensions
            and shape ``(..., 3, 3)``. Their symmetric part is used.
        gradient_tolerance: Points with gradient norm no larger than this
            threshold are treated as non-regular. The threshold has the units
            of the supplied field gradient and must match its normalization.

    Returns:
        Curvature values with shape ``gradients.shape[:-1]`` and inverse-length
        units. Single-point inputs produce a zero-dimensional NumPy array.
        Non-finite derivatives and non-regular points produce ``NaN``.

    Notes:
        With ``g = grad(F)`` and ``B = Hessian(F)``, the formula is
        ``(|g|**2 * trace(B) - g.T @ B @ g) / (2 * |g|**3)``.
        Its equivalent normal-based form avoids unnecessary powers of ``|g|``.
    """
    gradient, hessian, norms, valid, shape = _prepare_derivatives(
        gradients, hessians, gradient_tolerance
    )
    result = np.full(len(norms), np.nan, dtype=np.float64)
    if np.any(valid):
        normal = gradient[valid] / norms[valid, None]
        second = hessian[valid]
        normal_second = np.einsum("ni,nij,nj->n", normal, second, normal)
        with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
            curvature = (np.trace(second, axis1=-2, axis2=-1) - normal_second) / (
                2.0 * norms[valid]
            )
        result[valid] = np.where(np.isfinite(curvature), curvature, np.nan)
    return result.reshape(shape)


def _principal_curvatures(
    gradient: FloatArray,
    hessian: FloatArray,
    norms: FloatArray,
) -> FloatArray:
    """Restrict the normal derivative to a stable orthonormal tangent frame."""
    normal = gradient / norms[:, None]
    reference = np.eye(3)[np.argmin(np.abs(normal), axis=-1)]
    tangent_one = np.cross(normal, reference)
    tangent_one /= np.linalg.norm(tangent_one, axis=-1)[:, None]
    tangent_two = np.cross(normal, tangent_one)
    tangent = np.stack((tangent_one, tangent_two), axis=-1)
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        shape_operator = np.einsum("nia,nij,njb->nab", tangent, hessian, tangent)
        shape_operator /= norms[:, None, None]
    finite = np.all(np.isfinite(shape_operator), axis=(-2, -1))
    result = np.full((len(norms), 2), np.nan, dtype=np.float64)
    if np.any(finite):
        result[finite] = np.linalg.eigvalsh(shape_operator[finite])
    return result


def principal_curvatures_from_derivatives(
    gradients: ArrayLike,
    hessians: ArrayLike,
    *,
    gradient_tolerance: float = 1e-12,
) -> FloatArray:
    """Return ordered principal curvatures ``k1 <= k2`` with shape ``(..., 2)``.

    Inputs, gradient threshold, invalid values, units, and normal convention
    follow :func:`mean_curvature_from_derivatives`. The two curvatures are the
    eigenvalues of the Hessian restricted to the tangent plane, divided by
    the gradient norm. Under normal reversal their order also reverses:
    ``(k1, k2) -> (-k2, -k1)``.
    """
    gradient, hessian, norms, valid, shape = _prepare_derivatives(
        gradients, hessians, gradient_tolerance
    )
    result = np.full((len(norms), 2), np.nan, dtype=np.float64)
    if np.any(valid):
        result[valid] = _principal_curvatures(gradient[valid], hessian[valid], norms[valid])
    return result.reshape(shape + (2,))


def gaussian_curvature_from_derivatives(
    gradients: ArrayLike,
    hessians: ArrayLike,
    *,
    gradient_tolerance: float = 1e-12,
) -> FloatArray:
    """Compute Gaussian curvature ``K = k1*k2`` with inverse-length-squared units.

    Inputs and invalid-value handling follow
    :func:`mean_curvature_from_derivatives`. Gaussian curvature is independent
    of the orientation of the normal.
    """
    principal = principal_curvatures_from_derivatives(
        gradients, hessians, gradient_tolerance=gradient_tolerance
    )
    with np.errstate(over="ignore", invalid="ignore"):
        result = principal[..., 0] * principal[..., 1]
    return np.where(np.isfinite(result), result, np.nan)
