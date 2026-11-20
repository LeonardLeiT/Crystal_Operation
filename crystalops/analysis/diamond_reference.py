"""Exact Schwarz Diamond field and derivatives in a rectangular period cell."""

from __future__ import annotations

import numpy as np


def exact_diamond_field(points, lengths, *, origin=None, phase=None, order=0):
    """Evaluate the exact D zero set for one conventional cell in ``lengths``.

    The physical Cartesian points have shape (N, 3). ``origin`` is the box
    corner and ``phase`` is a dimensionless translation in fractions of its
    lengths. Return F, (F, gradient), or (F, gradient, Hessian) for orders
    0, 1, or 2. The latter two differentiate with respect to physical lengths.
    The homogeneous Jacobi equation clears the poles of sn*dn/cn.
    """
    from scipy.special import ellipj, ellipk

    points = np.asarray(points, dtype=float)
    lengths = np.asarray(lengths, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (N, 3)")
    if lengths.shape != (3,) or not np.isfinite(lengths).all() or np.any(lengths <= 0):
        raise ValueError("lengths must contain three finite positive period lengths")
    if order not in (0, 1, 2):
        raise ValueError("order must be 0, 1, or 2")

    relative = (points - (0 if origin is None else origin)) / lengths
    if phase is not None:
        relative -= phase
    m = 0.25
    period = 4 * ellipk(m)
    sn, cn, dn, _ = ellipj(period * (relative % 1), m)
    s = sn * dn
    value = s.prod(axis=1) - cn.prod(axis=1)
    if order == 0:
        return value

    slope = period / lengths
    ds = cn * (1 - 2 * m * sn**2) * slope
    dc = -s * slope
    gradient = np.empty_like(points)
    for axis in range(3):
        i, j = [other for other in range(3) if other != axis]
        gradient[:, axis] = ds[:, axis] * s[:, i] * s[:, j] - dc[:, axis] * cn[:, i] * cn[:, j]
    if order == 1:
        return value, gradient

    dds = -s * (1 + 4 * m - 6 * m * sn**2) * slope**2
    ddc = -ds * slope
    hessian = np.empty((len(points), 3, 3))
    for axis in range(3):
        i, j = [other for other in range(3) if other != axis]
        hessian[:, axis, axis] = dds[:, axis] * s[:, i] * s[:, j] - ddc[:, axis] * cn[:, i] * cn[:, j]
        hessian[:, i, j] = hessian[:, j, i] = ds[:, i] * ds[:, j] * s[:, axis] - dc[:, i] * dc[:, j] * cn[:, axis]
    return value, gradient, hessian
