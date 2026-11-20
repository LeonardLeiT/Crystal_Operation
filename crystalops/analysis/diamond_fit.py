"""Phase registration of a fixed exact Schwarz D template to atom coordinates."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .diamond_reference import exact_diamond_field


def _minimum_image(delta, lengths):
    return delta - lengths * np.rint(delta / lengths)


def _project_level_set(points, lengths):
    points = np.asarray(points, dtype=float).copy()
    for _ in range(5):
        value, gradient = exact_diamond_field(points, lengths, order=1)
        points -= value[:, None] * gradient / np.sum(gradient**2, axis=1)[:, None]
    return points


@dataclass(slots=True)
class DiamondTemplate:
    """Exact level-set samples, triangle quadrature and periodic sample index."""

    lengths: np.ndarray
    vertices: np.ndarray
    faces: np.ndarray
    centers: np.ndarray
    areas: np.ndarray
    tree: object
    resolution: int


def sample_diamond_template(lengths, *, resolution=64):
    """Sample one full conventional cell for registration and area integration."""
    from scipy.spatial import cKDTree
    from skimage.measure import marching_cubes

    lengths = np.asarray(lengths, dtype=float)
    if lengths.shape != (3,) or not np.isfinite(lengths).all() or np.any(lengths <= 0):
        raise ValueError("lengths must be three finite positive periods")
    if not isinstance(resolution, int) or resolution < 16:
        raise ValueError("resolution must be an integer of at least 16")
    offset = np.array([0.137271, 0.294832, 0.081117])
    axis_points = (np.arange(resolution)[:, None] / resolution + offset) * lengths
    grid = np.stack(np.meshgrid(*axis_points.T, indexing="ij"), axis=-1)
    field = exact_diamond_field(grid.reshape(-1, 3), lengths).reshape(grid.shape[:-1])
    vertices, faces, _, _ = marching_cubes(np.pad(field, (0, 1), mode="wrap"), level=0)
    vertices, inverse = np.unique(vertices % resolution, axis=0, return_inverse=True)
    vertices = _project_level_set((vertices.astype(float) / resolution + offset) * lengths, lengths)
    vertices = np.mod(vertices, lengths)
    faces = inverse.ravel()[faces]
    triangles = vertices[faces]
    first = _minimum_image(triangles[:, 1] - triangles[:, 0], lengths)
    second = _minimum_image(triangles[:, 2] - triangles[:, 0], lengths)
    areas = np.linalg.norm(np.cross(first, second), axis=1) / 2
    centers = _project_level_set(triangles[:, 0] + (first + second) / 3, lengths)
    return DiamondTemplate(lengths, vertices, faces, centers, areas, cKDTree(vertices, boxsize=lengths), resolution)


def _phase_candidates(points, template, width, count=8):
    """Coarse periodic cross-correlation of an atom histogram and template proximity."""
    resolution = 24
    half_lengths = template.lengths / 2
    step = half_lengths / resolution
    indices = np.rint(points / step).astype(int) % resolution
    counts = np.zeros((resolution,) * 3)
    np.add.at(counts, tuple(indices.T), 1)
    axes = np.arange(resolution)[:, None] * step
    grid = np.stack(np.meshgrid(*axes.T, indexing="ij"), axis=-1)
    distance, _ = template.tree.query(grid.reshape(-1, 3))
    proximity = np.exp(-0.5 * (distance / width)**2).reshape(counts.shape)
    correlation = np.fft.ifftn(np.fft.fftn(counts) * np.fft.fftn(proximity).conj()).real
    selected = []
    separation = 2.5 * max(step)
    for index in np.argsort(correlation, axis=None)[::-1]:
        candidate = np.asarray(np.unravel_index(index, counts.shape)) * step
        if all(np.linalg.norm(_minimum_image(candidate - previous, half_lengths)) >= separation for previous in selected):
            selected.append(candidate)
        if len(selected) == count:
            break
    return selected


def fit_diamond_phase(points, template, *, width=2.5, seed=23):
    """Estimate translation only, holding the exact shape and all three periods fixed.

    The fit loss is a bounded Gaussian of F/|grad F|, a near-surface distance
    proxy. Evaluate physical distances separately with ``diamond_distances``.
    """
    from scipy.optimize import minimize

    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 20 or not np.isfinite(points).all():
        raise ValueError("phase fitting needs at least 20 finite Cartesian points")
    if not np.isfinite(width) or width <= 0:
        raise ValueError("width must be a finite positive distance")
    lengths = template.lengths
    candidates = _phase_candidates(points, template, width)
    rng = np.random.default_rng(seed)
    training = points[rng.choice(len(points), min(len(points), 2500), replace=False)]

    def objective(translation):
        f, g, h = exact_diamond_field(training - translation, lengths, order=2)
        normal_length = np.maximum(np.linalg.norm(g, axis=1), 1e-8)
        d = f / normal_length
        weights = np.exp(-0.5 * (d / width)**2)
        hg = np.einsum("nij,nj->ni", h, g)
        grad_d = g / normal_length[:, None] - f[:, None] * hg / normal_length[:, None]**3
        jacobian = -np.mean((weights * d / width**2)[:, None] * grad_d, axis=0)
        return float(1 - weights.mean()), jacobian

    solutions = [minimize(objective, x, jac=True, method="L-BFGS-B", options={"maxiter": 60}) for x in candidates]
    optimum = min(solutions, key=lambda result: result.fun)
    return (optimum.x / lengths) % 0.5, {"proxy_loss": float(optimum.fun), "coarse_start_count": len(solutions)}


def diamond_distances(points, template, phase):
    """Estimate periodic closest-point distances on the exact zero set.

    Refine several nearby exact surface samples by solving the constrained
    closest-point equations. Keep the nearest certified level-set sample as
    a fallback. This local multistart estimate is not a global proof.
    """
    lengths = template.lengths
    translation = np.asarray(phase) * lengths
    query = (np.asarray(points) - translation) % lengths
    distances, selected = template.tree.query(query, k=3)
    candidate = template.vertices[selected]
    unwrapped = query[:, None, :] + _minimum_image(candidate - query[:, None, :], lengths)
    best = unwrapped[:, 0].copy()
    best_distance = distances[:, 0].copy()
    targets = np.repeat(query, 3, axis=0)
    fitted = unwrapped.reshape(-1, 3).copy()
    _, gradient = exact_diamond_field(fitted, lengths, order=1)
    multiplier = np.sum((targets - fitted) * gradient, axis=1) / np.sum(gradient**2, axis=1)
    for _ in range(7):
        f, g, h = exact_diamond_field(fitted, lengths, order=2)
        jac = np.zeros((len(fitted), 4, 4))
        jac[:, :3, :3] = np.eye(3) + multiplier[:, None, None] * h
        jac[:, :3, 3] = jac[:, 3, :3] = g
        residual = np.column_stack((fitted - targets + multiplier[:, None] * g, f))
        active = (np.linalg.norm(residual, axis=1) > 1e-9) & (np.abs(np.linalg.det(jac)) > 1e-12)
        if not active.any():
            break
        step = np.linalg.solve(jac[active], residual[active, :, None])[:, :, 0]
        factor = np.minimum(1, 2 / np.maximum(np.linalg.norm(step[:, :3], axis=1), 2))
        fitted[active] -= step[:, :3] * factor[:, None]
        multiplier[active] -= step[:, 3] * factor
    f, g = exact_diamond_field(fitted, lengths, order=1)
    stationary = np.linalg.norm(fitted - targets + multiplier[:, None] * g, axis=1) < 1e-5
    exact = (np.abs(f) < 1e-9) & stationary
    delta = _minimum_image(fitted - targets, lengths).reshape(len(query), 3, 3)
    refined = np.where(exact.reshape(len(query), 3), np.linalg.norm(delta, axis=2), np.inf)
    chosen = np.argmin(refined, axis=1)
    best_refined = refined[np.arange(len(query)), chosen]
    improved = best_refined < best_distance
    best_distance[improved] = best_refined[improved]
    best[improved] = query[improved] + delta[np.arange(len(query)), chosen][improved]
    return best_distance, (best + translation) % lengths, exact.reshape(len(query), 3).any(axis=1)
