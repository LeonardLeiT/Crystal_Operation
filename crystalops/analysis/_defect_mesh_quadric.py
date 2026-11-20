"""Local curvature fits on a periodic defect mesh's valid surface graph."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from .defect_mesh import DefectMeshResult


MIN_FIT_POINTS = 12
FIT_RCOND = 1e-6


def _fit_quadric(points, factors, reference_normal):
    """Fit a height graph in a weighted PCA frame using normalized lengths."""
    weights = factors**2
    centered = points - np.average(points, axis=0, weights=weights)
    weighted_points = centered * factors[:, None]
    _, basis = np.linalg.eigh(weighted_points.T @ weighted_points)
    normal_alignment = np.dot(basis[:, 0], reference_normal)
    if normal_alignment == 0:
        return None
    if normal_alignment < 0:
        basis[:, 0] *= -1

    z, x, y = (points @ basis).T
    terms = np.column_stack((np.ones(len(x)), x, y, x**2 / 2, x * y, y**2 / 2))
    coeff, _, rank, _ = np.linalg.lstsq(terms * factors[:, None], z * factors, rcond=FIT_RCOND)
    if rank < len(coeff):
        return None

    _, zx, zy, zxx, zxy, zyy = coeff
    metric = 1 + zx**2 + zy**2
    # Evaluate div(n)/2 at the central vertex's tangent-plane coordinates.
    curvature = -((1 + zy**2) * zxx - 2 * zx * zy * zxy + (1 + zx**2) * zyy) / (2 * metric**1.5)
    residual = np.sqrt(np.average((terms @ coeff - z)**2, weights=weights))
    return curvature, residual


def fit_local_quadric_curvature(
    mesh: DefectMeshResult,
    faces: np.ndarray,
    edges: np.ndarray,
    double_area: np.ndarray,
    face_normals: np.ndarray,
    radius: float,
):
    """Fit only through faces whose three vertex fans passed topology checks."""
    vertex_count = len(mesh.vertices)
    curvature = np.full(vertex_count, np.nan)
    residual = np.full(vertex_count, np.nan)
    diagnostics = {
        "fit_eligible_vertex_count": 0,
        "fit_insufficient_support_vertex_count": 0,
        "fit_failed_vertex_count": 0,
        "fit_median_residual": None,
    }
    if not len(faces):
        return curvature, diagnostics

    try:
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import dijkstra
    except ImportError as exc:
        raise ImportError(
            "local quadric DXA curvature requires SciPy; install crystalops[analysis]"
        ) from exc

    start = faces.ravel()
    end = np.roll(faces, -1, axis=1).ravel()
    pairs, first_edge = np.unique(
        np.column_stack((np.minimum(start, end), np.maximum(start, end))),
        axis=0, return_index=True,
    )
    edge_length = np.linalg.norm(edges.reshape(-1, 3)[first_edge], axis=1)
    graph = coo_matrix((edge_length, (pairs[:, 0], pairs[:, 1])), shape=(vertex_count, vertex_count)).tocsr()

    support_area = np.zeros(vertex_count)
    support_normal = np.zeros_like(mesh.vertices)
    np.add.at(support_area, start, np.repeat(double_area / 6, 3))
    np.add.at(support_normal, start, np.repeat(face_normals, 3, axis=0))
    eligible = np.flatnonzero(support_area > 0)
    diagnostics["fit_eligible_vertex_count"] = int(len(eligible))

    periodic = np.asarray(mesh.pbc, dtype=bool)
    fractional = (mesh.vertices - mesh.origin) @ np.linalg.inv(mesh.cell)
    # Limit the size of the temporary dense distance matrix.
    block_size = 64
    for start_index in range(0, len(eligible), block_size):
        block = eligible[start_index:start_index + block_size]
        distances = dijkstra(graph, directed=False, indices=block, limit=radius)
        for vertex, distance in zip(block, distances):
            indices = np.flatnonzero(distance < radius)
            if len(indices) < MIN_FIT_POINTS:
                diagnostics["fit_insufficient_support_vertex_count"] += 1
                continue
            delta = fractional[indices] - fractional[vertex]
            delta[:, periodic] -= np.rint(delta[:, periodic])
            points = (delta @ mesh.cell) / radius
            kernel = 1 - (distance[indices] / radius)**2
            factors = np.sqrt(support_area[indices]) * kernel
            fit = _fit_quadric(points, factors, support_normal[vertex])
            if fit is None:
                continue
            curvature[vertex] = fit[0] / radius
            residual[vertex] = fit[1] * radius

    succeeded = np.isfinite(curvature) & np.isfinite(residual)
    curvature[~succeeded] = np.nan
    diagnostics["fit_failed_vertex_count"] = int(len(eligible) - succeeded.sum())
    if succeeded.any():
        diagnostics["fit_median_residual"] = float(np.median(residual[succeeded]))
    return curvature, diagnostics
