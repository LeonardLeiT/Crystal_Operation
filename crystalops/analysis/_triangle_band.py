"""Signed distances in a narrow band around orthorhombic periodic triangles."""

from __future__ import annotations

from itertools import product

import numpy as np


def _row_dot(first, second):
    return np.einsum("ij,ij->i", first, second)


def _closest_triangle_offsets(query, triangle):
    """Return point minus closest point for corresponding point/triangle pairs."""
    a, b, c = triangle.transpose(1, 0, 2)
    ab = b - a
    ac = c - a
    face_normal = np.cross(ab, ac)
    normal_sq = _row_dot(face_normal, face_normal)
    ap = query - a
    height = np.divide(
        _row_dot(ap, face_normal), normal_sq, out=np.zeros(len(query)), where=normal_sq > 0,
    )
    plane_offset = height[:, None] * face_normal
    planar_ap = ap - plane_offset
    coef_b = _row_dot(np.cross(planar_ap, ac), face_normal)
    coef_c = _row_dot(np.cross(ab, planar_ap), face_normal)
    within_face = (normal_sq > 0) & (coef_b >= 0) & (coef_c >= 0) & (coef_b + coef_c <= normal_sq)

    best = plane_offset
    best_sq = np.where(within_face, _row_dot(best, best), np.inf)
    for start, end in ((a, b), (b, c), (c, a)):
        edge = end - start
        from_start = query - start
        edge_sq = _row_dot(edge, edge)
        t = np.divide(
            _row_dot(from_start, edge), edge_sq,
            out=np.zeros(len(query)), where=edge_sq > 0,
        )
        offset = from_start - np.clip(t, 0, 1)[:, None] * edge
        distance_sq = _row_dot(offset, offset)
        better = distance_sq < best_sq
        best[better] = offset[better]
        best_sq[better] = distance_sq[better]
    return best, best_sq


def nearest_triangle_band(points, triangles, normals, lengths, max_distance):
    """Return signed nearest-triangle distances, or ``inf`` outside the band.

    The orthorhombic cell has periods ``lengths`` in x, y, z. Each triangle
    must use three vertices in a single unwrapped image, but different
    triangles and points can be in arbitrary images. Use consistently
    oriented, nonzero ``normals`` to define which side is positive. At a
    distance tie the sign is taken from the first candidate visited. A
    displacement tangent to a triangle at its edge has positive sign.
    """
    points = np.asarray(points, dtype=float)
    triangles = np.asarray(triangles, dtype=float)
    normals = np.asarray(normals, dtype=float)
    lengths = np.asarray(lengths, dtype=float)
    if points.ndim != 2 or points.shape[1:] != (3,) or not np.isfinite(points).all():
        raise ValueError("points must be a finite array with shape (N, 3)")
    if triangles.ndim != 3 or triangles.shape[1:] != (3, 3) or not np.isfinite(triangles).all():
        raise ValueError("triangles must be a finite array with shape (M, 3, 3)")
    if normals.shape != (len(triangles), 3) or not np.isfinite(normals).all() or np.any(~normals.any(axis=1)):
        raise ValueError("normals must contain one finite nonzero vector per triangle")
    if lengths.shape != (3,) or not np.isfinite(lengths).all() or np.any(lengths <= 0):
        raise ValueError("lengths must contain three finite positive cell lengths")
    if not np.isfinite(max_distance) or max_distance < 0:
        raise ValueError("max_distance must be finite and nonnegative")

    result = np.full(len(points), np.inf)
    if not len(points) or not len(triangles):
        return result

    from scipy.spatial import cKDTree

    # Any triangle within the band has its centroid within this radius.
    centroid = triangles.mean(axis=1)
    local_triangles = triangles - centroid[:, None, :]
    triangle_radius = np.linalg.norm(local_triangles, axis=2).max(axis=1)
    radius = max_distance + triangle_radius.max()
    centroid %= lengths
    np.minimum(centroid, np.nextafter(lengths, 0), out=centroid)
    points = np.minimum(points % lengths, np.nextafter(lengths, 0))
    tree = cKDTree(centroid, boxsize=lengths)
    chunk_size = 512
    band_sq = max_distance**2
    # Usually only the nearest centroid image fits. For wider bands or
    # larger triangles, another image can contain the closest point.
    image_counts = np.floor(radius / lengths + 0.5).astype(int)
    images = [lengths * image for image in product(*(range(-n, n + 1) for n in image_counts))]
    for start in range(0, len(points), chunk_size):
        chunk = points[start:start + chunk_size]
        nearby = tree.query_ball_point(chunk, np.nextafter(radius, np.inf), return_sorted=True)
        counts = np.fromiter(map(len, nearby), dtype=np.intp, count=len(chunk))
        if not counts.any():
            continue
        point_ids = np.repeat(np.arange(len(chunk)), counts)
        triangle_ids = np.concatenate(nearby).astype(np.intp, copy=False)
        offsets = chunk[point_ids] - centroid[triangle_ids]
        offsets -= lengths * np.rint(offsets / lengths)
        best_sq = np.full(len(chunk), np.inf)
        best_signed = np.full(len(chunk), np.inf)
        reach = np.nextafter(max_distance + triangle_radius[triangle_ids], np.inf)

        for image in images:
            query = offsets + image
            inside = _row_dot(query, query) <= reach**2
            if not inside.any():
                continue
            p = point_ids[inside]
            t = triangle_ids[inside]
            displacement, distance_sq = _closest_triangle_offsets(query[inside], local_triangles[t])
            nearest_sq = np.full(len(chunk), np.inf)
            np.minimum.at(nearest_sq, p, distance_sq)
            update = (distance_sq == nearest_sq[p]) & (distance_sq < best_sq[p]) & (distance_sq <= band_sq)
            # Keep the first triangle per point among exact distance ties.
            positions = np.flatnonzero(update)
            if not len(positions):
                continue
            _, first = np.unique(p[positions], return_index=True)
            positions = positions[first]
            q = p[positions]
            d = distance_sq[positions]
            sign_dot = _row_dot(displacement[positions], normals[t[positions]])
            best_sq[q] = d
            best_signed[q] = np.where(sign_dot < 0, -np.sqrt(d), np.sqrt(d))
        result[start:start + len(chunk)] = best_signed
    return result
