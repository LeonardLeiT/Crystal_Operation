"""Signed mean curvature on the native periodic triangles of a DXA mesh."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np
from numpy.typing import NDArray

from ._defect_mesh_quadric import FIT_RCOND, MIN_FIT_POINTS, fit_local_quadric_curvature

if TYPE_CHECKING:
    from .defect_mesh import DefectMeshResult


@dataclass(slots=True)
class DefectMeshCurvature:
    """Curvature and area arrays in the ordering of the native DXA mesh.

    Vertex areas are one third of each incident usable triangle's area.
    Invalid vertex curvature is NaN. A face is valid only when all three
    vertices are valid. Summary statistics integrate over those valid faces.
    """

    vertex_mean_curvature: NDArray[np.float64]
    vertex_area: NDArray[np.float64]
    vertex_valid: NDArray[np.bool_]
    face_area: NDArray[np.float64]
    face_valid: NDArray[np.bool_]
    summary: dict[str, Any]
    parameters: dict[str, Any]


def _edge_vectors(mesh: DefectMeshResult) -> NDArray[np.float64]:
    triangles = mesh.vertices[mesh.faces]
    edges = np.roll(triangles, -1, axis=1) - triangles
    if not any(mesh.pbc) or not len(edges):
        return edges

    # Componentwise image selection is valid when the periodic lattice
    # vectors are orthogonal to every other cell vector, even in a rotated box.
    cell = mesh.cell
    lengths = np.linalg.norm(cell, axis=1)
    directions = cell / lengths[:, None]
    gram = directions @ directions.T - np.eye(3)
    periodic = np.asarray(mesh.pbc)
    if np.any(np.abs(gram[periodic]) > 1e-12):
        raise ValueError("DXA curvature requires periodic cell vectors orthogonal to the other cell vectors")

    fractional_periodic = (edges @ directions[periodic].T) / lengths[periodic]
    image = np.rint(fractional_periodic)
    wrapped = fractional_periodic - image
    ambiguous = np.any(np.abs(np.abs(wrapped) - 0.5) <= 1e-12, axis=(1, 2))
    if ambiguous.any():
        first = int(np.flatnonzero(ambiguous)[0])
        raise ValueError(f"native face {first} has an ambiguous half-cell periodic edge")
    nonclosing = np.any(image.sum(axis=1) != 0, axis=1)
    if nonclosing.any():
        first = int(np.flatnonzero(nonclosing)[0])
        raise ValueError(f"minimum-image edges of native face {first} do not form a closed triangle")
    return edges - image @ cell[periodic]


def _fan_counts(faces: np.ndarray, edge_pairs: np.ndarray, vertex_count: int) -> np.ndarray:
    """Count connected fans of incident corners via consistently paired edges."""
    corners = np.arange(faces.size).reshape(faces.shape)
    next_corner = np.roll(corners, -1, axis=1).ravel()
    parent = list(range(faces.size))

    def root(corner):
        while parent[corner] != corner:
            parent[corner] = parent[parent[corner]]
            corner = parent[corner]
        return corner

    for start, end in zip(edge_pairs.ravel(), next_corner[edge_pairs[:, ::-1]].ravel()):
        first, second = root(int(start)), root(int(end))
        if first != second:
            parent[first] = second
    is_root = np.asarray(parent) == corners.ravel()
    return np.bincount(faces.ravel()[is_root], minlength=vertex_count)


def _topology_support(faces: np.ndarray, vertex_count: int):
    """Find vertices with incomplete, non-manifold, or inconsistently wound fans."""
    invalid = np.zeros(vertex_count, dtype=bool)
    counts = dict.fromkeys((
        "open_edge_count", "nonmanifold_edge_count", "inconsistent_winding_edge_count",
        "disconnected_fan_vertex_count",
    ), 0)
    if not len(faces):
        return invalid, counts

    start = faces.ravel()
    end = np.roll(faces, -1, axis=1).ravel()
    pairs = np.column_stack((np.minimum(start, end), np.maximum(start, end)))
    order = np.lexsort((pairs[:, 1], pairs[:, 0]))
    ordered_pairs = pairs[order]
    groups = np.flatnonzero(np.r_[True, np.any(ordered_pairs[1:] != ordered_pairs[:-1], axis=1)])
    multiplicity = np.diff(np.r_[groups, len(order)])
    open_edges = multiplicity == 1
    nonmanifold_edges = multiplicity > 2
    internal_pairs = order[groups[multiplicity == 2, None] + np.arange(2)]
    consistent = start[internal_pairs[:, 0]] == end[internal_pairs[:, 1]]
    bad_edges = np.concatenate((
        ordered_pairs[groups[open_edges | nonmanifold_edges]],
        pairs[internal_pairs[~consistent, 0]],
    ))
    invalid[bad_edges.ravel()] = True
    fans = _fan_counts(faces, internal_pairs[consistent], vertex_count)
    disconnected = (fans > 1) & ~invalid
    invalid |= disconnected
    counts.update({
        "open_edge_count": int(open_edges.sum()),
        "nonmanifold_edge_count": int(nonmanifold_edges.sum()),
        "inconsistent_winding_edge_count": int((~consistent).sum()),
        "disconnected_fan_vertex_count": int(disconnected.sum()),
    })
    return invalid, counts


def _validate_fit_radius(mesh: DefectMeshResult, radius: float) -> float:
    try:
        if isinstance(radius, (bool, np.bool_)):
            raise TypeError
        radius = float(radius)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("fit_radius must be a finite positive number") from exc
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("fit_radius must be a finite positive number")
    periodic = np.asarray(mesh.pbc, dtype=bool)
    if periodic.any() and radius >= 0.5 * np.linalg.norm(mesh.cell[periodic], axis=1).min():
        raise ValueError("fit_radius must be less than half the shortest periodic cell length")
    return radius


def compute_defect_mesh_curvature(
    mesh: DefectMeshResult, *, fit_radius: float | None = None,
) -> DefectMeshCurvature:
    """Estimate mean curvature H=(k1+k2)/2 on native DXA triangles.

    By default the cotangent Laplacian of position uses barycentric vertex area and
    minimum-image edge vectors. The signed scalar is -0.5 times its projection
    onto the unit vertex normal, computed from the source face winding with
    area weights. An outward-oriented sphere has positive H. The absolute
    mean and RMS use this signed scalar, rather than the Laplacian magnitude.

    Set ``fit_radius`` in the mesh length unit to estimate curvature of a
    local quadratic surface over that graph-geodesic radius. The fit uses
    periodic offsets and triangle area and distance weights. The radius must
    be positive and less than half the shortest periodic cell length.

    Degenerate or duplicate triangles and open, non-manifold, or inconsistently
    wound support invalidate the affected vertices. The area statistics use
    only faces with three valid vertices and apply vertex quadrature separately
    to H, |H|, and H squared. The total area is the area of all native faces.
    """
    mesh.validate()
    if fit_radius is not None:
        fit_radius = _validate_fit_radius(mesh, fit_radius)
    faces = mesh.faces
    vertex_count = len(mesh.vertices)
    edges = _edge_vectors(mesh)
    cross = np.cross(edges[:, 0], -edges[:, 2])
    double_area = np.linalg.norm(cross, axis=1)
    longest_edge_sq = np.einsum("fij,fij->fi", edges, edges).max(axis=1)
    if not np.isfinite(double_area).all() or not np.isfinite(longest_edge_sq).all():
        raise ValueError("native triangles have non-finite edge lengths or areas")
    tolerance = 64 * np.finfo(float).eps
    degenerate = double_area <= tolerance * longest_edge_sq
    _, inverse, multiplicity = np.unique(
        np.sort(faces, axis=1), axis=0, return_inverse=True, return_counts=True,
    )
    duplicate = multiplicity[inverse.ravel()] > 1
    usable = ~(degenerate | duplicate)
    good_faces = faces[usable]
    good_edges = edges[usable]
    invalid, topology = _topology_support(good_faces, vertex_count)
    invalid[faces[~usable].ravel()] = True

    vertex_area = np.zeros(vertex_count, dtype=float)
    vertex_normal = np.zeros((vertex_count, 3), dtype=float)
    vertices = good_faces.ravel()
    np.add.at(vertex_area, vertices, np.repeat(double_area[usable] / 6, 3))
    np.add.at(vertex_normal, vertices, np.repeat(cross[usable], 3, axis=0))

    normal_length = np.linalg.norm(vertex_normal, axis=1)
    positive_area = vertex_area > 0
    unresolved_normal = positive_area & (normal_length <= 6 * tolerance * vertex_area)
    vertex_valid = ~invalid & positive_area & ~unresolved_normal
    fit_diagnostics = {}
    if fit_radius is None:
        laplace_sum = np.zeros_like(vertex_normal)
        # Along each directed edge use the cotangent at the opposite corner.
        inward_edges = -np.roll(good_edges, 1, axis=1)
        cotangent = np.einsum("fij,fij->fi", good_edges, inward_edges) / double_area[usable, None]
        edge_terms = (good_edges * np.roll(cotangent, 1, axis=1)[..., None]).reshape(-1, 3)
        np.add.at(laplace_sum, vertices, edge_terms)
        np.add.at(laplace_sum, np.roll(good_faces, -1, axis=1).ravel(), -edge_terms)
        vertex_curvature = np.full(vertex_count, np.nan)
        normals = vertex_normal[vertex_valid] / normal_length[vertex_valid, None]
        laplacian = laplace_sum[vertex_valid] / (2 * vertex_area[vertex_valid, None])
        vertex_curvature[vertex_valid] = -0.5 * np.einsum("ij,ij->i", normals, laplacian)
    else:
        support = usable & vertex_valid[faces].all(axis=1)
        vertex_curvature, fit_diagnostics = fit_local_quadric_curvature(
            mesh, faces[support], edges[support], double_area[support], cross[support], fit_radius,
        )
    vertex_valid &= np.isfinite(vertex_curvature)
    vertex_curvature[~vertex_valid] = np.nan
    face_valid = usable & vertex_valid[faces].all(axis=1)

    face_area = double_area / 2
    total_area = float(face_area.sum())
    valid_area = float(face_area[face_valid].sum())
    mean_signed = mean_abs = rms = None
    if valid_area > 0:
        face_curvature = vertex_curvature[faces[face_valid]]
        weights = face_area[face_valid, None] / (3 * valid_area)
        mean_signed = float(np.sum(weights * face_curvature))
        mean_abs = float(np.sum(weights * np.abs(face_curvature)))
        peak = float(np.max(np.abs(face_curvature)))
        rms = peak * float(np.sqrt(np.sum(weights * (face_curvature / peak) ** 2))) if peak else 0.0

    length_unit = mesh.parameters.get("length_unit", "input_length_unit")
    summary = {
        "mean_signed_H": mean_signed,
        "mean_abs_H": mean_abs,
        "rms_H": rms,
        "total_area": total_area,
        "valid_area": valid_area,
        "valid_area_fraction": valid_area / total_area if total_area > 0 else 0.0,
        "native_vertex_count": vertex_count,
        "native_face_count": int(len(faces)),
        "valid_vertex_count": int(vertex_valid.sum()),
        "valid_face_count": int(face_valid.sum()),
        "degenerate_face_count": int(degenerate.sum()),
        "duplicate_face_count": int(duplicate.sum()),
        "unresolved_normal_vertex_count": int(unresolved_normal.sum()),
        **topology,
        **fit_diagnostics,
    }
    parameters = {
        "method": "cotangent_laplacian",
        "vertex_area": "barycentric",
        "curvature_definition": "H=(k1+k2)/2; -0.5 times the normal component of the position Laplacian",
        "normal_orientation": "source face winding with area-weighted vertex normals",
        "area_integration": "vertex quadrature on faces with three valid vertices",
        "length_unit": length_unit,
        "curvature_unit": f"1/{length_unit}",
        "area_unit": f"{length_unit}^2",
        "surface_smoothing": mesh.parameters.get("surface_smoothing"),
    }
    if fit_radius is not None:
        parameters.update({
            "method": "local_quadric_fit",
            "fit_radius": fit_radius,
            "fit_radius_unit": length_unit,
            "fit_distance": "shortest path along native valid mesh edges",
            "fit_weights": "barycentric support area times (1 - (distance / radius)^2)^2",
            "fit_minimum_points": MIN_FIT_POINTS,
            "fit_lstsq_rcond": FIT_RCOND,
            "fit_residual_unit": length_unit,
            "curvature_definition": "H=(k1+k2)/2 on the fitted local height graph at the central vertex",
            "normal_orientation": "weighted PCA normal oriented by the winding of valid support faces",
        })
    return DefectMeshCurvature(
        vertex_curvature, vertex_area, vertex_valid, face_area, face_valid, summary, parameters,
    )
