"""Periodic grain-boundary reconstruction and scale-dependent curvature.

Grain labels define Gaussian occupancy fields. Curvature is evaluated from
their analytic derivatives, not from atomic steps or a curvature-flow fit.
An auxiliary grid locates surfaces and supplies area quadrature. Pairwise
meshes are clipped by local grain dominance; junction bands are excluded
from curvature statistics. They are not watertight multi-region meshes for
topological measurements. This implementation requires an orthogonal,
fully periodic simulation cell and complete current grain labels.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from numpy.typing import ArrayLike, NDArray

from crystalops.core import CrystalStructure
from .curvature import (
    gaussian_curvature_from_derivatives,
    mean_curvature_from_derivatives,
)


FloatArray = NDArray[np.float64]


@dataclass(slots=True)
class GrainBoundarySelection:
    grain_ids: NDArray[np.int64]
    boundary_mask: NDArray[np.bool_]
    junction_mask: NDArray[np.bool_]
    pair_atom_indices: dict[tuple[int, int], NDArray[np.int64]]
    neighbor_cutoff: float


@dataclass(slots=True)
class GrainBoundarySurface:
    grain_pair: tuple[int, int]
    vertices: FloatArray
    faces: NDArray[np.int64]
    face_centers: FloatArray
    face_area: FloatArray
    face_mean_curvature: FloatArray
    face_gaussian_curvature: FloatArray
    face_valid: NDArray[np.bool_]
    face_junction: NDArray[np.bool_]
    face_projection_residual: FloatArray
    face_gradient_strength: FloatArray
    face_projection_distance: FloatArray


@dataclass(slots=True)
class GrainBoundaryReconstruction:
    surfaces: list[GrainBoundarySurface]
    parameters: dict[str, Any]


@dataclass(slots=True)
class GrainBoundaryCurvatureResult:
    surfaces: list[GrainBoundarySurface]
    atom_mean_curvature: FloatArray
    atom_valid: NDArray[np.bool_]
    atom_grain_pairs: NDArray[np.int64]
    atom_boundary_mask: NDArray[np.bool_]
    atom_junction_mask: NDArray[np.bool_]
    atom_surface_distance: FloatArray
    summary: dict[str, Any]
    parameters: dict[str, Any]


def _box(structure: CrystalStructure) -> FloatArray:
    if not all(structure.pbc):
        raise NotImplementedError("grain-boundary analysis currently requires xyz periodicity")
    lengths = np.diag(structure.cell).copy()
    if np.any(lengths <= 0) or not np.allclose(structure.cell, np.diag(lengths), atol=1e-10):
        raise NotImplementedError("grain-boundary analysis currently requires an orthogonal cell")
    if not np.all(np.isfinite(structure.positions)):
        raise ValueError("atom positions must be finite")
    return lengths


def _labels(structure: CrystalStructure, grain_ids: ArrayLike | None) -> NDArray[np.int64]:
    values = structure.properties.get("grain_ids") if grain_ids is None else grain_ids
    if values is None:
        raise ValueError("grain_ids are required; load labels or run grain segmentation first")
    values = np.asarray(values)
    if values.shape != (structure.number_of_atoms,):
        raise ValueError("grain_ids must contain one label per input atom")
    if not np.issubdtype(values.dtype, np.number) or not np.all(np.isfinite(values)):
        raise ValueError("grain_ids must be finite integers")
    if np.any(values != np.floor(values)) or np.any(values < 0):
        raise ValueError("grain_ids must be nonnegative integers with no unassigned atoms")
    labels = values.astype(np.int64, copy=True)
    if len(np.unique(labels)) < 2:
        raise ValueError("at least two grains are required")
    return labels


def identify_grain_boundary_atoms(
    structure: CrystalStructure,
    *,
    neighbor_cutoff: float,
    grain_ids: ArrayLike | None = None,
) -> GrainBoundarySelection:
    """Identify unlike-grain neighbor contacts, retaining the complete input.

    A junction atom has contacts involving at least three distinct grains.
    This discrete mask is supplemented by a continuous junction mask during
    reconstruction. Element types are never interpreted as grain labels.
    """
    from scipy.spatial import cKDTree

    lengths = _box(structure)
    labels = _labels(structure, grain_ids)
    if not np.isfinite(neighbor_cutoff) or not 0 < neighbor_cutoff < 0.5 * min(lengths):
        raise ValueError("neighbor_cutoff must be positive and smaller than half the shortest box edge")
    positions = (structure.positions - structure.origin) % lengths
    tree = cKDTree(positions, boxsize=lengths)
    edges = tree.query_pairs(neighbor_cutoff, output_type="ndarray")
    edges = edges[labels[edges[:, 0]] != labels[edges[:, 1]]]
    if not len(edges):
        raise ValueError("no unlike-grain contacts found; check grain labels and neighbor_cutoff")
    pairs = np.sort(labels[edges], axis=1)
    pair_atoms = {}
    for pair in np.unique(pairs, axis=0):
        match = np.all(pairs == pair, axis=1)
        pair_atoms[tuple(int(x) for x in pair)] = np.unique(edges[match])
    boundary = np.zeros(len(labels), dtype=bool)
    boundary[np.unique(edges)] = True
    touched: dict[int, set[int]] = {}
    for a, b in edges:
        touched.setdefault(int(a), {int(labels[a])}).add(int(labels[b]))
        touched.setdefault(int(b), {int(labels[b])}).add(int(labels[a]))
    junction = np.zeros(len(labels), dtype=bool)
    for atom, neighbors in touched.items():
        junction[atom] = len(neighbors) >= 3
    return GrainBoundarySelection(labels, boundary, junction, pair_atoms, float(neighbor_cutoff))


class _PeriodicGaussianField:
    """Chunked direct Gaussian sums with an explicit finite support radius."""

    def __init__(self, structure, labels, sigma, kernel_cutoff=5.0, chunk_size=1024):
        from scipy.spatial import cKDTree

        self.lengths = _box(structure)
        self.origin = structure.origin.copy()
        self.sigma = float(sigma)
        self.kernel_cutoff = float(kernel_cutoff)
        self.chunk_size = int(chunk_size)
        if not np.isfinite(sigma) or sigma <= 0:
            raise ValueError("smoothing_length must be finite and positive")
        if not np.isfinite(kernel_cutoff) or kernel_cutoff < 3:
            raise ValueError("kernel_cutoff must be finite and at least 3 sigma")
        self.radius = self.sigma * self.kernel_cutoff
        if self.radius >= 0.5 * min(self.lengths):
            raise ValueError("kernel_cutoff * smoothing_length must be smaller than half the box edge")
        self.positions = (structure.positions - self.origin) % self.lengths
        self.grains, self.label_indices = np.unique(labels, return_inverse=True)
        self.grain_index = {int(g): i for i, g in enumerate(self.grains)}
        self.tree = cKDTree(self.positions, boxsize=self.lengths)

    def _chunks(self, points):
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        for start in range(0, len(points), self.chunk_size):
            p = (points[start:start + self.chunk_size] - self.origin) % self.lengths
            neighbors = self.tree.query_ball_point(p, self.radius)
            counts = np.fromiter((len(x) for x in neighbors), dtype=np.int64, count=len(p))
            rows = np.repeat(np.arange(len(p)), counts)
            if not len(rows):
                yield start, len(p), rows, np.empty(0, dtype=int), np.empty((0, 3)), np.empty(0)
                continue
            atoms = np.concatenate(neighbors).astype(np.int64, copy=False)
            displacement = p[rows] - self.positions[atoms]
            displacement -= self.lengths * np.rint(displacement / self.lengths)
            weights = np.exp(-np.einsum("ij,ij->i", displacement, displacement) / (2 * self.sigma**2))
            yield start, len(p), rows, atoms, displacement, weights

    def density(self, points):
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        result = np.zeros((len(points), len(self.grains)))
        for start, size, rows, atoms, _, weights in self._chunks(points):
            flat = rows * len(self.grains) + self.label_indices[atoms]
            result[start:start + size] = np.bincount(
                flat, weights=weights, minlength=size * len(self.grains)
            ).reshape(size, len(self.grains))
        return result

    def pair_derivatives(self, points, grain_pair, *, order=2):
        """F=q_j-q_i; its normal points from the lower ID to the higher ID."""
        points = np.asarray(points, dtype=float).reshape(-1, 3)
        values = np.zeros(len(points))
        gradients = np.zeros((len(points), 3))
        hessians = np.zeros((len(points), 3, 3)) if order == 2 else None
        i, j = (self.grain_index[g] for g in grain_pair)
        for start, size, rows, atoms, displacement, weights in self._chunks(points):
            labels = self.label_indices[atoms]
            signed = weights * ((labels == j).astype(float) - (labels == i))
            values[start:start + size] = np.bincount(rows, weights=signed, minlength=size)
            for axis in range(3):
                gradients[start:start + size, axis] = np.bincount(
                    rows, weights=-signed * displacement[:, axis] / self.sigma**2, minlength=size
                )
            if order == 2:
                for a in range(3):
                    for b in range(a, 3):
                        term = displacement[:, a] * displacement[:, b] / self.sigma**4
                        if a == b:
                            term -= 1 / self.sigma**2
                        value = np.bincount(rows, weights=signed * term, minlength=size)
                        hessians[start:start + size, a, b] = value
                        hessians[start:start + size, b, a] = value
        return values, gradients, hessians


def reconstruct_grain_boundary_surfaces(
    structure: CrystalStructure,
    boundary_atoms: GrainBoundarySelection,
    *,
    smoothing_length: float,
    sampling_spacing: float | None = None,
    junction_threshold: float = 0.1,
    kernel_cutoff: float = 5.0,
    max_grid_points: int = 2_000_000,
) -> GrainBoundaryReconstruction:
    """Reconstruct pairwise interfaces and evaluate analytic face curvatures.

    The junction score is the total occupancy of all grains other than the
    interface pair, divided by total occupancy. Curvature is invalid when
    this score exceeds ``junction_threshold``. The mesh is a centroid-clipped
    area approximation at junctions, not a conforming multi-region tessellation.
    """
    from skimage.measure import marching_cubes

    if not 0 <= junction_threshold < 1 / 3:
        raise ValueError("junction_threshold must lie in [0, 1/3)")
    labels = _labels(structure, boundary_atoms.grain_ids)
    field = _PeriodicGaussianField(structure, labels, smoothing_length, kernel_cutoff)
    requested = smoothing_length / 3 if sampling_spacing is None else float(sampling_spacing)
    if not np.isfinite(requested) or requested <= 0:
        raise ValueError("sampling_spacing must be finite and positive")
    shape = np.maximum(4, np.ceil(field.lengths / requested).astype(int))
    size = int(np.prod(shape))
    if size > max_grid_points or size * len(field.grains) * 4 > 256_000_000:
        raise ValueError("auxiliary sampling grid exceeds budget; increase sampling_spacing")
    spacing = field.lengths / shape
    axes = [np.arange(n) * step + origin for n, step, origin in zip(shape, spacing, field.origin)]
    points = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1).reshape(-1, 3)
    rho = field.density(points).astype(np.float32).reshape(*shape, len(field.grains))
    del points
    rho = np.pad(rho, ((0, 1), (0, 1), (0, 1), (0, 0)), mode="wrap")
    total_grid = np.sum(rho, axis=-1)
    surfaces = []
    for pair in sorted(boundary_atoms.pair_atom_indices):
        i, j = (field.grain_index[g] for g in pair)
        values = rho[..., j] - rho[..., i]
        # Permissive support mask; exact pair dominance is checked at quadrature points.
        support = (rho[..., i] + rho[..., j]) >= 0.25 * total_grid
        if not np.any(values[support] < 0) or not np.any(values[support] > 0):
            continue
        try:
            vertices, faces, _, _ = marching_cubes(
                values, level=0.0, spacing=tuple(spacing), mask=support,
                allow_degenerate=False,
            )
        except RuntimeError as error:
            if "No surface" in str(error):
                continue
            raise
        vertices = vertices.astype(float) + field.origin
        centers = vertices[faces].mean(axis=1)
        area = np.linalg.norm(np.cross(
            vertices[faces[:, 1]] - vertices[faces[:, 0]],
            vertices[faces[:, 2]] - vertices[faces[:, 0]],
        ), axis=1) / 2
        # Project quadrature points to the actual continuous field's zero surface.
        projected = centers.copy()
        for _ in range(4):
            f, g, _ = field.pair_derivatives(projected, pair, order=1)
            norm2 = np.einsum("ij,ij->i", g, g)
            displacement = f[:, None] * g / np.maximum(norm2[:, None], 1e-24)
            length = np.linalg.norm(displacement, axis=1)
            displacement *= np.minimum(1.0, min(spacing) / np.maximum(length, 1e-24))[:, None]
            projected -= displacement
        occupancy = field.density(projected)
        total = occupancy.sum(axis=1)
        pair_mean = (occupancy[:, i] + occupancy[:, j]) / 2
        others = [k for k in range(len(field.grains)) if k not in (i, j)]
        other_max = np.max(occupancy[:, others], axis=1) if others else np.zeros(len(centers))
        keep = (pair_mean >= other_max) & (total > 1e-10) & (area > 1e-12)
        # A face exactly on an upper periodic box plane belongs to the lower copy.
        for axis in range(3):
            upper = field.origin[axis] + field.lengths[axis]
            keep &= ~np.all(np.isclose(vertices[faces, axis], upper, atol=1e-7, rtol=0), axis=1)
        if not np.any(keep):
            continue
        faces, area, projected = faces[keep], area[keep], projected[keep]
        occupancy, total = occupancy[keep], total[keep]
        used, inverse = np.unique(faces, return_inverse=True)
        vertices, faces = vertices[used], inverse.reshape(-1, 3)
        f, g, hessian = field.pair_derivatives(projected, pair)
        mean = mean_curvature_from_derivatives(g, hessian)
        gaussian = gaussian_curvature_from_derivatives(g, hessian)
        pair_sum = occupancy[:, i] + occupancy[:, j]
        residual = np.abs(f) / np.maximum(pair_sum, 1e-24)
        projection_distance = np.linalg.norm(projected - centers[keep], axis=1)
        strength = np.linalg.norm(g, axis=1) * smoothing_length / np.maximum(total, 1e-24)
        junction = (1 - pair_sum / total) > junction_threshold
        valid = (~junction & (residual < 1e-4) & (strength > 1e-3)
                 & (projection_distance < 0.5 * min(spacing)) & np.isfinite(mean))
        mean[~valid] = np.nan
        gaussian[~valid] = np.nan
        # Orient each triangle consistently with the implicit normal.
        geometric = np.cross(vertices[faces[:, 1]] - vertices[faces[:, 0]],
                             vertices[faces[:, 2]] - vertices[faces[:, 0]])
        flip = np.einsum("ij,ij->i", geometric, g) < 0
        faces[flip] = faces[flip][:, [0, 2, 1]]
        surfaces.append(GrainBoundarySurface(
            pair, vertices, faces, projected, area, mean, gaussian, valid,
            junction, residual, strength, projection_distance,
        ))
    parameters = {
        "method": "periodic Gaussian grain occupancy; analytic implicit curvature",
        "smoothing_length": float(smoothing_length),
        "requested_sampling_spacing": requested,
        "actual_sampling_spacing": spacing.tolist(),
        "grid_cells": shape.tolist(),
        "kernel_cutoff_sigma": float(kernel_cutoff),
        "junction_other_occupancy_threshold": float(junction_threshold),
        "normal_convention": "lower grain ID to higher grain ID; sphere outward H=+1/R",
        "length_unit": structure.metadata.get("length_unit", "input_length_unit"),
        "curvature_unit": "1/" + structure.metadata.get("length_unit", "input_length_unit"),
        "mesh_junction_treatment": "pairwise centroid-dominance clipping; not watertight at junctions",
        "quadrature": "face centroid projected onto continuous F=0",
        "maximum_valid_projection_distance": float(0.5 * min(spacing)),
        "periodic_cell": structure.cell.tolist(),
        "periodic_origin": structure.origin.tolist(),
        "periodic_directions": list(structure.pbc),
    }
    return GrainBoundaryReconstruction(surfaces, parameters)


def _statistics(surfaces) -> dict[str, Any]:
    total_area = sum(float(s.face_area.sum()) for s in surfaces)
    areas = [s.face_area[s.face_valid] for s in surfaces]
    curvatures = [s.face_mean_curvature[s.face_valid] for s in surfaces]
    weights = np.concatenate(areas) if areas else np.empty(0)
    values = np.concatenate(curvatures) if curvatures else np.empty(0)
    valid_area = float(weights.sum())
    result = {
        "total_interface_area": total_area,
        "valid_interface_area": valid_area,
        "excluded_interface_area": total_area - valid_area,
        "valid_area_fraction": valid_area / total_area if total_area else 0.0,
        "triangle_count": sum(len(s.faces) for s in surfaces),
        "valid_triangle_count": int(len(values)),
        "junction_excluded_area": sum(float(s.face_area[s.face_junction].sum()) for s in surfaces),
        "numerically_excluded_area": sum(float(s.face_area[~s.face_valid & ~s.face_junction].sum()) for s in surfaces),
        "maximum_projection_distance": max((float(s.face_projection_distance.max()) for s in surfaces), default=0.0),
        "mean_signed_H": None,
        "mean_absolute_H": None,
        "rms_H": None,
        "median_absolute_H": None,
        "p95_absolute_H": None,
        "signed_cancellation_fraction": None,
        "exact_zero_H_valid_face_count": int(np.count_nonzero(values == 0)),
        "positive_H_area_fraction": None,
        "negative_H_area_fraction": None,
    }
    if valid_area:
        absolute = np.abs(values)
        order = np.argsort(absolute)
        cumulative = (np.cumsum(weights[order]) - 0.5 * weights[order]) / valid_area
        result.update(
            mean_signed_H=float(np.average(values, weights=weights)),
            mean_absolute_H=float(np.average(absolute, weights=weights)),
            rms_H=float(np.sqrt(np.average(values**2, weights=weights))),
            median_absolute_H=float(np.interp(0.5, cumulative, absolute[order])),
            p95_absolute_H=float(np.interp(0.95, cumulative, absolute[order])),
            signed_cancellation_fraction=(
                float(np.clip(1 - abs(np.sum(weights * values)) / np.sum(weights * absolute), 0, 1))
                if np.sum(weights * absolute) > 0 else None
            ),
            positive_H_area_fraction=float(weights[values > 0].sum() / valid_area),
            negative_H_area_fraction=float(weights[values < 0].sum() / valid_area),
        )
    return result


def analyze_grain_boundary_curvature(
    structure: CrystalStructure,
    *,
    neighbor_cutoff: float,
    smoothing_length: float,
    grain_ids: ArrayLike | None = None,
    sampling_spacing: float | None = None,
    junction_threshold: float = 0.1,
    kernel_cutoff: float = 5.0,
    grain_id_source: str | None = None,
) -> GrainBoundaryCurvatureResult:
    """Run boundary identification, reconstruction, area statistics, and mapping.

    Atomic values are nearest same-pair surface quadrature values, restricted
    to the detected boundary band. Invalid/junction values remain NaN.
    """
    from scipy.spatial import cKDTree

    selection = identify_grain_boundary_atoms(
        structure, neighbor_cutoff=neighbor_cutoff, grain_ids=grain_ids,
    )
    reconstruction = reconstruct_grain_boundary_surfaces(
        structure, selection, smoothing_length=smoothing_length,
        sampling_spacing=sampling_spacing, junction_threshold=junction_threshold,
        kernel_cutoff=kernel_cutoff,
    )
    surfaces = reconstruction.surfaces
    if not surfaces:
        raise ValueError("no shared grain-boundary surfaces reconstructed")
    lengths = _box(structure)
    atom_mean = np.full(structure.number_of_atoms, np.nan)
    atom_distance = np.full(structure.number_of_atoms, np.nan)
    atom_pairs = np.full((structure.number_of_atoms, 2), -1, dtype=np.int64)
    atom_valid = np.zeros(structure.number_of_atoms, dtype=bool)
    atom_junction = selection.junction_mask.copy()
    for surface in surfaces:
        candidates = selection.pair_atom_indices[surface.grain_pair]
        candidates = candidates[~selection.junction_mask[candidates]]
        if not len(candidates):
            continue
        tree = cKDTree((surface.face_centers - structure.origin) % lengths, boxsize=lengths)
        distances, indices = tree.query((structure.positions[candidates] - structure.origin) % lengths)
        atom_pairs[candidates] = surface.grain_pair
        atom_distance[candidates] = distances
        atom_junction[candidates] |= surface.face_junction[indices]
        good = surface.face_valid[indices] & (distances <= neighbor_cutoff)
        atom_mean[candidates[good]] = surface.face_mean_curvature[indices[good]]
        atom_valid[candidates[good]] = True
    summary = _statistics(surfaces)
    summary.update(
        atom_count=structure.number_of_atoms,
        grain_count=int(len(np.unique(selection.grain_ids))),
        boundary_atom_count=int(selection.boundary_mask.sum()),
        junction_atom_count=int(atom_junction.sum()),
        mapped_valid_atom_count=int(atom_valid.sum()),
        interface_pair_count=len(surfaces),
        per_pair=[dict(grain_pair=list(s.grain_pair), **_statistics([s])) for s in surfaces],
    )
    parameters = reconstruction.parameters
    parameters.update(
        neighbor_cutoff=float(neighbor_cutoff),
        grain_id_source=grain_id_source or structure.metadata.get("grain_id_source", "provided grain_ids"),
        source=structure.metadata.get("source"),
        atom_mapping="nearest same-pair face quadrature point; boundary band only",
        interpretation="geometric curvature of the labeled occupancy interface at the specified smoothing scale",
        signed_mean_interpretation=(
            "mean_signed_H is a network-wide signed area average, depends on grain-ID normal conventions, "
            "and is not a test that local mean curvature vanishes"
        ),
        reconstruction_validation_scope=(
            "grid refinement checks quadrature, not accuracy of the reconstructed physical interface; "
            "validate smoothing and labeling against a matched control before interpreting minimality or strain sensitivity"
        ),
        minimal_surface_assessment="not performed; no automatic threshold/classification",
        physical_interface_accuracy_status="not calibrated or validated for this input",
    )
    return GrainBoundaryCurvatureResult(
        surfaces, atom_mean, atom_valid, atom_pairs, selection.boundary_mask,
        atom_junction, atom_distance, summary, parameters,
    )
