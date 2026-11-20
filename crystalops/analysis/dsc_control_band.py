"""Select atoms within physical distances of the DSC builder's curved cuts."""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

from crystalops.core import CrystalStructure
from crystalops.modeling.builders.schwarz import (
    _DSC_TWIN_ANGLES,
    _DSC_TWIN_NORMALS,
    _DSC_TWIN_SIGNS,
    _d_schwarz_field,
    _dsc_coherent_geometry,
    _dsc_transition_planes,
    _dsc_translations,
)
from crystalops.modeling.operations import select_atoms

from ._triangle_band import nearest_triangle_band


def _clip_triangles(
    triangles: NDArray[np.float64], normal: NDArray[np.float64], offset: float
) -> NDArray[np.float64]:
    """Clip to normal·x >= offset without adding a surface on the cut."""
    side = triangles @ normal - offset
    counts = (side >= 0.0).sum(axis=1)
    pieces = [triangles[counts == 3]]
    for count in (1, 2):
        vertices = triangles[counts == count]
        values = side[counts == count]
        if not len(vertices):
            continue
        # Rotate each triangle so the unique inside or outside vertex is first.
        pivot = np.argmax(values, axis=1) if count == 1 else np.argmin(values, axis=1)
        order = (pivot[:, None] + np.arange(3)) % 3
        vertices = vertices[np.arange(len(vertices))[:, None], order]
        values = values[np.arange(len(values))[:, None], order]
        fraction = values[:, :1] / (values[:, :1] - values[:, 1:])
        intersections = vertices[:, :1] + fraction[..., None] * (vertices[:, 1:] - vertices[:, :1])
        if count == 1:
            pieces.append(np.concatenate((vertices[:, :1], intersections), axis=1))
        else:
            first, second = intersections.transpose(1, 0, 2)
            _, b, c = vertices.transpose(1, 0, 2)
            pieces.extend((np.stack((first, b, c), axis=1), np.stack((first, c, second), axis=1)))
    return np.concatenate(pieces)


def _control_mesh(lengths, lattice_constant, spacing):
    from skimage.measure import marching_cubes

    size = lengths / lattice_constant
    cubic = np.allclose(size, size[0], atol=1e-10, rtol=0.0)
    reference = np.full(3, size[0])
    aspect = size / reference
    translations = _dsc_translations(size)
    if not cubic:
        geometry, _ = _dsc_coherent_geometry(size, reference, translations)

    counts = np.maximum(2, np.ceil(lengths / spacing).astype(int))
    step = lengths / counts
    axes = [np.linspace(0.0, length, n + 1) for length, n in zip(lengths, counts)]
    grid = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)
    reference_grid = grid.reshape(-1, 3) / (lattice_constant * aspect)
    patches = []
    normals = []
    for index, (angle, sign) in enumerate(zip(_DSC_TWIN_ANGLES, _DSC_TWIN_SIGNS)):
        domain = index + 3
        shift = lattice_constant * translations[domain]
        field = sign * _d_schwarz_field(reference_grid, reference, angle=angle)
        vertices, faces, _, _ = marching_cubes(
            field.reshape(grid.shape[:-1]), level=0.0, spacing=step, allow_degenerate=False
        )
        triangles = vertices[faces]
        if cubic:
            for point, normal, sense in _dsc_transition_planes(size, domain):
                normal = np.asarray(normal, dtype=float)
                if sense[0] == "l":
                    normal = -normal
                triangles = _clip_triangles(triangles, normal, lattice_constant * (normal @ point))
        else:
            normal = _DSC_TWIN_NORMALS[index]
            low = lattice_constant * geometry.low[index] - normal @ shift
            high = lattice_constant * geometry.high[index] - normal @ shift
            triangles = _clip_triangles(triangles, normal, low)
            triangles = _clip_triangles(triangles, -normal, -high)

        edges = triangles[:, 1:] - triangles[:, :1]
        cross = np.cross(edges[:, 0], edges[:, 1])
        area = np.linalg.norm(cross, axis=1)
        valid = area > 1e-12
        triangles = triangles[valid]
        cross = cross[valid] / area[valid, None]
        center = triangles.mean(axis=1)
        delta = 1e-4 * cross
        plus = _d_schwarz_field((center + delta) / (lattice_constant * aspect), reference, angle=angle)
        minus = _d_schwarz_field((center - delta) / (lattice_constant * aspect), reference, angle=angle)
        cross[sign * (plus - minus) < 0] *= -1
        patches.append(triangles + shift)
        normals.append(cross)
    return np.concatenate(patches), np.concatenate(normals)


def select_dsc_control_band(
    structure: CrystalStructure,
    *,
    behind: float = 2.0,
    ahead: float = 2.0,
    lattice_constant: float = 3.61,
    mesh_spacing: float = 1.0,
) -> CrystalStructure:
    """Mark atoms near the DSC construction's four curved zero-set patches.

    Length arguments are in Å. The positive side of each patch faces into its
    DSC transition domain. Distance is measured to the nearest periodic
    triangle, including patch edges. Planar twin faces and ownership changes
    from boundary refinement are outside the scope of this reference surface.
    For a relaxed structure it is a control geometry band in the current cell.
    """
    distances = np.asarray([behind, ahead], dtype=float)
    if not np.isfinite(distances).all() or np.any(distances < 0) or not distances.any():
        raise ValueError("behind and ahead must be nonnegative finite distances, at least one positive")
    resolution = np.asarray([lattice_constant, mesh_spacing], dtype=float)
    if not np.isfinite(resolution).all() or np.any(resolution <= 0):
        raise ValueError("lattice_constant and mesh_spacing must be finite and positive")
    lengths = np.diag(structure.cell)
    if (
        not all(structure.pbc)
        or not np.array_equal(structure.cell, np.diag(lengths))
        or not np.isfinite(lengths).all()
        or np.any(lengths <= 0)
    ):
        raise ValueError("the DSC control band requires an orthorhombic cell periodic on all axes")

    triangles, normals = _control_mesh(lengths, lattice_constant, mesh_spacing)

    def region(points):
        signed = nearest_triangle_band(
            points - structure.origin, triangles, normals, lengths, float(distances.max())
        )
        return (-behind <= signed) & (signed <= ahead)

    return select_atoms(structure, region)
