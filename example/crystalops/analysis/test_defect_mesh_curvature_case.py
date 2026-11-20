"""Curvature sign, physical scale, and periodic seam checks."""

import numpy as np

from crystalops.analysis import DefectMeshResult, compute_defect_mesh_curvature


def _periodic_cylinder(radius):
    # An extruded regular polygon has area proportional to its radius.
    # Vertex symmetry makes its projected discrete H exactly 1/(2R).
    rings = 6
    sides = 24
    cell = np.diag([8 * radius, 8 * radius, 6 * radius])
    theta = np.arange(sides) * (2 * np.pi / sides)
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    vertices = np.array([(a, b, ring * radius) for ring in range(rings) for a, b in zip(x, y)])

    faces = []
    for ring in range(rings):
        for side in range(sides):
            next_ring = (ring + 1) % rings
            next_side = (side + 1) % sides
            a = ring * sides + side
            b = ring * sides + next_side
            c = next_ring * sides + next_side
            d = next_ring * sides + side
            faces.extend([(a, b, c), (a, c, d)])
    faces = np.asarray(faces)
    # Reposition a subset into another periodic image without changing geometry.
    vertices[::3, 2] += cell[2, 2]
    return DefectMeshResult(
        vertices=vertices, faces=faces, clipped_vertices=np.empty((0, 3)),
        clipped_faces=np.empty((0, 3), dtype=int), clipped_face_sources=np.empty(0, dtype=int),
        cell=cell, origin=np.array([-4 * radius, -4 * radius, 0.0]), pbc=(True, True, True),
        parameters={"length_unit": "angstrom"}, diagnostics={},
    )


def test_cylinder_curvature_through_periodic_seam():
    radius = 2.5
    mesh = _periodic_cylinder(radius)
    result = compute_defect_mesh_curvature(mesh)
    assert np.all(result.vertex_valid)
    np.testing.assert_allclose(result.vertex_mean_curvature, 1 / (2 * radius), atol=1e-13)
    np.testing.assert_allclose(result.summary["valid_area_fraction"], 1.0)
    np.testing.assert_allclose(result.summary["rms_H"], 1 / (2 * radius), atol=1e-13)


def test_open_edge_invalidates_estimate_without_zero_fill():
    mesh = _periodic_cylinder(2.5)
    removed = mesh.faces[-1].copy()
    mesh.faces = mesh.faces[:-1]
    result = compute_defect_mesh_curvature(mesh)
    assert not result.vertex_valid[removed].any()
    assert np.isnan(result.vertex_mean_curvature[removed]).all()
    assert 0 < result.summary["valid_area_fraction"] < 1


def test_quadratic_fit_uses_physical_radius_and_periodic_images():
    original = compute_defect_mesh_curvature(_periodic_cylinder(2.5), fit_radius=3.8)
    enlarged = compute_defect_mesh_curvature(_periodic_cylinder(5), fit_radius=7.6)
    assert original.vertex_valid.all() and enlarged.vertex_valid.all()
    assert np.min(original.vertex_mean_curvature) > 0
    np.testing.assert_allclose(enlarged.vertex_mean_curvature * 2, original.vertex_mean_curvature, atol=1e-12)
