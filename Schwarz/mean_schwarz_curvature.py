"""Compute theoretical S_D = 0 curvature in the boxes of the listed data files."""

from pathlib import Path

import numpy as np
from ovito.io import import_file
from skimage.measure import marching_cubes

folder = Path(__file__).resolve().parent
data_files = ["1-1.6.data", "HCP-0.data"]
grid_intervals = 128

# S_D is half the sum of these four sine waves in fractional coordinates.
modes = np.pi * np.array([[1, 1, -1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]])
axis = np.linspace(0, 1, grid_intervals + 1)
u, v, w = np.meshgrid(axis, axis, axis, indexing="ij", sparse=True)
field = np.zeros((axis.size,) * 3)
for ku, kv, kw in modes:
    field += 0.5 * np.sin(ku * u + kv * v + kw * w)

vertices, faces, _, _ = marching_cubes(field, level=0, allow_degenerate=False)
fractional_triangles = vertices[faces].astype(float) / grid_intervals
phase = fractional_triangles.mean(axis=1) @ modes.T

for filename in data_files:
    data = import_file(str(folder / filename)).compute()
    cell = data.cell[:, :3]
    cartesian_modes = np.linalg.solve(cell.T, modes.T).T
    gradient = 0.5 * np.cos(phase) @ cartesian_modes
    hessian = -0.5 * np.einsum("nm,mi,mj->nij", np.sin(phase), cartesian_modes, cartesian_modes)

    # H = (trace(Hessian) - normal.T @ Hessian @ normal) / (2 * |gradient|).
    gradient_norm = np.linalg.norm(gradient, axis=1)
    normal = gradient / gradient_norm[:, None]
    normal_second = np.einsum("ni,nij,nj->n", normal, hessian, normal)
    curvature = (np.trace(hessian, axis1=1, axis2=2) - normal_second) / (2 * gradient_norm)

    p0, p1, p2 = (fractional_triangles @ cell.T).transpose(1, 0, 2)
    areas = 0.5 * np.linalg.norm(np.cross(p1 - p0, p2 - p0), axis=1)
    mean_signed = np.average(curvature, weights=areas)
    mean_absolute = np.average(np.abs(curvature), weights=areas)
    normalized = np.cbrt(abs(np.linalg.det(cell))) * mean_absolute

    lengths = np.linalg.norm(cell, axis=0)
    a, b, c = (cell / lengths).T
    cosines = [b @ c, c @ a, a @ b]
    angles = np.degrees(np.arccos(np.clip(cosines, -1, 1)))
    print(f"\n{filename}")
    print("Edge ratio (a:b:c):", np.round(lengths / lengths[0], 6))
    print("Angles (alpha, beta, gamma) in degrees:", np.round(angles, 6))
    print("Orthogonal:", np.allclose(cosines, 0, atol=1e-10, rtol=0))
    print(f"Area-weighted H: {mean_signed:.8g} / length")
    print(f"Area-weighted |H|: {mean_absolute:.8g} / length")
    print(f"V^(1/3) * mean |H|: {normalized:.8g}")
