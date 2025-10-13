import os
import sys
import numpy as np
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from crystalops.config.voronoi import compute_voronoi_face_centers, plot_voronoi_polyhedron, plot_voronoi_polycombine
from crystalops.config.lattice import *


shape_lattice = BCC() 
vertices, face = compute_voronoi_face_centers(            
            shape_lattice.point, 
            shape_lattice.axis, 
            center_index=0,
            center_surface=False)

face_center = np.array([[0.5,0.5,0.5],
                        [-0.5,-0.5,-0.5],
                        [-0.5,0.5,0.5],
                        [0.5,-0.5,-0.5],
                        [0.5,-0.5,0.5],
                        [-0.5,0.5,-0.5],
                        [0.5,0.5,-0.5],
                        [-0.5,-0.5,0.5]
])

# print(face_center*2)

# print(face)
# print(vertices+np.array([0.5,0.5,0.5]))
# plot_voronoi_polyhedron(vertices+np.array([0,0,0]), face, save_path='FCO')



plot_voronoi_polycombine(vertices, face, face_center, save_path='FCO')