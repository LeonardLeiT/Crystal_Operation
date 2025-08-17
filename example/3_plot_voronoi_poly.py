import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from crystalops.config.voronoi import compute_voronoi_face_centers, plot_voronoi_polyhedron
from crystalops.config.lattice import *


shape_lattice = RHL(60) 
vertices, face = compute_voronoi_face_centers(            
            shape_lattice.point, 
            shape_lattice.axis, 
            center_index=0,
            center_surface=False)

print(face)
plot_voronoi_polyhedron(vertices, face, save_path='FCO')