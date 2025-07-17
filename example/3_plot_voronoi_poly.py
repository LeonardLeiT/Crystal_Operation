import os
import sys
from pathlib import Path

# 添加父目录到 sys.path，以便能够导入 crystal_config 包
sys.path.append(str(Path(__file__).resolve().parent.parent))

from crystal_operation.crystal_config import compute_voronoi_face_centers, plot_voronoi_polyhedron
from crystal_operation.crystal_config.lattice import *


shape_lattice = HCP()
vertices, face = compute_voronoi_face_centers(            
            shape_lattice.point, 
            shape_lattice.axis, 
            center_index=1,
            center_surface=False)

print(face)
plot_voronoi_polyhedron(vertices, face, save_path='HCP')