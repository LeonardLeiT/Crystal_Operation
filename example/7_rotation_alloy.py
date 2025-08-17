import numpy as np
import random

from crystalops.config.base_operation import *
from crystalops.config.lattice import *
from crystalops.config.voronoi import *
from crystalops.config.write import *
from crystalops.users.lattice_shape import *

orientations = [
    # Y
    [0, 0, 0],
    [0 , 180 ,0],
    # Z
    [90, 0, 0],
    [-90, 0, 0],
    # X
    [0, 0, 90],
    [-180, 0, 270]
]
base_lattice = {}
base_lattice[0] = FCC_X_1()
base_lattice[1] = FCC_X_2()
shape_lattice = BCC()
size = 14
# 1. Build crystals with different orientations
range_box = [size, size, size]
cube = {}
for i in range(len(base_lattice)):
    # R = rotation_matrix(ori[0], ori[1], ori[2])
    lattice_cube = Grow_Crystal(range_box, base_lattice[i])
    cube[i] = lattice_cube.point @ lattice_cube.axis

# 2. Compute Voronoi face centers for the shape lattice
face_slices = []
for i in range(len(shape_lattice.point)):
    face_centers, _ = compute_voronoi_face_centers(
        shape_lattice.point, 
        shape_lattice.axis * size, 
        center_index=i
    )
    face_slices.append(face_centers)

norms = np.linalg.norm(face_centers, axis=1, keepdims=True)
unit_vectors = face_centers / norms
print(unit_vectors)

# 3. Slice the base crystals into Voronoi-shaped domains
# inv_shape_axis = np.linalg.inv(base_lattice.axis @ R)
voronoi_crystals = {}

for i in range(len(base_lattice)):
    real_points = cube[i]
    methods = np.zeros(len(face_centers), dtype=int)
    
    # # Y
    # if i == 0 or i == 1: 
    #     target = np.array([0, 1, 0])
    #     matches = np.all(unit_vectors == target, axis=1)
    #     methods[matches] = 0
    # # Z
    # if i == 2 or i == 3:   
    #     target = np.array([0, 0, 1])
    #     matches = np.all(unit_vectors == target, axis=1)
    #     methods[matches] = 1
    # # X
    # if i == 4 or i == 5:   
    #     target = np.array([1, 0, 0])
    #     matches = np.all(unit_vectors == target, axis=1)
    #     methods[matches] = 1
    
    real_points = slice_center_points(real_points, [0, 0, 0], face_centers, methods=methods)
    
    write_lammps_data(
        real_points = real_points,
        axis = shape_lattice.axis,
        Size = [size, size, size],
        point_types = np.zeros(len(real_points), dtype=int),
        lattice_constant = 3.61,
        Name = f'rotation_{i}',
        element = ['Cu']
    )
    voronoi_crystals[i] = real_points
        
# 4. Combine the Voronoi-shaped crysystal
move_positions = Duplicate(shape_lattice, 2, 2, 2)
combined_crystals = {}
print(move_positions.point)
# Y
# combined_crystals[0] = Move(voronoi_crystals[0], move_positions.point[8], range_box)
# combined_crystals[1] = Move(voronoi_crystals[1], move_positions.point[12], range_box)
# # Z
# combined_crystals[2] = Move(voronoi_crystals[2], move_positions.point[4], range_box)
# combined_crystals[3] = Move(voronoi_crystals[3], move_positions.point[6], range_box)
# X
combined_crystals[0] = Move(voronoi_crystals[0], move_positions.point[2], range_box)
combined_crystals[1] = Move(voronoi_crystals[1], move_positions.point[10], range_box)

# exist_positions = [8, 12, 4, 6, 2, 10]
# next_key = 6
# for i in range(len(move_positions.point)):
#     if i in exist_positions:
#         continue
    
#     face_slice_idx = 1 if i % 2 else 0
    
#     orientation = rotation_matrix(random.randint(0, 180), random.randint(0, 180), random.randint(0, 180))
#     lattice_cube = Grow_Crystal(range_box, base_lattice, orientation)
#     real_points = lattice_cube.point @ lattice_cube.axis
#     real_points = slice_center_points(real_points, [0, 0, 0], face_slices[face_slice_idx])
#     combined_crystals[next_key] = Move(real_points, move_positions.point[i], range_box)
#     next_key += 1

all_crystals = Merge_Group(combined_crystals)
print(f'Output kelvin shape, the atoms num is {len(all_crystals)}')
# 6. Write the structure to a LAMMPS data file
write_lammps_data(
    real_points = all_crystals,
    axis = shape_lattice.axis,
    Size = [size * 2, size * 2, size * 2],
    point_types = np.zeros(len(all_crystals), dtype=int),
    lattice_constant = 3.61,
    Name = f'kelvin',
    element = ['Cu']
)