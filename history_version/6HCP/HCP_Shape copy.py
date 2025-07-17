import numpy as np
import Crystal_Operation as Lt
import random

# np.sqrt(2)
siz = 6
fcc = Lt.FCC()
Lattice = Lt.HCP(np.sqrt(siz)/3)
Size = 16
Range = [Size, Size, Size]
a = 3.61
M = 63.546
Cfg_name = f'HCP'
rotation = Lt.rotation_matrix(0, 0, 0)
print(rotation)
orientation  = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]],          # O 
                [[-1, 2, 2], [2, -1, 2], [2, 2, -1]],       # A
                [[-1, -2, -2], [-2, -1, 2], [-2, 2, -1]],   # B
                [[-1, -2, 2], [-2, -1, -2], [2, -2, -1]],   # C
                [[-1, 2, -2], [2, -1, -2], [-2, -2, -1]],   # D
                ]
print(rotation@orientation[0])
for i in range(len(orientation)):
    orientation[i] = rotation@orientation[i]
# 1. Bulit crystal with orientation
Cube = {}
for i in range(len(orientation)):
    Cube[i] = Lt.Grow_Crystal(Range, fcc.point, orientation[i])
    # Lt.Cfg(Cube[i], fcc.axis, Range, M, a, f'{i}')
# 2. Computer Voronoi Center point
face_slice = []
for i in range(len(Lattice.point)):
    face_centers, _ = Lt.compute_voronoi_face_centers(Lattice.point, Lattice.axis*Size, center_index=i)
    face_slice.append(face_centers)
# 3. Slice crystal to Voronoi shape
Inv = np.linalg.inv(Lattice.axis)
Voronoi_crytsal_1 = {}
for i in range(len(orientation)):
    real_points = Cube[i] @ fcc.axis
    real_points = Lt.slice_center_points(real_points, [0, 0, 0], face_slice[0])
    Voronoi_crytsal_1[i] = real_points @ Inv
Voronoi_crytsal_2 = {}
for i in range(len(orientation)):
    real_points = Cube[i] @ fcc.axis
    real_points = Lt.slice_center_points(real_points, [0, 0, 0], face_slice[1])
    Voronoi_crytsal_2[i] = real_points @ Inv
Lt.Cfg(Voronoi_crytsal_1[0], Lattice.axis,  Range, M, a, f'slice0')
Lt.Cfg(Voronoi_crytsal_2[0], Lattice.axis,  Range, M, a, f'slice1')
# 4. Move crystal to position
move_position = Lt.Duplicate(Lattice, 1, 1, 2)
print(len(move_position.point))
print(move_position.point)
combine_crystal = {}
# O
combine_crystal[0] = Lt.Move(Voronoi_crytsal_1[0], move_position.point[0], Range)
# combine_crystal[1] = Lt.Move(Voronoi_crytsal_1[0], move_position.point[14], Range)
# A
combine_crystal[2] = Lt.Move(Voronoi_crytsal_2[1], move_position.point[1], Range)
# combine_crystal[3] = Lt.Move(Voronoi_crytsal_2[1], move_position.point[3], Range)
# B
# combine_crystal[4] = Lt.Move(Voronoi_crytsal_1[2], move_position.point[2], Range)
# combine_crystal[5] = Lt.Move(Voronoi_crytsal_2[4], move_position.point[15], Range)
# C
# combine_crystal[6] = Lt.Move(Voronoi_crytsal_2[3], move_position.point[3], Range)
# combine_crystal[7] = Lt.Move(Voronoi_crytsal_2[3], move_position.point[7], Range)
# D
# combine_crystal[8] = Lt.Move(Voronoi_crytsal_2[4], move_position.point[15], Range)
# combine_crystal[9] = Lt.Move(Voronoi_crytsal_1[4], move_position.point[12], Range)

# j = 9
# exist_position = [12, 14, 13, 3, 5, 1, 7, 15]
# # exist_position = [1, 3, 0, 14, 8, 2, 10, 12]
# for i in range(len(move_position.point)):
#     if i%2:
#         face_centers = face_slice[1]
#     else:
#         face_centers = face_slice[0]
#     if i in exist_position:
#         continue
#     else:
#         j += 1
#         orientation = Lt.rotation_matrix(random.randint(0, 180), random.randint(0, 180), random.randint(0, 180))
#         # print(orientation)
#         Cube = Lt.Grow_Crystal(Range, fcc.point, orientation)
#         real_points = Cube @ fcc.axis
#         real_points = Lt.slice_center_points(real_points, [0, 0, 0], face_centers)
#         Voronoi_crytsal = real_points @ Inv
#         combine_crystal[j] = Lt.Move(Voronoi_crytsal, move_position.point[i], Range)
# 5. Combine the all crystal to fill the simulation BOX
all = Lt.Merge_Group(combine_crystal)
# 6. Output Cfg
Lt.Cfg(all, Lattice.axis, [Size, Size, Size*2], M, a, Cfg_name)