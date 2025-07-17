import numpy as np
import random
import Crystal_Operation as Lt

fcc = Lt.FCC()
Lattice = Lt.BCC()

Size = 8.30
Range = [Size, Size, Size]
a = 3.61
M = 63.546
Cfg_name = 'BCC'

orientation  = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]],          # O 
                [[-1, 2, 2], [2, -1, 2], [2, 2, -1]],       # A
                [[-1, -2, -2], [-2, -1, 2], [-2, 2, -1]],   # B
                [[-1, -2, 2], [-2, -1, -2], [2, -2, -1]],   # C
                [[-1, 2, -2], [2, -1, -2], [-2, -2, -1]],   # D
                # [[-1, 8, 4], [8, -1, 4], [4, 4, -7]],
                # [[-7, -4, 4], [-4, -1, -8], [4, -8, -1]],
                # [[1, 8, 4], [8, 1, -4], [-4, 4, -7]],
                ]

# 1. Bulit crystal with orientation
Cube = {}
for i in range(len(orientation)):
    Cube[i] = Lt.Grow_Crystal(Range, fcc.point, orientation[i])
    # Lt.Cfg(Cube[i], fcc.axis, Range, M, a, f'{i}')

# 2. Computer Voronoi Center point
face_centers, _ = Lt.compute_voronoi_face_centers(Lattice.point, Lattice.axis*Size)
# 3. Slice crystal to Voronoi shape
Inv = np.linalg.inv(Lattice.axis)
Voronoi_crytsal = {}
for i in range(len(orientation)):
    real_points = Cube[i] @ fcc.axis
    real_points = Lt.slice_center_points(real_points, [0, 0, 0], face_centers)
    Voronoi_crytsal[i] = real_points @ Inv
    # Lt.Cfg(Voronoi_crytsal[i], fcc.axis,  Range, M, a, f'slice{i}')
Lt.Cfg(Voronoi_crytsal[0], Lattice.axis,  Range, M, a, f'slice')
# 4. Move crystal to position
move_position = Lt.Duplicate(Lattice, 2, 2, 2)
# print(move_position.point)
combine_crystal = {}
# O
combine_crystal[0] = Lt.Move(Voronoi_crytsal[0], move_position.point[1], Range)
combine_crystal[1] = Lt.Move(Voronoi_crytsal[0], move_position.point[15], Range)
# A
# combine_crystal[2] = Lt.Move(Voronoi_crytsal[1], move_position.point[0], Range)
combine_crystal[3] = Lt.Move(Voronoi_crytsal[1], move_position.point[14], Range)
# B
combine_crystal[4] = Lt.Move(Voronoi_crytsal[2], move_position.point[8], Range)
# combine_crystal[5] = Lt.Move(Voronoi_crytsal[2], move_position.point[6], Range)
# C
# combine_crystal[6] = Lt.Move(Voronoi_crytsal[3], move_position.point[10], Range)
combine_crystal[7] = Lt.Move(Voronoi_crytsal[3], move_position.point[4], Range)
# D
combine_crystal[8] = Lt.Move(Voronoi_crytsal[4], move_position.point[2], Range)
# combine_crystal[9] = Lt.Move(Voronoi_crytsal[4], move_position.point[12], Range)

j = 9
exist_position = [1, 15, 14, 8, 4, 2]
for i in range(len(move_position.point)):
    if i in exist_position:
        continue
    else:
        j += 1
        orientation = Lt.rotation_matrix(random.randint(0, 180), random.randint(0, 180), random.randint(0, 180))
        # print(orientation)
        Cube = Lt.Grow_Crystal(Range, fcc.point, orientation)
        real_points = Cube @ fcc.axis
        real_points = Lt.slice_center_points(real_points, [0, 0, 0], face_centers)
        Voronoi_crytsal = real_points @ Inv
        combine_crystal[j] = Lt.Move(Voronoi_crytsal, move_position.point[i], Range)
# # E
# combine_crystal[10] = Lt.Move(Voronoi_crytsal[5], move_position.point[5], Range)
# combine_crystal[11] = Lt.Move(Voronoi_crytsal[5], move_position.point[11], Range)
# # F
# combine_crystal[12] = Lt.Move(Voronoi_crytsal[6], move_position.point[3], Range)
# combine_crystal[13] = Lt.Move(Voronoi_crytsal[6], move_position.point[13], Range)
# # I
# combine_crystal[14] = Lt.Move(Voronoi_crytsal[7], move_position.point[7], Range)
# combine_crystal[15] = Lt.Move(Voronoi_crytsal[7], move_position.point[9], Range)
# 5. Combine the all crystal to fill the simulation BOX
all = Lt.Merge_Group(combine_crystal)
# 6. Output Cfg
Lt.Cfg(all, Lattice.axis, np.dot(Range, 2), M, a, Cfg_name)