import numpy as np
import Crystal_Operation as Lt

fcc = Lt.FCC()
Lattice = Lt.BCO(1, 1.1, 1.2)

Size = 13
Range = [Size, Size, Size]
a = 3.61
M = 63.546
Cfg_name = 'BCO'

orientation  = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]],          # O 
                [[-1, 2, 2], [2, -1, 2], [2, 2, -1]],       # A
                [[-1, -2, -2], [-2, -1, 2], [-2, 2, -1]],   # B
                [[-1, -2, 2], [-2, -1, -2], [2, -2, -1]],   # C
                [[-1, 2, -2], [2, -1, -2], [-2, -2, -1]],   # D
                [[-1, 8, 4], [8, -1, 4], [4, 4, -7]],
                [[-7, -4, 4], [-4, -1, -8], [4, -8, -1]],
                [[1, 8, 4], [8, 1, -4], [-4, 4, -7]],
                ]

# 1. Bulit crystal with orientation
Cube = {}
for i in range(len(orientation)):
    Cube[i] = Lt.Grow_Crystal(np.dot(Range, 2), fcc.point, orientation[i])
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
Lt.Cfg(Voronoi_crytsal[0], Lattice.axis,  Range, M, a, f'slice')
# 4. Move crystal to position
move_position = Lattice
# print(move_position.point)
combine_crystal = {}
# O
combine_crystal[0] = Lt.Move(Voronoi_crytsal[0], move_position.point[0], Range)
combine_crystal[1] = Lt.Move(Voronoi_crytsal[1], move_position.point[1], Range)
# 5. Combine the all crystal to fill the simulation BOX
all = Lt.Merge_Group(combine_crystal)
# 6. Output Cfg
Lt.Cfg(all, Lattice.axis, Range, M, a, Cfg_name)