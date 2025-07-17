import numpy as np
import Crystal_Operation as Lt

fcc = Lt.FCC()
Lattice = Lt.SC()

Size = 10
Range = [Size, Size, Size]
a = 3.61
M = 63.546
Cfg_name = 'SC'

orientation  = [[[1, 0, 0], [0, 1, 0], [0, 0, 1]],          # O 
                [[-1, 2, 2], [2, -1, 2], [2, 2, -1]],       # A
                [[-1, -2, -2], [-2, -1, 2], [-2, 2, -1]],   # B
                [[-1, -2, 2], [-2, -1, -2], [2, -2, -1]],   # C
                [[-1, 2, -2], [2, -1, -2], [-2, -2, -1]],   # D
                ]

# 1. Bulit crystal with orientation
Cube = {}
for i in range(len(orientation)):
    Cube[i] = Lt.Grow_Crystal(Range, fcc.point, orientation[i])
    Lt.Cfg(Cube[i], fcc.axis, Range, M, a, f'{i}')
# 2. Computer Voronoi Center point
face_centers, _ = Lt.compute_voronoi_face_centers(Lattice.point, Lattice.axis*Size)
# 3. Slice crystal to Voronoi shape
Inv = np.linalg.inv(Lattice.axis)
Voronoi_crytsal = {}
for i in range(len(orientation)):
    real_points = Cube[i] @ fcc.axis
    real_points = Lt.slice_center_points(real_points, [0, 0, 0], face_centers)
    Voronoi_crytsal[i] = real_points @ Inv
# 4. Move crystal to position
move_position = Lt.Duplicate(Lattice, 2, 2, 2)
print(move_position.point)
combine_crystal = {}
combine_crystal[0] = Lt.Move(Voronoi_crytsal[0], move_position.point[0], Range)
combine_crystal[1] = Lt.Move(Voronoi_crytsal[3], move_position.point[1], Range)
combine_crystal[2] = Lt.Move(Voronoi_crytsal[2], move_position.point[2], Range)
combine_crystal[3] = Lt.Move(Voronoi_crytsal[1], move_position.point[3], Range)
combine_crystal[4] = Lt.Move(Voronoi_crytsal[2], move_position.point[4], Range)
combine_crystal[5] = Lt.Move(Voronoi_crytsal[1], move_position.point[5], Range)
combine_crystal[6] = Lt.Move(Voronoi_crytsal[3], move_position.point[6], Range)
combine_crystal[7] = Lt.Move(Voronoi_crytsal[0], move_position.point[7], Range)
combine_crystal[6] = Lt.Move(Voronoi_crytsal[3], move_position.point[6], Range)
combine_crystal[7] = Lt.Move(Voronoi_crytsal[0], move_position.point[7], Range)

# 5. Combine the all crystal to fill the simulation BOX
all = Lt.Merge_Group(combine_crystal)
# 6. Output Cfg
Lt.Cfg(all, Lattice.axis, np.dot(Range, 2), M, a, Cfg_name)