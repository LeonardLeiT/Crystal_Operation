import numpy as np
import random
import Crystal_Operation as Lt

fcc = Lt.FCC()
Lattice = Lt.HCP()

Size = 28
Range = [Size, Size, Size]
a = 3.61
M = 63.546
Cfg_name = 'HCP_Cell'

face_slice = []
# 1. Computer Voronoi Center point
for i in range(len(Lattice.point)):
    face_centers, _ = Lt.compute_voronoi_face_centers(Lattice.point, Lattice.axis*Size, center_index=i)
    face_slice.append(face_centers)
# 2. Slice crystal to Voronoi shape
Inv = np.linalg.inv(Lattice.axis)
# 3. Move crystal to position
move_position = Lattice.point
combine_crystal={}
for i in range(len(move_position)):
    orientation = Lt.rotation_matrix(random.randint(0, 180), random.randint(0, 180), random.randint(0, 180))
    # print(orientation)
    Cube = Lt.Grow_Crystal(np.dot(Range, 3), fcc.point, orientation)
    real_points = Cube @ fcc.axis
    real_points = Lt.slice_center_points(real_points, [0, 0, 0], face_slice[i])
    Voronoi_crytsal = real_points @ Inv
    combine_crystal[i] = Lt.Move(Voronoi_crytsal, move_position[i], Range)
    Lt.Cfg(Voronoi_crytsal, Lattice.axis,  Range, M, a, f'slice{i}')
# 4. Combine the all crystal to fill the simulation BOX
all = Lt.Merge_Group(combine_crystal)
# 5. Output Cfg
Lt.Cfg(all, Lattice.axis, Range, M, a, Cfg_name)