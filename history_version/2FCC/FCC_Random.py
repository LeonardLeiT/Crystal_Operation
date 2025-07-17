import numpy as np
import random
import Crystal_Operation as Lt

fcc = Lt.FCC()
Lattice = Lt.FCC()

Size = 27
Range = [Size, Size, Size]
a = 3.61
M = 63.546
Cfg_name = 'FCC_Cell'

# 1. Computer Voronoi Center point
face_centers, _ = Lt.compute_voronoi_face_centers(Lattice.point, Lattice.axis*Size)
# 2. Slice crystal to Voronoi shape
Inv = np.linalg.inv(Lattice.axis)
# 3. Move crystal to position
move_position = Lattice.point
combine_crystal={}
for i in range(len(move_position)):
    orientation = Lt.rotation_matrix(random.randint(0, 180), random.randint(0, 180), random.randint(0, 180))
    # print(orientation)
    Cube = Lt.Grow_Crystal(Range, fcc.point, orientation)
    real_points = Cube @ fcc.axis
    real_points = Lt.slice_center_points(real_points, [0, 0, 0], face_centers)
    Voronoi_crytsal = real_points @ Inv
    combine_crystal[i] = Lt.Move(Voronoi_crytsal, move_position[i], Range)

# 4. Combine the all crystal to fill the simulation BOX
all = Lt.Merge_Group(combine_crystal)
# 5. Output Cfg
Lt.Cfg(all, Lattice.axis, Range, M, a, Cfg_name)