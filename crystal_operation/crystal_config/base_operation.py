import numpy as np
from .lattice import Lattice, Duplicate

def slice_center_point(alloy, center, point):
    structure = np.copy(alloy)
    normal_vector = point - center
    
    # Calculate which points are on the desired side of the plane
    # Using vectorized operations for better performance
    dot_products = np.sum(normal_vector * (structure - point), axis=1)
    mask = dot_products <= 0
    return structure[mask]

def slice_center_points(alloy, center, points):
    slice_points = alloy.copy()
    for face_center in points:
        slice_points = slice_center_point(slice_points, center, face_center)
    return slice_points

def Move_lattice(lattice, vector):
    shift = np.array(vector) * lattice.size
    lattice.point = lattice.point + shift
    return lattice

def Move(points, vector, range):
    shift = np.array(vector) * np.array(range)
    Move = points + shift
    return Move

def Orhtogonality(axis):
    lenth = len(axis)
    for i in range(lenth-1):
        for j in range(i + 1, lenth):
            if np.dot(axis[i], axis[j]) >= 1e-14:
                return False
    return True

def Unitization(Vector):
    Vec = np.copy(Vector)
    square = np.square(Vec)
    return Vec / np.sqrt(np.sum(square))

def Grow_Crystal(Range, crystallographic,
                 axis=[[1, 0, 0], [0, 1, 0], [0, 0, 1]]):

    if not Orhtogonality(axis):
        print('axis is not orthogonal, please change!')
        return None

    # Normalize axis vectors
    axis = np.array([Unitization(a) for a in axis])
    
    # Transform crystallographic basis atoms
    basis_atoms = crystallographic @ axis.T
    basis_atoms = np.dot(crystallographic, axis)
    # Estimate replication range
    R = int(np.ceil(np.linalg.norm(Range)))
    
    # Generate lattice translations
    shifts = np.array([[i, j, k] for i in range(R) for j in range(R) for k in range(R)])
    displacements = shifts @ axis  # (N, 3)

    # Generate full structure
    structure = (basis_atoms[:, None, :] + displacements[None, :, :]).reshape(-1, 3)

    # Center the structure around the desired Center point (default [0,0,0])
    centroid = structure.mean(axis=0)
    structure -= centroid  # move to center at [0, 0, 0]
    return structure

#  晶体合并
def Merge(Crys1, Crys2):
    Crys3 = []
    for s in Crys1:
        Crys3.append(s)
    for s in Crys2:
        Crys3.append(s)
    return Crys3

# 对于多个晶体组合并
def Merge_Group(Group):
    return np.vstack(list(Group.values()))

# 晶体三维旋转矩阵
def rotation_matrix(x_angle, y_angle, z_angle):
    x_angle = np.radians(x_angle)  # Rotation about x-axis in degrees
    y_angle = np.radians(y_angle)  # Rotation about y-axis in degrees
    z_angle = np.radians(z_angle)  # Rotation about z-axis in degrees
    # Rotation matrix around x-axis
    R_x = np.array([
        [1, 0, 0],
        [0, np.cos(x_angle), -np.sin(x_angle)],
        [0, np.sin(x_angle), np.cos(x_angle)]
    ])

    # Rotation matrix around y-axis
    R_y = np.array([
        [np.cos(y_angle), 0, np.sin(y_angle)],
        [0, 1, 0],
        [-np.sin(y_angle), 0, np.cos(y_angle)]
    ])

    # Rotation matrix around z-axis
    R_z = np.array([
        [np.cos(z_angle), -np.sin(z_angle), 0],
        [np.sin(z_angle), np.cos(z_angle), 0],
        [0, 0, 1]
    ])

    # Combined rotation matrix: R = Rz * Ry * Rx
    R = R_z @ R_y @ R_x
    return R

    