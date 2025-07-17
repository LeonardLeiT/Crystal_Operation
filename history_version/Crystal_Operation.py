import numpy as np
import math
import random
import subprocess
import os
from scipy.spatial import Voronoi
##### 基础函数 ######
# 定义晶体结构
class Lattice:
    def __init__(self, point, axis, size=None):
        self.point = np.array(point)
        self.axis = np.array(axis)
        self.size = np.array(size) if size is not None else np.array([1, 1, 1])

# 立方
def SC():  # Simple Cubic
    return Lattice([[0, 0, 0]], np.eye(3))

def BCC():  # Body-Centered Cubic
    return Lattice([[0, 0, 0], [0.5, 0.5, 0.5]], np.eye(3))

def FCC():  # Face-Centered Cubic
    return Lattice([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]], np.eye(3))

# 六方 / 三方
def HCP(c_a_ratio=np.sqrt(8/3)):  # Hexagonal Close-Packed
    axis = [[1, 0, 0], 
            [-0.5, np.sqrt(3)/2, 0], 
            [0, 0, c_a_ratio]]
    return Lattice([[0, 0, 0], [2/3, 1/3, 0.5]], axis)

# 四方
def ST(z = 1.5):  # Simple Tetragonal
    axis = [[1, 0, 0], [0, 1, 0], [0, 0, z]]
    return Lattice([[0, 0, 0]], axis)

def BCT(z = 1.5):  # Body-Centered Tetragonal
    axis = [[1, 0, 0], [0, 1, 0], [0, 0, z]]
    return Lattice([[0, 0, 0], [0.5, 0.5, 0.5]], axis)

# 正交
def SO(x = 1, y = 2, z = 3):  # Simple Orthorhombic
    axis = [[x, 0, 0], [0, y, 0], [0, 0, z]]
    return Lattice([[0, 0, 0]], axis)

def BCO(x = 1, y = 2, z = 3):  # Base-Centered Orthorhombic
    axis = [[x, 0, 0], [0, y, 0], [0, 0, z]]
    return Lattice([[0, 0, 0], [0.5, 0.5, 0]], axis)

def FCO(x = 1, y = 2, z = 3):  # Face-Centered Orthorhombic
    axis = [[x, 0, 0], [0, y, 0], [0, 0, z]]
    return Lattice([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]], axis)

def ICO(ratio_yx = 2, ratio_zx = 3):  # Body-Centered Orthorhombic
    axis = [[1, 0, 0], [0, ratio_yx, 0], [0, 0, ratio_zx]]
    return Lattice([[0, 0, 0], [0.5, 0.5, 0.5]], axis)

# 单斜
def SM(x = 1, y = 2, z = 3, alpha_deg=100):  # Simple Monoclinic
    alpha = np.radians(alpha_deg)
    axis = [[x, 0, 0], [0, y, 0], [np.cos(alpha) * z, 0, np.sin(alpha) * z]]
    return Lattice([[0, 0, 0]], axis)

def BCM(x = 1, y = 2, z = 3, alpha_deg=100):  # Base-Centered Monoclinic
    alpha = np.radians(alpha_deg)
    axis = [[x, 0, 0], [0, y, 0], [np.cos(alpha) * z, 0, np.sin(alpha) * z]]
    return Lattice([[0, 0, 0], [0.5, 0.5, 0]], axis)

# 三斜
def TRI(length = [1, 2, 3], angle = [80, 85, 75]):  # Triclinic
    alpha = np.radians(angle[0])
    beta = np.radians(angle[1])
    gamma = np.radians(angle[2])
    a = length[0]
    b = length[1]
    c = length[2]
    ax = [a, 0, 0]
    ay = [b * np.cos(gamma), b * np.sin(gamma), 0]
    cz_x = c * np.cos(beta)
    cz_y = c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz_z = np.sqrt(c**2 - cz_x**2 - cz_y**2)
    az = [cz_x, cz_y, cz_z]
    return Lattice([[0, 0, 0]], [ax, ay, az])

# 菱方
def RHL(alpha_deg=75):  # Rhombohedral (Trigonal)
    alpha = np.radians(alpha_deg)
    a1 = [1, 0, 0]
    a2 = [np.cos(alpha), np.sin(alpha), 0]
    a3 = [
        np.cos(alpha),
        (np.cos(alpha) - np.cos(alpha)**2)/(np.sin(alpha)),
        np.sqrt(1 - 2 * np.cos(alpha)**2 + np.cos(alpha)**3)/(np.sin(alpha))
    ]
    return Lattice([[0, 0, 0]], [a1, a2, a3])
    
# 定义晶体结构复制
def Duplicate(lattice, nx=1, ny=1, nz=1):
    """
    Replicate a unit cell in the x, y, z directions.

    Parameters:
    - lattice: Lattice object with attributes .point (N,3), .axis (3,3), and .size (3,)
    - nx, ny, nz: Positive integers specifying the number of replications along each axis.

    Returns:
    - A new Lattice object with replicated atomic positions and updated size.
    """
    if not all(isinstance(n, int) and n > 0 for n in [nx, ny, nz]):
        raise ValueError("nx, ny, nz must be positive integers.")

    # Generate translation vectors in lattice units
    shifts = np.array([[i, j, k] for i in range(nx) for j in range(ny) for k in range(nz)])
    shifts = shifts * lattice.size  # convert to Cartesian translation in terms of unit cell

    # Add translations to each atom in the original lattice
    all_points = shifts[:, None, :] + lattice.point[None, :, :]
    all_points = all_points.reshape(-1, 3)

    # Update size
    new_size = lattice.size * np.array([nx, ny, nz])

    return Lattice(all_points, lattice.axis, new_size)

##### 写出为cfg文件 ######
def Cfg(points, axis, Size, M, lattice_constant=1, Name='structure', element='Cu'):
    """
    Write a CFG file from atomic positions and lattice vectors.

    Parameters:
    - Str: (N, 3) array-like, Cartesian coordinates of atoms.
    - axis: (3, 3) array-like, lattice vectors (as rows).
    - M: integer, number of atoms per type (repeats for multiple types if needed).
    - Name: string, output file name prefix (default 'structure').
    - Size: scaling factor for lattice size (currently unused but reserved).
    - element: string, chemical element symbol (default 'Cu').
    """

    # Scale lattice vectors
    real_points = np.array(points) @ np.array(axis) 
    axis_scaled = np.array(axis) * np.array(Size)
    # print(axis_scaled)
    # Convert to fractional coordinates (unitless)
    frac_coords = np.linalg.solve(axis_scaled.T, real_points.T).T
    axis_scaled = axis_scaled.astype('float64') * lattice_constant

    Natoms = len(frac_coords)
    filename = f"{Name}.cfg"

    with open(filename, 'w') as file:
        file.write(f"Number of particles = {Natoms}\n")
        file.write("# General structure with arbitrary cell vectors\n")
        file.write("A = 1.000000000 Angstrom (basic length-scale)\n")

        for i in range(3):
            for j in range(3):
                file.write(f"H0({i+1},{j+1}) =       {axis_scaled[i, j]:.8f}\n")

        file.write(".NO_VELOCITY.\n")
        file.write("entry_count = 3\n")
        file.write(f"  {M}\n")
        file.write(f"{element}\n")

        # Write atomic fractional coordinates
        for atom in frac_coords:
            file.write("{:.8f}    {:.8f}    {:.8f}\n".format(*atom))

def Cfg_lattice(lattice, M, lattice_constant=1, Name='structure', element='Cu'):
    Cfg(lattice.point, lattice.axis, lattice.size,M, lattice_constant, Name, element)

def compute_voronoi_face_centers(points, box_vectors = np.eye(3), center_index=0, layers=1, center_surface = 0):
    real_points = points @ box_vectors
    shifts = np.arange(-layers, layers+1)
    translations = np.array([[i, j, k] for i in shifts for j in shifts for k in shifts])
    real_translations = translations @ box_vectors
    all_points = (real_points[:, None, :] + real_translations[None, :, :]).reshape(-1, 3)
    vor = Voronoi(all_points)
    
    num_images = len(translations)
    zero_image = np.where((translations == np.array([0, 0, 0])).all(axis=1))[0][0]
    center_idx = center_index * num_images + zero_image
    
    # Identify the region of the center point
    region_index = vor.point_region[center_idx]
    region = vor.regions[region_index]
    
    if -1 in region:
        raise ValueError("Voronoi region is unbounded; consider increasing layers.")
    
    # List to store the face centers
    face_centers = []
    
    # Iterate over the ridge points to find the face centers
    for ridge_points, ridge_vertices in zip(vor.ridge_points, vor.ridge_vertices):
        if center_idx in ridge_points and all(v >= 0 for v in ridge_vertices):
            if center_surface == 0:
                p1 = all_points[ridge_points[0]]
                p2 = all_points[ridge_points[1]]
                center = (p1 + p2) / 2
            else:                   
                verts = vor.vertices[ridge_vertices]
                center = np.mean(verts, axis=0)  # Center of the face
            face_centers.append(center)
    # print(real_points[center_index])
    face_centers = np.array(face_centers) - real_points[center_index]
    inv_box = np.linalg.inv(box_vectors) 
    return face_centers, face_centers @ inv_box

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
    Crys = []
    for i in range(len(Group)):
        for s in Group[i]:
            Crys.append(s)
    return Crys

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



# 防止原子溢出边界
def Box_Slice(Str, Range, Center = [0.5 , 0.5, 0.5], Sign = 0):
    Other = []
    Structure = np.copy(Str).tolist()
    HRange = np.copy(Range) / 2
    TCenter = np.copy(Center) * Range
    UpSlice = TCenter + HRange
    DownSlice = TCenter - HRange
    delta = -0.1
    if Sign == 0:
        for s in Structure:
            a, b, c = s[0], s[1], s[2]
            B1 = (a >= DownSlice[0]) and (a <= UpSlice[0])
            B2 = (b >= DownSlice[1]) and (b <= UpSlice[1])
            B3 = (c >= DownSlice[2]) and (c <= UpSlice[2])
            if (B1 and B2 and B3):
                Other.append(s)
    else:
        for s in Structure:
            a, b, c = s[0], s[1], s[2]
            B1 = (a >= DownSlice[0]) and (a < UpSlice[0])
            B2 = (b >= DownSlice[1]) and (b < UpSlice[1])
            B3 = (c >= DownSlice[2]) and (c < UpSlice[2])
            if (B1 and B2 and B3):
                Other.append(s)
    return Other

# 镜面构造孪晶
def Mirror(s, Vector, Point):
    a, b, c = s[0], s[1], s[2]
    Down = Point
    Dis = (Vector[0]*(a-Down[0])+Vector[1]*(b-Down[1])+Vector[2]*(c-Down[2])) / \
          (math.pow(Vector[0], 2) + math.pow(Vector[1], 2) + math.pow(Vector[2], 2))
    am = a - 2 * Dis * Vector[0]
    bm = b - 2 * Dis * Vector[1]
    cm = c - 2 * Dis * Vector[2]
    return (Dis != 0), [am, bm, cm]

# 晶体绕轴旋转
def Rotate(s, Vector, angle, Point):
    s = s + [1]
    square = np.square(Vector)
    Vector = Vector/np.sqrt(np.sum(square))
    angle = math.radians(angle)
    K = 1 - math.cos(angle)
    M = np.dot(Vector, Point)
    Sin_An = math.sin(angle)
    rot = np.zeros((4, 4))
    # 构造旋转矩阵
    for i in range(3):
        rot[i, i] = math.cos(angle) + math.pow(Vector[i], 2) * (1 - math.cos(angle))
    rot[3, 3] = 1
    rot[0, 1] = Vector[0] * Vector[1] * K - Vector[2] * Sin_An
    rot[1, 0] = Vector[0] * Vector[1] * K + Vector[2] * Sin_An
    rot[0, 2] = Vector[0] * Vector[2] * K + Vector[1] * Sin_An
    rot[2, 0] = Vector[0] * Vector[2] * K - Vector[1] * Sin_An
    rot[1, 2] = Vector[1] * Vector[2] * K - Vector[0] * Sin_An
    rot[2, 1] = Vector[1] * Vector[2] * K + Vector[0] * Sin_An
    rot[0, 3] = (Point[0] - Vector[0] * M) * K + (Vector[2] * Point[1] - Vector[1] * Point[2]) * Sin_An
    rot[1, 3] = (Point[1] - Vector[1] * M) * K + (Vector[0] * Point[2] - Vector[2] * Point[0]) * Sin_An
    rot[2, 3] = (Point[2] - Vector[2] * M) * K + (Vector[1] * Point[0] - Vector[0] * Point[1]) * Sin_An
    return np.dot(rot, s)[:3]

# 晶体移动至目标位置
def Move_to_origin(Str, Range):
    Other = []
    Structure = np.copy(Str).tolist()
    X, Y, Z = Range[0], Range[1], Range[2]
    for s in Structure:
        Single = [s[0] - X, s[1] - Y, s[2] - Z]
        Other.append(Single)
    return Other

# 检验原子重合问题，将在阈值内的多原子视为重叠原子，只保留一个
def Coinciednt_atom(Alloy, Dis = 0.3):
    Structure = np.copy(Alloy).tolist()
    num = len(Structure)
    SDis = math.pow(Dis, 2)
    Str = []
    for i in range(num - 1):
        Target = Alloy[i]
        for j in range(i + 1, num):
            for s in Str:
                if Target[0] == s[0] and Target[1] == s[1] and Target[2] == s[2]:
                    break
            Other = Alloy[j]
            Dx = abs(Target[0] - Other[0])
            if Dx > Dis: continue
            Dy = abs(Target[1] - Other[1])
            if Dy > Dis: continue
            Dz = abs(Target[2] - Other[2])
            if Dz > Dis: continue
            square = math.pow(Dx, 2) + math.pow(Dx, 2) + math.pow(Dx, 2)
            if (square - SDis) <= 0:
                Str.append(Other)
                Structure.remove(Other)
        print(i)
    return Structure, Str

def Coinciednt_atom2(Alloy, Dis=0.3):
    Structure = np.copy(Alloy).tolist()
    num = len(Structure)
    SDis = math.pow(Dis, 2)
    Str = []
    index = []
    for i in range(num - 1):
        Target = Alloy[i]
        for j in range(i + 1, num):
            if i in index: break
            if j in index: continue
            Other = Alloy[j]
            Dx = abs(Target[0] - Other[0])
            if Dx > Dis: continue
            Dy = abs(Target[1] - Other[1])
            if Dy > Dis: continue
            Dz = abs(Target[2] - Other[2])
            if Dz > Dis: continue
            square = math.pow(Dx, 2) + math.pow(Dx, 2) + math.pow(Dx, 2)
            if (square - SDis) <= 0:
                Str.append(Other)
                Structure.remove(Other)
                index.append(j)
        if i % 1000 == 0:
            print(i)
    return Structure, Str

# 查找最小最大坐标，并返回值，方便移动
def coordinate_min(Alloy):
    Structure = np.copy(Alloy).tolist()
    num = len(Structure)
    Xmin, Ymin, Zmin= Alloy[0][0], Alloy[0][1], Alloy[0][2]
    for i in range(num):
        Target = Alloy[i]
        if Target[0] < Xmin: Xmin = Target[0]
        if Target[1] < Ymin: Ymin = Target[1]
        if Target[2] < Zmin: Zmin = Target[2]
    return [Xmin, Ymin, Zmin]

##### Schwarz-D 晶体 ######相关函数
def DSchwarz(Alloy, Range, Center = [1/2, 1/2, 1/2], Index = 1, eq = 2, angle = 0):
    Other = []
    X, Y, Z = Range[0], Range[1], Range[2]
    Structure = np.copy(Alloy).tolist()
    biasX, biasY, biasZ =  Center[0] - 1/2, Center[1] - 1/2, Center[2] - 1/2
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        sc =  rotation_mat @ np.array(s)
        if eq == 1:
            size = 2
            a, b, c = (sc[0] - biasX) * size, (sc[1] - biasY) * size, (sc[2]- biasZ) * size
            DSC = np.cos(a * np.pi / X) * np.cos(b * np.pi / Y) * np.cos(c * np.pi / Z) -\
              np.sin(a * np.pi / X) * np.sin(b * np.pi / Y) * np.sin(c * np.pi / Z)
        else:
            size = 1
            a, b, c = (sc[0] - biasX) * size, (sc[1] - biasY) * size, (sc[2]- biasZ) * size
            DSC = (np.sin(a * np.pi / X) * np.sin(b * np.pi / Y) * np.sin(c * np.pi / Z) +
                   np.sin(a * np.pi / X) * np.cos(b * np.pi / Y) * np.cos(c * np.pi / Z) +
                   np.cos(a * np.pi / X) * np.sin(b * np.pi / Y) * np.cos(c * np.pi / Z) +
                   np.cos(a * np.pi / X) * np.cos(b * np.pi / Y) * np.sin(c * np.pi / Z))
        if Index == 1:
            if DSC < 0:
                Other.append(s)
        elif Index ==2:
            if (DSC <= 0.1) and (DSC >= -0.1):
                Other.append(s)
        else:
            if DSC > 0:
                Other.append(s)
    return Other

##### Schwarz-G 晶体 ######相关函数
def GSchwarz(Alloy, Range, Center = [1/2, 1/2, 1/2], Index = 1):
    Other = []
    X, Y, Z = Range[0], Range[1], Range[2]
    Structure = np.copy(Alloy).tolist()
    biasX, biasY, biasZ =  Center[0] - 1/2, Center[1] - 1/2, Center[2] - 1/2
    for s in Structure:
        a, b, c = s[0] - biasX, s[1] - biasY, s[2]- biasZ
        size = 2
        GSC = (np.sin(size * a*np.pi/X)*np.cos(size * b*np.pi/Y) +
               np.sin(size * c*np.pi/Z)*np.cos(size * a*np.pi/X) +
               np.sin(size * b*np.pi/Y)*np.cos(size * c*np.pi/Z))
        if Index == 1:
            if GSC <= 0:
                Other.append(s)
        elif Index ==2:
            if (GSC <= 0.2) and (GSC >= -0.2):
                Other.append(s)
        else:
            if GSC > 0:
                Other.append(s)
    return Other


##### Schwarz-p 晶体 ######相关函数
def PSchwarz(Alloy, Range, Center = [1/2, 1/2, 1/2], Index = 1):
    Other = []
    X, Y, Z = Range[0], Range[1], Range[2]
    Structure = np.copy(Alloy).tolist()
    biasX, biasY, biasZ = Center[0] - 1/2, Center[1] - 1/2, Center[2] - 1/2
    for s in Structure:
        size = 2
        a, b, c = s[0] - biasX, s[1] - biasY, s[2]- biasZ
        PSC = np.cos(size*a*np.pi/X) + np.cos(size*b*np.pi/Y) + np.cos(size*c*np.pi/Z)
        if Index == 1:
            if PSC <= 0:
                Other.append(s)
        elif Index == 2:
            if (PSC <= 0.2) and (PSC >= -0.2):
                Other.append(s)
        else:
            if PSC > 0:
                Other.append(s)
    return Other

##### Schwarz-IWP 晶体 ######相关函数
def IWPSchwarz(Alloy, Range, Center = [1/2, 1/2, 1/2], Index = 1):
    Other = []
    X, Y, Z = Range[0], Range[1], Range[2]
    Structure = np.copy(Alloy).tolist()
    biasX, biasY, biasZ = Center[0] - 1/2, Center[1] - 1/2, Center[2] - 1/2
    for s in Structure:
        size = 2
        a, b, c = (s[0] - biasX) * size, (s[1] - biasY) * size, (s[2] - biasZ) * size
        ISC = 2*(np.cos(a*np.pi/X) * np.cos(b*np.pi/Y) +
                np.cos(b*np.pi/Y) * np.cos(c*np.pi/Z) +
                np.cos(c*np.pi/Z) * np.cos(a*np.pi/X)) - \
              (np.cos(2*a*np.pi/X) + np.cos(2*b*np.pi/Y) + np.cos(2*c*np.pi/Z))
        if Index == 1:
            if ISC <= 0:
                Other.append(s)
        elif Index == 2:
            if (ISC <= 0.2) and (ISC >= -0.2):
                Other.append(s)
        else:
            if ISC > 0:
                Other.append(s)
    return Other

##### Kelvin晶体 ######相关函数
# 定义Kelvin晶体形状
def Kelvin(Alloy, Range, Center = [0, 0, 0], Index = 1):
    delt = 0.0005
    Other = []
    Structure = np.copy(Alloy).tolist()
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        Top = [Center[0], Center[1], Center[2] + Range[2] / 2]
        Down = [Center[0], Center[1], Center[2] - Range[2] / 2]
        # 八面体
        # [111]
        Vector = [1, 1, 1]
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        B2 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= -delt
        # [-111]
        Vector = [-1, 1, 1]
        B3 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        B4 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= -delt
        # [1-11]
        Vector = [1, -1, 1]
        B5 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        B6 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= -delt
        # [11-1]
        Vector = [1, 1, -1]
        B8 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) <= delt
        B7 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= -delt
        # 截角
        # [100]
        Top = [Center[0] + Range[0] / 3, Center[1], Center[2]]
        Down = [Center[0] - Range[0] / 3, Center[1], Center[2]]
        Vector = [1, 0, 0]
        B9 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        B10 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= -delt
        # [010]
        Top = [Center[0], Center[1] + Range[0] / 3, Center[2]]
        Down = [Center[0], Center[1] - Range[0] / 3, Center[2]]
        Vector = [0, 1, 0]
        B11 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        B12 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= -delt
        # [001]
        Top = [Center[0], Center[1], Center[2] + Range[0] / 3]
        Down = [Center[0], Center[1], Center[2] - Range[0] / 3]
        Vector = [0, 0, 1]
        B13 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        B14 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= -delt
        if (B9 and B10 and B11 and B12 and B13 and B14):
            if Index == 1:
                if (B1 and B2 and B3 and B4 and B5 and B6 and B7 and B8):
                    Other.append(s)
            else:
                if (B1 and B2 and B3 and B4 and B5 and B6 and B7 and B8) == 0:
                    Other.append(s)
    return Other

# Kelvin晶体配套函数
# 找出原子最外层的八面体，用来决定是否要复制
def Sphere(s, Center, Range):
    a, b, c = s[0], s[1], s[2]
    Top = [Center[0], Center[1], Center[2] + Range[2]/2]
    Down = [Center[0], Center[1], Center[2] - Range[2]/2]
    # 八面体
    # [111]
    Vector = [1, 1, 1]
    B1 = (Vector[0]*(a-Top[0])+Vector[1]*(b-Top[1])+Vector[2]*(c-Top[2])) < 0
    B2 = (Vector[0]*(a-Down[0])+Vector[1]*(b-Down[1])+Vector[2]*(c-Down[2])) > 0
    # [-111]
    Vector = [-1, 1, 1]
    B3 = (Vector[0]*(a-Top[0])+Vector[1]*(b-Top[1])+Vector[2]*(c-Top[2])) < 0
    B4 = (Vector[0]*(a-Down[0])+Vector[1]*(b-Down[1])+Vector[2]*(c-Down[2])) > 0
    # [1-11]
    Vector = [1, -1, 1]
    B5 = (Vector[0]*(a-Top[0])+Vector[1]*(b-Top[1])+Vector[2]*(c-Top[2])) < 0
    B6 = (Vector[0]*(a-Down[0])+Vector[1]*(b-Down[1])+Vector[2]*(c-Down[2])) > 0
    # [11-1]
    Vector = [1, 1, -1]
    B7 = (Vector[0]*(a-Top[0])+Vector[1]*(b-Top[1])+Vector[2]*(c-Top[2])) > 0
    B8 = (Vector[0]*(a-Down[0])+Vector[1]*(b-Down[1])+Vector[2]*(c-Down[2])) < 0
    if (B1 and B2 and B3 and B4 and B5 and B6 and B7 and B8):
        return True

def Kelvin1(Alloy, Range, Center = [0, 0, 0]):
    delta = 0.05
    sigma = -0.05
    delta1 = 0
    sigma1 = 0
    Other = []
    Structure = np.copy(Alloy).tolist()
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        Top = [Center[0], Center[1], Center[2] + Range[2] / 2]
        Down = [Center[0], Center[1], Center[2] - Range[2] / 2]
        # 八面体
        # [111]
        Vector = [1, 1, 1]
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < sigma
        B2 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) > delta
        # [-111]
        Vector = [-1, 1, 1]
        B3 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < sigma
        B4 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) > delta
        # [1-11]
        Vector = [1, -1, 1]
        B5 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < sigma
        B6 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) > delta
        # [11-1]
        Vector = [1, 1, -1]
        B7 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) > delta
        B8 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) < sigma
        # 截角
        # [100]
        Top = [Center[0] + Range[0] / 3, Center[1], Center[2]]
        Down = [Center[0] - Range[0] / 3, Center[1], Center[2]]
        Vector = [1, 0, 0]
        B9 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < sigma1
        B10 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= delta1
        # [010]
        Top = [Center[0], Center[1] + Range[0] / 3, Center[2]]
        Down = [Center[0], Center[1] - Range[0] / 3, Center[2]]
        Vector = [0, 1, 0]
        B11 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < sigma1
        B12 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= delta1
        # [001]
        Top = [Center[0], Center[1], Center[2] + Range[0] / 3]
        Down = [Center[0], Center[1], Center[2] - Range[0] / 3]
        Vector = [0, 0, 1]
        B13 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < sigma1
        B14 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= delta1
        if (B1 and B2 and B3 and B4 and B5 and B6 and B7 and B8 and B9 and B10 and B11 and B12 and B13 and B14):
            Other.append(s)
    return Other

# 直接建立Kelvin构型，可以直接加热生成D-SC I型结构
def Kelvin_DSC_I(Size=16):
    X, Y, Z = Size, Size, Size
    Fcc = FCC()
    Range = [X, Y, Z]
    # Center = [Range[0] / 2, Range[1] / 2, Range[2] / 2]
    octahedron = {}
    # 用中心生成坐标点，设置晶相来建立晶体
    # 建立中心Kelvin晶体
    axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    Cube = Grow_Crystal(Range, Fcc, axis)
    # octahedron[0] = Kelvin1(Cube, Range)
    octahedron[0] = Kelvin(Cube, Range)
    # 建立111Kelvin晶体
    axis = [[-1, 2, 2], [2, -1, 2], [2, 2, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis)
    octahedron[1] = Kelvin1(Cube, Range)
    # 建立1-1-1Kelvin晶体
    axis = [[-1, -2, -2], [-2, -1, 2], [-2, 2, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis)
    octahedron[2] = Kelvin1(Cube, Range)
    # 建立1-11Kelvin晶体
    axis = [[-1, -2, 2], [-2, -1, -2], [2, -2, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis)
    octahedron[3] = Kelvin1(Cube, Range)
    # 建立-1-11Kelvin晶体
    axis = [[-1, 2, -2], [2, -1, -2], [-2, -2, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis)
    octahedron[4] = Kelvin1(Cube, Range)
    # 建立010Kelvin晶体
    axis = [[-1, 8, 4], [8, -1, 4], [4, 4, -7]]
    Cube = Grow_Crystal(Range, Fcc, axis)
    octahedron[5] = Kelvin1(Cube, Range)
    # 建立001Kelvin晶体
    axis = [[-7, -4, 4], [-4, -1, -8], [4, -8, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis)
    octahedron[6] = Kelvin1(Cube, Range)
    # 建立100Kelvin晶体
    axis = [[1, 8, 4], [8, 1, -4], [-4, 4, -7]]
    Cube = Grow_Crystal(Range, Fcc, axis)
    octahedron[7] = Kelvin1(Cube, Range)

    Schwarz = {}
    i = 0
    # 画中心Kelvin晶体-O
    KO = [
        [1 / 3, 1 / 3, 1 / 3],
        [1, 1, 1]
    ]
    for m in KO:
        Schwarz[i] = Move(octahedron[0], m, Range)
        i += 1
    # 画111方向Kelvin晶体-A
    KA = [
        [2 / 3, 2 / 3, 2 / 3],
    ]
    for m in KA:
        Schwarz[i] = Move(octahedron[1], m, Range)
        i += 1

    # 画-111方向Kelvin晶体-B
    KB = [
        [2 / 3, 0, 0],
    ]
    for m in KB:
        Schwarz[i] = Move(octahedron[2], m, Range)
        i += 1

    # 画1-11方向Kelvin晶体-C
    KC = [
        [0, 2 / 3, 0],
    ]
    for m in KC:
        Schwarz[i] = Move(octahedron[3], m, Range)
        i += 1

    # 画11-1方向Kelvin晶体-D
    KD = [
        [0, 0, 2 / 3],
    ]
    for m in KD:
        Schwarz[i] = Move(octahedron[4], m, Range)
        i += 1

    # 画010方向Kelvin晶体-E
    KE = [
        [1 / 3, 1, 1 / 3],
        [1, 1 / 3, 1]
    ]
    for m in KE:
        Schwarz[i] = Move(octahedron[5], m, Range)
        i += 1

    # 画001方向Kelvin晶体-F
    KF = [
        [1 / 3, 1 / 3, 1],
        [1, 1, 1 / 3]
    ]
    for m in KF:
        Schwarz[i] = Move(octahedron[6], m, Range)
        i += 1

    # 画100方向Kelvin晶体-I
    KI = [
        [1, 1 / 3, 1 / 3],
        [1 / 3, 1, 1]
    ]
    for m in KI:
        Schwarz[i] = Move(octahedron[7], m, Range)
        i += 1

    # 建立其他晶粒
    KR = [
        [0, 0, 0],
        [2 / 3, 0, 2 / 3],
        [2 / 3, 2 / 3, 0],
        [0, 2 / 3, 2 / 3],
    ]
    for m in KR:
        axis = rotation_matrix()
        Cube = Grow_Crystal(Range, Fcc, axis)
        Schwarz[i] = Move(Kelvin(Cube, Range), m, Range)
        i += 1
    All = Merge_Group(Schwarz)
    # All = Lt.Box_Slice(All, Range, Sign = 1)
    N = len(All)
    print("Ave:", N / 16)
    N = round(N / 1e4, 1)
    print("All:", N)
    Range = np.dot(Range, 4/3)
    return All, Range

  

##### fcc空间分布晶体 ######相关函数 20面体有望出现新结构
def FCC_S(Alloy, Range, Center = [0, 0, 0], Index = 1, angle = 0):
    delt = 0.0000
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    # rotation_mat = rotation_matrix(0, 0, 0)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        Top = [Center[0], Center[1], Center[2] + Range[2] / 2]
        Top = rotation_mat @ np.array(Top)
        Down = [Center[0], Center[1], Center[2] - Range[2] / 2]
        Down = rotation_mat @ np.array(Down)
        # 八面体
        # [101]
        Vector = [1, 0, 1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        B2 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= -delt
        # [10-1]
        Vector = [1, 0, -1]
        Vector = rotation_mat @ np.array(Vector)
        B3 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= delt
        B4 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) <= -delt
        # [011]
        Vector = [0, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B5 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        B6 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= -delt
        # [01-1]
        Vector = [0, 1, -1]
        Vector = rotation_mat @ np.array(Vector)
        B8 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) <= delt
        B7 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= -delt
        # 截角
        # [110]
        Top = [Center[0] + Range[0] / 2, Center[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        Down = [Center[0] - Range[0] / 2, Center[1], Center[2]]
        Down = rotation_mat @ np.array(Down)
        Vector = [1, 1, 0]
        Vector = rotation_mat @ np.array(Vector)
        B9 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        B10 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= -delt
        # [1-10]
        Top = [Center[0], Center[1] + Range[0] / 2, Center[2]]
        Top = rotation_mat @ np.array(Top)
        Down = [Center[0], Center[1] - Range[0] / 2, Center[2]]
        Down = rotation_mat @ np.array(Down)
        Vector = [1, -1, 0]
        Vector = rotation_mat @ np.array(Vector)
        B11 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= delt
        B12 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) <= -delt
        if (B9 and B10 and B11 and B12):
            if Index == 1:
                if (B1 and B2 and B3 and B4 and B5 and B6 and B7 and B8):
                # if (B1):
                    Other.append(s)
            else:
                if (B1 and B2) == 0:
                    Other.append(s)
    return Other

##### 八面体空间分布晶体 ######相关函数  无法填满周期性空间
def Octahedron(Alloy, Range, Center = [0, 0, 0], Index = 1):
    delt = 0.0005
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]] 
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        Top = [Center[0], Center[1], Center[2] + Range[2] / 2]
        Down = [Center[0], Center[1], Center[2] - Range[2] / 2]
        # 八面体
        # [111]
        Vector = [1, 1, 1]
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        B2 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= -delt
        # [-111]
        Vector = [-1, 1, 1]
        B3 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        B4 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= -delt
        # [1-11]
        Vector = [1, -1, 1]
        B5 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        B6 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) >= -delt
        # [11-1]
        Vector = [1, 1, -1]
        B8 = (Vector[0] * (a - Down[0]) + Vector[1] * (b - Down[1]) + Vector[2] * (c - Down[2])) <= delt
        B7 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= -delt
        if Index == 1:
            if (B1 and B2 and B3 and B4 and B5 and B6 and B7 and B8):
                Other.append(s)
        else:
            if (B1 and B2 and B3 and B4 and B5 and B6 and B7 and B8) == 0:
                Other.append(s)
    return Other

##### Schwarz D 空间分布晶体 ######相关函数
def Math_Schwarz_M1(Alloy, Range, Center = [0, 0, 0], Index = 1, angle = 0):
    delt = 0.0000
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        # [11-1]
        Top = [Center[0], Center[1], Center[2] + Range[2] / 8]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, 1, -1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= delt
        Top = [Center[0] + Range[2], Center[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B2 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= -delt
        # [-111]
        Top = [Center[0] + Range[2] / 8, Center[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [-1, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B3 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= delt
        Top = [Center[0], Center[1] + Range[2], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B4 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= -delt
        # [1-11]
        Top = [Center[0], Center[1] + Range[2] / 8, Center[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, -1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B5 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= delt
        Top = [Center[0] + Range[2], Center[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B6 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= -delt
        # [111]
        Top = [Center[0], Center[1] + Range[2], Center[2] + Range[2] / 8]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B7 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        Top = [Center[0], Center[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B8 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= -delt
        if (B1 and B2 and B3 and B4 and B5 and B6 and B7 and B8):
            Other.append(s)
    return Other

def Math_Schwarz_M2(Alloy, Range, Center = [0, 0, 0], Index = 1, angle = 0):
    delt = 0.0000
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        # [11-1]
        Top = [Center[0] + Range[0], Center[1] + Range[1], Center[2] + Range[2] / 8 * 7]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, 1, -1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        Top = [Center[0] + Range[0], Center[1], Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        B2 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= -delt
        # [-111]
        Top = [Center[0] + Range[0] / 8 * 7, Center[1] + Range[1], Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [-1, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B3 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        Top = [Center[0] + Range[0], Center[1] + Range[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B4 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= -delt
        # [1-11]
        Top = [Center[0] + Range[0], Center[1] + Range[1] / 8 * 7, Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, -1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B5 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= delt
        Top = [Center[0] + Range[0], Center[1] + Range[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B6 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= -delt
        # [111]
        Top = [Center[0] + Range[0], Center[1] + Range[1] / 8 * 7, Center[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B7 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) >= delt
        Top = [Center[0] + Range[0], Center[1] + Range[1], Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        B8 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) <= -delt
        if (B1 and B2 and B3 and B4 and B5 and B6 and B7 and B8):
            Other.append(s)
    return Other

def Math_Schwarz_T1(Alloy, Range, Center = [0, 0, 0], Index = 1, angle = 0):
    delt = 0
    lamb = 0
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        # [111]
        Top = [Center[0] + Range[0], Center[1] + Range[1] / 8 * 7, Center[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < delt
        Top = [Center[0], Center[1] + Range[2], Center[2] + Range[2] / 8]
        Top = rotation_mat @ np.array(Top)
        B2 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) > lamb
        if (B1 and B2):
            Other.append(s)
    Other = Box_Slice(Other, Range, Sign=1)
    return DSchwarz(Other, Range, Index=3, eq=2, Center=Center)

def Math_Schwarz_T2(Alloy, Range, Center = [0, 0, 0], angle = 0):
    delt = 0.0
    lamb = 0
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        # [11-1]
        Top = [Center[0] + Range[0], Center[1] + Range[1] / 8 * 7, Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, 1, -1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < delt
        Top = [Center[0] + Range[2] / 8, Center[1], Center[2]]
        Top = rotation_mat @ np.array(Top)
        B2 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) > lamb
        if (B1 and B2):
            Other.append(s)
    Other = Box_Slice(Other, Range, Sign = 1)
    return DSchwarz(Other, Range, Index=3, eq=2, Center=Center, angle = 180)

def Math_Schwarz_T3(Alloy, Range, Center = [0, 0, 0], angle = 0):
    delt = 0
    lamb = 0.0
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        # [-111]
        Top = [Center[0] + Range[0] / 8 * 7, Center[1], Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [-1, 1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) > lamb
        Top = [Center[0] + Range[1] / 8, Center[1], Center[2] + Range[2]]
        Top = rotation_mat @ np.array(Top)
        B2 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < delt
        if (B1 and B2):
            Other.append(s)
    Other = Box_Slice(Other, Range, Sign = 1)
    return DSchwarz(Other, Range, Index=1, eq=2, Center=Center, angle = -90)

def Math_Schwarz_T4(Alloy, Range, Center = [0, 0, 0], angle = 0):
    delt = 0
    lamb = 0.0
    Other = []
    Structure = np.copy(Alloy).tolist()
    Center = [Center[0] * Range[0], Center[1] * Range[1], Center[2] * Range[2]]
    # 生成旋转矩阵
    rotation_mat = rotation_matrix(0, 0, angle)
    for s in Structure:
        a, b, c = s[0], s[1], s[2]
        # [-111]
        Top = [Center[0] + Range[0], Center[1] + Range[1] / 8 * 1, Center[2]]
        Top = rotation_mat @ np.array(Top)
        Vector = [1, -1, 1]
        Vector = rotation_mat @ np.array(Vector)
        B1 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) < delt
        Top = [Center[0] + Range[0], Center[1] + Range[1] / 8 * 7, Center[2]]
        Top = rotation_mat @ np.array(Top)
        B2 = (Vector[0] * (a - Top[0]) + Vector[1] * (b - Top[1]) + Vector[2] * (c - Top[2])) > lamb
        if (B1 and B2):
            Other.append(s)
    Other = Box_Slice(Other, Range, Sign = 1)
    return DSchwarz(Other, Range, Index=1, eq=2, Center=Center, angle = 90)

def Math_D_SCI(Range=None, a=3.61, M = 63.546, SeedModel=0):
    if Range is None:
        value = random.random()  # 生成0到1之间的随机浮点数
        print("Range Not Input, Random value is:", value)
        Range = [value, value, value]  
    Schwarz = {}
    Fcc = FCC()
    center = [1 / 2, 1 / 2, 1 / 2]
    # 用中心生成坐标点，设置晶相来建立晶体
    axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = center)
    Schwarz[0] = Math_Schwarz_M1(Cube, Range)
    if SeedModel == 1: Cfg(Schwarz[0], Range, a, M, 'D-SCI-1')
        
    axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = center)
    Schwarz[1] = Math_Schwarz_M2(Cube, Range)
    if SeedModel == 1: Cfg(Schwarz[1], Range, a, M, 'D_SCI-2')

    # 建立111Kelvin晶体
    axis = [[-1, 2, 2], [2, -1, 2], [2, 2, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = center)
    Schwarz[2] = Math_Schwarz_T1(Cube, Range)
    if SeedModel == 1: Cfg(Schwarz[2], Range, a, M, 'D_SCI-3')

    # 建立11-1Kelvin晶体
    axis = [[-1, 2, -2], [2, -1, -2], [-2, -2, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = center)
    Schwarz[3] = Math_Schwarz_T2(Cube, Range)
    Schwarz[3] = Move(Schwarz[3], [-1 / 2, -1 / 2, 0], Range)
    if SeedModel == 1: Cfg(Schwarz[3], Range, a, M, 'D_SCI-4')

    # 建立-111Kelvin晶体
    axis = [[-1, -2, -2], [-2, -1, 2], [-2, 2, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = center)
    Schwarz[4] = Math_Schwarz_T3(Cube, Range)
    Schwarz[4] = Move(Schwarz[4], [0, -1 / 2, -1 / 2], Range)
    if SeedModel == 1: Cfg(Schwarz[4], Range, a, M, 'D_SCI-5')

    # 建立1-11Kelvin晶体
    axis = [[-1, -2, 2], [-2, -1, -2], [2, -2, -1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = center)
    Schwarz[5] = Math_Schwarz_T4(Cube, Range)
    Schwarz[5] = Move(Schwarz[5], [-1 / 2, 0, -1 / 2], Range)
    if SeedModel == 1: Cfg(Schwarz[5], Range, a, M, 'D_SCI-6')

    All = Merge_Group(Schwarz)
    # Range = np.dot(Range)
    # All = Lt.Box_Slice(All, Range, Sign = 1)
    N = len(All)
    print("Num of Atoms:", N)
    Size = round(Range[0] * a / 2 * math.sqrt(2) / 10, 2)
    print("D-SC Ds:", Size)
    return All

def Math_D_SCI_NT(Range=None, a=3.61, M = 63.546, pitch=None, head=None, roll=None, SeedModel=0):
    if Range is None:
        value = random.random()  # 生成0到1之间的随机浮点数
        print("Range Not Input, Random value is:", value)
        Range = [value, value, value]  
    Fcc = FCC()
    center = [1 / 2, 1 / 2, 1 / 2]
    DSC = {}
    axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = [1/2, 1/2, 1/2])
    DSC[0] = DSchwarz(Cube, Range, Index = 3, eq = 1, Center = center)
    axis = rotation_matrix(pitch, head, roll)
    # axis = rotation_matrix(30, 60, 27)
    # axis = [[6, -5, 0], [5, 6, 0], [0, 0, 1]]
    Cube = Grow_Crystal(Range, Fcc, axis, Center = [1/2, 1/2, 1/2])
    DSC[1] = DSchwarz(Cube, Range, Index = 1, eq = 1, Center = center)
    All = Merge_Group(DSC)
    All = Box_Slice(All, Range, Sign = 1)
    N = len(All)
    print("All:", N)
    print("Num of Atoms:", N)
    Size = round(Range[0] * a / 2 * math.sqrt(2) / 10, 2)
    print("D-SC Ds:", Size)
    return All

def Model_to_lmp(file_Str=None, file_lmp=None):
    if file_Str==None:
        print("Unfind the file of Sturcture")    
    if file_lmp==None:
        print("Don't kown how to name file of lmp, Name it model.lmp")
        file_lmp='model.lmp'
    command = (f" atomsk {file_Str} {file_lmp}")
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"成功生成lammps结构文件: {file_lmp}")
    except subprocess.CalledProcessError as e:
        print(f"生成lammps结构失败: {e}")