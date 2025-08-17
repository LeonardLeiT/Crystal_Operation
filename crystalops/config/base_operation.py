import numpy as np
from .lattice import Lattice, Duplicate

def slice_center_point(alloy, center, point, method=0):
    """
    Slice the structure by a plane passing through a given point and perpendicular 
    to the vector from center to point.

    Removes all particles on the positive side of the plane defined by the vector 
    from `center` to `point`.

    Args:
        alloy (ndarray): (n_particles, 3) array of particle positions.
        center (ndarray): (3,) array representing the reference center.
        point (ndarray): (3,) array defining the slicing plane.

    Returns:
        ndarray: Filtered subset of `alloy`.
    """
    normal_vector = point - center
    dot_products = np.sum(normal_vector * (alloy - point), axis=1)
    mask = dot_products <= -method
    return alloy[mask]

def slice_center_points(alloy, center, points, methods=None):
    """
    Apply multiple slicing planes to the structure.

    Each plane is defined by a point and normal vector from `center` to that point.

    Args:
        alloy (ndarray): (n_particles, 3) array of particle positions.
        center (ndarray): (3,) array representing the reference center.
        points (ndarray): (n_planes, 3) array, each row is a slicing point.


    Returns:
        ndarray: Filtered alloy points after all slicing operations.
    """
    if methods is None:
        methods = np.zeros(len(points), dtype=int)
    elif len(methods) != len(points):
        raise ValueError("Length of methods must match number of points.")

    sliced = alloy.copy()
    for point, method in zip(points, methods):
        sliced = slice_center_point(sliced, center, point, method)
    return sliced

def Move_lattice(lattice, vector):
    """
    Translate a lattice structure by a scaled vector.

    Computes the shift vector by scaling `vector` with the lattice size, 
    and applies it to the lattice points.

    Args:
        lattice: An object with attributes `point` (positions) and `size` .
        vector (array-like): Translation vector (3,).

    Returns:
        The updated `lattice` object with shifted points.
    """
    shift = np.array(vector) * lattice.size
    lattice.point = lattice.point + shift
    return lattice

def Move(points, vector, range, axis = np.eye(3)):
    """
    Translate a set of points by a scaled vector.

    Computes the shift vector by element-wise multiplication of `vector` and `range`,
    and applies it to all points.

    Args:
        points (ndarray): Array of shape (npoints, 3) representing positions.
        vector (array-like): Translation vector (3,).
        range (array-like): Scaling factors (3,).

    Returns:
        ndarray: Translated points with same shape as `points`.
    """
    axis_scaled = np.array(axis) * np.array(range)
    shift = np.array(vector) @ axis_scaled
    Move = points + shift
    return Move

def Orhtogonality(axis):
    """
    Check orthogonality of a set of vectors.

    Tests whether all pairs of vectors in `axis` are mutually orthogonal
    within a numerical tolerance.

    Args:
        axis (array-like): Sequence of vectors, each of shape (3, 3).

    Returns:
        bool: True if all vectors are orthogonal to each other, False otherwise.
    """
    lenth = len(axis)
    for i in range(lenth-1):
        for j in range(i + 1, lenth):
            if np.dot(axis[i], axis[j]) >= 1e-14:
                return False
    return True

def Unitization(Vector):
    """
    Normalize a vector to unit length.

    Scales the input vector so that its Euclidean norm becomes 1.

    Args:
        Vector (array-like): Input vector of shape (3, 3).

    Returns:
        ndarray: Normalized vector of same shape as `Vector`.
    """
    Vec = np.copy(Vector)
    square = np.square(Vec)
    return Vec / np.sqrt(np.sum(square))

def Grow_Crystal(Range, lattice,
                 orientation=[[1, 0, 0], [0, 1, 0], [0, 0, 1]]):
    """
    Generate a crystal structure by replicating a crystallographic basis within a specified range.

    Builds a 3D crystal by placing copies of the `crystallographic` basis atoms on a lattice
    defined by the orthogonal vectors `axis`. The structure is centered at the origin.

    Args:
        Range (ndarray): Array of shape (3, ), Approximate extent of the crystal in lattice units.
        crystallographic (ndarray): Array of shape (n_atoms_basis, 3) specifying the fractional positions of atoms in the unit cell.
        axis (array-like, optional): List of three orthogonal vectors defining the lattice directions.
            Default is Cartesian axes [[1,0,0], [0,1,0], [0,0,1]].

    Returns:
        ndarray or None: Array of shape (n_atoms, 3) with the positions of all atoms in the generated crystal.
            Returns None if `axis` is not orthogonal.
    """

    if not Orhtogonality(orientation):
        print('orientation is not orthogonal, please change!')
        return None

    # Normalize axis vectors
    orientation = np.array([Unitization(a) for a in orientation])
    # Transform crystallographic basis atoms
    basis_atoms = lattice.point
    # Estimate replication range
    R = int(np.ceil(np.linalg.norm(Range)))
    # Generate lattice translations
    shifts = np.array([[i, j, k] for i in range(R) for j in range(R) for k in range(R)])
    displacements = shifts

    # Generate full structure
    structure = (basis_atoms[:, None, :] + displacements[None, :, :]).reshape(-1, 3)

    # Center the structure around the desired Center point (default [0,0,0])
    centroid = structure.mean(axis=0)
    structure -= centroid
    return Lattice(lattice.name, structure, lattice.axis @ orientation, R)

def Merge(Crys1, Crys2):
    """
    Merge two crystal structures into a single list.

    Concatenates all elements of `Crys1` and `Crys2` into a new list.

    Args:
        Crys1 (list): First crystal structure represented as a list.
        Crys2 (list): Second crystal structure represented as a list.

    Returns:
        list: Combined list containing all elements of `Crys1` and `Crys2`.
    """
    Crys3 = []
    for s in Crys1:
        Crys3.append(s)
    for s in Crys2:
        Crys3.append(s)
    return Crys3

def Merge_Group(Group):
    """
    Merge multiple crystal structures from a dictionary into a single array.

    Takes all values from `Group` and stacks them vertically into a single numpy array.

    Args:
        Group (dict): Dictionary where each value is an array-like structure representing a crystal.

    Returns:
        ndarray: Combined array containing all structures stacked along the first axis.
    """
    return np.vstack(list(Group.values()))

def rotation_matrix(x_angle, y_angle, z_angle):
    """
    Compute a 3D crystal rotation matrix for given Euler angles.

    Constructs the combined rotation matrix by sequentially rotating about the x, y, and z axes
    by the specified angles (in degrees). The rotations are applied in the order: Rx → Ry → Rz.

    Args:
        x_angle (float): Rotation angle around the x-axis, in degrees.
        y_angle (float): Rotation angle around the y-axis, in degrees.
        z_angle (float): Rotation angle around the z-axis, in degrees.

    Returns:
        ndarray: A (3, 3) rotation matrix representing the combined rotation.
    """
    x_angle = np.radians(x_angle)
    y_angle = np.radians(y_angle)
    z_angle = np.radians(z_angle)
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

    