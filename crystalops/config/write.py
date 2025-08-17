import numpy as np
from .element import get_element
from .lattice import *

##### 写出为cfg文件 ######
def Cfg(points, axis, Size, lattice_constant=1, Name='structure', element='Cu'):
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

    element_info = get_element(element)
    # Scale lattice vectors
    real_points = np.array(points) @ np.array(axis) 
    axis_scaled = np.array(axis) * np.array(Size)
    # print(axis_scaled* lattice_constant)
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
        file.write(f"  {element_info.atomic_mass}\n")
        file.write(f"{element_info.symbol}\n")

        # Write atomic fractional coordinates
        for atom in frac_coords:
            file.write("{:.8f}    {:.8f}    {:.8f}\n".format(*atom))

def Cfg_lattice(lattice, lattice_constant=1, Name='structure', element='Cu'):
    Cfg(lattice.point, lattice.axis, lattice.size, lattice_constant, Name, element)
    
##### lammps_data ######
def write_lammps_data(real_points, axis, Size, point_types, lattice_constant=1, Name='structure', element=['Cu']):
    """
    Write a CFG file from atomic positions and lattice vectors.

    Parameters:
    - Str: (N, 3) array-like, Cartesian coordinates of atoms.
    - axis: (3, 3) array-like, lattice vectors (as rows).
    - point_type (N, 1) arrary-like
    - Size: scaling factor for lattice size (currently unused but reserved).
    - Name: string, output file name prefix (default 'structure').
    - element: string, chemical element symbol (default 'Cu').
    """
    unique_type = len(np.unique(point_types))
    if len(element) != unique_type:
        print(element)
        print(unique_type)
        raise ValueError("The num element not match atoms type")
    Natoms = len(real_points)
    if len(point_types) != Natoms:
        raise ValueError("The num atoms not match atoms type")      
    
    # Scale lattice vectors
    real_points = np.array(real_points) * lattice_constant
    axis_scaled = np.array(axis) * np.array(Size) * lattice_constant
    
    xlo, xhi = 0.0, axis_scaled[0][0]
    ylo, yhi = 0.0, axis_scaled[1][1]
    zlo, zhi = 0.0, axis_scaled[2][2]
    xy, xz, yz = axis_scaled[1][0], axis_scaled[2][0], axis_scaled[2][1]

    filename = f"{Name}.data"
    

    with open(filename, 'w') as file:
        file.write(f"# LAMMPS data file via write_lammps_data\n\n")
        file.write(f"{Natoms} atoms\n")
        file.write(f"{unique_type} atom types\n\n")

        file.write(f"{xlo:.12f} {xhi:.12f} xlo xhi\n")
        file.write(f"{ylo:.12f} {yhi:.12f} ylo yhi\n")
        file.write(f"{zlo:.12f} {zhi:.12f} zlo zhi\n")
        file.write(f"{xy:.12f} {xz:.12f} {yz:.12f} xy xz yz\n\n")

        file.write("Masses\n\n")
        for i, elem in enumerate(element):
            element_info = get_element(elem)
            file.write(f"{i+1} {element_info.atomic_mass:.4f} # {element_info.symbol}\n")

        file.write("\nAtoms  # atomic\n\n")
        for i, (pos, t) in enumerate(zip(real_points, point_types)):
            file.write(f"{i+1} {t+1} {pos[0]:.8f} {pos[1]:.8f} {pos[2]:.8f}\n")
            
def write_lammps_data_lattice(lattice, lattice_constant=1, Name='structure', element=['Cu']):
    write_lammps_data(
        real_points = lattice.point @ lattice.axis,
        axis = lattice.axis,
        Size = lattice.size,
        point_types = lattice.point_types,
        lattice_constant = lattice_constant,
        Name = Name,
        element = element)
