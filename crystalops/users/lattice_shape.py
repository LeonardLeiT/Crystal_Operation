import numpy as np
import random
from ..config.lattice import *
from ..config.base_operation import *
from ..config.write import *
from ..config.voronoi import *
def generate_lattice_kelvin_shape(
    base_lattice=FCC(),
    shape_lattice=BCC(), 
    size=11.63, 
    lattice_constant=3.61,
    file_name='lattice_shape',
    element=['Cu'],
    twin_type = 'DSC1'
):
    """
    Generate a crystal structure with a specified shape and orientations.
    
    Parameters:
        base_lattice: The base crystal structure. Default is FCC.
        shape_lattice: The lattice defining the target shape. Default is BCM.
        size: The size of the simulation box.
        lattice_constant: The lattice constant.
        file_name: Name of the output CFG/LAMMPS data file.
        element: List of element symbols (default: ['Cu']).
        twin_type: Preset type for position selection. 
              Options: 'full', 'DSC1', 'DSC2'. Default is 'DSC1'.
    
    Returns:
        A Lattice object representing the generated crystal structure.
    """
    allowed_types = ['BCC', 'BCT', 'BCO', 'ACO', 'HCP', 'BCM']
    shape_lattice_name = shape_lattice.name
    if shape_lattice_name not in allowed_types:
        raise ValueError(f"shape_lattice must be one of: {', '.join(allowed_types)}")
    
    # Define position presets
    type_map = {
        'full': [1, 15, 0, 14, 8, 6, 10, 4, 2, 12],
        'DSC1': [1, 15, 14, 8, 4, 2],
        'DSC2': [1, 15, 0, 6, 10, 12],
        
        'FCC': [1, 15, 6, 10, 12],
        'anti_FCC': [1, 15, 8, 4, 2],
        
        'NB': [1, 15, 14, 4, 2],
        'anti_NB': [1, 15, 0, 10, 12],
        
        'NC': [1, 15, 14, 8, 2],
        'anti_NC': [1, 15, 0, 6, 12],
        
        'ND': [1, 15, 14, 8, 4],
        'anti_ND': [1, 15, 0, 6, 10],
        
        'NAB': [1, 15, 4, 2],
        'anti_NAB': [1, 15, 10, 12],
    }
    if twin_type not in type_map:
        raise ValueError(f"type must be one of: {', '.join(type_map.keys())}")
    exist_positions = type_map[twin_type]
        
    # Define default crystal orientations
    orientations = [
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],          # O 
        [[-1, 2, 2], [2, -1, 2], [2, 2, -1]],       # A
        [[-1, -2, -2], [-2, -1, 2], [-2, 2, -1]],   # B
        [[-1, -2, 2], [-2, -1, -2], [2, -2, -1]],   # C
        [[-1, 2, -2], [2, -1, -2], [-2, -2, -1]],   # D
    ]
    
    # 1. Build crystals with different orientations
    range_box = [size, size, size]
    cube = {}
    for i, ori in enumerate(orientations):
        lattice_cube = Grow_Crystal(range_box, base_lattice, ori)
        cube[i] = lattice_cube.point @ lattice_cube.axis
    
    # 2. Compute Voronoi face centers for the shape lattice
    face_slices = []
    for i in range(len(shape_lattice.point)):
        face_centers, _ = compute_voronoi_face_centers(
            shape_lattice.point, 
            shape_lattice.axis * size, 
            center_index=i
        )
        face_slices.append(face_centers)
    
    # 3. Slice the base crystals into Voronoi-shaped domains
    voronoi_crystals = []
    
    for face_center in face_slices:
        crystals = {}
        for i in range(len(orientations)):
            crystals[i] = slice_center_points(cube[i], [0, 0, 0], face_center)
        voronoi_crystals.append(crystals)
    
    # 4. Place crystals at specific positions
    move_positions = Duplicate(shape_lattice, 2, 2, 2)
    combined_crystals = {}
    
    # Assign specific orientations to specific positions
    position_map = {
        0: (0, 1),   # Orientation O at position 1
        1: (0, 15),  # Orientation O at position 15
        2: (1, 0),   # Orientation A at position 0
        3: (1, 14),  # Orientation A at position 14
        4: (2, 8),   # Orientation B at position 8
        5: (2, 6),   # Orientation B at position 6
        6: (3, 10),  # Orientation C at position 10
        7: (3, 4),   # Orientation C at position 4
        8: (4, 2),   # Orientation D at position 2
        9: (4, 12),  # Orientation D at position 12
    }
    
    for key, (ori_idx, pos_idx) in position_map.items():
        if pos_idx in exist_positions:
            face_slice_idx = 1 if pos_idx % 2 else 0
            combined_crystals[key] = Move(
                voronoi_crystals[face_slice_idx][ori_idx], 
                move_positions.point[pos_idx], 
                range_box,
                shape_lattice.axis
            )
    
    # Assign random orientations to remaining positions
    next_key = max(position_map.keys()) + 1
    for i in range(len(move_positions.point)):
        if i in exist_positions:
            continue
        
        face_slice_idx = 1 if i % 2 else 0
        
        orientation = rotation_matrix(random.randint(0, 180), random.randint(0, 180), random.randint(0, 180))
        lattice_cube = Grow_Crystal(range_box, base_lattice, orientation)
        real_points = lattice_cube.point @ lattice_cube.axis
        real_points = slice_center_points(real_points, [0, 0, 0], face_slices[face_slice_idx])
        combined_crystals[next_key] = Move(real_points, move_positions.point[i], range_box, shape_lattice.axis)
        next_key += 1
    
    # 5. Merge all crystals into a single structure
    all_crystals = Merge_Group(combined_crystals)
    print(f'Output kelvin shape, the atoms num is {len(all_crystals)}')
    # 6. Write the structure to a LAMMPS data file
    write_lammps_data(
        real_points = all_crystals,
        axis = shape_lattice.axis,
        Size = [size * 2, size * 2, size * 2],
        point_types = np.zeros(len(all_crystals), dtype=int),
        lattice_constant = lattice_constant,
        Name = file_name,
        element = element
    )
    
    return Lattice(shape_lattice.name, all_crystals, shape_lattice.axis, [size*2, size*2, size*2])

def generate_lattice_face_shape(
    base_lattice=FCC(),
    shape_lattice=FCO(1, 1.3, 0.8),
    size=24.38,
    lattice_constant=3.61,
    file_name='FCO',
    element=['Cu']
):
    """
    Generate a crystal structure with a specified shape and orientations based on a given shape Face-Center lattice.
    
    Parameters:
    base_lattice: The base crystal structure. Default is FCC.
    shape_lattice: The lattice defining the target shape. Default is  FCO(1, 1.3, 0.8)
    size: The size of the simulation box.
    lattice_constant: The lattice constant.
    file_name: Name of the output CFG/LAMMPS data file.
    element: List of element symbols (default: ['Cu']).
    
    Returns:
        A Lattice object representing the generated crystal structure.
    """
    
    allowed_types = ['FCC', 'FCO']
    shape_lattice_name = shape_lattice.name
    
    if shape_lattice_name not in allowed_types:
        raise ValueError(f"shape_lattice must be one of: {', '.join(allowed_types)}")
    
    # Define default crystal orientations
    orientations = [
        [[1, 0, 0], [0, 1, 0], [0, 0, 1]],          # O 
        [[-1, 2, 2], [2, -1, 2], [2, 2, -1]],       # A
        [[-1, -2, -2], [-2, -1, 2], [-2, 2, -1]],   # B
        [[-1, -2, 2], [-2, -1, -2], [2, -2, -1]],   # C
        [[-1, 2, -2], [2, -1, -2], [-2, -2, -1]],   # D
    ]
    
    # 1. Build crystals with different orientations
    range_box = [size, size, size]
    cube = {}
    for i, ori in enumerate(orientations):
        lattice_cube = Grow_Crystal(range_box, base_lattice, ori)
        cube[i] = lattice_cube.point @ lattice_cube.axis
    
    # 2. Compute Voronoi face centers for the shape lattice
    face_centers, _ = compute_voronoi_face_centers(
        shape_lattice.point, 
        shape_lattice.axis * size
    )
    
    # 3. Slice the base crystals into Voronoi-shaped domains
    inv_shape_axis = np.linalg.inv(shape_lattice.axis)
    voronoi_crystals = {}
    
    for i in range(len(orientations)):
        voronoi_crystals[i] = slice_center_points(cube[i], [0, 0, 0], face_centers)
    
    # 4. Place crystals at specific positions
    move_positions = Duplicate(shape_lattice, 1, 1, 1)
    
    # 组合晶体    
    combined_crystals = {}
    # Assign specific orientations to specific positions
    combined_crystals[0] = Move(voronoi_crystals[0], move_positions.point[0], range_box, shape_lattice.axis)  # O
    combined_crystals[1] = Move(voronoi_crystals[1], move_positions.point[1], range_box, shape_lattice.axis)  # A
    combined_crystals[2] = Move(voronoi_crystals[2], move_positions.point[2], range_box, shape_lattice.axis)  # B
    combined_crystals[3] = Move(voronoi_crystals[3], move_positions.point[3], range_box, shape_lattice.axis)  # C
    
    # 5. Merge all crystals into a single structure
    all_crystals = Merge_Group(combined_crystals)
    print(f'Output {shape_lattice.name} shape, the atoms num is {len(all_crystals)}')   
    
    write_lammps_data(
        real_points = all_crystals,
        axis = shape_lattice.axis,
        Size = [size, size, size],
        point_types = np.zeros(len(all_crystals), dtype=int),
        lattice_constant = lattice_constant,
        Name = file_name,
        element = element
    )
    
    return Lattice(shape_lattice.name, all_crystals, shape_lattice.axis, [size, size, size])

def generate_lattice_cell_shape(
    base_lattice=FCC(),
    shape_lattice=FCO(1, 1.3, 0.8),
    size=24.38,
    lattice_constant=3.61,
    file_name='FCO',
    element=['Cu']
):
    """
    
    参数:
    base_lattice: 基础晶体结构，默认为FCC
    shape_lattice: 决定形状的晶格结构，默认为FCO(1, 1.3, 0.8)
    size: 模拟箱大小
    lattice_constant: 晶格常数
    atomic_mass: 原子质量
    cfg_name: 输出CFG文件名
    
    返回:
    生成的晶体结构
    """

    range_box = [size, size, size]
    cube = {}
    for i in range(len(shape_lattice.point)):
        orientation = rotation_matrix(random.randint(0, 180), random.randint(0, 180), random.randint(0, 180))
        # print(orientation)
        cube[i] = Grow_Crystal(range_box, base_lattice.point, orientation)
    print(len(cube))
    # 1. Computer Voronoi Center point
    face_slices = []
    for i in range(len(shape_lattice.point)):
        face_centers, _ = compute_voronoi_face_centers(
            shape_lattice.point, 
            shape_lattice.axis * size, 
            center_index=i
        )
        face_slices.append(face_centers)

    # 2. Slice crystal to Voronoi shape
    inv_shape_axis = np.linalg.inv(shape_lattice.axis)
    voronoi_crystals = {}

    for i in range(len(shape_lattice.point)):
        face_center = face_slices[i]      
        real_points = cube[i] @ base_lattice.axis
        real_points = slice_center_points(real_points, [0, 0, 0], face_center)
        voronoi_crystals[i]= real_points @ inv_shape_axis
    
    # 3. Move crystal to position
    combined_crystals = {}
    for i in range(len(shape_lattice.point)):
        combined_crystals[i] = Move(
                voronoi_crystals[i], 
                shape_lattice.point[i], 
                range_box
            )
    
    all_crystals = Merge_Group(combined_crystals)
    # 6. 输出为lammps_data文件
    # Cfg(
    #     all_crystals, 
    #     shape_lattice.axis, 
    #     [size*2, size*2, size*2], 
    #     lattice_constant, 
    #     file_name,
    #     element
    
    write_lammps_data(
        points = all_crystals,
        axis = shape_lattice.axis,
        Size = [size, size, size],
        point_types = np.zeros(len(all_crystals), dtype=int),
        lattice_constant = lattice_constant,
        Name = file_name,
        element = element
    )
    
    return Lattice(shape_lattice.name, all_crystals, shape_lattice.axis, [size, size, size])

def generate_lattice_kelvin_111_shape(
    base_lattice=FCC_111(),
    shape_lattice=BCC(), 
    size=11.63, 
    lattice_constant=3.61,
    file_name='kelvin_111',
    element=['Cu'],
):

    allowed_types = ['BCC', 'BCT', 'BCO', 'ABO', 'HCP', 'BCM']
    shape_lattice_name = shape_lattice.name
    
    if shape_lattice_name not in allowed_types:
        raise ValueError(f"shape_lattice must be one of: {', '.join(allowed_types)}")
        
    # Define default crystal orientations
    orientations = [
        # Y
        [0, 0, 0],
        [0 , 180 ,0],
        # Z
        [90, 0, 0],
        [-90, 0, 0],
        # X
        [0, 0, 90],
        [-180, 0, 270]
    ]
    
    # 1. Build crystals with different orientations
    range_box = [size, size, size]
    cube = {}
    for i, ori in enumerate(orientations):
        R = rotation_matrix(ori[0], ori[1], ori[2])
        lattice_cube = Grow_Crystal(range_box, base_lattice, R)
        cube[i] = lattice_cube.point @ lattice_cube.axis
    
    # 2. Compute Voronoi face centers for the shape lattice
    face_slices = []
    for i in range(len(shape_lattice.point)):
        face_centers, _ = compute_voronoi_face_centers(
            shape_lattice.point, 
            shape_lattice.axis * size, 
            center_index=i
        )
        face_slices.append(face_centers)
        
    # 3. Slice the base crystals into Voronoi-shaped domains
    voronoi_crystals = []
    
    for face_center in face_slices:
        crystals = {}
        norms = np.linalg.norm(face_center, axis=1, keepdims=True)
        unit_vectors = face_center / norms
        for i in range(len(orientations)):
            methods = np.zeros(len(face_centers), dtype=int)
            # Y
            if i == 0 or i == 1: 
                target = np.array([0, 1, 0])
                matches = np.all(unit_vectors == target, axis=1)
                methods[matches] = 0
            # Z
            if i == 2 or i == 3:   
                target = np.array([0, 0, 1])
                matches = np.all(unit_vectors == target, axis=1)
                methods[matches] = 1
            # X
            if i == 4 or i == 5:   
                target = np.array([1, 0, 0])
                matches = np.all(unit_vectors == target, axis=1)
                methods[matches] = 1
            crystals[i] = slice_center_points(cube[i], [0, 0, 0], face_center, methods)
            # write_lammps_data(
            #     real_points = crystals[i],
            #     axis = shape_lattice.axis,
            #     Size = [size, size, size],
            #     point_types = np.zeros(len(crystals[i]), dtype=int),
            #     lattice_constant = 3.61,
            #     Name = f'rotation_{i}',
            #     element = ['Cu']
            # )
        voronoi_crystals.append(crystals)
    
    # 4. Place crystals at specific positions
    move_positions = Duplicate(shape_lattice, 2, 2, 2)
    combined_crystals = {}
    
    # Assign specific orientations to specific positions
    position_map = {
        0: (0, 8),   # Orientation Y at position 8
        1: (1, 12),  # Orientation Y at position 12
        2: (2, 4),   # Orientation Z at position 4
        3: (3, 6),  # Orientation Z at position 6
        4: (4, 2),   # Orientation X at position 2
        5: (5, 10),   # Orientation X at position 10
    } 
    
    exist_positions = [8, 12, 4, 6, 2, 10]
    for key, (ori_idx, pos_idx) in position_map.items():
        if pos_idx in exist_positions:
            face_slice_idx = 1 if pos_idx % 2 else 0
            combined_crystals[key] = Move(
                voronoi_crystals[face_slice_idx][ori_idx], 
                move_positions.point[pos_idx], 
                range_box,
                shape_lattice.axis
            )
    
    # Assign random orientations to remaining positions
    # next_key = max(position_map.keys()) + 1
    # for i in range(len(move_positions.point)):
    #     if i in exist_positions:
    #         continue
        
    #     face_slice_idx = 1 if i % 2 else 0
        
    #     orientation = rotation_matrix(random.randint(0, 180), random.randint(0, 180), random.randint(0, 180))
    #     lattice_cube = Grow_Crystal(range_box, base_lattice, orientation)
    #     real_points = lattice_cube.point @ lattice_cube.axis
    #     real_points = slice_center_points(real_points, [0, 0, 0], face_slices[face_slice_idx])
    #     combined_crystals[next_key] = Move(real_points, move_positions.point[i], range_box, shape_lattice.axis)
    #     next_key += 1
    
    # 5. Merge all crystals into a single structure
    all_crystals = Merge_Group(combined_crystals)
    print(f'Output kelvin shape, the atoms num is {len(all_crystals)}')
    # 6. Write the structure to a LAMMPS data file
    write_lammps_data(
        real_points = all_crystals,
        axis = shape_lattice.axis,
        Size = [size * 2, size * 2, size * 2],
        point_types = np.zeros(len(all_crystals), dtype=int),
        lattice_constant = lattice_constant,
        Name = file_name,
        element = element
    )
    
    return Lattice(shape_lattice.name, all_crystals, shape_lattice.axis, [size*2, size*2, size*2])