from crystalops.config.lattice import *
from crystalops.config.element import get_element
from crystalops.users.lattice_shape import *
from crystalops.config.select_region import select_sphere_region
from crystalops.config.write import write_lammps_data_lattice

print(get_element('Cu'))
# Test1: generate_lattice_kelvin_111_shape
# for idx, i in enumerate(np.linspace(0, 1, 1), start=1):
#     print(idx)
#     kelvin_111 = generate_lattice_kelvin_111_shape(    
#         base_lattice=FCC_111(), 
#         shape_lattice=BCC(), 
#         size=13+i, 
#         lattice_constant=3.61,
#         file_name=f'kelvin_111_{idx}',  # Now names are kelvin_111_1, ..., kelvin_111_50
#         element=['Cu'],
#     )

# Test2: generate_lattice_kelvin_shape
for idx, i in enumerate(np.linspace(0, 0.2, 10), start=1):
    print(idx)
    kelvin = generate_lattice_kelvin_shape(    
        base_lattice=FCC(), 
        shape_lattice=BCM(1, 1.1, 1.2, 60),    
        size=11.7+i, 
        lattice_constant=3.61,
        file_name=f'kelvin_{BCM().name}_{idx}',
        element=['Cu'],
        twin_type = 'DSC1'
    )

# Test3: generate_lattice_face_shape
# for idx, i in enumerate(np.linspace(0, 0.5, 10), start=1):
#     print(idx)
#     kelvin = generate_lattice_face_shape(    
#         base_lattice=FCC(), 
#         shape_lattice=FCC(), 
#         size=25+i, 
#         lattice_constant=3.61,
#         file_name=f'kelvin_{FCO().name}_{idx}',
#         element=['Cu'],
#     )

# choose_type = 'anti_NAB'
# FCC_DSC =  generate_lattice_kelvin_shape(    
#     base_lattice=FCC(), 
#     shape_lattice=BCC(), 
#     size=13.6, 
#     file_name = choose_type + '_Kelvin',
#     element=['Cu'],
#     twin_type = choose_type
# )

# move_positions = Duplicate(BCT(1.7), 2, 2, 2)
# # print(move_positions.point)
# exist_positions = [1, 15, 14, 8, 4, 2]

# for i in range(len(move_positions.point)):
#     # if i not in exist_positions:
#     #     continue
#     point = move_positions.point[i]
#     DSC_index = select_sphere_region(FCC_DSC.point, point * FCC_DSC.size / 2, 4)
#     FCC_DSC = setting_type(FCC_DSC, DSC_index, 1)

# write_lammps_data_lattice(FCC_DSC, 3.61, "sphere", ['Cu', 'Cu'])