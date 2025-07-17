from crystalops.config.lattice import *
from crystalops.config.element import get_element
from crystalops.users.lattice_shape import generate_lattice_kelvin_shape
from crystalops.config.select_region import select_sphere_region
from crystalops.config.crystal_write import write_lammps_data_lattice

print(get_element('Cu'))

choose_type = 'anti_FCC'
FCC_DSC =  generate_lattice_kelvin_shape(    
    base_lattice=FCC(), 
    shape_lattice=BCC(), 
    size=13.6, 
    file_name = choose_type + '_Kelvin',
    element=['Cu'],
    twin_type = choose_type
)

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