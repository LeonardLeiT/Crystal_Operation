from crystal_operation.crystal_config.lattice_shape import generate_lattice_kelvin_shape
from crystal_operation.crystal_config.element import get_element
from crystal_operation.crystal_config.select_region import select_sphere_region
from crystal_operation.crystal_config.lattice import *
from crystal_operation.crystal_config.crystal_write import write_lammps_data_lattice

print(get_element('Cu'))

FCC_DSC =  generate_lattice_kelvin_shape(    
    base_lattice=FCC(), 
    shape_lattice=BCT(1.7), 
    size=20, 
    file_name=BCC.__name__ + '_Kelvin',
    element=['Cu']
)

move_positions = Duplicate(BCT(1.7), 2, 2, 2)
# print(move_positions.point)
exist_positions = [1, 15, 14, 8, 4, 2]

for i in range(len(move_positions.point)):
    # if i not in exist_positions:
    #     continue
    point = move_positions.point[i]
    DSC_index = select_sphere_region(FCC_DSC.point, point * FCC_DSC.size / 2, 4)
    FCC_DSC = setting_type(FCC_DSC, DSC_index, 1)

write_lammps_data_lattice(FCC_DSC, 3.61, "sphere", ['Cu', 'Cu'])