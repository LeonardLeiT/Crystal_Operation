import os
import sys
from pathlib import Path

# 添加父目录到 sys.path，以便能够导入 crystal_config 包
sys.path.append(str(Path(__file__).resolve().parent.parent))

from crystal_operation.crystal_config import generate_lattice_kelvin_shape, generate_lattice_face_shape, get_element, write_lammps_data, select_sphere_region, setting_type, write_lammps_data_lattice
from crystal_operation.crystal_config.lattice import *

# FCC_DSC =  generate_lattice_face_shape(    
#     base_lattice=FCC(), 
#     shape_lattice=FCC(), 
#     size=5, 
#     cfg_name='Test'
# )

print(get_element('Cu'))

FCC_DSC =  generate_lattice_kelvin_shape(    
    base_lattice=FCC(), 
    shape_lattice=BCT(1.7), 
    size=20, 
    file_name=BCC.__name__ + '_Kelvin',
    element=['Cu']
)

move_positions = Duplicate(BCC(), 2, 2, 2)
# print(move_positions.point)
exist_positions = [1, 15, 14, 8, 4, 2]

for i in range(len(move_positions.point)):
    # if i not in exist_positions:
    #     continue
    point = move_positions.point[i]
    DSC_index = select_sphere_region(FCC_DSC.point, point * FCC_DSC.size / 2, 4)
    FCC_DSC = setting_type(FCC_DSC, DSC_index, 1)

# print(FCC_DSC.point)

write_lammps_data_lattice(FCC_DSC, 3.61, "sphere", ['Cu', 'Cu'])