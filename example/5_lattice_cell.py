import os
import sys
from pathlib import Path

# 添加父目录到 sys.path，以便能够导入 crystal_config 包
sys.path.append(str(Path(__file__).resolve().parent.parent))

from crystal_operation.crystal_config import generate_lattice_cell_shape, get_element
from crystal_operation.crystal_config.lattice import *

# FCC_DSC =  generate_lattice_face_shape(    
#     base_lattice=FCC(), 
#     shape_lattice=FCC(), 
#     size=5, 
#     cfg_name='Test'
# )

print(get_element('Cu'))

FCC_cell =  generate_lattice_cell_shape(    
    base_lattice=FCC(), 
    shape_lattice=BCC(), 
    size=20, 
    file_name=BCC.__name__,
    element=['Cu']
)