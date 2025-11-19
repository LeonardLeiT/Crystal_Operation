import os
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parent.parent))

from crystalops.users.lattice_shape import generate_lattice_cell_shape
from crystalops.config.element import get_element
from crystalops.config.lattice import *

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