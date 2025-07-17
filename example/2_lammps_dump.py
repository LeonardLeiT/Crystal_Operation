import numpy as np
import pandas as pd
import sys
from pathlib import Path

# 添加父目录到 sys.path，以便能够导入 crystal_config 包
sys.path.append(str(Path(__file__).resolve().parent.parent))
from crystal_operation.crystal_config.lattice import *
from crystal_operation.lammps_read import *

filename = 'test.dump'
filename = 'filename'

dump_reader = read_dump_lammps(filename)
print(dump_reader)