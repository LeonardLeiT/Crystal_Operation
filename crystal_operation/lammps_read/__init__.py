# 导入所有模块
from . import read_lammps_dump

# 将模块中的函数添加到包的命名空间
from .read_lammps_dump import *

# 可选：定义 __all__ 变量，指定哪些函数是公共的
__all__ = []
__all__.extend(read_lammps_dump.__all__)
