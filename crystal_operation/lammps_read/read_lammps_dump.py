__all__ = ['read_dump_lammps']

import numpy as np
import pandas as pd
from ..crystal_config.lattice import Lattice

class LammpsDumpReader:
    """LAMMPS Dump文件读取器类"""
    
    def __init__(self, filename: str):
        """
        初始化读取器
        
        参数:
        - filename: LAMMPS dump文件路径
        """
        self.filename = filename
        self.file_num = 0
        self.lattices = []  # 存储每一帧的Lattice对象
    
    def read_file(self) -> None:
        """读取并解析整个dump文件，创建Lattice对象列表"""
        with open(self.filename, 'r') as f:
            lines = f.readlines()
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            
            # 检查是否是时间步标记
            if line.startswith("ITEM: TIMESTEP"):
                i += 1
                timestep = int(lines[i].strip())  
                
                i += 2
                n_atoms = int(lines[i].strip())
                
                i += 1
                box_line = lines[i].strip()
                boundary_types = box_line.split()[3:]

                box_bounds = []
                size = np.array([1, 1, 1])
                if len(boundary_types) == 3:
                    for _ in range(3):
                        axis = np.array([
                            [1, 0, 0],
                            [0, 1, 0],
                            [0, 0, 1]
                        ])
                        
                        i += 1
                        bounds_line = lines[i].strip()
                        min_val, max_val = map(float, bounds_line.split())
                        box_bounds.append((min_val, max_val)) 
                    xlo, xhi = box_bounds[0]
                    ylo, yhi = box_bounds[1]
                    zlo, zhi = box_bounds[2]
                    size = np.array([xhi - xlo, yhi - ylo, zhi - zlo])

                i += 1
                header_line = lines[i].strip()
                headers = header_line.split()[2:]
                atoms = []
                for _ in range(n_atoms):
                    i += 1
                    atom_line = lines[i].strip()
                    atom_data = atom_line.split()
                    atoms.append(atom_data)
                df = pd.DataFrame(atoms, columns=headers)
                
                if len(boundary_types) == 3:
                    point_columns = ['x', 'y', 'z']
                    if not all(col in df.columns for col in point_columns):
                        raise ValueError("Dump文件中缺少坐标列")
                    points = df[point_columns].values
                    other_columns = [col for col in df.columns if col not in ['id', 'type', 'x', 'y', 'z']]
                    point_properties = df[other_columns] if other_columns else None
                point_types = df['type'].values if 'type' in df.columns else None
                
                lattice_name = f"timestep_{timestep}"
                lattice = Lattice(
                    name=lattice_name,
                    point=points,
                    axis=axis,
                    size=size,
                    point_types=point_types,
                    point_properties=point_properties
                )
                self.lattices.append(lattice)  
            i += 1
        self.file_num = len(self.lattices)
    
    def get_lattices(self) -> list:
        """获取所有帧的Lattice对象列表"""
        return self.lattices

def read_dump_lammps(filename: str) -> Lattice:
    """
    读取LAMMPS dump文件并返回第一个帧的Lattice对象
    
    参数:
    - filename: LAMMPS dump文件路径
    
    返回:
    - Lattice对象
    """
    reader = LammpsDumpReader(filename)
    reader.read_file()
    
    if not reader.get_lattices():
        raise ValueError("未能从dump文件中读取到有效数据")
    
    return reader.get_lattices()
# 示例使用
if __name__ == "__main__":
    # 读取dump文件并创建Lattice对象
    lattice = read_dump_lammps("example.dump")
    
    # 打印基本信息
    print(f"Lattice名称: {lattice.name}")
    print(f"原子数量: {len(lattice.point)}")
    print(f"盒子向量:\n{lattice.axis}")
    
    # 打印前10个原子的信息
    print("\n前10个原子的信息:")
    print("索引\t类型\t坐标")
    for i in range(min(10, len(lattice.point))):
        point_type = lattice.point_types[i]
        point = lattice.point[i]
        print(f"{i}\t{point_type}\t{point}")
    
    # 如果有结构类型信息，打印分布
    if lattice.point_properties and 'StructureType' in lattice.point_properties:
        struct_types = np.array(lattice.point_properties['StructureType'])
        unique_types, counts = np.unique(struct_types, return_counts=True)
        
        print("\n结构类型分布:")
        for t, c in zip(unique_types, counts):
            print(f"类型 {t}: {c}个原子 ({c/len(struct_types)*100:.2f}%)")
