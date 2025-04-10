import numpy as np
import subprocess
import math
import Crystal_Operation as Lt
import random
import os
from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier

# 定义FCC坐标
FCC = Lt.FCC()
# 定义晶格常数
a = 3.61
# 定义摩尔质量
M = 63.546
# 定义XYZ轴单位长度
# X, Y, Z = 5, 5, 5
# X, Y, Z = X * 3, Y * 3, Z * 3
# 循环生成不同的num,从5到20，
sequence = np.linspace(11.5, 13.5, 100)
for num in sequence:
    X, Y, Z = num, num, num
    Range = [X, Y, Z]
    D_SC = Lt.Math_D_SCI(Range)
    Atoms = len(D_SC) / 1000
    folder_path = "model_cfg"
    os.makedirs(folder_path, exist_ok=True)
    name = f'Schwarz-Math-{Atoms:.3f}k'
    file_path = os.path.join(folder_path, name)
    Lt.Cfg(D_SC, Range, a, M, file_path)
    file_name = file_path + '.cfg'
    try:
        # 加载当前文件
        pipeline = import_file(file_name)
        
        # 添加 HCP 结构检测
        cna = CommonNeighborAnalysisModifier()
        cna.structures[cna.Type.HCP].enabled = True
        pipeline.modifiers.append(cna)
        
        # 计算并获取结果
        data = pipeline.compute()
        structure_types = data.particles['Structure Type']
        
        # 统计 HCP 原子数量
        num_hcp = (structure_types == cna.Type.HCP).sum()
        total_atoms = len(structure_types)
        percentage = (num_hcp / total_atoms) * 100
        
        # 输出结果
        print(f"文件: {file_name}")
        print(f"  HCP原子数: {num_hcp}/{total_atoms}")
        print(f"  HCP百分比: {percentage:.2f}%")
        
        # 如果 HCP 百分比小于 5%，删除文件
        if percentage < 5:
            os.remove(file_name)
            print(f"  文件 {file_name} 已删除（HCP 百分比 < 5%）\n")
        else:
            print(f"  文件 {file_name} 保留\n")
            folder_lmp = "model_DSC"
            os.makedirs(folder_lmp, exist_ok=True)
            file_lmp = folder_lmp + "/" +name + '.lmp'
            command = (f" atomsk {file_name} {file_lmp}")
            try:
                subprocess.run(command, shell=True, check=True)
                print(f"成功生成lammps结构文件: {file_lmp}")
            except subprocess.CalledProcessError as e:
                print(f"生成lammps结构失败: {e}")                         
    except Exception as e:
        print(f"处理文件 {file_name} 时出错: {str(e)}\n")

# DSC = {}
# axis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
# Cube = Lt.Grow_Crystal(Range, FCC, axis, Center = [1/2, 1/2, 1/2])
# DSC[0] = Lt.DSchwarz(Cube, Range, Index = 3, eq = 1, Center = center)
# axis = Lt.rotation_matrix(30, 60, 27)
# # axis = [[6, -5, 0], [5, 6, 0], [0, 0, 1]]
# Cube = Lt.Grow_Crystal(Range, FCC, axis, Center = [1/2, 1/2, 1/2])
# DSC[1] = Lt.DSchwarz(Cube, Range, Index = 1, eq = 1, Center = center)
# All = Lt.Merge_Group(DSC)
# All = Lt.Box_Slice(All, Range, Sign = 1)
# N = len(All)
# print("All:", N)
# Size = round(X * a / 2 * math.sqrt(2) / 10, 2)
# print("kelvin-Size:", Size)
# name = 'Schwarz-Math-2'
# Lt.Cfg(All, Range, a, M, name)