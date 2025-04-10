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
sequence = np.linspace(38, 45, 10)
Nummin = 1000
Nummax = 0
for num in sequence:
    X, Y, Z = num, num, num
    Range = [X, Y, Z]
    attempt = 0
    while attempt < 3:  # 最多尝试三次
        D_SC = Lt.Math_D_SCI_NT(Range)
        Atoms = len(D_SC) / 1000
        folder_cfg = "model_cfg"
        os.makedirs(folder_cfg, exist_ok=True)
        name = f'DSC-NT-{Atoms:.0f}k'
        file_path = os.path.join(folder_cfg, name)
        Lt.Cfg(D_SC, Range, a, M, file_path)
        file_path = file_path + '.cfg'
        try:
            # 加载当前文件
            pipeline = import_file(file_path)
            
            # 添加 OTHER 结构检测
            cna = CommonNeighborAnalysisModifier()
            cna.structures[cna.Type.OTHER].enabled = True
            pipeline.modifiers.append(cna)
            
            # 计算并获取结果
            data = pipeline.compute()
            structure_types = data.particles['Structure Type']
            
            # 统计 OTHER 原子数量
            num_OTHER = (structure_types == cna.Type.OTHER).sum()
            total_atoms = len(structure_types)
            percentage = (num_OTHER / total_atoms) * 100
            
            # 输出结果
            print(f"文件: {file_path}")
            print(f"  OTHER原子数: {num_OTHER}/{total_atoms}")
            print(f"  OTHER百分比: {percentage:.2f}%")
            
            if percentage > 16 or percentage <10:
                os.remove(file_path)
                print(f"  文件 {file_path} 已删除（Other 百分比 < 10% and >17%）\n")
                attempt += 1
            else:
                print(f"  文件 {file_path} 保留\n")    
                if num < Nummin: Nummin = num
                if num > Nummax: Nummax = num 
                break               
        except Exception as e:
            print(f"处理文件 {file_path} 时出错: {str(e)}\n")
print("范围为：", Nummin, " to ", Nummax)