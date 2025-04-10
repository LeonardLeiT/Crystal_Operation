import subprocess
import math
import Crystal_Operation as Lt
import random
import os
from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier

import numpy as np
from scipy.optimize import minimize
F_Save = 100
Angle = 30
def Math_DSCI_NT_Min(x):
    global F_Save
    global Angle
    a = 3.61
    M = 63.546
    num = 36
    pitch= x
    head=0
    roll=0
    X, Y, Z = num, num, num
    Range = [X, Y, Z]
    D_SC = Lt.Math_D_SCI_NT(Range,a, M, pitch, head, roll)
    Atoms = len(D_SC) / 1000
    folder_cfg = "model_cfg"
    os.makedirs(folder_cfg, exist_ok=True)
    name = f'DSC-NT-{Atoms:.2f}k-{x:.2f}'
    file_path = os.path.join(folder_cfg, name)
    Lt.Cfg(D_SC, Range, a, M, file_path)
    file_path = file_path + '.cfg'
    try:
        # 加载当前文件
        pipeline = import_file(file_path)
            
        # 添加 OTHER 结构检测
        cna = CommonNeighborAnalysisModifier()
        cna.structures[cna.Type.OTHER].enabled = True
        cna.structures[cna.Type.HCP].enabled = True
        pipeline.modifiers.append(cna)
            
        # 计算并获取结果
        data = pipeline.compute()
        structure_types = data.particles['Structure Type']
            
        # 统计 OTHER 原子数量
        num_OTHER = (structure_types == cna.Type.OTHER).sum()
        num_HCP = (structure_types == cna.Type.HCP).sum()
        total_atoms = len(structure_types)
        percentage = (num_OTHER / total_atoms) * 100
        F = percentage-13.1
        percentage = (num_HCP / total_atoms) * 100
        print(f"文件: {file_path}")
        print(f"F_value:{F:.2f}%")
        print(f"HCP_Percentage:{percentage:.2f}%")
        # # 输出结果
        # print(f"文件: {file_path}")
        # print(f"  OTHER原子数: {num_OTHER}/{total_atoms}")
        # print(f"  OTHER百分比: {percentage:.2f}%")
            
        if F > 1 or F < 0 or percentage > 0.2:
            os.remove(file_path)
            # print(f"  文件 {file_path} 已删除（Other 百分比 < 10% and >17%）\n")
        else:
            # F_Save = F
            # Angle = x
            print(f"  文件 {file_path} 保留\n")             
    except Exception as e:
        print(f"处理文件 {file_path} 时出错: {str(e)}\n")
    return F

def find_minimum(bounds, initial_guess=None):
    if initial_guess is None:
        initial_guess = np.array([(low + high) / 2 for low, high in bounds])
    # 使用支持边界约束的 L-BFGS-B 方法
    result = minimize(Math_DSCI_NT_Min, initial_guess, method='CG', bounds=bounds)
    return result

# bounds = [(35, 43), (0, 180)]
# initial_guess = [39, 30]
# res = find_minimum(bounds, initial_guess)
# print("F 的最小值为：", res.fun)
# print("对应的变量值为：", res.x)

sequence = np.linspace(10, 80, 50)
for i in sequence:
    error = Math_DSCI_NT_Min(i)