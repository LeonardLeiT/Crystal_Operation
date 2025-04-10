import numpy as np
import os
import Crystal_Operation as Lt
from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier

# 定义晶格常数
a = 3.61
# 定义摩尔质量
M = 63.546
sequence = np.linspace(8, 10, 50)
for num in sequence:
    All, Range = Lt.Kelvin_DSC_I(num)
    Atoms = len(All) / 1000
    folder_cfg = "model_kelvin"
    os.makedirs(folder_cfg, exist_ok=True)
    name = f'kelvin-{Atoms:3f}-{num:3f}'
    file_path = os.path.join(folder_cfg, name)
    Lt.Cfg(All, Range, a, M, file_path)
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
        percentage = (num_HCP / total_atoms) * 100
        print(f"文件: {file_path}")
        print(f"HCP_Percentage:{percentage:.2f}%")
            
        if percentage < 1.0:
            os.remove(file_path)
        else:
            print(f"  文件 {file_path} 保留\n")             
    except Exception as e:
        print(f"处理文件 {file_path} 时出错: {str(e)}\n")
