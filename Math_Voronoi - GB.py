import numpy as np
import subprocess
import math
import Crystal_Operation as Lt
import random
import pandas as pd
import os
from ovito.io import import_file
from ovito.modifiers import CommonNeighborAnalysisModifier
import shutil
def remove_unwanted_folders(directory):
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        if os.path.isdir(item_path):
            # 如果是文件夹，删除它
            shutil.rmtree(item_path)
        elif os.path.isfile(item_path):
            # 如果是文件，检查是否是 polycrystal.txt 或以 .lmp 结尾
            if item == 'polycrystal.txt' or item.endswith('.lmp'):
                continue
            else:
                # 删除不需要的文件
                os.remove(item_path)


# 定义多晶结构的参数
material = "Cu"  
lattice_constant = 3.61
lattice_structure = "fcc"

folder_path = "model_Polycrystal"
shutil.rmtree(folder_path)
os.makedirs(folder_path, exist_ok=True)
unitcell_file = os.path.join(folder_path, "Cu.xsf")

Size_I_need = [910, 2088, 4192, 4962, 9640, 11207, 13685, 16710, 23819, 26500, 39202, 44164, 55708, 55882, 76948, 77584, 105880, 128851, 132180, 168018, 209294]
# Size_I_need = [910, 2088]
num_grains = 2

polysize = 22
result = pd.DataFrame(columns=['PolySize', 'Ave_num', '1', '2', '3', '4', '5'])
for i in range(len(Size_I_need)):
    print(f"Finding the same size of {Size_I_need[i]}")
    num_need = Size_I_need[i]
    num_error = -1
    run_time = 0
    while num_error < 0:
        run_time += 1
        if run_time>200:
            print(f'{num_need}超出计算极限，需要更换')
            break
        remove_unwanted_folders(folder_path)
        print("生成晶粒参数文件...")
        grain_file = f"{folder_path}/polycrystal.txt"
        with open(grain_file, "w") as f:
            f.write(f"box {polysize} {polysize} {polysize}\n")
            f.write(f"random {num_grains}\n")
        command = f"atomsk --create {lattice_structure} {lattice_constant} {material} {unitcell_file}"
        subprocess.run(command, shell=True, check=True)
        polycrystal_cfg  = f"{folder_path}/Voronoi_{polysize:.2f}.cfg"
        command = f"atomsk --polycrystal {unitcell_file} {grain_file} {polycrystal_cfg} -wrap"
        try:
            # 执行命令
            subprocess.run(command, shell=True, check=True)
            try:
                # 加载当前文件
                pipeline = import_file(polycrystal_cfg)
                    
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
                num_error = total_atoms-num_need
                percentage = (num_HCP / total_atoms) * 100
                print(f"文件: {polycrystal_cfg}")
                print(f"F_value:{num_error:.2f}%")
                print(f"HCP_Percentage:{percentage:.2f}%")
                all_atoms = total_atoms / 1000
                print(f"polysize:{polysize}%")
                if  num_error < 0:
                    polysize += 0.1
                else:
                    break
            except Exception as e:
                print(f"处理文件 {polycrystal_cfg} 时出错: {str(e)}\n")
        except subprocess.CalledProcessError as e:
            print(f"生成多晶结构失败: {e}")
    num_error = total_atoms-num_need
    if num_error < 0:
        print(f"Can't find size {num_need}, the current size is {total_atoms}, and curren ploysize {polysize}")
        break
    print("================================================================================")
    print(f"The size I need is {polysize}")
    Ave_num = 0
    OtherPercentage = []
    for j in range(5):
        remove_unwanted_folders(folder_path)
        print("生成晶粒参数文件...")
        grain_file = f"{folder_path}/polycrystal.txt"
        with open(grain_file, "w") as f:
            f.write(f"box {polysize} {polysize} {polysize}\n")
            f.write(f"random {num_grains}\n")
        command = f"atomsk --create {lattice_structure} {lattice_constant} {material} {unitcell_file}"
        subprocess.run(command, shell=True, check=True)
        polycrystal_cfg = f"{folder_path}/Voronoi_{polysize:.2f}.cfg"
        command = f"atomsk --polycrystal {unitcell_file} {grain_file} {polycrystal_cfg} -wrap"
        try:
            # 执行命令
            subprocess.run(command, shell=True, check=True)
            try:
                # 加载当前文件
                pipeline = import_file(polycrystal_cfg)
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
                print(f"文件: {polycrystal_cfg}")
                print(f"percentage:{percentage:.2f}%")
                OtherPercentage.append(percentage)
                Ave_num += total_atoms
                polycrystal_lmp = f"{folder_path}/Voronoi_{num_need}_{j}.lmp"
                command = f"atomsk {polycrystal_cfg} {polycrystal_lmp}"
                subprocess.run(command, shell=True, check=True)
            except Exception as e:
                print(f"处理文件 {polycrystal_cfg} 时出错: {str(e)}\n")
        except subprocess.CalledProcessError as e:
            print(f"生成多晶结构失败: {e}")
    Ave_num = Ave_num / len(OtherPercentage)
    new_row = {'PolySize': polysize, 'Target_num': num_need, 'Ave_num': Ave_num, '1': OtherPercentage[0], '2': OtherPercentage[1], '3': OtherPercentage[2], '4': OtherPercentage[3], '5': OtherPercentage[4]}
    result = pd.concat([result, pd.DataFrame([new_row])], ignore_index=True)
    result.to_excel("Voronoi_GB.xlsx", index=False)

