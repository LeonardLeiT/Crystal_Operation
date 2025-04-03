import numpy as np
import subprocess
import math
import Crystal_Operation as Lt
import random
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

folder_path = "model_Voronoi"
os.makedirs(folder_path, exist_ok=True)
unitcell_file = os.path.join(folder_path, "Cu.xsf")

num_I_need = [2843, 4122, 6789, 9640, 13685]
error_I_accpet=[200, 250, 300, 350, 400]
sizes = [32.7, 37, 43.64, 48.8, 54.55]
num_grains = 2

for i in range(len(num_I_need)):
    structure = 0
    run_time = 0
    while structure < 5:
        run_time += 1
        if run_time>50:
            print(f'{sizes[i]}超出计算极限，需要更换')
            break
        remove_unwanted_folders(folder_path)
        print("生成晶粒参数文件...")
        grain_file = f"{folder_path}/polycrystal.txt"
        with open(grain_file, "w") as f:
            f.write(f"box {sizes[i]} {sizes[i]} {sizes[i]}\n")
            f.write(f"random {num_grains}\n")
        command = f"atomsk --create {lattice_structure} {lattice_constant} {material} {unitcell_file}"
        subprocess.run(command, shell=True, check=True)
        polycrystal_cfg  = f"{folder_path}/Voronoi_{sizes[i]:.2f}.cfg"
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
                F = total_atoms-num_I_need[i]
                percentage = (num_HCP / total_atoms) * 100
                print(f"文件: {polycrystal_cfg}")
                print(f"F_value:{F:.2f}%")
                print(f"HCP_Percentage:{percentage:.2f}%")
                all_atoms = total_atoms / 1000
                if  F > error_I_accpet[i] or F < 0 or percentage > 0.2:
                    os.remove(polycrystal_cfg)
                else:
                    print(f"  文件 {polycrystal_cfg} 保留\n")
                    polycrystal_lmp  = f"{folder_path}/Voronoi_{num_I_need[i]}_{structure}.lmp"
                    command = f"atomsk {polycrystal_cfg} {polycrystal_lmp}" 
                    subprocess.run(command, shell=True, check=True)     
                    structure += 1             
            except Exception as e:
                print(f"处理文件 {polycrystal_cfg} 时出错: {str(e)}\n")
        except subprocess.CalledProcessError as e:
            print(f"生成多晶结构失败: {e}")
