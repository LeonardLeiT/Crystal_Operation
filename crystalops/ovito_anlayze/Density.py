import os
from ovito.io import import_file
# from ovito.modifiers import AtomicPropertyModifier
import numpy as np

# 设置原子质量（单位：g/mol）
atomic_masses = {
    1: 28.085,  # Si
    2: 26.982,  # Al
    3: 24.305,  # Mg
    4: 15.999,  # O
}

# 读取文件
pipeline = import_file("Equli_1023.15_13.data")  # 修改为你的文件路径

# 获取第一帧数据
data = pipeline.compute()

# 提取原子类型和坐标
types = data.particles['Particle Type']
positions = data.particles['Position']

# 统计每种type原子数
unique_types, counts = np.unique(types, return_counts=True)
type_counts = dict(zip(unique_types, counts))

# 获取盒子体积 (Å^3)
box = data.cell
volume_A3 = box.volume
# 转换为 cm^3
volume_cm3 = volume_A3 * 1e-24

# 计算总质量 (g)
NA = 6.02214076e23  # 阿伏伽德罗常数
total_mass_g = 0
for t, n in type_counts.items():
    mass_t = atomic_masses[t]  # g/mol
    total_mass_g += n * mass_t / NA  # 转换为 g

# 计算密度 g/cm³
density = total_mass_g / volume_cm3
print(f"Density: {density:.4f} g/cm³")
