import os
import pandas as pd
import matplotlib.pyplot as plt

time = 1
# 设置数据文件夹路径
data_dir = f"data"  # 数据文件夹路径
file_pattern = "enthalpy"  # 文件名前缀

# 获取所有符合条件的文件名
files = [f for f in os.listdir(data_dir) if f.startswith(file_pattern) and f.endswith("100.txt")]

# 检查是否找到文件
if not files:
    raise FileNotFoundError(f"未找到符合条件的文件。请检查文件夹路径和文件名格式: {data_dir}")

print(f"找到 {len(files)} 个文件: {files}")

# 初始化一个空的DataFrame来存储所有数据
all_data = pd.DataFrame()

# 读取每个文件并合并数据
for file in files:
    # 构建完整文件路径
    file_path = os.path.join(data_dir, file)

    # 检查文件是否为空
    if os.stat(file_path).st_size == 0:
        print(f"文件 {file} 是空的，跳过该文件。")
        continue

    # 读取文件数据，跳过以 '#' 开头的注释行
    try:
        df = pd.read_csv(file_path, delim_whitespace=True, header=None, names=["Temperature", "Potential", "Volume"], comment="#")
    except Exception as e:
        print(f"读取文件 {file} 时出错: {e}")
        continue  # 跳过出错的文件

    # 提取结构名称 {Str}
    structure = file.replace("enthalpy", "").replace("w-100.txt", "")  # 提取 {Str} 部分

    # 添加结构名称作为新列
    # df["Structure"] = float(structure)
    df["Structure"] = structure
    # 将当前文件的数据合并到总数据中
    all_data = pd.concat([all_data, df], ignore_index=True)

Potential = pd.DataFrame()
# 遍历所有不同的结构，并横向拼接数据
for structure, group in all_data.groupby("Structure"):
    # 重置索引以确保拼接时索引对齐
    group_data = group[["Temperature", "Potential"]].reset_index(drop=True)
    # 修改列名，避免不同结构的数据覆盖
    group_data.columns = [f"Temp_{structure}", f"Potential_{structure}"]
    # 横向拼接
    Potential = pd.concat([Potential, group_data], axis=1)
# 查看结果
Potential.to_excel(f'Potential-{time}.xlsx', index=False)

Volume = pd.DataFrame()
# 遍历所有不同的结构，并横向拼接数据
for structure, group in all_data.groupby("Structure"):
    # 重置索引以确保拼接时索引对齐
    group_data = group[["Temperature", "Volume"]].reset_index(drop=True)
    # 修改列名，避免不同结构的数据覆盖
    group_data.columns = [f"Temp_{structure}", f"Vol_{structure}"]
    # 横向拼接
    Volume = pd.concat([Volume, group_data], axis=1)
# 查看结果
Volume.to_excel(f'Volume-{time}.xlsx', index=False)