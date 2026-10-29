import os
import pandas as pd
import matplotlib.pyplot as plt

# 设置数据文件夹路径
data_dir = "data"  # 数据文件夹路径
file_pattern = "enthalpy"  # 文件名前缀

# 获取所有符合条件的文件名
files = [f for f in os.listdir(data_dir) if f.startswith(file_pattern) and f.endswith("5000.txt")]

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

# all_data = all_data.sort_values(by=["Structure", "Temperature"])
# 检查数据
print("汇总数据的前几行:")
print(all_data.head())

# 保存汇总数据为Excel文件
output_excel = "all_data_summary.xlsx"
all_data.to_excel(output_excel, index=False)
print(f"汇总数据已保存为: {output_excel}")

# 绘制温度和势能的汇总图
plt.figure(figsize=(10, 6))
for structure, group in all_data.groupby("Structure"):
    plt.plot(group["Temperature"], group["Potential"], label=f"{structure} w")
plt.xlabel("Temperature K")
plt.ylabel("Potential")
# plt.title("Temperature vs Potential")
plt.legend()
plt.grid(True)
plt.savefig("temperature_vs_potential.png")  # 保存图像
plt.show()

# 绘制温度和体积的汇总图
plt.figure(figsize=(10, 6))
for structure, group in all_data.groupby("Structure"):
    plt.plot(group["Temperature"], group["Volume"], label=f"{structure} w")
plt.xlabel("Temperature K")
plt.ylabel("Volume")
# plt.title("Temperature vs Volume")
plt.legend()
plt.grid(True)
plt.savefig("temperature_vs_volume.png")  # 保存图像
plt.show()

import pandas as pd

# Read the data from the Excel file
df = pd.read_excel('all_data_summary.xlsx')

# Create an empty list to hold the processed data
processed_data = []

# Get the unique values in the 'Structure' column
structures = df['Structure'].unique()

# For each structure, process the corresponding rows
for structure in structures:
    structure_data = df[df['Structure'] == structure]

    # Start with the header for the current structure
    header = ['Temp', f'Structure {structure}']

    # Append the temperature and structure values to the header
    temp_values = structure_data['Temperature'].values
    structure_values = structure_data['Structure'].values

    # Combine 'Temp' and 'Structure' values
    temp_structure_combined = list(zip(temp_values, structure_values))

    # Add 'Potential' values in the corresponding row
    potential_values = structure_data['Potential'].values

    # Add all this information to the processed data
    processed_data.append(header)
    processed_data.append(list(zip(temp_structure_combined, potential_values)))

# Convert the processed data into a DataFrame
processed_df = pd.DataFrame(processed_data)

# Save the result to a new Excel file
processed_df.to_excel('processed_data.xlsx', index=False, header=False)
