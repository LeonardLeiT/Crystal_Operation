import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

time = 1
# 读取数据 Potential Volume
# file_path = f'Volume-{time}.xlsx'
file_path = f'Potential-{time}.xlsx'
data = pd.read_excel(file_path)
# 计算数据的中间索引
# mid_index = len(data) // 10  * 1
# # 取后一半数据并重置索引
# data = data.iloc[mid_index:].reset_index(drop=True)
# 获取所有构型的列名
columns = data.columns
configs = list(set(col.split('_')[1] for col in columns))  # 提取唯一的构型编号
# configs.sort(key=lambda x: float(x))  # 按构型编号排序
configs.sort(key=lambda x: x)  # 按构型编号排序
# 定义平滑阶跃函数
def smooth_step(x, x0, a1, a2, k):
    return a1 + (a2 - a1) / (1 + np.exp(-k * (x - x0)))

# 存储拟合结果
fit_results = {}

# 遍历不同构型
for config in configs:
    temp_col = f"Temp_{config}"
    # potential_col = f"Vol_{config}"
    potential_col = f"Potential_{config}"
    if temp_col in data.columns and potential_col in data.columns:
        temp = data[temp_col].dropna().values
        potential = data[potential_col].dropna().values

        # 估计初始参数
        x0_guess = np.median(temp)  # 中位数作为初始阶跃点
        a1_guess, a2_guess = np.percentile(potential, [5, 95])  # 避免极端值影响
        k_guess = 1.0  # 斜率参数
        p0 = [x0_guess, a1_guess, a2_guess, k_guess]

        # 设置参数边界
        bounds = (
            [min(temp), min(potential), min(potential), 0],  # 下界
            [max(temp), max(potential), max(potential), 10]  # 上界
        )

        # 执行拟合
        try:
            popt, _ = curve_fit(smooth_step, temp, potential, p0=p0, bounds=bounds, maxfev=10000)
            x0_fit, a1_fit, a2_fit, k_fit = popt
            fit_results[config] = x0_fit  # 存储拟合结果

            # 生成拟合曲线
            temp_sorted = np.linspace(min(temp), max(temp), 500)
            potential_fit = smooth_step(temp_sorted, *popt)

            # 绘制拟合曲线
            plt.figure(figsize=(6, 4))
            plt.scatter(temp, potential, label=f'Data ({config})', color='blue', s=10)
            plt.plot(temp_sorted, potential_fit, label='Smooth Step Fit', color='red', linestyle='--')
            plt.axvline(x0_fit, color='green', linestyle=':', label=f'Step Center: {x0_fit:.3f}')

            plt.xlabel('Temperature (K)')
            plt.ylabel('Potential')
            plt.title(f'Fit for Config {config}')
            plt.legend()
            plt.grid()
            plt.show()

        except RuntimeError:
            print(f"拟合失败: 构型 {config}，请检查数据或调整初始参数。")

# 输出拟合结果
for config, x0_fit in fit_results.items():
    print(f"构型 {config} 拟合阶跃中心点: {x0_fit:.3f}")
