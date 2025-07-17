import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.cm as cm

# 设置中文字体
plt.rcParams["font.family"] = ["SimHei", "WenQuanYi Micro Hei", "Heiti TC"]
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# Define parameters
a = 1
b_values = np.linspace(0.4, 1.6, 50)  # y-axis
c_values = np.linspace(0.4, 1.6, 50)  # x-axis

# Create grid
B, C = np.meshgrid(b_values, c_values)

# Calculate ratios
ratio_x = C / a  # c/a
ratio_y = B / a  # b/a

# Initialize classification array (0: undefined, 1: BCC, 2: BCT, 3: ICO)
classification = np.zeros_like(ratio_x, dtype=int)
ico_product = np.zeros_like(ratio_x)  # 存储ICO构型的计算值

# Classify each point
for i in range(ratio_x.shape[0]):
    for j in range(ratio_x.shape[1]):
        rx = ratio_x[i, j]
        ry = ratio_y[i, j]
        
        if np.isclose(rx, 1, atol=1e-6) and np.isclose(ry, 1, atol=1e-6):
            classification[i, j] = 1  # BCC
        elif (np.isclose(ry, 1, atol=1e-6) and not np.isclose(rx, 1, atol=1e-6)) or \
             (np.isclose(rx, 1, atol=1e-6) and not np.isclose(ry, 1, atol=1e-6)):
            classification[i, j] = 2  # BCT
        elif not np.isclose(rx, ry, atol=1e-6) and not np.isclose(rx, 1, atol=1e-6) and not np.isclose(ry, 1, atol=1e-6):
            classification[i, j] = 3  # ICO
            
            # 计算三个值的最小值
            min_val = min(rx, ry, a)
            
            # 每个值除以最小值
            normalized_rx = rx / min_val
            normalized_ry = ry / min_val
            normalized_a = a / min_val
            
            # 取最大的两个值的乘积
            sorted_vals = sorted([normalized_rx, normalized_ry, normalized_a], reverse=True)
            ico_product[i, j] = sorted_vals[0] * sorted_vals[1]  # 最大的两个值的乘积

# 定义颜色映射
base_colors = ['white', '#FF9999', '#99CCFF']  # 白色(未定义), 红色(BCC), 蓝色(BCT)

# 创建ICO构型的绿色渐变
min_product = np.min(ico_product[classification == 3])
max_product = np.max(ico_product[classification == 3])
norm = plt.Normalize(min_product, max_product)
green_cmap = cm.Greens  # 使用绿色渐变

# 准备绘图数据
x_flat = ratio_x.flatten()
y_flat = ratio_y.flatten()
class_flat = classification.flatten()
product_flat = ico_product.flatten()

# 绘制构型分类图
plt.figure(figsize=(12, 10))

# 绘制未定义、BCC和BCT构型
for class_id, color in enumerate(base_colors, start=0):
    mask = class_flat == class_id
    if np.sum(mask) > 0:
        plt.scatter(x_flat[mask], y_flat[mask], c=color, s=100, edgecolors='k', alpha=0.8, 
                   label=['未定义', 'BCC', 'BCT'][class_id])

# 绘制ICO构型，颜色根据计算值变化
mask_ico = class_flat == 3
if np.sum(mask_ico) > 0:
    scatter_ico = plt.scatter(x_flat[mask_ico], y_flat[mask_ico], c=product_flat[mask_ico], 
                             cmap=green_cmap, norm=norm, s=100, edgecolors='k', alpha=0.8,
                             label='ICO (颜色表示 归一化后最大两值乘积)')
    # 添加ICO颜色条
    cbar_ico = plt.colorbar(scatter_ico)
    cbar_ico.set_label('归一化后最大两值乘积', fontsize=12)

# 添加网格线和参考线
plt.grid(True, linestyle='--', alpha=0.6)
plt.axhline(y=1, color='k', linestyle='-', alpha=0.5)
plt.axvline(x=1, color='k', linestyle='-', alpha=0.5)
plt.plot([0.4, 1.6], [0.4, 1.6], 'k--', alpha=0.5)

# 添加标签和标题
plt.xlabel('c/a 比值', fontsize=14)
plt.ylabel('b/a 比值', fontsize=14)
plt.title('基于 b/a 和 c/a 比值的晶体构型分类', fontsize=16)
plt.xlim(0.35, 1.65)
plt.ylim(0.35, 1.65)
# plt.legend(loc='upper right', fontsize=12)

plt.tight_layout()
plt.show()    