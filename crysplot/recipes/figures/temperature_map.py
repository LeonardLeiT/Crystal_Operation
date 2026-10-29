import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# 数据
x = [1, 1.2, 1.4]
y = [1, 1.3, 1.6]
z = np.array([
    [1319, 1304, 1296],
    [1305, 1309, 1300],
    [1302, 1305, 1302]
])

# 转为 DataFrame
df = pd.DataFrame(z, index=x, columns=y)

# 设置字体
plt.rcParams.update({
    "font.family": "Times New Roman",
    "font.size": 30,
    "axes.labelsize": 30,
    "axes.titlesize": 30,
    "xtick.labelsize": 30,
    "ytick.labelsize": 30,
    "legend.fontsize": 30,
})

fig, ax = plt.subplots(figsize=(12,6))

# 绘制热力图
sns.heatmap(
    df,
    annot=True,
    fmt=".0f",
    cmap="coolwarm",  # 也可尝试 "magma", "coolwarm"
    cbar_kws={'label': r'$T_{\mathrm{m}}$ (K)'},
    linewidths=0.5,
    linecolor='white',
    square=True,
    ax=ax
)

# 标签和标题
ax.set_xlabel("c")
ax.set_ylabel("b")

# 美化刻度
ax.set_xticklabels([f"{v:.1f}" for v in df.columns])
ax.set_yticklabels([f"{v:.1f}" for v in df.index], rotation=0)

plt.tight_layout()
fig.savefig("heatmap_temperature.svg", bbox_inches='tight')
plt.show()
