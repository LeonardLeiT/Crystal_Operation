import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline

# 1. Data preparation
# Cu
data1 = np.array([
    [104, 305.4, 361.5, 902.1, 1040.5, 960.5, 1259.4, 1264.7, 1306.4, 1315.1],
    [91.8, 139, 758.4, 769.3, 989.8, 1199, 1240, 1273.1, 1296.3, 1317.3],
    [99.3, 92, 818.5, 998.2, 962.1, 1098.3, 1248.3, 1285.1, 1308.3, 1314.7]
])
# Ni
data2 = np.array([
    [360.7, 822.3, 1086.4, 867.2, 1564, 1660.6, 1665.2, 1700, 1725.8, 1745.4],
    [162.5, 1024.9, 1271.9, 1056.9, 810.3, 1460.9, 1655.6, 1692.7, 1721.6, 1746.2]
])

# NiCoCr
data3 = np.array([
    [204.8, 458.3, 1034.5, 1280.1, 1200.4, 1330.3, 1354.2, 1380.5, 1394.2, 1408],
    [207, 406, 1132, 907.1, 1129.5, 1306, 1344.3, 1345.2, 1390.8, 1404.4],
    [314.7, 514.5, 1096.1, 947.4, 924.6, 1279.3, 1339, 1363.1, 1391.4, 1403.6]
])

# FeNiCrCoCu
data4 = np.array([
    [211.4, 667.7, 1785, 1788.4, 1875.4, 1788.1, 1930.4, 1941.4, 2007.8, 2029],
    [267.4, 229.6, 1519, 1769.8, 1767.7, 1918.4, 1861.1, 1931.7, 2011.8, 2024.7],
    [175.6, 845.4, 1364.4, 1841.5, 1562.9, 1924.5, 1900.3, 1968.1, 2005.1, 2030.7],
])

# Put all datasets in a list
all_data = [data1/1358, data2/1748, data3/1410, data4/2200]
all_data = [data1/1358, data2/1748]

x_values = np.array([
    [1.65, 2.07, 2.37, 2.62, 2.82, 3.08, 3.44, 4.13, 5.48, 6.16],
    [1.64, 2.00, 2.30, 2.53, 2.74, 2.99, 3.33, 4.00, 5.31, 5.97],
    [1.66, 2.04, 2.33, 2.57, 2.78, 3.03, 3.38, 4.06, 5.39, 6.06],
    [1.67, 2.04, 2.34, 2.58, 2.78, 3.04, 3.39, 4.06, 5.39, 6.06]
])

# 2. Plot Settings
plt.figure(figsize=(6.5, 5))

## 2.1 Color Setting
rgb = [(232, 68, 69), (85, 103, 117), (75, 101, 175), (128, 198, 109), (136, 135, 203)]
rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
markers = ['o', 's', '^', 'D', '*']
linestyles = ['-', '--', '-.', ':']
markerfacecolor = 'white'       # color: the color in the data point
alpha = 0.9                     # alpha: the transparency of the data point

## 2.2 line Size Setting
linewidth = 2                   # line: the line width
elinewidth = 2                  # error bar: the error bar line width
capsize = 10                    # error bar: the error bar head size
markersize = 13                 # marker: the data point size

## 2.3 font Setting 
ticksize = 20                   # ticks: the size of the tick labels
labelsize = 20                  # labels: the size of the label
legendsize = 20                 # legends: the size of the legend
plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.labelsize": labelsize,
    "xtick.labelsize": ticksize,
    "ytick.labelsize": ticksize,
    "legend.fontsize": legendsize,
})

## 2.4 labels setting
plt.ylabel(r'$T/T_{\mathrm{m}}$')
plt.xlabel(r'$D_S$ (nm)')
labels = ['Unary (Cu)', 'Unary (Ni)', 'Ternary (NiCoCr)', "Quinary (FeNiCrCoCu)"]
labels = ['Unary (Cu)', 'Unary (Ni)']

## 2.5 Plot Figure
for i, data in enumerate(all_data):
    mean_vals = np.mean(data, axis=0)
    stderr_vals = np.std(data, axis=0, ddof=1) / np.sqrt(data.shape[0])
    plt.errorbar(
        x_values[i % len(x_values)], mean_vals, yerr=stderr_vals,
        fmt=markers[i % len(markers)],
        color=rgb[i % len(rgb)], ecolor=rgb[i % len(rgb)],
        markerfacecolor=markerfacecolor, markeredgecolor=rgb[i % len(rgb)],
        elinewidth=elinewidth, capsize=capsize, 
        markersize=markersize, linewidth=linewidth, alpha = alpha, zorder=7,
        label=labels[i]
    )
    x_smooth = np.linspace(x_values[i % len(x_values)].min(), x_values[i % len(x_values)].max(), 200)
    spl = make_interp_spline(x_values[i % len(x_values)], mean_vals, k=1)
    y_smooth = spl(x_smooth)
    plt.plot(x_smooth, y_smooth, linewidth=linewidth, color=rgb[i % len(rgb)], linestyle=linestyles[i % len(linestyles)], zorder=3)
plt.axhline(y = 1, color=(115/255, 115/255, 115/255), linestyle=':', linewidth=linewidth, alpha = alpha)

# 2.6 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)

major_x_spacing = 1
major_y_spacing = 0.2
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.set_xlim(1.3, 6.5)
ax.set_ylim(0, 1.1)
ax.legend(frameon=True)


plt.tight_layout()
plt.savefig("thermal_stablity_Unary.svg", bbox_inches='tight')
plt.savefig("thermal_stablity_Unary.pdf", bbox_inches='tight')
plt.show()
