import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline

# 1. Data preparation
# Cu
data1 = np.array([
    [305.4, 361.5, 902.1, 1040.5, 960.5, 1259.4, 1264.7],
    [139, 758.4, 769.3, 989.8, 1199, 1240, 1273.1],
    [92, 818.5, 998.2, 962.1, 1098.3, 1248.3, 1285.1]
])
# CuTa 0.1
data2 = np.array([
    [203, 683.1, 913.4, 1048.4, 1155.6, 1208.1, 1262.4],
    [351, 698.6, 770.9, 1087.4, 991.1, 1220.6, 1264.1],
    [130.6, 782.2, 810.1, 1052.4, 914.7, 1239.5, 1280.2]
])
# CuTa 0.5
data3 = np.array([
    [193.558, 617.235, 1026.7, 1102.3, 1117.2, 1174.2, 1244.1],
    [381.769, 809.512, 829, 1145.4, 1157, 1233.2, 1239.6],
    [135.346, 899.764, 1000, 1029.9, 1152.8, 1251.1, 1251.2]
])
# CuTa 1.0
data4 = np.array([
    [478.7, 964.1, 1024.1, 1199.3, 1195.1, 1221.9, 1227],
    [659.868, 999.9, 1026.1, 1106.7, 1189.4, 1225.7, 1248],
    [383.1, 912.69, 1023.86, 1117.5, 1157.9, 1213.4, 1238.7]
])

# CuTa 2.0
data5 = np.array([
    [433.816, 989.342, 931.925, 956.053, 1101.04, 1220.05, 1210.37],
    [683.087, 763.684, 1122.8, 1087.59, 1218.54, 1208.08, 1197.91],
    [417.618, 818.5, 1063.833, 1105.33, 1141.18, 1188.65, 1204.11]
])
# Put all datasets in a list
all_data = [data1/1358, data2/1353.2, data3/1339, data4/1310, data5/1282]

x_values = np.array([
    [2.07, 2.37, 2.62, 2.82, 3.08, 3.44, 4.13],
])

# 2. Plot Settings
plt.figure(figsize=(6.5, 5))

## 2.1 Color Setting
rgb = [(232, 68, 69), (85, 103, 117), (75, 101, 175), (36, 135, 103), (128, 198, 109)]
# rgb = [(38, 70, 83), (42, 157, 142), (233, 196, 107), (243, 162, 97), (230, 111, 81)]
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
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman", 
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
    "axes.labelsize": labelsize,
    "xtick.labelsize": ticksize,
    "ytick.labelsize": ticksize,
    "legend.fontsize": legendsize,
})

## 2.4 labels setting
plt.ylabel(r'$T/T_{\mathrm{m}}$')
plt.xlabel(r'$D_S$ (nm)')
labels = ['Cu', 'Cu-0.1at.%Ta', 'Cu-0.5at.%Ta', 'Cu-1.0at.%Ta', 'Cu-2.0at.%Ta']

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
    plt.plot(x_smooth, y_smooth, linewidth=linewidth, color=rgb[i % len(rgb)], linestyle=linestyles[2 % len(linestyles)], zorder=3)
plt.axhline(y = 1, color=(115/255, 115/255, 115/255), linestyle=':', linewidth=linewidth, alpha = alpha)

# 2.6 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)

major_x_spacing = 0.5
major_y_spacing = 0.2
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.set_xlim(1.9, 4.2)
ax.set_ylim(0, 1.05)
ax.legend(frameon=False)


plt.tight_layout()
plt.savefig("Cu_CuTa.svg", bbox_inches='tight')
plt.savefig("Cu_CuTa.pdf", bbox_inches='tight')
plt.show()
