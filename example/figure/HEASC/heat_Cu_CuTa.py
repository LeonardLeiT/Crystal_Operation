import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline

# 1. Data preparation
# Cu
data1 = np.array([
    [305.4, 361.5, 902.1, 1040.5, 960.5, 1259.4],
    [139, 758.4, 769.3, 989.8, 1199, 1240],
    [92, 818.5, 998.2, 962.1, 1098.3, 1248.3]
])
# CuTa
data2 = np.array([
    [473.5, 953.1, 1013.1, 1187.5, 1188, 1127.4],
    [659.868, 980.457, 1015.68, 1096.04, 1178.15, 1223.41],
    [373.094, 912.69, 1023.86, 1111.65, 1145.65, 1208.04]
])

# Put all datasets in a list
all_data = [data1/1358, data2/1310]

x_values = np.array([
    [2.07, 2.37, 2.62, 2.82, 3.08, 3.44],
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
labels = ['Cu', 'Cu-1.0at.%Ta']

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

major_x_spacing = 0.5
major_y_spacing = 0.2
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.set_xlim(1.8, 3.8)
ax.set_ylim(0, 1.1)
ax.legend(frameon=False)


plt.tight_layout()
plt.savefig("Cu_CuTa.svg", bbox_inches='tight')
plt.savefig("Cu_CuTa.pdf", bbox_inches='tight')
plt.show()
