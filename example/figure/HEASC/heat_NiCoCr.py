import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline

# 1. Data preparation
# NiCoCr Random
data1 = np.array([
    [204.8, 458.3, 1034.5, 1280.1, 1200.4, 1330.3, 1354.2],
    [207, 406, 1132, 1107.1, 1129.5, 1306, 1344.3],
    [314.7, 514.5, 1096.1, 947.4, 924.6, 1279.3, 1339]
])
# NiCoCr MC
data2 = np.array([
    [424.292, 650.637, 1215.12, 1357.48, 1325.49, 1361.73, 1372.07],
    [372.474, 796.959, 1198.47, 1315.28, 1324.58, 1365.65, 1365.24],
    [341.1, 1036, 1236, 1310.3, 1337.3, 1367.5, 1366.9]
])

# Put all datasets in a list
all_data = [data1/1390, data2/1390]

x_values = np.array([
    [1.66, 2.04, 2.33, 2.57, 2.78, 3.03, 3.38],
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
labels = ['NiCoCr-Random', 'NiCoCr-MDMC']

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
ax.set_xlim(1.5, 3.8)
ax.set_ylim(0, 1.1)
ax.legend(frameon=False)


plt.tight_layout()
plt.savefig("NiCoCr.svg", bbox_inches='tight')
plt.savefig("NiCoCr.pdf", bbox_inches='tight')
plt.show()
