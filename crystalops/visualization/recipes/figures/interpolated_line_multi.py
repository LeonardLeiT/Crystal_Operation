import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline

# 1. Data preparation
# Cu
data1 = np.array(
    [104, 305.4, 361.5, 902.1, 1040.5, 960.5, 1259.4, 1264.7, 1306.4, 1315.1, 1329]   
)
# Ni
data2 = np.array(
    [360.7, 822.3, 1086.4, 867.2, 1590, 1600, 1644.7, 1700, 1710.3, 1720.6, 1730.9]
)
# CuTa
data3 = np.array(
    [105.1, 559.2, 463.8, 226.6, 887.9, 949.8, 947.8, 1226, 1265, 1289.4, 1288.2]
)
# NiCoCr
data4 = np.array(
    [204.8, 458.3, 1034.5, 1280.1, 1200.4, 1330.3, 1354.2, 1380.5, 1394.2, 1408, 1409.3]
)
# FeNiCrCoCu
data5 = np.array(
    [211.4, 667.7, 1785, 1788.4, 1875.4, 1788.1, 1930.4, 1941.4, 2007.8, 2029, 2050.2]
)

# Put all datasets in a list
all_data = [data1/1358, data2/1728, data3/1310, data4/1410, data5/2260]

x_values = np.array(
    [1.65, 2.07, 2.37, 2.62, 2.82, 3.08, 3.44, 4.13, 5.48, 6.16, 7.22]
)

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
markersize = 100                # marker: the data point size

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
label = ['Unary (Cu)', 'Unary (Ni)', 'Binary (CuTa)', 'Ternary (NiCoCr)', "Quinary (FeNiCrCoCu)"]


## 2.5 Plot Figure
for i, data in enumerate(all_data):
    # mean_vals = np.mean(data, axis=1)
    # stderr_vals = np.std(data, axis=1, ddof=1) / np.sqrt(data.shape[1])
    # plt.errorbar(
    #     x_values[i % len(x_values)], mean_vals, yerr=stderr_vals,
    #     fmt=markers[i % len(markers)],
    #     color=rgb[i % len(rgb)], ecolor=rgb[i % len(rgb)],
    #     markerfacecolor=markerfacecolor, markeredgecolor=rgb[i % len(rgb)],
    #     elinewidth=elinewidth, capsize=capsize,
    #     markersize=markersize, linewidth=linewidth, alpha = alpha, zorder=7
    # )
    x_smooth = np.linspace(x_values.min(), x_values.max(), 200)
    spl = make_interp_spline(x_values, data, k=2)
    y_smooth = spl(x_smooth)
    plt.plot(x_smooth, y_smooth, linewidth=linewidth, color=rgb[i % len(rgb)], linestyle=linestyles[i % len(linestyles)], zorder=3)
    plt.scatter(x_values, data, label=label[i], 
                marker=markers[i % len(markers)], s=markersize,
                facecolor=(1, 1, 1, 0.5), edgecolor=rgb[i % len(rgb)], 
                linewidth=linewidth, zorder=5)
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
ax.set_xlim(1.5, 6.5)
ax.set_ylim(0, 1.1)
ax.legend(frameon=False)


plt.tight_layout()
plt.savefig("thermal_stablity.svg", bbox_inches='tight')
plt.savefig("thermal_stablity.pdf", bbox_inches='tight')
plt.show()
