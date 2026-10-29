import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline

# 1. Data preparation
# FeNiCrCoCu Random
data1 = np.array([
    [211.4, 667.7, 1785, 1788.4, 1875.4, 1788.1, 1930.4, 1941.4],
    [267.4, 229.6, 1519, 1769.8, 1767.7, 1918.4, 1861.1, 1931.7],
    [175.6, 845.4, 1364.4, 1841.5, 1562.9, 1924.5, 1900.3, 1968.1]
])

# FeNiCrCoCu MC
data2 = np.array([
    [1015.55, 1291.7, 1763.02, 1876.75, 1881.33, 1903.89,1951.93, 1986.2],
    [1112.1, 1202.49, 1646.86, 1845.36, 1774.94, 1950.39, 1964.65, 1942.4],
    [950.4, 1311.2, 1778.6, 1848, 1932.5, 1799, 1969.6, 1976.1]
])

# Put all datasets in a list
all_data = [data1/2120, data2/2120]

x_values = np.array([
    [1.67, 2.04, 2.34, 2.58, 2.78, 3.04, 3.39, 4.06],
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
plt.ylabel(r'$\mathrm{T / T_m}$')
plt.xlabel(r'$\mathrm{D_S}$ (nm)')
labels = ['FeNiCrCoCu-Random', 'FeNiCrCoCu-MDMC']

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
ax.set_xlim(1.51, 4.2)
ax.set_ylim(0, 1.1)
ax.legend(frameon=False)


plt.tight_layout()
plt.savefig("heat-FeNiCrCoCu.svg", bbox_inches='tight')
plt.savefig("heat-FeNiCrCoCu.pdf", bbox_inches='tight')
plt.show()
