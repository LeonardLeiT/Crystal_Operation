import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 1. Data preparation
data1 = np.array([
    [1306.4,1305.4,1310.8,1308.7],
    [1301.0,1310.5,1305.3,1306],
    [1303.3,1300.9,1296.5,1305.2],
    [1306.1,1298.8,1297.9,1301.0],
])

data2 = np.array([
    [1308.0,1310.0,1324.7,1306.2],
    [1302.5,1308.6,1318.1,1310],
    [1313.2,1308.1,1305.4,1305.4],
    [1297.3,1308.1,1301.1,1304.2],
])

x = np.array([1.05, 1.10, 1.15, 1.20])

# Put all datasets in a list
all_data = [data1, data2]

# Calculate mean and standard error across columns (per row)
mean_vals1 = np.mean(data1, axis=1)
stderr_vals1 = np.std(data1, axis=1, ddof=1) / np.sqrt(data1.shape[1])
mean_vals2 = np.mean(data2, axis=1)
stderr_vals2 = np.std(data2, axis=1, ddof=1) / np.sqrt(data2.shape[1])

# Plot Settings
# 1. Input Data
x_values = x
mean_values1 = mean_vals1
y_error1 = stderr_vals1
mean_values2 = mean_vals2
y_error2 = stderr_vals2

# 2. Plot error bars
plt.figure(figsize=(6, 5))
## 2.1 Color Setting
rgb = [(232, 68, 69), (85, 103, 117), (75, 101, 175), (128, 198, 109)]
rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
markers = ['o', 's', '^', 'D', '*']
linestyles = ['-', '--', '-.', ':']
markerfacecolor = 'white'       # color: the color in the data point
alpha = 0.9                       # alpha: the transparency of the data point
## 2.2 line Size Setting
linewidth = 2                   # line: the line width
elinewidth = 2                  # error bar: the error bar line width
capsize = 8                     # error bar: the error bar head size
markersize = 20                 # marker: the data point size
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
plt.xlabel(f'c')
plt.ylabel(r'$T_{\mathrm{m}}$ (K)')
label = ['b = 0.8', 'b = 0.9']

for i, data in enumerate(all_data):
    mean_vals = np.mean(data, axis=1)
    stderr_vals = np.std(data, axis=1, ddof=1) / np.sqrt(data.shape[1])
    plt.errorbar(
        x, mean_vals, yerr=stderr_vals,
        fmt=markers[i % len(markers)],
        color=rgb[i % len(rgb)], ecolor=rgb[i % len(rgb)],
        markerfacecolor=markerfacecolor, markeredgecolor=rgb[i % len(rgb)],
        elinewidth=elinewidth, capsize=capsize,
        markersize=markersize, linewidth=linewidth, alpha = alpha, zorder=5
    )
    plt.plot(x, mean_vals, linewidth=linewidth, color=rgb[i % len(rgb)], linestyle=linestyles[i % len(linestyles)], label=label[i], zorder=3)

# 2.5 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)

major_x_spacing = 0.05
major_y_spacing = 5
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.set_xlim(1.025, 1.225)
ax.set_ylim(1297.5, 1317.5)
ax.legend(frameon=False)

plt.savefig('ave_stderr.svg', format='svg', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()
