import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 1. Data preparation
data1 = np.array([
    [5.345322, 5.246334, 5.840259],
    [4.685405, 5.411313, 4.619414],
    [2.870636, 3.068611, 2.969623],
    [1.029469, 1.544204, 1.425419],
])

data2 = np.array([
    [2.078736, 2.573673, 2.52418],
    [1.821369, 2.415293, 2.336104],
    [1.346229, 1.227444, 1.227444],
    [-2.454888, -1.583799, -2.256914],
])

data3 = np.array([
    [2.870636, 3.217092, 1.930255],
    [2.276711, 2.326205, 1.682786],
    [1.504609, 1.583799, 1.029469],
    [0.39595, 0.47514, 0.23757],
])

x = np.array([0.6, 0.7, 0.8, 0.9])

# Put all datasets in a list
all_data = [data1, data2, data3]

# Calculate mean and standard error across columns (per row)
mean_vals1 = np.mean(data1, axis=1)
stderr_vals1 = np.std(data1, axis=1, ddof=1) / np.sqrt(data1.shape[1])
mean_vals2 = np.mean(data2, axis=1)
stderr_vals2 = np.std(data2, axis=1, ddof=1) / np.sqrt(data2.shape[1])
mean_vals3 = np.mean(data3, axis=1)
stderr_vals3 = np.std(data3, axis=1, ddof=1) / np.sqrt(data2.shape[1])

# Plot Settings
# 1. Input Data
x_values = x
mean_values1 = mean_vals1
y_error1 = stderr_vals1
mean_values2 = mean_vals2
y_error2 = stderr_vals2
mean_values3 = mean_vals3
y_error3 = stderr_vals3

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
capsize = 10                     # error bar: the error bar head size
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
plt.xlabel(r'$T/T_{\mathrm{m}}$')
plt.ylabel(r'$U$ (m/s)')
label = ['Cube', 'Kelvin', "Schwarz"]

for i, data in enumerate(all_data):
    mean_vals = np.mean(data, axis=1)
    stderr_vals = np.std(data, axis=1, ddof=1) / np.sqrt(data.shape[1])
    plt.errorbar(
        x, mean_vals, yerr=stderr_vals,
        fmt=markers[i % len(markers)],
        color=rgb[i % len(rgb)], ecolor=rgb[i % len(rgb)],
        markerfacecolor=markerfacecolor, markeredgecolor=rgb[i % len(rgb)],
        elinewidth=elinewidth, capsize=capsize,
        markersize=markersize, linewidth=linewidth, alpha = alpha, zorder=7
    )
    plt.plot(x, mean_vals, linewidth=linewidth, color=rgb[i % len(rgb)], linestyle=linestyles[i % len(linestyles)], label=label[i], zorder=3)
plt.axhline(y = 0, color=(115/255, 115/255, 115/255), linestyle=':', linewidth=linewidth, alpha = alpha)

# 2.5 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)

major_x_spacing = 0.1
major_y_spacing = 2.0
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.set_xlim(0.55, 0.95)
ax.set_ylim(-3, 6)
ax.legend(frameon=True)

plt.savefig('growth.svg', format='svg', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()
