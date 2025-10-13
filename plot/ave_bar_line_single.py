import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 1. Data preparation
data = np.array([
    [1319.272, 1320.147, 1316.615],
    [1310.018, 1312.036, 1304.236],
    [1313.706, 1304.91, 1297.621],
    [1308.799, 1293.043, 1313.741],
    [1304.409, 1300.728, 1306.58],
    [1303.15, 1298.74, 1304.826],
    [1299.857, 1300.229, 1300.875],
    [1297.557, 1308.423, 1298.782],
    [1297.15, 1292.74, 1296.826]
])

x = np.array([0, 0.186339, 0.2, 0.243252, 0.3, 0.319438, 0.372678, 0.4, 0.6])

# Calculate mean and standard error across columns (per row)
mean_vals = np.mean(data, axis=1)
stderr_vals = np.std(data, axis=1, ddof=1) / np.sqrt(data.shape[1])


# Plot Settings
# 1. Input Data
x_values = x
mean_values = mean_vals
y_error = stderr_vals

# 2. Plot error bars
plt.figure(figsize=(6, 6))
## 2.1 Color Setting
rgb_color = (88/255, 89/255, 91/255)
color = rgb_color               # color: data point and line
ecolor = rgb_color              # color: error bar
markerfacecolor = 'white'       # color: the color in the data point
markeredgecolor = rgb_color     # color: the color of the edge of the data point
alpha = 1                       # alpha: the transparency of the data point
## 2.2 line Size Setting
linewidth = 1.5                 # line: the line width
elinewidth = 1.5                # error bar: the error bar line width
capsize = 8                     # error bar: the error bar head size
markersize = 15                 # marker: the data point size
## 2.3 font Setting 
ticksize = 20                   # ticks: the size of the tick labels
labelsize = 20                  # labels: the size of the label
plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.labelsize": labelsize,
    "xtick.labelsize": ticksize,
    "ytick.labelsize": ticksize,
})
## 2.4 labels setting
plt.xlabel(f'$\Delta$ lattice')
plt.ylabel(r'$T_{\mathrm{m}}$ (K)')


plt.errorbar(
    x = x_values, y = mean_vals, yerr=stderr_vals,
    fmt='o-', 
    color=color, ecolor=ecolor, 
    markerfacecolor=markerfacecolor, markeredgecolor=markeredgecolor, 
    elinewidth=elinewidth, capsize=capsize,
    markersize=markersize, linewidth=linewidth
)

# 2.5 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)

ax.tick_params(axis='both', which='major',
               length=8, width=1.5, direction='in',
               top=False, right=False, bottom=True, left=True)

ax.tick_params(axis='both', which='minor',
               length=4, width=1.5, direction='in',
               top=False, right=False, bottom=True, left=True)
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())

plt.savefig('ave_stderr.svg', format='svg', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()
