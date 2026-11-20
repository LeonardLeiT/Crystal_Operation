import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 1. Data preparation
data1 = np.array([
    [868.8,1079.9,808],
    [977.1,1025.5,970.4],
    [1055,1028,951],
    [1060.6,1092.4,1000.1],
])

x_labels = [1.52, 1.55, 1.60, 1.633]
x_pos = np.arange(len(x_labels))

# Calculate mean and standard error across columns (per row)
mean_vals1 = np.mean(data1, axis=1)
stderr_vals1 = np.std(data1, axis=1, ddof=1) / np.sqrt(data1.shape[1])

# Plot Settings
# 1. Input Data
mean_values1 = mean_vals1
y_error1 = stderr_vals1

# 2. Plot error bars
plt.figure(figsize=(8, 5.6))

## 2.1 Color Setting
rgb = [(232, 68, 69), (85, 103, 117), (75, 101, 175), (128, 198, 109)]
rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
alpha = 0.8                     # alpha: the transparency of the data point

## 2.2 line Size Setting
linewidth = 2.0                   # line: the line width
capsize = 10                     # error bar: the error bar head size

## 2.3 font Setting 
ticksize = 30                   # ticks: the size of the tick labels
labelsize = 30                  # labels: the size of the label
plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.labelsize": labelsize,
    "xtick.labelsize": ticksize,
    "ytick.labelsize": ticksize,
})

## 2.4 labels setting
plt.xlabel(f'c')
plt.ylabel(r'$T_{\mathrm{m}}$ (K)')

# dx = np.diff(np.sort(x)).min()      
# bar_w = 0.6 * dx                      width=bar_w,      
bars = plt.bar(x_pos, mean_vals1, yerr=stderr_vals1, 
               capsize=capsize, alpha=alpha, 
               color=plt.cm.viridis(np.linspace(0.3,0.7,len(x_pos))), 
               edgecolor="black", linewidth=linewidth)

# 2.5 tick setting
framewidth = 2.0                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)
    
plt.xticks(x_pos, x_labels)

major_y_spacing = 100
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.tick_params(axis='both', which='major', length=8, width=2, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=2, direction='in')
ax.set_ylim(700, 1100)
ax.legend(frameon=False)

plt.savefig('ave_stderr.svg', format='svg', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()
