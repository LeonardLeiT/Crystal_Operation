import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 1. Data preparation
data1 = np.array([
    [890.8,	800.4,676.2],
    [1305,1290.873,	1300.956]
])
data2 = np.array([
    [1240.4,1067.8,	1107],
    [1314,	1316,	1286]
])

datasets = [data1, data2]

x_labels = ['(1,0.9,1.1)', '(1,1.2,1.1)']
x_pos = np.arange(len(x_labels))


# Plot Settings
# 1. Plot error bars
plt.figure(figsize=(9, 6))

## 1.1 Color Setting
rgb = [(232, 68, 69), (85, 103, 117), (75, 101, 175), (128, 198, 109)]
rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
alpha = 0.8                     # alpha: the transparency of the data point

## 1.2 line Size Setting
linewidth = 2.0                   # line: the line width
capsize = 10                     # error bar: the error bar head size

## 1.3 font Setting 
ticksize = 30                   # ticks: the size of the tick labels
labelsize = 30                  # labels: the size of the label
legendsize = 30                 # legends: the size of the legend
plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.labelsize": labelsize,
    "xtick.labelsize": ticksize,
    "ytick.labelsize": ticksize,
    "legend.fontsize": legendsize,
})

## 2.4 labels setting
plt.xlabel(f'(a, b, c)')
plt.ylabel(r'$T_{\mathrm{m}}$ (K)')

label = ['60', '65']

width = 0.35   # bar width
offsets = np.linspace(-width/2, width/2, len(datasets))

bars = []
for i, data in enumerate(datasets):
    mean_vals = np.mean(data, axis=1)
    stderr_vals = np.std(data, axis=1, ddof=1) / np.sqrt(data.shape[1])

    b = plt.bar(
        x_pos + offsets[i], mean_vals, width=width/len(datasets)*2,
        yerr=stderr_vals, capsize=capsize, alpha=alpha,
        color=rgb[i], edgecolor="black", linewidth=linewidth,
        label=label[i]
    )
    bars.append(b)

# 2.5 tick setting
framewidth = 2.0                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)
    
plt.xticks(x_pos, x_labels)

major_y_spacing = 200
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.tick_params(axis='both', which='major', length=8, width=2, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=2, direction='in')
ax.set_ylim(600, 1400)
ax.legend(frameon=False)

plt.savefig('ave_stderr.svg', format='svg', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()
