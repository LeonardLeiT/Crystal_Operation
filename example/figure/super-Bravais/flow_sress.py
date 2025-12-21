import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline

# 1. Data preparation
data = {}
axis = 'Y'
time=3
# data[0] = np.loadtxt(f'1-1/{time}/data/strain_{axis}_300_0.txt')
# data[1] = np.loadtxt(f'1.2-1/{time}/data/strain_{axis}_300_0.txt')
# data[2] = np.loadtxt(f'1-1.3/{time}/data/strain_{axis}_300_0.txt')
# data[3] = np.loadtxt(f'1.4-1/{time}/data/strain_{axis}_300_0.txt')
# data[4] = np.loadtxt(f'1-1.6/{time}/data/strain_{axis}_300_0.txt')

data[0] = np.loadtxt(f'1-0.7/1/data/strain_Z_300_0.txt')
data[1] = np.loadtxt(f'1-1/2/data/strain_Z_300_0.txt')
data[2] = np.loadtxt(f'1.2-1/1/data/strain_Z_300_0.txt')
data[3] = np.loadtxt(f'1.4-1/1/data/strain_X_300_0.txt')
data[4] = np.loadtxt(f'1-1.6/3/data/strain_X_300_0.txt')

# 2. Plot Settings
plt.figure(figsize=(6.5, 5))

## 2.1 Color Setting
rgb = [(85, 103, 117), (232, 68, 69), (75, 101, 175), (128, 198, 109), (136, 135, 203), (251, 132, 2)]
rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
markers = ['o', 's', '^', 'D', '*']
linestyles = ['-', '--', '-.', ':']
markerfacecolor = 'white'       # color: the color in the data point
alpha = 0.9                     # alpha: the transparency of the data point

## 2.2 line Size Setting
linewidth = 2                   # line: the line width
elinewidth = 2                  # error bar: the error bar line width
capsize = 10                    # error bar: the error bar head size
markersize = 15                 # marker: the data point size

## 2.3 font Setting 
ticksize = 20
labelsize = 20
legendsize = 20
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
plt.xlabel(r'BCT Z-axis $c$')
plt.ylabel('Flow stress (GPa)')
# label = ["BCC (Perfect SC)", "BCT-1.2", "BCT-1.3", "BCT-1.4", "BCT-1.6"]
labels = ["BCT-0.7", "BCC (Perfect SC)", "BCT-1.2", "BCT-1.4", "BCT-1.6"]
strain_ranges = [
    [2.2, 5.5],
    [3.6, 5.5],
    [2.9, 5.5],
    [2.4, 5.5],
    [2.4, 5.5],
         ]
x_values = np.array([0.7, 1, 1.2, 1.4, 1.6])
mean_values = []
for i in range(len(data)):
    strain = data[i][:, 0] * 100
    stress = data[i][:, 1]
    mask = (strain > strain_ranges[i][0]) & (strain < strain_ranges[i][1])
    filtered_stress = stress[mask]
    filtered_strain = strain[mask]
    mean_vals = np.mean(filtered_stress)
    mean_values.append(mean_vals)
    stderr_vals = np.std(filtered_stress, ddof=1) / np.sqrt(len(filtered_stress))
    stderr_vals = np.std(filtered_stress, ddof=1)
    plt.errorbar(
        x_values[i % len(x_values)], mean_vals, yerr=stderr_vals,
        fmt=markers[i % len(markers)],
        color=rgb[i % len(rgb)], ecolor=rgb[i % len(rgb)],
        markerfacecolor=markerfacecolor, markeredgecolor=rgb[i % len(rgb)],
        elinewidth=elinewidth, capsize=capsize, 
        markersize=markersize, linewidth=linewidth, alpha = alpha, zorder=7,
        label=labels[i], 
    )

x_smooth = np.linspace(x_values.min()-0.05, x_values.max()+0.05, 200)
spl = make_interp_spline(x_values, mean_values, k=1)
y_smooth = spl(x_smooth)
plt.plot(x_smooth, y_smooth, linewidth=linewidth, color='black', linestyle=linestyles[3 % len(linestyles)], zorder=3)


# for i, data in enumerate(all_data):
#     mean_vals = np.mean(data, axis=0)
#     stderr_vals = np.std(data, axis=0, ddof=1) / np.sqrt(data.shape[0])
#     plt.errorbar(
#         x_values[i % len(x_values)], mean_vals, yerr=stderr_vals,
#         fmt=markers[i % len(markers)],
#         color=rgb[i % len(rgb)], ecolor=rgb[i % len(rgb)],
#         markerfacecolor=markerfacecolor, markeredgecolor=rgb[i % len(rgb)],
#         elinewidth=elinewidth, capsize=capsize, 
#         markersize=markersize, linewidth=linewidth, alpha = alpha, zorder=7,
#         label=labels[i]
#     )

# plt.axhline(y = 1, color=(115/255, 115/255, 115/255), linestyle=':', linewidth=linewidth, alpha = alpha)


# 2.5 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)
    
major_x_spacing = 0.2
major_y_spacing = 0.2
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.set_xlim(0.6001, 1.7)
ax.set_ylim(1.6, 2.7)
ax.legend(frameon=False)

plt.savefig('flow_stress_SC.svg', format='svg', dpi=300, bbox_inches='tight')
plt.savefig('flow_stress_SC.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()