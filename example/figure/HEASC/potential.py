import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 1. Data preparation
data = {}
data[0] = np.loadtxt(f'data/enthalpy-6966-500.txt')
data[1] = np.loadtxt(f'data/enthalpy-6966-600.txt')
# data[2] = np.loadtxt(f'cube/data/strain_Para_Z_0.txt')

# 2. Plot Settings
plt.figure(figsize=(6, 5))

## 2.1 Color Setting
rgb = [(102, 168, 69), (85, 103, 117), (75, 101, 175), (128, 198, 109)]
rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
markers = ['o', 's', '^', 'D', '*']
linestyles = ['-', '--', '-.', ':']

## 2.2 line Size Setting
linewidth = 2.0
marker_size = 200

## 2.3 font Setting 
ticksize = 20
labelsize = 20
legendsize = 20
plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.labelsize": labelsize,
    "xtick.labelsize": ticksize,
    "ytick.labelsize": ticksize,
    "legend.fontsize": legendsize,
})

## 2.4 labels setting
plt.xlabel('Time (ps)')
plt.ylabel('Pe (ev/atom)')
label = ["500", "600"]

for i in range(len(data)):
    time = data[i][:, 0]
    pe = data[i][:, 1]
    plt.plot(time, pe, label=label[i], linewidth=linewidth, color=rgb[i % len(rgb)])

# 2.5 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)
    
major_x_spacing = 10
major_y_spacing = 0.01
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.set_xlim(1, 50)
ax.set_ylim(-4.28, -4.23)
ax.legend(frameon=False)

plt.savefig('MDMC_pe.svg', format='svg', dpi=300, bbox_inches='tight')
plt.savefig('MDMC_pe.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()