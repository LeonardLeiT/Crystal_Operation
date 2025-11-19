import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 1. Data preparation
data = {}
data[0] = np.loadtxt(f'xtal/data/strain_XZ_300_0.txt')
data[1] = np.loadtxt(f'paracrystal-3/data/strain_XZ_300_0.txt')
# data[1] = np.loadtxt(f'polycrystal-5/data/strain_Z_300_0.txt')
data[2] = np.loadtxt(f'amorphous-1/data/strain_XZ_300_0.txt')

# 2. Plot Settings
plt.figure(figsize=(6.5, 5))

## 2.1 Color Setting
rgb = [(232, 68, 69), (85, 103, 117), (75, 101, 175), (128, 198, 109), (136, 135, 203), (251, 132, 2)]
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
plt.xlabel('Strain (%)')
plt.ylabel('Stress (GPa)')
label = ["Xtal", "Paracrystalline", "Glassy"]

for i in range(len(data)):
    strain = data[i][:, 0] * 100
    stress = data[i][:, 1]
    plt.plot(strain, stress, label=label[i], linewidth=linewidth, color=rgb[i % len(rgb)])

# 2.5 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)
    
major_x_spacing = 3
major_y_spacing = 1
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.set_xlim(0, 15)
ax.set_ylim(0.01, 6)
ax.legend(frameon=True)

plt.savefig('shear-XZ.svg', format='svg', dpi=300, bbox_inches='tight')
plt.savefig('shear-XZ.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()