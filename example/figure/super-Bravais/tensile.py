import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 1. Data preparation
data = {}
axis = 'Z'
time=1
# data[0] = np.loadtxt(f'1-1/{time}/data/strain_{axis}_300_0.txt')
# data[1] = np.loadtxt(f'1-0.7/{time}/data/strain_{axis}_300_0.txt')
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
# label = ["BCC (Perfect SC)", "BCT-1.2", "BCT-1.3", "BCT-1.4", "BCT-1.6"]
label = ["BCT-0.7", "BCC (Perfect SC)", "BCT-1.2", "BCT-1.4", "BCT-1.6"]
for i in range(len(data)):
    strain = data[i][:, 0] * 100
    stress = data[i][:, 1]
    plt.plot(strain, stress, label=label[i], linewidth=linewidth, color=rgb[i % len(rgb)])

# 2.5 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)
    
major_x_spacing = 1
major_y_spacing = 0.5
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.set_xlim(0.0, 5.5)
ax.set_ylim(0.001, 3)
ax.legend(frameon=False)

plt.savefig('tensile_SC.svg', format='svg', dpi=300, bbox_inches='tight')
plt.savefig('tensile_SC.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()