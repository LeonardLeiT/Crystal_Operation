import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 1. Data preparation
data = {}
data[0] = np.loadtxt(f'xtal/data/strain_Z_300_0.txt')
data[1] = np.loadtxt(f'amorphous-1/data/strain_Z_300_0.txt')
data[2] = np.loadtxt(f'polycrystal-10/data/strain_Z_300_0.txt')
data[3] = np.loadtxt(f'polycrystal-10-3.34-8/data/strain_Z_300_0.txt')
data[4] = np.loadtxt(f'polycrystal-10-3.34-8/data/strain_Z_300_8.3.txt')
data[5] = np.loadtxt(f'polycrystal-10-3.56-16/data/strain_Z_300_0.txt')
data[6] = np.loadtxt(f'polycrystal-10-3.56-16/data/strain_Z_300_15.9.txt')

# 2. Plot Settings
plt.figure(figsize=(6.5, 5))

## 2.1 Color Setting
rgb = [(232, 68, 69), (85, 103, 117), (75, 101, 175), (128, 198, 109), (136, 135, 203), (251, 132, 2), (254, 183, 5)]
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
label = ["Xtal", "Glassy", "Poly", "Para-3.3-0", "Para-3.3-8", "Para-3.5-0", "Para-3.5-16"]
# label = ["Paracrystalline-0", "Paracrystalline-15", "Paracrystalline"]

for i in range(len(data)):
    strain = data[i][:, 0] * -1
    stress = data[i][:, 1] * -1
    strain = strain - strain[0]
    stress = stress - stress[0]
    # print(len(strain), len(stress))
    # print(f"Strain range: {strain.min():.4f} to {strain.max():.4f}")
    # ---- Select elastic region cutoff: strain < 0.20 ----
    mask = strain < 0.03
    strain_elastic = strain[mask]
    stress_elastic = stress[mask]
    # print(len(strain_elastic), len(stress_elastic))
    # Linear fit y = kx (Young's modulus)
    k, _ = np.polyfit(strain_elastic, stress_elastic, 1)
    print(f"{label[i]}: Young's modulus E = {k:.3f} GPa")

    # Fitted line for plotting
    strain_fit = np.linspace(0, 0.03, 200)
    stress_fit = k * strain_fit

    # Raw curve
    plt.plot(strain * 100, stress, label=label[i], linewidth=linewidth, color=rgb[i % len(rgb)])
    # Fitted line
    plt.plot(strain_fit * 100, stress_fit,
             linestyle='--',
             linewidth=linewidth,
             color=rgb[i % len(rgb)])

# 2.5 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)
    
major_x_spacing = 3
major_y_spacing = 3
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.set_xlim(0, 18)
ax.set_ylim(0.01, 18)
# ax.set_ylim(-5, 14)
ax.legend(frameon=True)

plt.savefig('compress.svg', format='svg', dpi=300, bbox_inches='tight')
plt.savefig('compress.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()