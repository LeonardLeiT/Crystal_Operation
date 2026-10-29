import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 1. Data preparation
data = {}

data[0] = np.loadtxt(f'glass_300_2.751_1/data/strain_300_2.751_1.txt')
data[1] = np.loadtxt(f'para_300_3.346_5_1/data/strain_300_3.346_5.txt')
data[2] = np.loadtxt(f'para_300_3.346_5_2/data/strain_300_3.346_6.txt')
data[3] = np.loadtxt(f'para_300_3.346_5_3/data/strain_300_3.346_6.txt')
data[4] = np.loadtxt(f'para_300_3.346_5_4/data/strain_300_3.346_6.txt')
data[5] = np.loadtxt(f'para_300_3.346_5_5/data/strain_300_3.346_6.txt')

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
label = ["Glassy", "Para-1", "Para-2", "Para-3", "Para-4", "Para-5"]

for i in range(len(data)):
    strain = data[i][:, 0]
    stress = data[i][:, 1]
    strain = strain - strain.min()
    stress = stress - stress.min()
    print(f"Strain: {strain[502]}, Stress: {stress[502]}")
    imax = np.argmax(stress)
    print(f"Strain: {strain[imax]}, Index: {imax}")
    print(f"Stress range: {stress.min():.4f} to {stress.max():.4f}")
    # ---- Select elastic region cutoff: strain < 0.20 ----
    mask = strain < 0.01
    strain_elastic = strain[mask]
    stress_elastic = stress[mask]
    # print(len(strain_elastic), len(stress_elastic))
    # Linear fit y = kx (Young's modulus)
    k, _ = np.polyfit(strain_elastic, stress_elastic, 1)
    print(f"{label[i]}: Young's modulus E = {k:.3f} GPa")

    # Fitted line for plotting
    strain_fit = np.linspace(0, 0.1, 200)
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
major_y_spacing = 2
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.set_xlim(0, 19)
ax.set_ylim(0.01, 8)
ax.legend(frameon=True)

plt.savefig('crack-stress.svg', format='svg', dpi=300, bbox_inches='tight')
plt.savefig('crack-stress.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()
