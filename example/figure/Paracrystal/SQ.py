import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# 1. Data preparation
data = {}
# data[0] = np.loadtxt(f'Xtal_SQ.txt', encoding='utf-8')
data[0] = np.loadtxt(f'Equli_300_0_SQ.txt', encoding='utf-8')
data[1] = np.loadtxt(f'min_700_0_SQ.txt', encoding='utf-8')
data[2] = np.loadtxt(f'amorphous-1_SQ.txt', encoding='utf-8')
# data[0] = np.loadtxt(f'single_SQ.txt', encoding='utf-8')
# data[1] = np.loadtxt(f'polycrystal-2_SQ.txt', encoding='utf-8')
# data[2] = np.loadtxt(f'polycrystal-5_SQ.txt', encoding='utf-8')
# data[3] = np.loadtxt(f'polycrystal-8_SQ.txt', encoding='utf-8')

# 2. Plot Settings
plt.figure(figsize=(6.5, 5))

## 2.1 Color Setting
rgb = [(85, 103, 117), (232, 68, 69), (75, 101, 175), (128, 198, 109), (136, 135, 203), (251, 132, 2)]
rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
markers = ['o', 's', '^', 'D', '*']
linestyles = ['--', '-', '-.', ':']

# rgb = [(232, 68, 69), (75, 101, 175), (128, 198, 109), (136, 135, 203), (251, 132, 2)]
# rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
# linestyles = ['-', '-.', ':']

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
plt.xlabel("Q (Å⁻¹)")
plt.ylabel("S (Q)")
label = ["Xtal", "Paracrystalline", "Glassy"]
# label = ["Paracrystalline", "Glassy"]

for i in range(len(data)):
    r = data[i][:, 0]
    g_r = data[i][:, 1]
    plt.plot(r, g_r, label=label[i], linewidth=linewidth, color=rgb[i % len(rgb)], linestyle=linestyles[i % len(linestyles)])

# 2.5 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)
    
major_x_spacing = 3
major_y_spacing = 0.5
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.set_xlim(0.5, 15)
ax.set_ylim(0, 2.5)
ax.legend(frameon=False)

plt.savefig('SQ-all.svg', format='svg', dpi=300, bbox_inches='tight')
plt.savefig('SQ-all.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()