import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import interpolate

# 1. Data Setting
Temp = 800
size = 5344
# 2. Plot Settings
plt.figure(figsize=(6.5, 5))

## 2.1 Color Setting
rgb = [(181, 206, 78), (151, 208, 197), (247, 172, 83), (145, 204, 192), (236, 110, 102), 
       (189, 119, 149), (124, 121, 121), (150, 59, 121), (127, 171, 209), (243, 152, 101),
       (82, 170, 220), (199, 193, 222), (238, 182 ,212), (200, 151, 54), (45, 136, 117)]
rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
markers = ['o', 's', '^', 'D', '*']
linestyles = [
    '-', 
    '--',  
    '-.',
    ':',
    (0, (3, 1, 1, 1, 1, 1)),
    (0, (5, 5)), 
    (0, (5, 1)), 
    (0, (3, 1, 1, 1)), 
    (0, (1, 1)),  
    (0, (5, 2, 1, 2)), 
]

## 2.2 line Size Setting
linewidth = 2.5
marker_size = 200

## 2.3 font Setting 
ticksize = 20
labelsize = 20
legendsize = 20
fontsize = 20
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
plt.xlabel('Time (ps)')
plt.ylabel(r'$\mathrm{MSD}~(\mathrm{\AA^2})$')
label = ["Fe", 'Ni', "Cr", "Co", "Cu"]

data = {}
# time = np.linspace(0, 1000, 11)
# print(time)

data[0] = np.loadtxt(f'Random/{size}/data/msd_{size}_{Temp}.txt')
data[1] = np.loadtxt(f'MCMD/{size}/data/msd_{size}_{Temp}.txt')
for i in range(len(label)):
    msd = data[0][:, i+2]
    time = data[0][:, 0] / 1000
    x_smooth = np.linspace(time.min(), time.max(), 50)
    tck = interpolate.splrep(time, msd, s=2)
    y_smooth = interpolate.splev(time, tck, der=0)
    plt.plot(time, y_smooth-y_smooth[0], linewidth=linewidth, color=rgb[i+2], linestyle = linestyles[i])
    # plt.plot(time, msd-msd[1], linewidth=linewidth, color=rgb[i], linestyle = linestyles[i])
    msd = data[1][:, i+2]
    time = data[1][:, 0] / 1000
    x_smooth = np.linspace(time.min(), time.max(), 100)
    tck = interpolate.splrep(time, msd, s=1)
    y_smooth = interpolate.splev(time, tck, der=0)
    plt.plot(time, y_smooth-y_smooth[0], linewidth=linewidth, color=rgb[i+2], linestyle = linestyles[i], label=label[i])
    # plt.plot(time, msd-msd[1], linewidth=linewidth, color=rgb[i], linestyle = linestyles[i], label=label[i])


plt.text(830, 0.9, 'Random', fontsize=fontsize)
plt.text(830, 0.25, 'MDMC', fontsize=fontsize)   
# plt.text(150, 0.5, '700 K', fontsize=fontsize)
# plt.text(150, 1.3, '900 K', fontsize=fontsize) 
# plt.text(50, 0.83, 'Atoms = 5344', fontsize=fontsize)  
# plt.text(50, 0.71, r'$\mathrm{D_S}$ = 2.78 nm', fontsize=fontsize)  

# 2.5 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)
    
major_x_spacing = 200
major_y_spacing = 0.3
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
# ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
# ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
ax.set_xlim(0, 1050)
ax.set_ylim(0.0, 1.35)
ax.legend(frameon=False)

plt.savefig('MSD-FeNiCrCoCu-element.svg', format='svg', dpi=300, bbox_inches='tight')
plt.savefig('MSD-FeNiCrCoCu-element.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()