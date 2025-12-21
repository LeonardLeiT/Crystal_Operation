import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import interpolate

# 1. Data Setting
Temp = [800, 1000]
size = 5344
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
label = ["FeNiCrCoCu-Random", 'FeNiCrCoCu-MCMD']

data = {}
time = np.linspace(0, 1000, 11)
# print(time)
for T in Temp:
    data[0] = np.loadtxt(f'Random/{size}/data/msd_{size}_{T}.txt')
    data[1] = np.loadtxt(f'MCMD/{size}/data/msd_{size}_{T}.txt')
    for i in range(len(label)):
        time = data[i][:, 0] / 1000
        msd = data[i][:, 1]
        # # msd = msd - msd[1]
        # step = 10
        # time_sub = time[::step]
        # print(time_sub)
        # msd_sub = msd[::step]
        # print(msd_sub)
        x_smooth = np.linspace(time.min(), time.max(), 100)
        
        # Use time as x and msd as y (not the other way around)
        tck = interpolate.splrep(time, msd, s=4)
        y_smooth = interpolate.splev(x_smooth, tck, der=0)
        
        if T == Temp[0]:
            # tck = interpolate.splrep(time, data[i], s=0)
            # y_smooth = interpolate.splev(x_smooth, tck, der=0)
            plt.plot(x_smooth, y_smooth-y_smooth[0], linewidth=linewidth, color=rgb[i % 2], linestyle = linestyles[i % 2])
            # plt.plot(time, msd-msd[1], linewidth=linewidth, color=rgb[i % 2], linestyle = linestyles[i % 2])
            # plt.plot(time, data[i], 'o')
        else:
            # tck = interpolate.splrep(time, data[i], s=0)
            # y_smooth = interpolate.splev(x_smooth, tck, der=0)
            # coeffs = np.polyfit(time, data[i+2], 3)
            # poly = np.poly1d(coeffs)
            # y_smooth = poly(x_smooth)
            # y_smooth = y_smooth - y_smooth[0]
            plt.plot(x_smooth, y_smooth- y_smooth[0], label=label[i % 2], linewidth=linewidth, color=rgb[i % 2], linestyle = linestyles[i % 2])
            # plt.plot(time, msd-msd[1], label=label[i % 2], linewidth=linewidth, color=rgb[i % 2], linestyle = linestyles[i % 2])
            # plt.plot(time, data[i+2], 'o')
plt.text(700, 0.08, '800 K', fontsize=fontsize)
plt.text(700, 0.5, '1000 K', fontsize=fontsize)   
plt.text(350, 0.53, '800 K', fontsize=fontsize)
plt.text(350, 0.95, '1000 K', fontsize=fontsize) 
# plt.text(430, 1.06, r'Atoms = 5344, $\mathrm{D_S}$ = 2.78', fontsize=fontsize)  

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
ax.set_ylim(0.0, 1.5)
ax.legend(frameon=True)

plt.savefig('MSD-FeNiCrCoCu.svg', format='svg', dpi=300, bbox_inches='tight')
plt.savefig('MSD-FeNiCrCoCu.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()