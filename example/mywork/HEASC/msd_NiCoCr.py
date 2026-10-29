import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy import interpolate

# 1. Data Setting
Temp = [700, 900]
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
label = ["NiCoCr-Random", 'NiCoCr-MCMD']

data = {}
data=[[0.156109642945214, 0.418007141050415, 0.526886153392257, 0.643048679950303, 0.727411516333147, 
       0.824540196788257, 0.906199342526078, 0.99314232078384, 1.08767327076622, 1.15626078218249, 1.22387509884051],
      [0.126113041079261, 0.133975551591252, 0.138962349868493, 0.141663271713949, 0.14477654989258, 
       0.14400015030914, 0.147685099350079, 0.155730205086094, 0.153521853655357, 0.169337078517438, 0.1724576955942],
      [0.216538350287761, 1.31472379078327, 1.97250401219108, 2.5691448744911, 2.96175285924021, 3.58619237104039, 4.03107568214751, 4.52073081022375, 4.86408448907885, 5.16566842045225, 5.42164676644732],
      [0.175292506141337, 0.253294688857122, 0.34884764181344, 0.418852790824976, 0.493972981819248, 0.576980114281083, 0.688167873654442, 0.783994101580208, 0.916281046057887, 1.02619103923783, 1.17421403134943]]
time = np.linspace(0, 1000, 11)
# print(time)
for T in Temp:
    # data[0] = np.loadtxt(f'Random/{size}/data/msd_{size}_{T}.txt')
    # data[1] = np.loadtxt(f'MCMD/{size}/data/msd_{size}_{T}.txt')
    for i in range(len(label)):
        # time = data[i][:, 0] / 1000
        # msd = data[i][:, 1]
        # # msd = msd - msd[1]
        # step = 10
        # time_sub = time[::step]
        # print(time_sub)
        # msd_sub = msd[::step]
        # print(msd_sub)
        x_smooth = np.linspace(time.min(), time.max(), 100)
        
        # Use time as x and msd as y (not the other way around)
        # tck = interpolate.splrep(time_sub, msd_sub, s=4)
        # y_smooth = interpolate.splev(x_smooth, tck, der=0)
        
        if T == Temp[0]:
            tck = interpolate.splrep(time, data[i], s=0)
            y_smooth = interpolate.splev(x_smooth, tck, der=0)
            plt.plot(x_smooth, y_smooth-y_smooth[0], linewidth=linewidth, color=rgb[i % 2], linestyle = linestyles[i % 2])
            # plt.plot(time, msd-msd[1], linewidth=linewidth, color=rgb[i % 2], linestyle = linestyles[i % 2])
            # plt.plot(time, data[i], 'o')
        else:
            tck = interpolate.splrep(time, data[i+2], s=0)
            y_smooth = interpolate.splev(x_smooth, tck, der=0)
            # coeffs = np.polyfit(time, data[i+2], 3)
            # poly = np.poly1d(coeffs)
            # y_smooth = poly(x_smooth)
            # y_smooth = y_smooth - y_smooth[0]
            plt.plot(x_smooth, y_smooth- y_smooth[0], label=label[i % 2], linewidth=linewidth, color=rgb[i % 2], linestyle = linestyles[i % 2])
            # plt.plot(time, msd-msd[1], label=label[i % 2], linewidth=linewidth, color=rgb[i % 2], linestyle = linestyles[i % 2])
            # plt.plot(time, data[i+2], 'o')
plt.text(700, 0.05, '700 K', fontsize=fontsize)
plt.text(700, 0.5, '900 K', fontsize=fontsize)   
plt.text(150, 0.5, '700 K', fontsize=fontsize)
plt.text(150, 1.3, '900 K', fontsize=fontsize) 
plt.text(430, 1.06, r'Atoms = 5344, $\mathrm{D_S}$ = 2.78', fontsize=fontsize)  

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
ax.legend(frameon=False, loc='upper right')

plt.savefig('MSD-NiCoCr.svg', format='svg', dpi=300, bbox_inches='tight')
plt.savefig('MSD-NiCoCr.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()