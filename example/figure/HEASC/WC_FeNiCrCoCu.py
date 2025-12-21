import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline
import pandas as pd

# 1. Data preparation
data = np.array([
    [-0.09200114, -0.115549946, -0.105960516, -0.10641376, -0.125361793, -0.088809165, -0.060473678, -0.116019805, -0.055306604, -0.125488153, -0.094182949],
    [0.036723642, -0.118959542, -0.113773662, -0.099402848, -0.073629736, -0.05826757, -0.116181348, -0.083409993, -0.085837862, -0.102805464,-0.100423488],
    [-0.045729396, -0.011931007, 0.021951396, 0.009106529, -0.023863491, -0.067960141, -0.016733912, -0.036059243, -0.042839374, -0.02141561, -0.046213533],
    [0.048059175, 0.053665486, 0.064661659, 0.014649747, 0.062679849, 0.040972988, 0.011878376, 0.066660525, 0.008634458, 0.05598083, 0.033189054],
    [0.038023193, 0.175664752, 0.116674677, 0.17526762, 0.147897347, 0.157441249, 0.167893017, 0.149037368, 0.165734774, 0.173079598, 0.190336593],
    [-0.000212677, -0.073331537, -0.084974792, -0.115886042, -0.084229095, -0.073742841, -0.105227449, -0.094810008, -0.045705173, -0.092824288, -0.038145679],
    [0.009408083, -0.002386992, 0.033216007, 0.004206212, 0.030990002, 0.019630855, 0.038712651, -0.018475725, -0.008246856, 0.035720515, 0.019139518],
    [0.00024072, 0.017855887, -0.001650096, -0.001975882, 0.000210231, -0.033106919, 0.014262701, 0.02681016, -0.019035479, -0.009676078, -0.005854232],
    [-0.023359674, 0.192567877, 0.178520452, 0.217837438, 0.141543641, 0.159015206, 0.180270932, 0.17262916, 0.170136441, 0.180877699, 0.140886049],
    [0.057416283, 0.016011048, 0.011789769, -0.009120817, -0.01210908, 0.016582204, -0.009167959, 0.028565706, 0.053627676, -0.028438408, 0.051029537],
    [-0.031048247, -0.057250196, -0.078642828, -0.030072846, 0.000853952, -0.007451713, -0.032141319, -0.035586815, -0.040978026, -0.028465203, -0.005701483],
    [0.023266407, 0.061902803, 0.031663691, 0.038048478, 0.022371875, 0.059222753, 0.026373112, 0.077829539, 0.046325387, 0.051804848, 0.003194542],
    [0.003472721, -0.017410264, -0.031323236, -0.017850892, -0.077009378, -0.024902138, -0.022926172, -0.061950157, 0.022260507, -0.030217989, -0.053762967],
    [-0.042153647, -0.01111636, 0.02501966, 0.019979119, -0.006474642, -0.00205658, 0.014565095, -0.001577376, 0.017137474, 0.006859714, 0.014131285],
    [0.004413756, -0.409778315, -0.344864828, -0.446051612, -0.306495163, -0.364070888, -0.380037191, -0.391514162, -0.396966568, -0.406965835, -0.350348927],
])

data = pd.read_excel('FiNiCrCoCu_wc_results.xlsx')

x_values = np.array(
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
)

x_values = data.iloc[:, 0] * 1e-4

# 2. Plot Settings
plt.figure(figsize=(9, 6))

## 2.1 Color Setting
# rgb = [(232, 68, 69), (85, 103, 117), (75, 101, 175), (128, 198, 109), (136, 135, 203), (249, 159, 35)]
rgb = [(181, 206, 78), (151, 208, 197), (247, 172, 83), (145, 204, 192), (236, 110, 102), 
       (189, 119, 149), (124, 121, 121), (150, 59, 121), (127, 171, 209), (243, 152, 101),
       (82, 170, 220), (199, 193, 222), (238, 182 ,212), (200, 151, 54), (45, 136, 117)]
rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
markers = ['o', 's', '^', 'D', '*', 'v', '>']
linestyles = ['-', '--', '-.', ':']

markerfacecolor = 'white'       # color: the color in the data point
alpha = 0.9                     # alpha: the transparency of the data point

## 2.2 line Size Setting
linewidth = 2                   # line: the line width
markeredgewidth  = 2            # markeredgewidth: line width of markeredgewidth
markersize = 10                 # marker: the data point size

## 2.3 font Setting 
ticksize = 20                   # ticks: the size of the tick labels
labelsize = 20                  # labels: the size of the label
legendsize = 18                 # legends: the size of the legend
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
plt.ylabel(r'Warren-Cowley parameter $\alpha$')
plt.xlabel(r'MC steps ($\times 10^4$)')
labels = ['Fe-Fe', 'Fe-Ni', 'Fe-Cr', 'Fe-Co', 'Fe-Cu',
          'Ni-Ni', 'Ni-Cr', 'Ni-Co', 'Ni-Cu', 'Cr-Cr',
          'Cr-Co','Cr-Cu', 'Co-Co', 'Co-Cu', 'Cu-Cu']

{1: "Co", 2: "Cr", 3: "Cu", 4: "Fe", 5: "Ni"}
labels = ['Co-Co', 'Co-Cr', 'Co-Cu', 'Co-Fe', 'Co-Ni',
          'Cr-Cr', 'Cr-Cu', 'Cr-Fe', 'Cr-Ni', 'Cu-Cu',
          'Cu-Fe','Cu-Ni', 'Fe-Fe', 'Fe-Ni', 'Ni-Ni']
## 2.5 Plot Figure
for i in range(len(labels)):
    plt.plot(x_values, data.iloc[:, i+3], linewidth=linewidth, color=rgb[i % len(rgb)], linestyle=linestyles[i % len(linestyles)], marker=markers[i % len(markers)], label=labels[i], markersize= markersize)
plt.axhline(y = 0, color=(115/255, 115/255, 115/255), linestyle=':', linewidth=linewidth, alpha = alpha)

# 2.6 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)

major_x_spacing = 3
major_y_spacing = 0.2
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.set_xlim(0, 12)
ax.set_ylim(-0.65, 0.25)
ax.legend(frameon=True, loc="lower right", ncol=5, columnspacing=0.4)


plt.tight_layout()
plt.savefig("WC_FeNiCrCoCu.svg", bbox_inches='tight')
plt.savefig("WC_FeNiCrCoCu.pdf", bbox_inches='tight')
plt.show()

