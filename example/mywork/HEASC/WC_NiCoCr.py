import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import make_interp_spline

# 1. Data preparation
data = np.array([
    [-0.02626596, -0.523263751, -0.612855729, -0.659476427, -0.671902966, -0.701904134, -0.705717139, -0.697348366, -0.708176244, -0.70319855, -0.713181999],
    [0.013829124, 0.317727045, 0.359158377, 0.376143393, 0.380817021, 0.395367708, 0.390825729, 0.399514035, 0.391782608, 0.398861298, 0.406609289],
    [0.0182234, 0.216388821, 0.261173246, 0.293425106, 0.300601128, 0.314628767, 0.324645366, 0.307906259, 0.327103663, 0.313183214, 0.316492616],
    [0.102978584, 0.046965699, 0.02763238, 0.033463605, 0.036399951, 0.032108412, 0.036844034, 0.021940517, 0.037496956, 0.029293463, 0.027046813],
    [-0.109491541, -0.356624186, -0.376916706, -0.399411186, -0.40628779, -0.415910844, -0.417819213, -0.413771587, -0.420750034, -0.417294334, -0.424000748],
    [0.078144135, 0.121281427, 0.098365886, 0.085664352, 0.085208936, 0.081593608, 0.073536686, 0.088079488, 0.074373493, 0.084372972, 0.087900056]
])

x_values = np.array(
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
)

x_values = x_values

# 2. Plot Settings
plt.figure(figsize=(6.5, 5))

## 2.1 Color Setting
rgb = [(232, 68, 69), (85, 103, 117), (75, 101, 175), (128, 198, 109), (136, 135, 203), (249, 159, 35)]
rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
markers = ['o', 's', '^', 'D', '*', 'v', '>']
# linestyles = ['-', '--', '-.', ':']
linestyles = ['-', '--', '-.']
markerfacecolor = 'white'       # color: the color in the data point
alpha = 0.9                     # alpha: the transparency of the data point

## 2.2 line Size Setting
linewidth = 2                   # line: the line width
markeredgewidth  = 2            # markeredgewidth: line width of markeredgewidth
markersize = 10                 # marker: the data point size

## 2.3 font Setting 
ticksize = 20                   # ticks: the size of the tick labels
labelsize = 20                  # labels: the size of the label
legendsize = 20                 # legends: the size of the legend
plt.rcParams.update({
    "font.family": "Times New Roman",
    "axes.labelsize": labelsize,
    "xtick.labelsize": ticksize,
    "ytick.labelsize": ticksize,
    "legend.fontsize": legendsize,
})

## 2.4 labels setting
plt.ylabel(r'Warren-Cowley parameter $\alpha$')
plt.xlabel(r'MC steps ($\times 10^4$)')
labels = ['Ni-Ni', 'Ni-Co', 'Ni-Cr', 'Co-Co', 'Co-Cr', 'Cr-Cr']

## 2.5 Plot Figure
for i in range(len(data)):
    plt.plot(x_values, data[i], linewidth=linewidth, color=rgb[i % len(rgb)], linestyle=linestyles[i % len(linestyles)], marker=markers[i % len(markers)], label=labels[i], markersize= markersize)
plt.axhline(y = 0, color=(115/255, 115/255, 115/255), linestyle=':', linewidth=linewidth, alpha = alpha)

# 2.6 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)

major_x_spacing = 3
major_y_spacing = 0.3
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.set_xlim(0, 12)
ax.set_ylim(-0.9, 0.6)
ax.legend(frameon=True, loc="lower right")


plt.tight_layout()
plt.savefig("WC_NiCoCr.svg", bbox_inches='tight')
plt.savefig("WC_NiCoCr.pdf", bbox_inches='tight')
plt.show()

