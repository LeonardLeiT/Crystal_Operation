import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import matplotlib.patches as patches

# 1. Data preparation
Stru_name = 'FCC'
Stru_index = 2 + 6
data = pd.read_excel('FiNiCrCoCu_structure_fraction.xlsx')
print(data.info)

time = data.iloc[:, 1] * 1e-6

# 2. Plot Settings
# plt.figure(figsize=(22, 5))

## 2.1 Color Setting
rgb = [(251, 132, 2), (75, 101, 175), (232, 68, 69), (128, 198, 109), (254, 179, 174), (136, 135, 203), (85, 103, 117)]
rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
markers = ['o', 's', '^', 'D', '*']
linestyles = ['-', '--', '-.', ':']

## 2.2 line Size Setting
linewidth = 2.0
marker_size = 200

## 2.3 font Setting 
ticksize = 26
labelsize = 26
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
fig, ax1 = plt.subplots(figsize=(27, 4))

# Left y-axis for Structure Fraction
ax1.set_xlabel(r'$t$ (ns)')
ax1.set_ylabel(f'{Stru_name}\nFraction (%)')
ax1.tick_params(axis='y')
ax1.set_zorder(1)
ax1.patch.set_visible(False)

# Right y-axis for Elements
ax2 = ax1.twinx()
ax2.set_ylabel('Element Fraction (%)')
ax2.tick_params(axis='y')
ax2.patch.set_visible(False)
ax2.set_zorder(2)

label = ["Structure Fraction", "Co", "Cr", "Cu", "Fe", "Ni"]

# Plot Elements on right axis
for i in range(5):
    element = data.iloc[:, Stru_index + i + 1] * 100
    ax2.plot(time, element, label=label[i + 1], linewidth=linewidth, color=rgb[(i) % len(rgb)])

# Plot Structure Fraction on left axis
Struc_Frac = data.iloc[:, Stru_index] * 100
ax1.plot(time, Struc_Frac, label=label[0], linewidth=linewidth+2, color="black")

## 2.5 tick setting
framewidth = 2
for spine in ax1.spines.values():
    spine.set_linewidth(framewidth)
for spine in ax2.spines.values():
    spine.set_linewidth(framewidth)

# X-axis ticks (shared)
major_x_spacing = 10
ax1.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))

# Left y-axis ticks
major_y_spacing_left = 20
ax1.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing_left))
ax1.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing_left / 2))
ax1.set_ylim(1, 90)

# Right y-axis ticks (adjust range as needed for elements)
major_y_spacing_right = 10  # or adjust based on your element data range
ax2.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing_right))
ax2.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing_right / 2))
ax2.set_ylim(14, 45)  # Set appropriate range for elements

# Tick parameters
ax1.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax1.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax2.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax2.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')

ax1.set_xlim(0, 100)

ax1.legend(loc="upper left")
ax2.legend(bbox_to_anchor=(0, 0.8), loc="upper left", ncol=3)

x_min1, x_max1 = 0, 5
y_min1, y_max1 = 15, 25

rect = patches.Rectangle((x_min1, y_min1), 
                        x_max1 - x_min1, 
                        y_max1 - y_min1,
                        linewidth=linewidth, 
                        edgecolor=(1, 0, 0, 0.8), 
                        facecolor=(1, 0, 0, 0.1),
                        linestyle='--',
                        zorder=10)
ax2.add_patch(rect)

x_min2, x_max2 = 27, 31
y_min2, y_max2 = 15, 25

rect = patches.Rectangle((x_min2, y_min2), 
                        x_max2 - x_min2, 
                        y_max2 - y_min2,
                        linewidth=linewidth, 
                        edgecolor=(1, 0, 0, 0.8), 
                        facecolor=(1, 0, 0, 0.1),
                        linestyle='--',
                        zorder=10)
ax2.add_patch(rect)

x_min3, x_max3 = 86, 94
y_min3, y_max3 = 15, 25

rect = patches.Rectangle((x_min3, y_min3), 
                        x_max3 - x_min3, 
                        y_max3 - y_min3,
                        linewidth=linewidth, 
                        edgecolor=(1, 0, 0, 0.8), 
                        facecolor=(1, 0, 0, 0.1),
                        linestyle='--',
                        zorder=10)
ax2.add_patch(rect)

# 3. add subgraph
left = 0.39
bottom = 0.38  # Other 0.48 HCP 0.45 FCC 0.38
# Setting Position [left, bottom, width, height] 
subax_left = plt.axes([left, bottom, 0.08, 0.35])
subax_left.set_zorder(2)
subax_left.patch.set_visible(False)

subax_right = subax_left.twinx()
subax_right.set_zorder(1)

# Plot Elements on right axis
for i in range(5):
    element = data.iloc[:, Stru_index + i + 1] * 100
    subax_right.plot(time, element, linewidth=linewidth, 
                     color=rgb[i % len(rgb)], zorder=1)

# Plot Structure Fraction on left axis
Struc_Frac = data.iloc[:, Stru_index] * 100
subax_left.plot(time, Struc_Frac, linewidth=linewidth+2, 
                color="black", zorder=5)

# tick setting
major_x_spacing = 2
major_y_spacing_left = 50
major_y_spacing_right = 5

subax_left.set_xlim(x_min1, x_max1)
subax_left.set_ylim(0, 100)
subax_right.set_ylim(y_min1, y_max1)

subax_left.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
subax_left.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing_left))
subax_right.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing_right))

subax_left.tick_params(axis='both', which='major', length=4, width=1, direction='in')
subax_right.tick_params(axis='y', which='major', length=4, width=1, direction='in')

for spine in subax_left.spines.values():
    spine.set_linewidth(framewidth)
for spine in subax_right.spines.values():
    spine.set_linewidth(framewidth)

# 4. add subgraph
# Setting Position [left, bottom, width, height] 
subax_left1 = plt.axes([left + 0.15, bottom, 0.08, 0.35])
subax_left1.set_zorder(2)
subax_left1.patch.set_visible(False)

subax_right1 = subax_left1.twinx()
subax_right1.set_zorder(1)

# Plot Elements on right axis
for i in range(5):
    element = data.iloc[:, Stru_index + i + 1] * 100
    subax_right1.plot(time, element, linewidth=linewidth, 
                     color=rgb[i % len(rgb)], zorder=1)

# Plot Structure Fraction on left axis
Struc_Frac = data.iloc[:, Stru_index] * 100
subax_left1.plot(time, Struc_Frac, linewidth=linewidth+2, 
                color="black", zorder=5)

# tick setting
major_x_spacing = 2
major_y_spacing_left = 50
major_y_spacing_right = 5

subax_left1.set_xlim(x_min2, x_max2)
subax_left1.set_ylim(0, 100)
subax_right1.set_ylim(y_min2, y_max2)

subax_left1.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
subax_left1.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing_left))
subax_right1.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing_right))

subax_left1.tick_params(axis='both', which='major', length=4, width=1, direction='in')
subax_right1.tick_params(axis='y', which='major', length=4, width=1, direction='in')

for spine in subax_left1.spines.values():
    spine.set_linewidth(framewidth)
for spine in subax_right1.spines.values():
    spine.set_linewidth(framewidth)

# 5. add subgraph
# Setting Position [left, bottom, width, height] 
subax_left2 = plt.axes([left + 0.3, bottom, 0.08, 0.35])
subax_left2.set_zorder(2)
subax_left2.patch.set_visible(False)

subax_right2 = subax_left2.twinx()
subax_right2.set_zorder(1)

# Plot Elements on right axis
for i in range(5):
    element = data.iloc[:, Stru_index + i + 1] * 100
    subax_right2.plot(time, element, linewidth=linewidth, 
                     color=rgb[i % len(rgb)], zorder=1)

# Plot Structure Fraction on left axis
Struc_Frac = data.iloc[:, Stru_index] * 100
subax_left2.plot(time, Struc_Frac, linewidth=linewidth+2, 
                color="black", zorder=5)

# tick setting
major_x_spacing = 4
major_y_spacing_left = 50
major_y_spacing_right = 5

subax_left2.set_xlim(x_min3, x_max3)
subax_left2.set_ylim(0, 100)
subax_right2.set_ylim(y_min3, y_max3)

subax_left2.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
subax_left2.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing_left))
subax_right2.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing_right))

subax_left2.tick_params(axis='both', which='major', length=4, width=1, direction='in')
subax_right2.tick_params(axis='y', which='major', length=4, width=1, direction='in')

for spine in subax_left2.spines.values():
    spine.set_linewidth(framewidth)
for spine in subax_right2.spines.values():
    spine.set_linewidth(framewidth)

plt.savefig(f'{Stru_name}.svg', format='svg', dpi=300, bbox_inches='tight')
plt.savefig(f'{Stru_name}.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()