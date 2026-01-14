import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.animation import FuncAnimation, PillowWriter, FFMpegWriter

# 1. Data preparation
data = {}

# data[0] = np.loadtxt(f'glass_300_2.751_1/data/strain_300_2.751_1.txt')
# data[1] = np.loadtxt(f'para_300_3.346_5_1/data/strain_300_3.346_5.txt')
# data[2] = np.loadtxt(f'para_300_3.346_5_2/data/strain_300_3.346_6.txt')
# data[3] = np.loadtxt(f'para_300_3.346_5_3/data/strain_300_3.346_6.txt')
data[0] = np.loadtxt(f'para_300_3.346_5_4/data/strain_300_3.346_6.txt')
# data[5] = np.loadtxt(f'para_300_3.346_5_5/data/strain_300_3.346_6.txt')

# 2. Plot Settings
## 2.1 Color Setting
rgb = [(232, 68, 69), (75, 101, 175), (85, 103, 117), (128, 198, 109), (136, 135, 203), (251, 132, 2), (254, 183, 5)]
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
    "font.family": ["Times New Roman"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman", 
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
    "axes.labelsize": labelsize,
    "xtick.labelsize": ticksize,
    "ytick.labelsize": ticksize,
    "legend.fontsize": legendsize,
})
fig, ax = plt.subplots(figsize=(6.5, 5))

## 2.4 labels setting
ax.set_xlabel('Strain (%)')
ax.set_ylabel('Stress (GPa)')
label = ["Glassy", "Para-1", "Para-2", "Para-3", "Para-4", "Para-5"]
label = ["Paracrystalline"]

# 2.5 tick setting
framewidth = 1.5                # frame: the frame line width
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

# 3. Animation
## 3.1 Initialize lines
lines = []
for i in range(len(data)):
    line, = ax.plot([0], [0], label=label[i], linewidth=linewidth, color=rgb[i % len(rgb)])
    lines.append(line)
fig.tight_layout()
## 3.2 Animation function
max_len = max(len(d[:, 0]) for d in data.values())

def animate(frame):
    for i, d in enumerate(data.values()):
        strain = d[:, 0] - d[:, 0].min()
        stress = d[:, 1] - d[:, 1].min()
        if frame < len(strain):
            lines[i].set_data(strain[:frame]*100, stress[:frame])
        else:
            lines[i].set_data(strain*100, stress)
    return lines
ax.legend(frameon=False)

anim = FuncAnimation(fig, animate, frames=max_len, interval=100, blit=True)

# 5. Save as GIF or MOV
FPS_Speed = 80
anim.save(
    'crack-stress-para.gif',
    writer=PillowWriter(fps=FPS_Speed),
    dpi=300
)

anim.save(
    'crack-stress-para.mov',
    writer=FFMpegWriter(fps=FPS_Speed, bitrate=3000),
    dpi=300
)

plt.show()
