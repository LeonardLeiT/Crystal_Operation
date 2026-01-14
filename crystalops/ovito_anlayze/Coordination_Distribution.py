import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.optimize import curve_fit
import pandas as pd

# 1. Data preparation
filenames = ["Cool_500_0.data", "Cool_1000_0.data"]
coord_vals = {}
rates = {}
for i, fname in enumerate(filenames):
    df = pd.read_excel(f'{fname}_Coordination_Distribution.xlsx')
    coord_vals[i] = df['coord_vals'].values
    rates[i] = df['rate'].values

# 2. Plot Settings
plt.figure(figsize=(7, 5))

## 2.1 Color Setting
rgb = [(232, 68, 69), (85, 103, 117), (75, 101, 175), (128, 198, 109), (136, 135, 203), (251, 132, 2)]
rgb = [(r/255, g/255, b/255) for r, g, b in rgb]
markers = ['s', '*', '^', 'D', '*']
linestyles = ['-', '--', '-.', ':']

## 2.2 line Size Setting
linewidth = 3.0
marker_size = 15

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
plt.xlabel("Coordination Number")
plt.ylabel("Rate (%)")
label = ["Paracrystalline", "Glassy"]

## 2.5 tick setting
framewidth = 1.5                # frame: the frame line width
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(framewidth)
    
major_x_spacing = 2
major_y_spacing = 20
ax.xaxis.set_major_locator(ticker.MultipleLocator(major_x_spacing))
ax.yaxis.set_major_locator(ticker.MultipleLocator(major_y_spacing))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(major_x_spacing / 2))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(major_y_spacing / 2))
ax.tick_params(axis='both', which='major', length=8, width=1.5, direction='in')
ax.tick_params(axis='both', which='minor', length=4, width=1.5, direction='in')
ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))

x_min = 1.5
x_max = 9.5
ax.set_xlim(x_min, x_max)
ax.set_ylim(0, 100)

# ---------- Plot ----------
for i in range(len(filenames)):
    # Compute coordination distribution
    coord_val = coord_vals[i]
    rate = rates[i] * 100

    # μ = coordination number mode
    mu = coord_val[np.argmax(rate)]

    # Gaussian function with fixed μ
    def gaussian_fixed_mu(x, a, sigma):
        return a * np.exp(-(x - mu)**2 / (2 * sigma**2))

    # Initial guess: a=max(rate), sigma=std deviation
    p0 = [max(rate), np.std(coord_val)]

    try:
        popt, _ = curve_fit(gaussian_fixed_mu, coord_val, rate, p0=p0)
        a_fit, sigma_fit = popt
        x_smooth = np.linspace(x_min, x_max, 400)
        fit_smooth = gaussian_fixed_mu(x_smooth, a_fit, sigma_fit)
    except RuntimeError:
        popt = None
        x_smooth = None
        fit_smooth = None

    plt.plot(coord_val, rate, marker=markers[i % len(markers)],
             markeredgecolor=rgb[i % len(rgb)],
             markerfacecolor=rgb[i % len(rgb)],
             markersize=marker_size,
             linestyle='None',  # points only
             label=label[i])

    if fit_smooth is not None:
        plt.plot(x_smooth, fit_smooth, color=rgb[i % len(rgb)],
                #  label=f"Gaussian Fit ({label[i]})",
                 label=f"Gaussian Fit",
                 linewidth=linewidth,
                 linestyle=linestyles[i % len(linestyles)])


ax.legend(frameon=False)
plt.savefig('coord_rate_1.svg', format='svg', dpi=300, bbox_inches='tight')
plt.savefig('coord_rate_1.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.tight_layout()
plt.show()