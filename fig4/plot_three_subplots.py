import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import FuncFormatter

plt.rcParams.update({'font.size': 12})
# Define fitting functions
def func2(params, R, xi, phi):
    w1, w2 = params
    R_xi = R / xi
    return w1 * R_xi**(3 - 1/0.588) + w2 * R_xi**2 - np.log(1 - phi)

def saw(vf, c):
    nu = 0.588
    return c * 0.6 * vf ** (-nu / (3 * nu - 1))

def scientific_formatter(x, p):
    if x == 0:
        return '0'
    exp = np.floor(np.log10(x))
    coef = x / 10**exp

    return f'${coef:.1f}\\times10^{{{int(exp)}}}$'

fig = plt.figure(figsize=(7, 11))

gs = fig.add_gridspec(2, 2, height_ratios=[2, 1], wspace=0.3, hspace=0.25)

ax1 = fig.add_subplot(gs[0, :])

ax2 = fig.add_subplot(gs[1, 0])

ax3 = fig.add_subplot(gs[1, 1])

for ax in [ax1, ax2, ax3]:
    ax.set_xticks([])
    ax.set_yticks([])

# Plot 1: Original fitting plot
filenames = [
    'partitioning_for_vf29.dat', 'partitioning_for_vf33.dat', 'partitioning_for_vf37.dat',
    'partitioning_for_vf41.dat', 'partitioning_for_vf46.dat', 'partitioning_for_vf51.dat',
    'partitioning_for_vf57.dat', 'partitioning_for_vf64.dat', 'partitioning_for_vf71.dat',
    'partitioning_for_vf79.dat'
][::-1]
phis = np.exp(np.linspace(np.log(0.03), np.log(0.08), 10))[::-1]
press = np.loadtxt('pressures.txt')[::-1]
nu = 0.588
xis = (4.14 / (3 * nu - 1) / press) ** (1/3)

all_R = []
all_F = []
all_xis = []
all_phis = []


for filename, xi, phi in zip(filenames, xis, phis):
    data = np.loadtxt(filename, delimiter=' ')
    R = data[:18, 0]
    y = data[:18, 1]
    F = -np.log(y) - (4 * np.pi / 3) / (3 * 0.588 - 1) * (R / xi)**3 
    all_R.extend(R)
    all_F.extend(F)
    all_xis.extend([xi] * len(R))
    all_phis.extend([phi] * len(R))

all_R = np.array(all_R)
all_F = np.array(all_F)
all_xis = np.array(all_xis)
all_phis = np.array(all_phis)

def curve_func(R, w1, w2):
    return func2((w1, w2), R, all_xis, all_phis)

params, _ = curve_fit(curve_func, all_R, all_F)
print(params, " a1 and a2 parameters")
# Plot in ax1
colors = plt.cm.viridis(np.linspace(0, 1, len(filenames)))
fitted_line_label_added = False

for i, (filename, xi, phi) in enumerate(zip(filenames, xis, phis)):
    data = np.loadtxt(filename, delimiter=' ')
    R = data[:18, 0]
    y = data[:18, 1]
    F = -np.log(y) - (4 * np.pi / 3) / (3 * 0.588 - 1) * (R / xi)**3
    
    R_fit = np.linspace(min(R), max(R), 500)
    F_fit = func2(params, R_fit, xi, phi)
    
    if not fitted_line_label_added:
        ax1.plot(R_fit, F_fit, color='black', linestyle='-', linewidth=1, label='Fitted curves')
        fitted_line_label_added = True
    else:
        ax1.plot(R_fit, F_fit, color='black', linestyle='-', linewidth=1, zorder=-1)
    
    ax1.scatter(R, F, color=colors[i], marker='o', s=80, label=f'ξ {xi:.2f} nm')
    ax1.set_xticks([])
    ax1.set_yticks([])

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.legend(title='Correlation lengths', fontsize=8, title_fontsize=11)

print(data[:18, 0])

ax1.minorticks_off()
log_yticks = np.linspace(np.log10(0.1), np.log10(10), 11)  # log values from -1 to 1
ytick_values = 10**log_yticks  # Convert to actual values
ytick_labels = [f'$10^{0}$' if abs(val) < 1e-1 else f'$10^{{{val:.1f}}}$' for val in log_yticks]  # Format as 10^x
ax1.set_yticks(ytick_values)
ax1.set_yticklabels(ytick_labels)
ax1.set_ylim([0.13, 10])

log_xticks = np.linspace(np.log10(0.5), np.log10(2.4843), 8)  # Adjust range as needed
xtick_values = 10**log_xticks
xtick_labels = [f'$10^{0}$' if abs(val) < 1e-2 else f'$10^{{{val:.1f}}}$' for val in log_xticks]
ax1.set_xticks(xtick_values)
ax1.set_xticklabels(xtick_labels)

# Plot 2: Correlation length fitting


phis = np.exp(np.linspace(np.log(0.03), np.log(0.08), 10))
press = np.loadtxt('pressures.txt')
nu = 0.588
xis = (4.14 / (3 * nu - 1) / press) ** (1/3)

popt2, _ = curve_fit(saw, phis, xis)
c = popt2[0]
print(c, ' c for correllation length')

vf_fine = np.linspace(0.027, 0.09, 100)
xi_line = saw(vf_fine, c)

colors = plt.cm.viridis(np.linspace(1, 0, len(phis)))
for i, (vf_value, xi_value) in enumerate(zip(phis, xis)):
    ax2.plot(vf_value, xi_value, 'o', markersize=10, color=colors[i])

ax2.plot(vf_fine, xi_line, '-', c='black')
ax2.set_xlabel('phi', fontsize=12)
ax2.set_ylabel('xi', fontsize=12)

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.minorticks_off()
ytick2_values = xis
ytick2_labels = [f'{val:.2f}' for val in ytick2_values]  # Format labels as you prefer
ax2.set_yticks(ytick2_values)  # Set the positions of the y-ticks
ax2.set_yticklabels(ytick2_labels)  # Set the labels of the y-ticks

xtick2_values = phis
xtick2_values = xtick2_values[::3]
xtick2_labels = [f'{val:.2f}' for val in xtick2_values]  # Format labels as you prefer
ax2.set_xticks(xtick2_values)  # Set the positions of the y-ticks
ax2.set_xticklabels(xtick2_labels)  # Set the labels of the y-ticks


# Plot 3: Partitioning coefficient prediction
R_pred = np.linspace(0.01, 3, 100)
for i, phi in enumerate(phis):
    xi_pred = saw(phi, c)
    F_pred = func2(params, R_pred, xi_pred, phi)
    y_pred = np.exp(-F_pred)
    print(np.exp(- func2(params, 2, xi_pred, phi)))
    ax3.plot(R_pred, y_pred, color=colors[i], alpha = 0.4)

BOND_LENGTH = 0.38
BEAD_RADIUS = 0.3
R = BEAD_RADIUS
H = R - BOND_LENGTH / 2
V1 = np.pi / 3 * (3 * R - H) * H**2
V2 = 4 * np.pi / 3 * R**3 - 2 * V1
phi = 0.4 * 250 / 110* 0.6 * V2
xi_pred = saw(phi, c)
print(xi_pred, 'xi_pred')
F_pred = func2(params, R_pred, xi_pred, phi) 
y_pred = np.exp(-F_pred)
print(np.exp(- func2(params, 2, xi_pred, phi)), func2(params, 2, xi_pred, phi),'exclusion')
ax3.plot(R_pred, y_pred, color='black', label='Predicted')

ax3.set_xlabel('R (nm)', fontsize=12)
ax3.set_ylabel('Partitioning Coefficient', fontsize=12)
ax3.minorticks_off()
ytick3_values = np.linspace(0,1,6)
ytick3_labels = [f'{val:.1f}' for val in ytick3_values]  # Format labels as you prefer
ax3.set_yticks(ytick3_values)  # Set the positions of the y-ticks
ax3.set_yticklabels(ytick3_labels)  # Set the labels of the y-ticks

xtick3_values = np.linspace(0,3,7)
xtick3_labels = [f'{val:.1f}' for val in xtick3_values]  # Format labels as you prefer
ax3.set_xticks(xtick3_values)  # Set the positions of the y-ticks
ax3.set_xticklabels(xtick3_labels)  # Set the labels of the y-ticks

plt.tight_layout()
plt.savefig('three_subplots_fig4', dpi = 1000)
plt.show()
