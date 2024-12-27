import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

plt.rcdefaults()
def scientific_formatter(x, p):
    if x == 0:
        return '0'
    exp = np.floor(np.log10(x))
    coef = x / 10**exp

    return f'${coef:.1f}\\times10^{{{int(exp)}}}$'


plt.rcParams.update({'font.size': 12})
def func(R_xi, w1, w2):
    nu = 0.588
    return w1 * (R_xi)**(3 - 1/nu) + w2 * (R_xi)**2

filenames = [
    'partitioning_for_vf29.dat', 'partitioning_for_vf33.dat', 'partitioning_for_vf37.dat',
    'partitioning_for_vf41.dat', 'partitioning_for_vf46.dat', 'partitioning_for_vf51.dat',
    'partitioning_for_vf57.dat', 'partitioning_for_vf64.dat', 'partitioning_for_vf71.dat',
    'partitioning_for_vf79.dat'
]

phis = np.exp(np.linspace(np.log(0.03), np.log(0.08), 10))

press_file = 'pressures.txt'
if not os.path.isfile(press_file):
    raise FileNotFoundError(f"Pressure file not found: {press_file}")

press = np.loadtxt(press_file)
nu = 0.588
xis = (4.14 / (3 * nu - 1) / press) ** (1/3)
all_R_xi = []
all_F = []

for filename, xi, phi in zip(filenames, xis, phis):
    file_path = os.path.join('.', filename)
    if not os.path.isfile(file_path):
        print(f"File not found: {file_path}. Skipping.")
        continue
    data = np.loadtxt(file_path, delimiter=' ')
    if data.shape[1] < 2:
        print(f"Unexpected data format in {filename}. Skipping.")
        continue
    x = data[:18, 0]
    y = data[:18, 1]

    R_xi = x / xi
    F = -np.log(y) - (4 * np.pi / 3) / (3 * nu - 1) * (R_xi)**3 + np.log(1 - phi) 

    positive_indices = (R_xi > 0) & (F > 0)
    R_xi = R_xi[positive_indices]
    F = F[positive_indices]

    all_R_xi.extend(R_xi)
    all_F.extend(F)

all_R_xi = np.array(all_R_xi)
all_F = np.array(all_F)

if len(all_R_xi) == 0 or len(all_F) == 0:
    raise ValueError("No valid data available for plotting after filtering.")

def curve_func(R_xi, w1, w2):
    return func(R_xi, w1, w2)

params, covariance = curve_fit(curve_func, all_R_xi, all_F)
w1, w2 = params  # Extract fitted parameters

fig, ax = plt.subplots(figsize=(3, 2))

fig.patch.set_facecolor('white')  # Figure background
ax.set_facecolor('white')          # Axes background

for spine in ax.spines.values():
    spine.set_edgecolor('black')  # Spine color
    spine.set_linewidth(1)        # Spine width

ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

colors = plt.cm.viridis(np.linspace(0, 1, len(filenames)))
for i, (filename, xi, phi) in enumerate(zip(filenames[::-1], xis[::-1], phis[::-1])):
    file_path = os.path.join('.', filename)
    if not os.path.isfile(file_path):
        continue  # Already printed a message earlier
    data = np.loadtxt(file_path, delimiter=' ')
    if data.shape[1] < 2:
        continue  # Already printed a message earlier
    x = data[:18, 0]
    y = data[:18, 1]

    R_xi = x / xi
    F = -np.log(y) + np.log(1 - phi) - (4 * np.pi / 3) / (3 * nu - 1) * (R_xi)**3

    positive_indices = (R_xi > 0) & (F > 0)
    R_xi = R_xi[positive_indices]
    F = F[positive_indices]

    if len(R_xi) == 0 or len(F) == 0:
        continue  # No valid data to plot for this file

    ax.scatter((R_xi), (F), color=colors[i], marker='o', s=80, label=f'xi {xi:.1f} nm')

R_xi_fit = np.linspace(min(all_R_xi), max(all_R_xi), 500)
F_fit = func(R_xi_fit, w1, w2)
ax.plot((R_xi_fit), (F_fit), 'black', linewidth=2, label=f'Fitted Curve\n a1 = {w1:.2f}\n a2 = {w2:.2f}')
ax.set_ylim([0.11, 15])

ax.set_axisbelow(True)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.savefig('fitting_collapse', dpi = 500)
plt.show()
