import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

# Function to fit
def func(R_xi, w1, w2, w3):
    return w1 * (R_xi)**(3 - 1/0.588) + w2 * (R_xi)**2 + w3 * (R_xi)**3

filenames = ['partitioning_for_vf29.dat', 'partitioning_for_vf33.dat', 'partitioning_for_vf37.dat', 'partitioning_for_vf41.dat', 'partitioning_for_vf46.dat',
             'partitioning_for_vf51.dat', 'partitioning_for_vf57.dat', 'partitioning_for_vf64.dat', 'partitioning_for_vf71.dat', 'partitioning_for_vf79.dat']

phis = np.exp(np.linspace(np.log(0.03), np.log(0.08), 10))

press = np.loadtxt('/home/vgrigor2/scratch4-yzhan567/vgrigor/idr/figure3_boxes/multiple_boxes_3to8/pressures.txt')
nu = 0.588
press = np.array(press)
xis = (4.14 / (3 * nu - 1) / press) ** (1/3)


# Prepare for collecting all data
all_R_xi = []
all_F = []

# Read data from each file and store
for filename, xi, phi in zip(filenames, xis, phis):
    file_path = os.path.join('.', filename)
    data = np.loadtxt(file_path, delimiter=' ')
    x = data[:25, 0]
    y = data[:25, 1]

    R_xi = x / xi
    F = -np.log(y) + np.log(1-phi)

    all_R_xi.extend(R_xi)
    all_F.extend(F)

# Fit the model to all data points collectively
params, _ = curve_fit(func, all_R_xi, all_F)

# Create a plot
fig, ax = plt.subplots(figsize=(12, 6))
colors = plt.cm.tab10(np.linspace(0, 1, len(filenames)))

# Plot data from each file with fitting
for i, (filename, xi, phi) in enumerate(zip(filenames, xis, phis)):
    file_path = os.path.join('.', filename)
    data = np.loadtxt(file_path, delimiter=' ')
    x = data[:25, 0]
    y = data[:25, 1]

    R_xi = x / xi
    F = -np.log(y) + np.log(1-phi)

    ax.plot(R_xi, F, color=colors[i], marker='o', label=f'xi {xi:.1f} nm')

# Plot the fitted curve over a range
R_xi_fit = np.linspace(min(all_R_xi), max(all_R_xi), 500)
F_fit = func(R_xi_fit, *params)
ax.loglog(R_xi_fit, F_fit, 'black', linewidth=2, label=f'Fitted Curve,\n a1 = {params[0]:.2f}\n a2 = {params[1]:.2f},\n a3 = {params[2]:.2f}')

# Finalizing the plot
ax.set_xlabel('R, nm', fontsize=16)
ax.set_ylabel('F, kT', fontsize=16)
ax.legend(title='Correlation lengths', fontsize=12)
plt.show()
