import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

# Define the fitting function, now including xi as a parameter
def func(params, R, xi):
    w1, w2, w3 = params
    R_xi = R / xi  # Calculating R_xi within the function
    return w1 * R_xi**(3 - 1/0.588) + w2 * R_xi**2 + w3 * R_xi**3

filenames = ['partitioning_for_vf29.dat', 'partitioning_for_vf33.dat', 'partitioning_for_vf37.dat', 'partitioning_for_vf41.dat', 'partitioning_for_vf46.dat',
             'partitioning_for_vf51.dat', 'partitioning_for_vf57.dat', 'partitioning_for_vf64.dat', 'partitioning_for_vf71.dat', 'partitioning_for_vf79.dat']

phis = np.exp(np.linspace(np.log(0.03), np.log(0.08), 10))

press = np.loadtxt('/home/vgrigor2/scratch4-yzhan567/vgrigor/idr/figure3_boxes/multiple_boxes_3to8/pressures.txt')
nu = 0.588
press = np.array(press)
xis = (4.14 / (3 * nu - 1) / press) ** (1/3)


# Prepare for collecting all data
all_R = []
all_F = []
all_xis = []

# Read data from each file and store
for filename, xi, phi in zip(filenames, xis, phis):
    file_path = os.path.join('.', filename)
    data = np.loadtxt(file_path, delimiter=' ')
    R = data[:25, 0]
    y = data[:25, 1]
    F = -np.log(y) + np.log(1-phi)

    all_R.extend(R)
    all_F.extend(F)
    all_xis.extend([xi] * len(R))

# Prepare data for curve fitting
all_R = np.array(all_R)
all_F = np.array(all_F)
all_xis = np.array(all_xis)

# Function for the curve_fit to call, must only take R, F as parameters
def curve_func(R, w1, w2, w3):
    return func((w1, w2, w3), R, all_xis)

# Fit the model to all data points collectively with the same weights
params, _ = curve_fit(curve_func, all_R, all_F)

# Create a plot
fig, ax = plt.subplots(figsize=(12, 6))
colors = plt.cm.tab10(np.linspace(0, 1, len(filenames)))

# Plot data and fitted curve for each file
for i, (filename, xi, phi) in enumerate(zip(filenames, xis, phis)):
    file_path = os.path.join('.', filename)
    data = np.loadtxt(file_path, delimiter=' ')
    R = data[:25, 0]
    y = data[:25, 1]
    F = -np.log(y) + np.log(1-phi)

    ax.scatter(R, F, color=colors[i], marker='o', label=f'xi {xi:.1f} nm')

    # Plot the fitted curve
    R_fit = np.linspace(min(R), max(R), 500)
    F_fit = func(params, R_fit, xi)
    ax.loglog(R_fit, F_fit, color=colors[i], linewidth=1, linestyle='--')

# Finalizing the plot
ax.set_title('Fitted curves with shared weights')
ax.set_xlabel('R (nm)', fontsize=16)
ax.set_ylabel('F (kT)', fontsize=16)
ax.legend(title='Correlation lengths', fontsize=12)
plt.show()
