import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

# Define the fitting function, now including xi as a parameter
def func2(params, R, xi):
    w1, w2 = params
    
    R_xi = R / xi  # Calculating R_xi within the function
    return w1 * R_xi**(3 - 1/0.588) + w2 * R_xi**2

def func3(params, R, xi):
    w1, w2 = params
    
    R_xi = R / xi  # Calculating R_xi within the function
    return w1 * R_xi**(3 - 1/0.588) + w2 * R_xi**2  + (4 * np.pi / 3) / (3 * 0.588 - 1) * R_xi**3





# List of filenames
filenames = [
    'partitioning_for_vf29.dat', 'partitioning_for_vf33.dat', 'partitioning_for_vf37.dat',
    'partitioning_for_vf41.dat', 'partitioning_for_vf46.dat', 'partitioning_for_vf51.dat',
    'partitioning_for_vf57.dat', 'partitioning_for_vf64.dat', 'partitioning_for_vf71.dat',
    'partitioning_for_vf79.dat'
]

# Calculate phis using exponential spacing
phis = np.exp(np.linspace(np.log(0.03), np.log(0.08), 10))

# Load pressures and calculate xis
press = np.loadtxt('pressures.txt')
nu = 0.588
xis = (4.14 / (3 * nu - 1) / press) ** (1/3)

# Prepare lists to collect all data
all_R = []
all_F = []
all_xis = []

# Read data from each file and store
for filename, xi, phi in zip(filenames, xis, phis):
    file_path = os.path.join('.', filename)
    data = np.loadtxt(file_path, delimiter=' ')
    R = data[:18, 0]
    y = data[:18, 1]
    F = -np.log(y) - (4 * np.pi / 3) / (3 * 0.588 - 1) * (R / xi)**3 + np.log(1 - phi) # If we only want to keep first two terms
    #F = -np.log(y) + np.log(1 - phi) # If we are keeping three terms

    all_R.extend(R)
    all_F.extend(F )
    all_xis.extend([xi] * len(R))

# Convert lists to numpy arrays for processing
all_R = np.array(all_R)
all_F = np.array(all_F)
all_xis = np.array(all_xis)

# Define the curve fitting function that only takes R and F as inputs
def curve_func(R, w1, w2):
    return func2((w1, w2), R, all_xis) # If we only want to keep first two terms
    #return func3((w1, w2), R, all_xis) # If we are keeping three terms

# Perform curve fitting
params, _ = curve_fit(curve_func, all_R, all_F)
print(params)
# Transform data to log scale

plt.rcdefaults()

# Create a plot with white background
fig, ax = plt.subplots(figsize=(8, 6))

# Set figure and axis background to white
fig.patch.set_facecolor('white')
ax.set_facecolor('white')

# Define a colormap
colors = plt.cm.viridis(np.linspace(0, 1, len(filenames)))

# Initialize a flag to label the fitted line only once
fitted_line_label_added = False

# Plot data and fitted curve for each file
for i, (filename, xi, phi) in enumerate(zip(filenames, xis, phis)):
    file_path = os.path.join('.', filename)
    data = np.loadtxt(file_path, delimiter=' ')
    R = data[:18, 0]
    y = data[:18, 1]
    F = -np.log(y) + np.log(1 - phi) # If we only want to keep first two terms
    #F = -np.log(y) + np.log(1 - phi) # If we are keeping three terms
    

    
    # Plot the fitted curve first (black solid line)
    R_fit = np.linspace(min(R), max(R), 500)
    #F_fit = func2(params, R_fit, xi) # If we only want to keep first two terms
    F_fit = func3(params, R_fit, xi) # If we are keeping three terms
    log_R_fit = (R_fit)
    log_F_fit = (F_fit)
    
    if not fitted_line_label_added:
        ax.plot(log_R_fit, log_F_fit, color='black', linestyle='-', linewidth=1, label='Fitted Line')
        fitted_line_label_added = True
    else:
        ax.plot(log_R_fit, log_F_fit, color='black', linestyle='-', linewidth=1, zorder = -1)
    
    # Scatter plot for the transformed data on top of the fitted line
    ax.scatter(R, F, color=colors[i], marker='o',s=80, label=f'ξ {xi:.2f} nm')

# Finalizing the plot
ax.set_title(f'Fitted Curves with a1 = {params[0]:.3f}; a2 = {params[1]:.3f}', fontsize=14)
ax.set_xlabel('log(R) (log nm)', fontsize=12)
ax.set_ylabel('log(F) (log kT)', fontsize=12)

# Enable grid with specific settings
# ax.grid(True, linestyle='--', linewidth=0.5, color='grey', which='both', alpha=0.5)

# Ensure grid lines are behind the data points
ax.set_axisbelow(True)
plt.xscale('log')
plt.yscale('log')
# Add legend
ax.legend(title='Correlation Lengths', fontsize=8, title_fontsize=11)

# Adjust layout for better spacing
plt.tight_layout()
#plt.savefig('fitting_saw_3', dpi=1000)
# Display the plot
plt.show()