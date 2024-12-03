import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit

# Reset Matplotlib to default settings
plt.rcdefaults()
#plt.style.use('default')
def scientific_formatter(x, p):
    if x == 0:
        return '0'
    exp = np.floor(np.log10(x))
    coef = x / 10**exp

    return f'${coef:.1f}\\times10^{{{int(exp)}}}$'


plt.rcParams.update({'font.size': 12})
# Define the fitting function with two parameters
def func(R_xi, w1, w2):
    """
    Fitting function with two weights.

    Parameters:
    - R_xi: Dimensionless ratio R/xi
    - w1, w2: Fitting parameters

    Returns:
    - F: Calculated force value based on the function
    """
    nu = 0.588
    return w1 * (R_xi)**(3 - 1/nu) + w2 * (R_xi)**2 # + (4 * np.pi / 3) / (3 * nu - 1) * (R_xi)**3

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
press_file = 'pressures.txt'
if not os.path.isfile(press_file):
    raise FileNotFoundError(f"Pressure file not found: {press_file}")

press = np.loadtxt(press_file)
nu = 0.588
xis = (4.14 / (3 * nu - 1) / press) ** (1/3)

# Prepare for collecting all data
all_R_xi = []
all_F = []

# Read data from each file and store
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

    # Filter out non-positive values to avoid log scaling issues
    positive_indices = (R_xi > 0) & (F > 0)
    R_xi = R_xi[positive_indices]
    F = F[positive_indices]

    all_R_xi.extend(R_xi)
    all_F.extend(F)

# Convert lists to numpy arrays for processing
all_R_xi = np.array(all_R_xi)
all_F = np.array(all_F)

# Check if data arrays are not empty
if len(all_R_xi) == 0 or len(all_F) == 0:
    raise ValueError("No valid data available for plotting after filtering.")

# Define the curve fitting function that takes R and xi into account
def curve_func(R_xi, w1, w2):
    return func(R_xi, w1, w2)

# Perform curve fitting with two parameters
params, covariance = curve_fit(curve_func, all_R_xi, all_F)
w1, w2 = params  # Extract fitted parameters

# Create a plot with white background
fig, ax = plt.subplots(figsize=(3, 2))

# Set figure and axis background to white
fig.patch.set_facecolor('white')  # Figure background
ax.set_facecolor('white')          # Axes background

# Explicitly set spine colors to black and ensure they are visible
for spine in ax.spines.values():
    spine.set_edgecolor('black')  # Spine color
    spine.set_linewidth(1)        # Spine width

# Ensure all spines are visible
ax.spines['top'].set_visible(True)
ax.spines['right'].set_visible(True)
ax.spines['bottom'].set_visible(True)
ax.spines['left'].set_visible(True)

# Define a colormap
colors = plt.cm.viridis(np.linspace(0, 1, len(filenames)))
# 
# Plot data from each file with larger markers
for i, (filename, xi, phi) in enumerate(zip(filenames, xis, phis)):
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

    # Filter out non-positive values
    positive_indices = (R_xi > 0) & (F > 0)
    R_xi = R_xi[positive_indices]
    F = F[positive_indices]

    if len(R_xi) == 0 or len(F) == 0:
        continue  # No valid data to plot for this file

    # Scatter plot for the data with larger markers and black edges
    ax.scatter((R_xi), (F), color=colors[i], marker='o', s=80, label=f'xi {xi:.1f} nm')

# Plot the fitted curve over a range
R_xi_fit = np.linspace(min(all_R_xi), max(all_R_xi), 500)
F_fit = func(R_xi_fit, w1, w2)
ax.plot((R_xi_fit), (F_fit), 'black', linewidth=2, label=f'Fitted Curve\n a1 = {w1:.2f}\n a2 = {w2:.2f}')

# Finalizing the plot
#ax.set_title(f'Fitted Curves with a1 = {w1:.3f}; a2 = {w2:.3f}', fontsize=16)
#ax.set_xlabel('R/ξ (dimensionless)', fontsize=14)
#ax.set_ylabel('F (kT)', fontsize=14)

# Enable grid with specific settings
#ax.grid(True, linestyle='--', linewidth=0.5, color='grey', which='both', alpha=0.5)

#ax.minorticks_off()
#log_yticks = np.linspace(np.log10(0.11), np.log10(10), 11)  # log values from -1 to 1
#ytick_values = 10**log_yticks  # Convert to actual values
#ytick_labels = [f'$10^{0}$' if abs(val) < 1e-1 else f'$10^{{{val:.1f}}}$' for val in log_yticks]  # Format as 10^x

#ax.set_yticks(ytick_values)
#ax.set_yticklabels(ytick_labels)
ax.set_ylim([0.11, 15])



# Ensure grid lines are behind the data points
ax.set_axisbelow(True)
plt.xscale('log')
plt.yscale('log')
# Add legend with better placement and larger font size
#ax.legend(title='Correlation Lengths', fontsize=12, title_fontsize=14, loc='best')

#log_xticks = np.linspace(np.log10(0.1), np.log10(1), 3)  # Adjust range as needed
#xtick_values = 10**log_xticks
#xtick_labels = [f'$10^{0}$' if abs(val) < 1e-2 else f'$10^{{{val:.1f}}}$' for val in log_xticks]
#ax.set_xticks(xtick_values)
#ax.set_xticklabels(xtick_labels)

# Adjust layout for better spacing
plt.tight_layout()
plt.savefig('fitting_collapse', dpi = 500)
# Display the plot
plt.show()