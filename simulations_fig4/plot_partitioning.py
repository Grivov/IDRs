import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# Set the color for the plots
orng = '#FF8000'

# Initialize lists to store the data
x_values = []
log_y_values = []

# Iterate through all files matching the pattern
for file_path in glob.glob('partitioning_for_vf*'):
    # Load the data from the file
    with open(file_path, 'r') as f:
        data = np.loadtxt(f, delimiter=' ')
    x = data[:, 0]
    y = data[:, 1]

    # Compute the log of partitioning values
    log_y = -np.log(y)

    # Append the data to the lists
    x_values.append(x)
    log_y_values.append(log_y)

# Create the plot
fig = plt.figure(figsize=(12, 6))

# Plot for the logarithmic transformation
for i, log_y in enumerate(log_y_values):
    plt.loglog(x_values[i], log_y, color=orng, marker='o', markersize=12, linewidth=3, label=f'IDP {i+1}')
plt.xlabel('R, nm', fontsize=16)
plt.ylabel('F, kT', fontsize=16)
#plt.tick_params(axis='both', labelsize=16)
plt.legend(fontsize=16)

#plt.savefig('all_partitioning_log_plots.png', dpi=600)
plt.show()
