import numpy as np
import matplotlib.colorbar as colorbar
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.cm as cm
import os

# Constants
RADIUS = 1
V2 = 4 * np.pi / 3 * RADIUS**3

def plot_profile(replica_data, ax):
    """Fit curve to averaged data and plot results."""
    y_data_avg = np.mean(replica_data, axis=0)
    print(y_data_avg)
    x_values = np.arange(25)
    # Plot raw data
    ax.plot(x_values, y_data_avg, linewidth=7)
    ax.scatter(x_values, y_data_avg, s=200, zorder = 5)

    # Fit curve and plot
    popt, _ = curve_fit(tanh_func, np.arange(25), y_data_avg, p0=[0.1, 0.1, 0.1, 12])
    A, B, C, D = popt
    #print(D, 1/B)

    fitted_curve = tanh_func(np.arange(25), A, B, C, D)
    ax.plot(fitted_curve, '--')

    edge1, edge2 = int(D - (2/B)), int(D + (2/B))
    print("edge1, edge2", edge1, edge2)
    #density1, density2 = np.mean(y_data_avg[:edge1]), np.mean(y_data_avg[edge2:25])
    density1, density2 = C-A, C+A
    print(density1, density2)

    return density1, density2

def tanh_func(x, A, B, C, D):
    """Hyperbolic tangent function for curve fitting."""
    return A * np.tanh(B * (x - D)) + C

def process_density_file(input_file):
    """Process a single density data file and return the density values."""
    densities = []
    
    with open(input_file, 'r') as file:
        for line in file:
            split_line = line.strip().split()
            if len(split_line) == 4:
                try:
                    density = float(split_line[3])
                    densities.append(density)
                except ValueError:
                    continue
    
    # Return None if file was empty or had no valid data
    return np.array(densities) if densities else None

def calculate_average_profile(densities):
    """Calculate average density profile from raw density data."""
    # Skip if densities is None (empty file)
    if densities is None:
        return None
        
    average_density = np.mean(densities.reshape(-1, 50), axis=0)
    average_vf = 0.5 * (average_density + average_density[::-1]) * V2
    return average_vf[:25]

def analyze_density_profiles():
    fig, ax = plt.subplots(figsize=(20, 5))
    den_densities = []
    dil_densities = []
    # Process data for each T
    temperatures = [290, 295, 305]
    for T in temperatures:
        replica_profiles = []

        for replica in range(1, 10):
            file_path = f"Out_profile_{T}/density_pols_T{T}_rep{replica}.dat"
            if os.path.exists(file_path):
                densities = process_density_file(file_path)
                # Only add profile if densities is not None
                if densities is not None:
                    profile = calculate_average_profile(densities)
                    if profile is not None:
                        replica_profiles.append(profile)

        # Only proceed with plot if we have valid profiles
        if replica_profiles:
            density1, density2 = plot_profile(
                replica_profiles,
                ax
            )
        den_densities.append(density2)
        dil_densities.append(density1)
    print('done')
    ax.set_xlabel("Position")
    ax.set_ylabel("Density")
    plt.show()
    fig, ax = plt.subplots(figsize=(10, 10))
    ax.scatter(den_densities, temperatures, s=300)
    ax.scatter(dil_densities, temperatures, s=300)
    plt.show()
    
    

analyze_density_profiles()
