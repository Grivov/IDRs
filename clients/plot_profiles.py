import numpy as np
import matplotlib.colorbar as colorbar
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib.cm as cm
import os 

# Constants
RADIUS = 1.5
V2 = 4 * np.pi / 3 * RADIUS**3 

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
    
    return np.array(densities)

def calculate_average_profile(densities):
    """Calculate average density profile from raw density data."""
    average_density = np.mean(densities.reshape(-1, 50), axis=0)
    average_vf = 0.5 * (average_density + average_density[::-1]) * V2
    return average_vf#[:25]

def plot_profile(epsilon, replica_data, color, ax):
    """Fit curve to averaged data and plot results."""
    y_data_avg = np.mean(replica_data, axis=0)
    x_values = np.linspace(-50, 50, 50)    
    # Plot raw data
    ax.plot(x_values, y_data_avg, color=color, linewidth=7, label=f"epsilon {epsilon}")
    ax.scatter(x_values, y_data_avg, color=color, s=200, zorder = 5)
    
    # Fit curve and plot
    #popt, _ = curve_fit(tanh_func, np.arange(25), y_data_avg, p0=[0.1, 0.1, 0.1, 12])
    #A, B, C, D = popt
    #fitted_curve = tanh_func(np.arange(25), A, B, C, D)
    #ax.plot(fitted_curve, '--', color=color, alpha=0.7)
    
    # Calculate and return metrics
    density_ratio = np.mean(y_data_avg[18:25]) / np.mean(y_data_avg[:7])
    #curve_ratio = (C + A) / (C - A)
    #print(C+A)
    
    return density_ratio

def analyze_density_profiles():
    """Main function to analyze density profiles across different epsilon values."""
    # Setup plotting
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots(figsize=(20, 5))
    cmap = cm.get_cmap('viridis', 10)
    
    # Process data for each epsilon value
    for idx, epsilon in enumerate(range(1,11)):
        replica_profiles = []
        color = cmap(idx)
        
        # Process all replicas for current epsilon
        for replica in range(1, 26):
            file_path = f"density_clients_eps{epsilon}_rep{replica}.dat"
            if os.path.exists(file_path):
                densities = process_density_file(file_path)
                profile = calculate_average_profile(densities)
                replica_profiles.append(profile)
    
        
        # Fit and plot averaged data
        #density_ratio, curve_ratio = fit_and_plot_profile(
        #    epsilon, 
        #    replica_profiles, 
        #    color, 
        #    ax
        #)
        density_ratio = plot_profile(
            epsilon, 
            replica_profiles, 
            color, 
            ax
        )
        
        #print(f"Epsilon {epsilon}:")
        print(f"  Density ratio: {density_ratio:.3f}")
        #print(f"  Curve ratio: {curve_ratio:.3f}")
    print('done') 
    # Finalize plot
    ax.set_xlabel("Position")
    ax.set_ylabel("Density")
    #cbar = plt.colorbar()
    #cbar.set_label('Epsilon', rotation=270, labelpad=15)

    #plt.savefig('clients_profiles.png', dpi=600, bbox_inches='tight')
    #plt.show()

def analyze_density_profiles():
    """Main function to analyze density profiles across different epsilon values."""
    # Setup plotting
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots(figsize=(25, 5))
    cmap = cm.get_cmap('viridis', 10)

    # Normalize the epsilon values for the color bar
    norm = plt.Normalize(vmin=1, vmax=10)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # Dummy array for colorbar

    # Process data for each epsilon value
    for idx, epsilon in enumerate(range(1, 11)):
        replica_profiles = []
        color = cmap(norm(epsilon))  # Map epsilon to color

        # Process all replicas for current epsilon
        for replica in range(1, 26):
            file_path = f"density_clients_eps{epsilon}_rep{replica}.dat"
            if os.path.exists(file_path):
                densities = process_density_file(file_path)
                profile = calculate_average_profile(densities)
                replica_profiles.append(profile)

        # Fit and plot averaged data
        density_ratio = plot_profile(
            epsilon,
            replica_profiles,
            color,
            ax
        )
        
        print(f"  Density ratio: {density_ratio:.3f}")
        
    # Finalize plot with labels and colorbar
    ax.set_xlabel("Position")
    ax.set_ylabel("Density")
    cbar = fig.colorbar(sm, ax=ax, orientation="vertical")

    # Save and display the plot
    #plt.savefig('clients_profiles.png', dpi=600, bbox_inches='tight')
    plt.show()


analyze_density_profiles()
