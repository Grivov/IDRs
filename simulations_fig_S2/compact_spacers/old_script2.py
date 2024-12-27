import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Constants
V2 = 4*np.pi/3 

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

    return np.array(densities) if densities else None

def calculate_average_profile(densities):
    """Calculate average density profile from raw density data."""
    if densities is None:
        return None

    average_density = np.mean(densities.reshape(-1, 50), axis=0)
    average_vf = 0.5 * (average_density + average_density[::-1]) * V2
    return average_vf[:25]

def plot_profile(replica_data, ax):
    """Fit curve to averaged data and plot results."""
    y_data_avg = np.mean(replica_data, axis=0)
    print(y_data_avg)
    x_values = np.arange(25)
    
    # Plot raw data
    ax.plot(x_values, y_data_avg, linewidth=7)
    ax.scatter(x_values, y_data_avg, s=200, zorder=5)

    # Fit curve and plot
    popt, _ = curve_fit(tanh_func, np.arange(25), y_data_avg, p0=[0.1, 0.1, 0.1, 12])
    A, B, C, D = popt

    fitted_curve = tanh_func(np.arange(25), A, B, C, D)
    ax.plot(fitted_curve, '--')

    edge1, edge2 = int(D - (2/B)), int(D + (2/B))
    print("edge1, edge2", edge1, edge2)
    density1, density2 = np.mean(y_data_avg[:2]), np.mean(y_data_avg[23:25])
    print(density1, density2)

    return density1, density2

def fit_connected_parabolas(x_data, y_data, x_split=None):
    """
    Fit two parabolas that connect at their maximum point.
    For parabolas: y = a*x^2 + b*x + c
    Connection point has zero derivative (2ax + b = 0)
    """
    # Sort data by x values
    sort_idx = np.argsort(x_data)
    x_data = x_data[sort_idx]
    y_data = y_data[sort_idx]
    
    # Split data into left and right sides
    middle_idx = len(x_data) // 2
    left_mask = np.arange(len(x_data)) < middle_idx
    right_mask = ~left_mask
    
    # Initial separate fits to get starting parameters
    def fit_side(x, y):
        popt, _ = curve_fit(lambda x, a, b, c: a*x**2 + b*x + c, x, y)
        return popt
    
    left_params = fit_side(x_data[left_mask], y_data[left_mask])
    right_params = fit_side(x_data[right_mask], y_data[right_mask])
    
    def two_parabolas(x, a1, b1, c1, a2, b2):
        """
        Two parabolas that meet at x_connect where both have zero derivative.
        Left parabola: y = a1*x^2 + b1*x + c1
        Right parabola: y = a2*x^2 + b2*x + c2
        """
        # Connection points (where derivatives are zero)
        x_connect_left = -b1/(2*a1)
        x_connect_right = -b2/(2*a2)
        
        # Calculate c2 to ensure continuity
        y_connect = a1*x_connect_left**2 + b1*x_connect_left + c1
        c2 = y_connect - (a2*x_connect_right**2 + b2*x_connect_right)
        
        # Return piecewise function
        y = np.zeros_like(x)
        left_side = x <= x_connect_left
        right_side = x > x_connect_left
        
        y[left_side] = a1*x[left_side]**2 + b1*x[left_side] + c1
        y[right_side] = a2*x[right_side]**2 + b2*x[right_side] + c2
        
        return y
    
    # Initial guess for parameters
    p0 = [left_params[0], left_params[1], left_params[2],
          right_params[0], right_params[1]]
    
    # Fit the connected parabolas
    popt, _ = curve_fit(two_parabolas, x_data, y_data, p0=p0)
    
    # Return the fitting function with optimal parameters
    def fitted_curve(x):
        return two_parabolas(x, *popt)
    
    # Calculate connection point
    x_connect = -popt[1]/(2*popt[0])
    
    return fitted_curve, x_connect, popt

def analyze_density_profiles():
    # First figure: density profiles
    fig, ax = plt.subplots(figsize=(20, 5))
    den_densities = []
    dil_densities = []
    
    temperatures = [290, 295, 305, 307, 308]
    for T in temperatures:
        replica_profiles = []

        for replica in range(1, 10):
            file_path = f"Out_profile_{T}/density_pols_T{T}_rep{replica}.dat"
            if os.path.exists(file_path):
                densities = process_density_file(file_path)
                if densities is not None:
                    profile = calculate_average_profile(densities)
                    if profile is not None:
                        replica_profiles.append(profile)

        if replica_profiles:
            density1, density2 = plot_profile(
                replica_profiles,
                ax
            )
            den_densities.append(density2)
            dil_densities.append(density1)
    
    ax.set_xlabel("Position")
    ax.set_ylabel("Density")
    plt.show()

    # Second figure: phase diagram with connected parabolas
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Convert to arrays for fitting
    den_densities = np.array(den_densities)
    dil_densities = np.array(dil_densities)
    temperatures = np.array(temperatures)
    
    # Combine data for fitting
    x_data = np.concatenate([dil_densities, den_densities])
    y_data = np.concatenate([temperatures, temperatures])
    
    # Fit connected parabolas
    fitted_curve, x_connect, params = fit_connected_parabolas(x_data, y_data)
    
    # Generate points for smooth curve
    x_smooth = np.linspace(min(x_data), max(x_data), 200)
    y_smooth = fitted_curve(x_smooth)
    
    # Plot experimental points
    ax.scatter(den_densities, temperatures, s=300, color='blue', 
              label='Dense phase')
    ax.scatter(dil_densities, temperatures, s=300, color='red',
              label='Dilute phase')
    
    # Plot fitted curve
    ax.plot(x_smooth, y_smooth, '--', color='black', label='Fitted curve')
    
    # Plot connection point
    y_connect = fitted_curve(x_connect)
    ax.scatter([x_connect], [y_connect], color='green', s=200,
              label='Maximum point')
    
    ax.set_xlabel('Density')
    ax.set_ylabel('Temperature (K)')
    ax.legend()
    plt.show()
    np.savetxt('compact_phase_diagram.dat', np.vstack([x_data,y_data])) 
    return fitted_curve, x_connect, params

if __name__ == "__main__":
    analyze_density_profiles()
