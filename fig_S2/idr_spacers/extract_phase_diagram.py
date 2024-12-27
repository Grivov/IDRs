import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, minimize
from scipy.interpolate import UnivariateSpline

    # Setting plot aesthetics
plt.rcParams['font.size'] = 20  # Increase font size for all text
plt.rcParams['axes.titlesize'] = 20  # Font size for titles
plt.rcParams['axes.labelsize'] = 20  # Font size for labels
plt.rcParams['xtick.labelsize'] = 20  # Font size for X tick labels
plt.rcParams['ytick.labelsize'] = 20  # Font size for Y tick labels
    

def process_volume_fraction(file_name, x):
    """
    Processes a volume fraction data file, symmetrizes it, and calculates averages.
    
    Parameters:
    - file_name: str, path to the input text file containing data.
    - x: int, threshold index for calculating c_dense and c_dilute.
    
    Returns:
    - c_dense: float, mean concentration in the dense phase.
    - c_dilute: float, mean concentration in the dilute phase.
    """
    volume_fraction = np.loadtxt(file_name)
    volume_fraction_sym = (volume_fraction + volume_fraction[::-1]) / 2
    c_dense = np.mean(volume_fraction_sym[int(25 - x):25])   # From index x to 25
    c_dilute = np.mean(volume_fraction_sym[:x])              # From start to index x
    return c_dense, c_dilute

def fit_critical_point(T_values, c_dense_values, c_dilute_values):
    """
    Fits the critical point (T_c, c_c) using the law of coexistence densities and rectilinear diameters.
    
    Parameters:
    - T_values: array-like, temperatures.
    - c_dense_values: array-like, dense phase concentrations.
    - c_dilute_values: array-like, dilute phase concentrations.
    
    Returns:
    - critical_params: tuple, containing T_c, c_c, d, and A.
    """
    # Convert concentrations to densities for clarity
    rho_high = np.array(c_dense_values)
    rho_low = np.array(c_dilute_values)
    T = np.array(T_values)

    # Define the coexistence density function (Eq. 6)
    def coexistence_func(T, T_c, d):
        return (rho_high - rho_low) ** 3.06 - d * (1 - T / T_c)

    # Define the rectilinear diameter function (Eq. 7)
    def rectilinear_func(T, T_c, c_c, A):
        return (rho_high + rho_low) / 2 - c_c - A * (T - T_c)

    # Initial guesses
    T_c_guess = np.mean(T)
    c_c_guess = np.mean((rho_high + rho_low) / 2)
    d_guess = 1.0
    A_guess = 0.1

    # Fit the coexistence density function
    popt_coexistence, _ = curve_fit(lambda T, T_c, d: coexistence_func(T, T_c, d), T, np.zeros_like(T),
                                    p0=[T_c_guess, d_guess])

    # Use the fitted T_c to fit the rectilinear diameter function
    T_c_fit, d_fit = popt_coexistence
    popt_rectilinear, _ = curve_fit(lambda T, c_c, A: rectilinear_func(T, T_c_fit, c_c, A), T, np.zeros_like(T),
                                    p0=[c_c_guess, A_guess])

    c_c_fit, A_fit = popt_rectilinear
    return T_c_fit, c_c_fit, d_fit, A_fit

def plot_phase_diagram(T_values, c_dense_values, c_dilute_values, T_c, c_c):
    """
    Plots the phase diagram using a UnivariateSpline for smoothing with enhanced visuals.

    Parameters:
    - T_values: array-like, temperatures.
    - c_dense_values: array-like, dense phase concentrations.
    - c_dilute_values: array-like, dilute phase concentrations.
    - T_c: float, critical temperature.
    - c_c: float, critical concentration.
    """
    # Combine dense, dilute, and critical points
    concentrations = np.concatenate([c_dense_values, c_dilute_values, [c_c]])
    temperatures = np.concatenate([T_values, T_values, [T_c]])

    # Sort points by concentration
    sort_idx = np.argsort(concentrations)
    x = concentrations[sort_idx]
    y = temperatures[sort_idx]

    # Create a UnivariateSpline
    spl = UnivariateSpline(x, y, k=2, s=1)

    # Generate points on the spline
    x_new = np.linspace(x.min(), x.max(), 500)
    y_new = spl(x_new)

    # Plotting
    plt.figure(figsize=(10, 10))

    # Plot dense phase points
    plt.scatter(
        c_dense_values, T_values, 
        color='lightcoral', s=200, label='Dense Phase'
    )
    
    # Plot dilute phase points
    plt.scatter(
        c_dilute_values, T_values, 
        color='lightskyblue', s=200, label='Dilute Phase'
    )
    
    # Plot critical point
    plt.scatter(
        [c_c], [T_c], 
        color='gold', s=200, label='Critical Point'
    )

    # Plot the spline curve
    plt.plot(x_new, y_new, 'k--', label='Phase Boundary (Spline)', linewidth=3)

    # Labels, legend, and aesthetics
    plt.xlabel('Volume fraction', fontsize=20)
    plt.ylabel('T (K)', fontsize=20)
    plt.title('Phase Diagram with Smoothed Spline', fontsize=20)
    plt.legend(fontsize=20)
    plt.tight_layout()
    plt.savefig('PD_IDR', dpi = 500)
    plt.show()

def main():
    temperatures = [300, 308, 310, 315]  # List of temperatures
    x = 5  # Threshold index for dilute and dense calculations
    dataset = []  # To store (T, c_dense) and (T, c_dilute) pairs

    # Loop through files for each temperature
    for T in temperatures:
        file_name = f'volume_frac_pol_T{T}_IDR.txt'
        print(f"Processing file: {file_name}")
        c_dense, c_dilute = process_volume_fraction(file_name, x)
        if c_dense is not None and c_dilute is not None:
            dataset.append((T, c_dense))   # Dense phase
            dataset.append((T, c_dilute)) # Dilute phase

    # Convert dataset to structured NumPy array
    dataset = np.array(dataset, dtype=[('Temperature', float), ('Concentration', float)])

    # Extract data for fitting
    T_values = dataset['Temperature'][::2]  # Unique temperatures
    c_dense_values = dataset['Concentration'][::2]  # Dense phase
    c_dilute_values = dataset['Concentration'][1::2]  # Dilute phase

    # Fit critical point
    T_c, c_c, d, A = fit_critical_point(T_values[1:], c_dense_values[1:], c_dilute_values[1:])

    # Print critical point results
    print(f"\nEstimated Critical Point:")
    print(f"T_c = {T_c:.2f} K")
    print(f"c_c = {c_c:.5f}")
    print(f"d = {d:.5f}")
    print(f"A = {A:.5f}")

    plot_phase_diagram(T_values, c_dense_values, c_dilute_values, T_c, c_c)

if __name__ == "__main__":
    main()
