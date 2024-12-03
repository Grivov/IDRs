import numpy as np
import matplotlib.pyplot as plt
import argparse

residues = ['C', 'M', 'F', 'I', 'L', 'V', 'W', 'Y', 'A', 'G',
            'T', 'S', 'N', 'Q', 'D', 'E', 'H', 'R', 'K', 'P']

sigma_array = np.array([5.48, 5.83, 5.92, 5.83, 5.83,
                        5.67, 6.13, 5.97, 5.26, 4.99,
                        5.75, 5.58, 5.33, 5.55, 5.53,
                        5.70, 5.78, 6.02, 5.92, 5.52,
                        6.18, 6.27, 6.18, 6.18, 6.02,
                        6.48, 6.32, 5.61, 5.34, 5.90,
                        5.68, 5.93, 6.10, 5.88, 6.05,
                        6.13, 6.37, 6.27, 5.87, 6.36,
                        6.27, 6.27, 6.11, 6.57, 6.41,
                        5.70, 5.43, 5.99, 5.77, 6.02,
                        6.19, 5.97, 6.14, 6.22, 6.46,
                        6.36, 5.96, 6.18, 6.18, 6.02,
                        6.48, 6.32, 5.61, 5.34, 5.90,
                        5.68, 5.93, 6.10, 5.88, 6.05,
                        6.13, 6.37, 6.27, 5.87, 6.18,
                        6.02, 6.48, 6.32, 5.61, 5.34,
                        5.90, 5.68, 5.93, 6.10, 5.88,
                        6.05, 6.13, 6.37, 6.27, 5.87,
                        5.86, 6.32, 6.16, 5.45, 5.18,
                        5.74, 5.52, 5.77, 5.94, 5.72,
                        5.89, 5.97, 6.21, 6.11, 5.71,
                        6.78, 6.62, 5.91, 5.64, 6.20,
                        5.98, 6.23, 6.40, 6.18, 6.35,
                        6.43, 6.67, 6.57, 6.17, 6.46,
                        5.75, 5.48, 6.04, 5.82, 6.07,
                        6.24, 6.02, 6.19, 6.27, 6.51,
                        6.41, 6.01, 5.04, 4.77, 5.33,
                        5.11, 5.36, 5.53, 5.31, 5.48,
                        5.56, 5.80, 5.70, 5.30, 4.50,
                        5.06, 4.84, 5.09, 5.26, 5.04,
                        5.21, 5.29, 5.53, 5.43, 5.03,
                        5.62, 5.40, 5.65, 5.82, 5.60,
                        5.77, 5.85, 6.09, 5.99, 5.59,
                        5.18, 5.43, 5.60, 5.38, 5.55,
                        5.63, 5.87, 5.77, 5.37, 5.68,
                        5.85, 5.63, 5.80, 5.88, 6.12,
                        6.02, 5.62, 6.02, 5.80, 5.97,
                        6.05, 6.29, 6.19, 5.79, 5.58,
                        5.75, 5.83, 6.07, 5.97, 5.57,
                        5.92, 6.00, 6.24, 6.14, 5.74,
                        6.08, 6.32, 6.22, 5.82, 6.56,
                        6.46, 6.06, 6.36, 5.96, 5.56])

qi_qj_array = np.array([[0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0],
                    [0, 0], [0, 0], [0, 0], [0, 0], [0, -1], [0, -1], [0, 0.25], [0, 1], [0, 1],
                    [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0],
                    [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, -1], [0, -1], [0, 0.25], [0, 1],
                    [0, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0],
                    [0, 0], [0, 0], [0, 0], [0, 0], [0, -1], [0, -1], [0, 0.25], [0, 1], [0, 1],
                    [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0],
                    [0, 0], [0, 0], [0, -1], [0, -1], [0, 0.25], [0, 1], [0, 1], [0, 0], [0, 0],
                    [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, -1],
                    [0, -1], [0, 0.25], [0, 1], [0, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0],
                    [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, -1], [0, -1], [0, 0.25], [0, 1],
                    [0, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0],
                    [0, -1], [0, -1], [0, 0.25], [0, 1], [0, 1], [0, 0], [0, 0], [0, 0], [0, 0],
                    [0, 0], [0, 0], [0, 0], [0, 0], [0, -1], [0, -1], [0, 0.25], [0, 1], [0, 1],
                    [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, -1], [0, -1],
                    [0, 0.25], [0, 1], [0, 1], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0], [0, 0],
                    [0, -1], [0, -1], [0, 0.25], [0, 1], [0, 1], [0, 0], [0, 0], [0, 0], [0, 0],
                    [0, 0], [0, -1], [0, -1], [0, 0.25], [0, 1], [0, 1], [0, 0], [0, 0], [0, 0],
                    [0, 0], [0, -1], [0, -1], [0, 0.25], [0, 1], [0, 1], [0, 0], [0, 0], [0, 0],
                    [0, -1], [0, -1], [0, 0.25], [0, 1], [0, 1], [0, 0], [0, 0], [0, -1], [0, -1],
                    [0, 0.25], [0, 1], [0, 1], [0, 0], [-1, -1], [-1, -1], [-1, 0.25], [-1, 1],
                    [-1, 1], [-1, 0], [-1, -1], [-1, 0.25], [-1, 1], [-1, 1], [-1, 0], [0.25, 0.25],
                    [0.25, 1], [0.25, 1], [0.25, 0], [1, 1], [1, 1], [1, 0], [1, 1], [1, 0], [0, 0]])

epsilon_IJ_array = np.array([5.40611759,  3.09808962,  0.51767336, -0.62764204,  3.88542973,
                 2.61443827,  0.51202621,  6.12446591,  1.37576539,  0.32210565,
                 5.06627117,  0.9647362 ,  4.32564458, -1.50980237,  1.46005216,
                 3.07880601,  2.04061218,  0.29545953,  2.83643323,  1.45560946,
                -0.30497317,  1.52988766,  3.70894349,  5.9557993 ,  5.03456066,
                 4.74468201, -0.01674577,  1.20937036, -2.22109131,  2.95550346,
                -0.35016995, -1.21752917,  6.73392849, -1.50708021,  4.62869075,
                -0.59512001, -0.55112695, -0.3656899 , -1.9658966 ,  2.68406431,
                 2.26592239,  5.70205144,  2.30349623,  7.3825577 ,  6.56432493,
                 1.80507238, -0.26861043,  3.54022756,  1.14694091, -0.19925119,
                -0.40476933, -0.51448164, -0.54374523,  2.82747419, -0.43225249,
                -0.04843008,  2.24079093,  3.14378417,  5.02962327,  3.3988123 ,
                -0.77766301,  6.53000233,  2.33972407,  0.60255702,  3.08127543,
                 2.24597019, -0.76336263,  3.44296293, -0.09434225,  2.8402411 ,
                 1.99283104,  4.73850077,  1.62716708,  1.18248019, -0.09998319,
                 0.6486047 ,  6.05758959, -0.65870019, -0.1109588 ,  0.3378409 ,
                -0.16326205, -0.17247464, -0.2934108 ,  1.31294028, -0.39801773,
                -0.20975644, -0.25145029,  1.86812694,  1.17256655, -0.18712225,
                 2.32733932,  3.58517465,  0.48914052,  1.31955928,  0.38585403,
                -0.13374552,  0.24067764, -0.38287574,  1.43475793, -0.15824449,
                -0.27761803,  2.91040309,  3.24258453,  0.65174409, -1.46007226,
                 5.88896039,  7.1430776 ,  3.88745143,  2.30743914,  2.48624402,
                 4.76125029,  3.37968047, -0.92619807, -0.52295748,  2.34658075,
                -0.34883458,  2.86150258,  8.44190868,  0.18479902,  5.90139303,
                 3.46739114, -0.65254193,  2.20577995,  1.57548556, -0.91784509,
                 4.23399264, -0.91084029,  1.96250241, -0.28071996,  6.50263792,
                -0.03736976,  0.85891016,  0.90442359, -0.43892065, -0.21006736,
                -0.13621232,  1.41282455,  0.02545168, -0.43625392, -0.24314636,
                -0.02268125,  0.02151196, -0.20146637, -0.17790325,  1.05866876,
                 0.85346621, -0.83219928, -0.16794813,  0.11123743, -1.07286937,
                 0.04053044, -0.00905679,  2.45017145,  0.5493021 , -0.12796812,
                 0.76363091,  1.02420736, -0.1558867 ,  1.01381073,  2.58471218,
                 1.63556123,  0.85508171, -0.41092507, -0.84079636, -0.4768993 ,
                 0.33689536, -0.31075561, -0.62846966, -0.25737436, -0.19516141,
                 0.43765623, -1.69732705,  0.28593947, -0.52365284, -0.42854419,
                 0.40018055,  1.10525538, -1.68142683, -0.50747077, -2.63998244,
                 1.08612754, -2.20352417, -0.7196425 , -0.38100844, -1.83995051,
                -0.1972894 ,  1.58274944,  2.39362824, -1.43086374,  0.96101211,
                 1.73861261,  0.39563111, -0.56566581, -0.99271574, -0.80221829,
                -1.76250634, -2.40922128, -0.91943366, -0.95607598, -0.7649912 ,
                -0.06799879,  2.0053653 ,  4.28148442,  2.33409089,  5.74551336,
                -1.75083667, -1.08100162,  0.03677466, -0.40148038, -1.1971242])

def potential_energy(r, qi, qj, sigma, epsilon_IJ):
    # Constants
    ionic_strength = 150e-3  # Ionic strength in M
    k = 1.380649e-23  # Boltzmann constant (J/K)
    T = 300  # Temperature (K)
    eps0 = 8.8541878128e-12  # Permittivity of free space (F/m)
    charge = 1.60217662e-19  # Charge of an electron (C)
    No = 6.02214076e23  # Avogadro's number (1/mol)
    conv = 2000  # Converts mol/L to mol/m^3 and accounts for ion pairs

    eta = 0.7  # Smoothness parameter
    r_cut = 8  # r_cut is the switching distance
    eps_w = 78.4  # Relative permittivity of water
    A = -8.5525
    kappa = 7.7839
    lam = 0.003627 # Debye length parameter (A)
    B = eps_w - A
    k_e = 332.0636 # (Kcal - Ang) / (q^2 - mole)(convert Coulomb eng to Kcal/mole)

    # Distance dependent dielectric constant
    eps_r = A + B / (1 + kappa * np.exp(-lam * B * r))
    # Debye length in Angstroms
    lambda_D = np.sqrt((eps0 * eps_r * k * T) / (conv * No * ionic_strength * charge ** 2))*1e10
    # Electrostatic potential energy
    electrostatic = k_e*(qi * qj) / (eps_r * r) * np.exp(-r / lambda_D)
    # Lennard-Jones potential (excluded volume)
    lj_potential = abs(epsilon_IJ) * (sigma / r) ** 12
    
    # Contact potential (smooth cutoff function)
    C_r = 0.5 * (1 + np.tanh(eta * (r_cut - r)))
    
    # Pairwise potential energy
    pairwise_potential = (lj_potential - epsilon_IJ * C_r) / 4.184  # Convert J to kcal/mol
    
    # Total potential energy
    total_potential = electrostatic + pairwise_potential
    
    #return total_potential / 4.184  # Convert J to kcal/mol
    return total_potential


potentials = {}

n_amino_acids = 20
pair_index = 1
n_amino_acids = 20
pair_index = 1
for i in range(n_amino_acids):
    for j in range(n_amino_acids):
        if i <= j:
            key = f"{residues[i]}-{residues[j]}"
            if (qi_qj_array[pair_index-1][0] * qi_qj_array[pair_index-1][1]) != 0:
    
                potentials[key] = {
                    "qi": qi_qj_array[pair_index-1][0],
                    "qj": qi_qj_array[pair_index-1][1],
                    "sigma": sigma_array[pair_index-1],
                    "epsilon_IJ": epsilon_IJ_array[pair_index-1],
                    "start": 0.1,
                    "cutoff": 20.0,
                    "npoints": 1000
                }
            else:
                potentials[key] = {
                    "qi": qi_qj_array[pair_index-1][0],
                    "qj": qi_qj_array[pair_index-1][1],
                    "sigma": sigma_array[pair_index-1],
                    "epsilon_IJ": epsilon_IJ_array[pair_index-1],
                    "start": 0.1,
                    "cutoff": 12.0,
                    "npoints": 1000
                }
            pair_index += 1

# Generate table data for each potential
tables = {}
for key, params in potentials.items():
    r_values = np.linspace(params["start"], params["cutoff"], params["npoints"])
    potential_values = potential_energy(r_values, params["qi"], params["qj"], params["sigma"], params["epsilon_IJ"])
    forces = -np.gradient(potential_values, r_values[1]-r_values[0])
    tables[key] = (r_values, potential_values, forces)


# Function to calculate the Lennard-Jones potential
def lj_potential(r, A):
    return A * ((7.0 / r)**12 - (7.0 / r)**6)

r_values = np.linspace(1.0, 8.0, 1000)
potential_values = lj_potential(r_values, 0.1)
forces = [-1.0 * 0.1 * ((-12*(7.0**12) / r**13) + (6*(7.0**6) / r**7)) for r in r_values]
tables["wall-aminoacid"] = (r_values, potential_values, forces)

def generate_table(output_filename):
    with open(output_filename, 'w') as output:
        # Header
        output.write("# UNITS: real\n")
        output.write("# Pairwise potential table\n\n")
        for key, (r_values, potential_values, forces) in tables.items():
            output.write(f"{key.upper()}\n")
            output.write(f"N {len(r_values)} R {r_values[0]} {r_values[-1]}\n\n")
            for i, (r, E, F) in enumerate(zip(r_values, potential_values, forces), 1):
                output.write(f"{i} {r} {E} {F}\n")
            output.write("\n")  # Separate tables with a newline for clarity

# Main function to handle command line arguments
def main():
    parser = argparse.ArgumentParser(description='Generate tables for pairwise potential')
    parser.add_argument('output_filename', type=str, help='The name of the output file')

    args = parser.parse_args()

    generate_table(args.output_filename)

if __name__ == '__main__':
    main()

