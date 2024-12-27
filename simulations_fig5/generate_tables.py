import numpy as np
import argparse

def cos_potential(r, epsilon, sigma):
    return np.where((r > sigma/2-0.3)&(r < sigma/2+0.3),
                  -0.5 * epsilon * (1 + np.cos(np.pi * (r-sigma/2 ) / (0.3))) ,
                  0)
tables = {}

epsilons = np.arange(1,11,1)*4.14

sigma = 3.0

for epsilon in epsilons:
    key = f"sig{int(sigma*10)}eps{int(epsilon/4.14)}"

    r_values = np.linspace(sigma/20, 0.3+sigma/2, 10000)
    potential_values = cos_potential(r_values, epsilon, sigma)
 
    forces = -np.gradient(potential_values, r_values[1]-r_values[0])
    tables[key] = (r_values, potential_values, forces)

def generate_table(output_filename):
    with open(output_filename, 'w') as output:
        # Header
        output.write("# UNITS: nano\n")
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
