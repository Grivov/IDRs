import numpy as np
import random
from scipy.optimize import minimize

def thomson_problem(N, radius=2):
    def energy(x):
        x = x.reshape(-1, 3)
        distances = np.sqrt(((x[:, np.newaxis, :] - x[np.newaxis, :, :])**2).sum(axis=-1))
        np.fill_diagonal(distances, np.inf)
        return (1 / distances).sum()

    def constraint(x):
        return np.sum(x.reshape(-1, 3)**2, axis=1) - radius**2

    initial_guess = np.random.randn(N * 3)
    initial_guess = initial_guess.reshape(-1, 3)
    initial_guess /= np.linalg.norm(initial_guess, axis=1)[:, np.newaxis]
    initial_guess = (initial_guess * radius).flatten()

    result = minimize(
        energy,
        initial_guess,
        method='SLSQP',
        constraints={'type': 'eq', 'fun': constraint},
        options={'ftol': 1e-5, 'maxiter': 1000}
    )

    if result.success:
        return result.x.reshape(-1, 3)
    else:
        raise ValueError(f"Optimization failed: {result.message}")

def create_molecule(N, radius=2):
    coordinates = []
    coordinates.append({'type': 4, 'pos': [0, 0, 0]})  # Center bead of type 4
    
    if N == 1:
        coordinates.extend([
            {'type': 5, 'pos': [0, 0, radius]}
        ])
    elif N == 2:
        coordinates.extend([
            {'type': 5, 'pos': [0, 0, radius]},
            {'type': 5, 'pos': [0, 0, -radius]}
        ])
    else:
        try:
            points = thomson_problem(N, radius)
            for point in points:
                coordinates.append({'type': 5, 'pos': point.tolist()})
        except Exception as e:
            print(f"Error in thomson_problem: {e}")
            # Fallback to a simple spherical distribution
            for i in range(N):
                phi = np.arccos(1 - 2 * (i + 0.5) / N)
                theta = np.pi * (1 + 5**0.5) * i
                x = radius * np.sin(phi) * np.cos(theta)
                y = radius * np.sin(phi) * np.sin(theta)
                z = radius * np.cos(phi)
                coordinates.append({'type': 5, 'pos': [x, y, z]})

    return coordinates

def generate_coms():
    coms = []
    while len(coms) < 25:
        com = np.array([random.uniform(-45, -35), random.uniform(-9.5, 9.5), random.uniform(-9.5, 9.5)])
        if not is_overlap(com, coms):
            coms.append(com)
    return coms

def is_overlap(com, all_coms, bead_radius=2):
    for other_com in all_coms:
        if np.linalg.norm(com - other_com) < 2 * bead_radius:
            return True
    return False

coms = generate_coms()

both_coms = np.concatenate((coms, coms + np.array([80, 0, 0])))

print(len(both_coms))

def write_lammps_data_file(filename, N):
    #N = 6
    molecule = create_molecule(N)
    with open(filename, 'w') as f:
        f.write("LAMMPS data file for client system\n\n")
        f.write(f"{(N+1)*len(both_coms)} atoms\n")
        f.write("0 bonds\n\n")
        f.write("5 atom types\n")
        f.write("0 bond types\n\n")
        f.write("-50.0 50.0 xlo xhi\n")
        f.write("-12.5 12.5 ylo yhi\n")
        f.write("-12.5 12.5 zlo zhi\n\n")
        #f.write("Masses\n\n")
        #f.write("1 1\n")     
        #f.write("2 1\n")
        #f.write("3 1\n\n")
        #f.write("4 1\n")
        #f.write("5 1\n\n")
        f.write("Atoms\n\n")
        atom_id = 1
        for i, com in enumerate(both_coms):
            for atom in molecule:
                x, y, z = atom['pos']
                f.write(f"{atom_id} {i+1} {atom['type']} {com[0]+x:.6f} {com[1]+y:.6f} {com[2]+z:.6f}\n")
                atom_id += 1

for N in range(1,16):
    write_lammps_data_file(f"client_system_{N}.data", N)



#for N in [6]:

#    molecule = create_molecule(N)
#    filename = f"client_{N}_stickers.txt"

    # Open file for writing
#    with open(filename, 'w') as file:
#        file.write(f"\n")
#        file.write(f"{N+1} atoms\n\n")
#        file.write("Coords\n\n")
#        for i, atom in enumerate(molecule):
#            x, y, z = atom['pos']
#            file.write(f"{i+1}\t{x:.5f}\t{y:.5f}\t{z:.5f}\n")

#        file.write("\nTypes\n\n")
#        for i, atom in enumerate(molecule):
#            file.write(f"{i+1}\t{atom['type']}\n")

#    print(f"Output written to {filename}")
