import numpy as np
import random
import os
# Constants
NUM_POLYMERS = 50
BEADS_PER_POLYMER = 450
BOND_LENGTH = 0.38
BEAD_RADIUS = 0.3

# Calculate box size
R = BEAD_RADIUS
H = R - BOND_LENGTH / 2
V1 = np.pi / 3 * (3 * R - H) * H**2
V2 = 4 * np.pi / 3 * R**3 - 2 * V1
single_pol_volume = BEADS_PER_POLYMER * V2

def generate_random_unit_vector():
    theta = 2 * np.pi * random.random()
    phi = np.arccos(1 - 2 * random.random())
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    return np.array([x, y, z])

def is_overlap(bead, all_beads):
    for other_bead in all_beads:
        distance = np.linalg.norm(bead - other_bead)
        if distance < 2 * BEAD_RADIUS:
            return True
    return False

def generate_polymer(box_side):
    polymer = []
    all_beads = []

    # Generate first bead
    bead = np.array([random.uniform(0, box_side) for _ in range(3)])
    polymer.append(bead)
    all_beads.append(bead)

    # Generate remaining beads
    for _ in range(1, BEADS_PER_POLYMER):
        attempts = 0
        while attempts < 300:
            direction = generate_random_unit_vector()
            new_bead = polymer[-1] + direction * BOND_LENGTH

            # Wrap around periodic boundary conditions
            new_bead = new_bead % box_side

            if not is_overlap(new_bead, all_beads[:-1]):
                polymer.append(new_bead)
                all_beads.append(new_bead)
                break
            attempts += 1

        if attempts == 300:
            return None

    return polymer

def generate_polymers(box_side):
    all_polymers = []
    all_beads = []

    while len(all_polymers) < NUM_POLYMERS:
        polymer = generate_polymer(box_side)
        if polymer is not None:
            all_polymers.append(polymer)
            all_beads.extend(polymer)

    return all_polymers, all_beads


def write_lammps_data_file(filename, polymers, all_beads, box_side):
    with open(filename, 'w') as f:
        f.write("LAMMPS data file for polymer system\n\n")
        f.write(f"{NUM_POLYMERS * BEADS_PER_POLYMER} atoms\n")
        f.write(f"{NUM_POLYMERS * (BEADS_PER_POLYMER - 1)} bonds\n\n")
        f.write("1 atom types\n")
        f.write("1 bond types\n\n")
        f.write(f"0.0 {box_side} xlo xhi\n")
        f.write(f"0.0 {box_side} ylo yhi\n")
        f.write(f"0.0 {box_side} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write("1 5.65487\n\n")
        f.write("Atoms\n\n")
        for i, bead in enumerate(all_beads):
            f.write(f"{i+1} 1 1 {bead[0]:.6f} {bead[1]:.6f} {bead[2]:.6f}\n")
        f.write("\nBonds\n\n")
        bond_id = 1
        for i in range(NUM_POLYMERS):
            for j in range(BEADS_PER_POLYMER - 1):
                atom1 = i * BEADS_PER_POLYMER + j + 1
                atom2 = atom1 + 1
                f.write(f"{bond_id} 1 {atom1} {atom2}\n")
                bond_id += 1

# Volume fractions to generate
VF_VALUES = [29, 33, 37, 41, 46, 51, 57, 64, 71, 79]

for vf_int in VF_VALUES:
    VOLUME_FRACTION = vf_int / 100.0
    path = f"./initial_vf{vf_int}"
    os.makedirs(path, exist_ok=True)

    box_side = (NUM_POLYMERS * single_pol_volume / VOLUME_FRACTION)**(1/3)
    polymers, all_beads = generate_polymers(box_side)

    # Write LAMMPS data file for rep1 only
    output_file = f"{path}/polymer_system_vf{vf_int}_rep1.data"
    write_lammps_data_file(output_file, polymers, all_beads, box_side)

    print(f"Generated polymer system for vf={VOLUME_FRACTION}")
    print(f"  Box side length: {box_side:.6f}")
    print(f"  Output: {output_file}")
