import os
import numpy as np
import argparse
from widoms_insertion import WidomsMethod

def process_directory(directory_path, box_size, bead_type, bead_radius, probing_radius, tmax, spacing=0):
    xyz_files = [f for f in os.listdir(directory_path) if f.endswith('.xyz')]
    if not xyz_files:
        print(f"No .xyz files found in {directory_path}")
        return None

    all_profiles = []
    for xyz_file in xyz_files:
        file_path = os.path.join(directory_path, xyz_file)
        print(f"Processing {xyz_file}...")
        widom = WidomsMethod(file_path, box_size, spacing, tmax)
        profile = widom.compute_profile(bead_type, bead_radius, probing_radius)
        all_profiles.append(profile)

    return np.mean(all_profiles, axis=0) 

def main():
    parser = argparse.ArgumentParser(description="Process XYZ files in a directory and output results.")
    parser.add_argument('directory_path', type=str, help="Directory containing XYZ files.")
    parser.add_argument('output_suffix', type=str, help="Suffix for the output file.")
    args = parser.parse_args()
    
    # Constants
    NUM_POLYMERS = 50 #60
    BEADS_PER_POLYMER = 450
    BOND_LENGTH = 0.38
    BEAD_RADIUS = 0.3

    # Calculate box size
    R = BEAD_RADIUS
    H = R - BOND_LENGTH / 2
    V1 = np.pi / 3 * (3 * R - H) * H**2
    V2 = 4 * np.pi / 3 * R**3 - 2 * V1
    single_pol_volume = BEADS_PER_POLYMER * V2

    VOLUME_FRACTION = float(args.output_suffix) / 1000
    box_side = (NUM_POLYMERS * single_pol_volume / VOLUME_FRACTION)**(1/3)
    box_size = [box_side-0.5, box_side-0.5, box_side-0.5]
    bead_type = 1
    bead_radius = 0.3
    tmax = 499 #1900 
    spacing = 0.2  # Meaning inserting only at 0
    average_profile_array = []
    
    r_values = np.arange(0, 3.5, 0.5) #r_values = np.exp(np.linspace(np.log(0.5), np.log(3),20))
    for probing_radius in r_values:
        average_profile = process_directory(args.directory_path, box_size, bead_type, bead_radius, probing_radius, tmax, spacing)
        if average_profile is not None:
            average_profile_array.append(average_profile)
    
    partitioning_values = np.array(average_profile_array)
    output_data = np.column_stack((r_values, partitioning_values))
    #output_file = f'partitioning_for_vf{args.output_suffix}.dat'
    output_file = f'partitioning_for_vf{args.output_suffix}.dat'
    np.savetxt(output_file, output_data)

if __name__ == "__main__":
    main()

