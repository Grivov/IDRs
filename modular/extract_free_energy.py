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
    #parser.add_argument('output_suffix', type=str, help="Suffix for the output file.")
    args = parser.parse_args()

    box_size = [100, 25, 25]
    bead_type = 1
    bead_radius = 1.0
    tmax = 100
    average_profile_array = []

    for probing_radius in np.arange(0.0, 4.5, 0.5):
        average_profile = process_directory(args.directory_path, box_size, bead_type, bead_radius, probing_radius, tmax, spacing)
        if average_profile is not None:
            average_profile_array.append(average_profile)

    average_profile_array = np.array(average_profile_array)
    output_file = f'P_array_large_spacers.dat'
    np.savetxt(output_file, average_profile_array)

if __name__ == "__main__":
    main()
