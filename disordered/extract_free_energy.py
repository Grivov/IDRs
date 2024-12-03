#import os
#import numpy as np
#from widoms_insertion import WidomsMethod  # Assuming the class is in a file named widoms_method.py

#def process_directory(directory_path, box_size, bead_type, bead_radius, probing_radius, tmax, spacing=0):
#    xyz_files = [f for f in os.listdir(directory_path) if f.endswith('.xyz')]
    
#    if not xyz_files:
#        print(f"No .xyz files found in {directory_path}")
#        return None

#    all_profiles = []

#    for xyz_file in xyz_files:
#        file_path = os.path.join(directory_path, xyz_file)
#        print(f"Processing {xyz_file}...")
        
#        widom = WidomsMethod(file_path, box_size, spacing, tmax)
#        profile = widom.compute_profile(bead_type, bead_radius, probing_radius)
#        all_profiles.append(profile)

    # Calculate the average profile
#    average_profile = np.mean(all_profiles, axis=0)
    
#    return average_profile

#if __name__ == "__main__":
    # Set your parameters here
    #directory_path = "/home/vgrigor2/scratch4-yzhan567/idr/figure2_matlab_scripts/StickerSpacer_Chain450/Out2"  # Replace with your directory path
#    directory_path = "/home/vgrigor2/scratch4-yzhan567/idr/figure3_matlab_scripts/StickerSpacer_Chain450/Out9"
#    box_size = [100, 25, 25] 
#    bead_type = 1 
#    bead_radius = 0.3 
#    probing_radius = 1.5 
#    tmax = 50 
#    spacing = 0.0 #meaninig inserting only at 0
#    average_profile_array = []

#    for probing_radius in [0.5, 1.0, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]:
#        average_profile = process_directory(directory_path, box_size, bead_type, bead_radius, probing_radius, tmax)
#        average_profile_array.append(average_profile)
    
#    average_profile_array = np.array(average_profile_array)
#    np.savetxt('par_prob_distr_sp9.txt', average_profile_array)

    #if average_profile is not None:
    #    print("Average profile calculated successfully.")
    #    print(average_profile)

        # If you want to plot the average profile
    #    widom = WidomsMethod(None, box_size, spacing)  # Create a dummy WidomsMethod object for plotting
    #    widom.plot_profile(average_profile, "Average y_coordinate")
    #    np.savetxt('par_prob_distr_sp9.txt', average_profile)

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
    bead_radius = 0.3
    tmax = 100
    spacing = 0.2  # Meaning inserting only at 0
    average_profile_array = []

    for probing_radius in np.arange(0.0, 4.5, 0.5):
        average_profile = process_directory(args.directory_path, box_size, bead_type, bead_radius, probing_radius, tmax, spacing)
        if average_profile is not None:
            average_profile_array.append(average_profile)

    average_profile_array = np.array(average_profile_array)
    output_file = f'P_array_small_spacers.dat'
    np.savetxt(output_file, average_profile_array)

if __name__ == "__main__":
    main()
