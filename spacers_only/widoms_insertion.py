import MDAnalysis as mda
import numpy as np
from scipy.spatial import cKDTree
import argparse
import matplotlib.pyplot as plt

class WidomsMethod:
    def __init__(self, input_file, box_size, spacing=0, tmax=0):
        self.input_file = input_file
        self.box_x = box_size[0]
        self.box_y = box_size[1]
        self.box_z = box_size[2]
        self.spacing = spacing
        self.tmax = tmax

    def compute_profile(self, bead_type, bead_radius, probing_radius):
        u = mda.Universe(self.input_file, format='XYZ')
        #for atom in u.atoms:
        #    atom.mass = 1.0 

        insertion_positions = []
        limit_y = self.box_y - probing_radius
        limit_z = self.box_z - probing_radius
        limit_x = self.box_x - probing_radius
        grid_x, grid_y, grid_z = np.meshgrid(np.arange(probing_radius, limit_x + limit_x/100, limit_x/100),
                                     np.arange(probing_radius, limit_y + limit_y/100, limit_y/100),
                                     np.arange(probing_radius, limit_z + limit_z/100, limit_z/100))
        insertion_positions = np.column_stack([grid_x.ravel(), grid_y.ravel(), grid_z.ravel()])
        #x = 0
        #insertion_positions_x = np.column_stack([np.full(grid_y.ravel().shape, x), grid_y.ravel(), grid_z.ravel()])
        #insertion_positions.append(insertion_positions_x)

        #slice_length = insertion_positions[0].shape[0]
        boltzmann_factors = []
        
        
        for ts in (u.trajectory[:self.tmax] if self.tmax else u.trajectory):
            beads = u.select_atoms(f"type {int(bead_type)}")
            #com = beads.center_of_mass()
            positions = beads.positions# - com 
            tree = cKDTree(positions)
            boltzmann_factors_timestep = []
                
            #all_insertion_positions = np.vstack(insertion_positions)
            distances_to_beads, _ = tree.query(insertion_positions)
            boltzmann_factors_all = (distances_to_beads >= bead_radius + probing_radius)
            #reshaped_boltzmann_factors = boltzmann_factors_all.reshape(-1, slice_length)
            #boltzmann_factors_timestep = np.mean(reshaped_boltzmann_factors, axis=1)
            boltzmann_factors_timestep = np.mean(boltzmann_factors_all) #averaging for all slices
            boltzmann_factors.append(boltzmann_factors_timestep)
                #for x_slice in insertion_positions:
                #    distances_to_beads, _ = tree.query(x_slice.reshape(-1, 3))
                #    boltzmann_factors_x = (distances_to_beads >= bead_radius + probing_radius)
                #    boltzmann_factors_timestep.append(np.mean(boltzmann_factors_x))

                #boltzmann_factors_timestep = np.array(boltzmann_factors_timestep)
                #boltzmann_factors.append(boltzmann_factors_timestep)

        boltzmann_factors = np.array(boltzmann_factors)
        return np.mean(boltzmann_factors, axis=0)
    
#if __name__ == "__main__":
#    parser = argparse.ArgumentParser(description="Compute volume fraction profile.")
#    parser.add_argument("-i", "--input_file", type=str, required=True, help="Path to the input file.")
#    parser.add_argument("-lo", "--ts_low", type=int, required=True, help="Lower timestep for analysis.")
#    parser.add_argument("-hi", "--ts_high", type=int, required=True, help="Upper timestep for analysis.")
#    parser.add_argument("-s", "--spacing", type=float, required=True, help="Spacing for grid points.")
#    parser.add_argument("-o", "--output_filename", type=str, required=True, help="Output filename for saving the results.")
#    parser.add_argument("-t", "--types", nargs='+', type=int, required=True, help="Bead types for analysis.")
#    parser.add_argument("-r", "--radii", nargs='+', type=float, required=True, help="Radii for each bead type.")
#    parser.add_argument("-p", "--probing_radius", type=float, required=True, help="Probing radius.")

#    args = parser.parse_args()

#    widom = WidomsMethod(args.input_file, args.ts_low, args.ts_high, args.spacing, args.output_filename)
#    profile = widom.compute_profile(args.types, args.radii, args.probing_radius)
#    print(profile)
#    np.savetxt(args.output_filename, profile)

