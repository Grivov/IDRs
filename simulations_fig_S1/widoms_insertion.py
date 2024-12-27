import MDAnalysis as mda
import numpy as np
from scipy.spatial import cKDTree
import argparse
import matplotlib.pyplot as plt

class WidomsMethod:
    def __init__(self, input_file, box_size, tmax=0):
        self.input_file = input_file
        self.box_x = box_size[0]
        self.box_y = box_size[1]
        self.box_z = box_size[2]
        self.tmax = tmax

    def compute_profile(self, probing_radius, shift):
        u = mda.Universe(self.input_file, format='XYZ')
        #for atom in u.atoms:
        #    atom.mass = 1.0 

        insertion_positions = []
        limit_y = self.box_y/2 - probing_radius
        limit_z = self.box_z/2 - probing_radius
        limit_x = self.box_x/10 - probing_radius
        grid_x, grid_y, grid_z = np.meshgrid(np.arange(-limit_x + shift, limit_x + shift + limit_x/30, limit_x/30),
                                     np.arange(-limit_y, limit_y + limit_y/30, limit_y/30),
                                     np.arange(-limit_z, limit_z + limit_z/30, limit_z/30))
        insertion_positions = np.column_stack([grid_x.ravel(), grid_y.ravel(), grid_z.ravel()])
        #x = 0
        #insertion_positions_x = np.column_stack([np.full(grid_y.ravel().shape, x), grid_y.ravel(), grid_z.ravel()])
        #insertion_positions.append(insertion_positions_x)

        #slice_length = insertion_positions[0].shape[0]
        boltzmann_factors = []


        for ts in (u.trajectory[:self.tmax] if self.tmax else u.trajectory):
            #beads = u.select_atoms(f"type {int(bead_type)}")
            beads1 = u.select_atoms(f"type 1")
            beads4 = u.select_atoms(f"type 4")
            #com = beads.center_of_mass()
            positions1 = beads1.positions# - com
            positions4 = beads4.positions 
            tree1 = cKDTree(positions1)
            tree4 = cKDTree(positions4)
            boltzmann_factors_timestep = []

            #all_insertion_positions = np.vstack(insertion_positions)
            distances_to_beads1, _ = tree1.query(insertion_positions)
            distances_to_beads4, _ = tree4.query(insertion_positions)
            boltzmann_factors_all = (distances_to_beads1 >= 0.3 + probing_radius) & (distances_to_beads4 >= 2.0 + probing_radius)
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
        return np.mean(boltzmann_factors)

