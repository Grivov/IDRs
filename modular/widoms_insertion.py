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
        limit_y = self.box_y/2 - probing_radius
        limit_z = self.box_z/2 - probing_radius
        limit_x = self.box_x/10 - probing_radius
        grid_x, grid_y, grid_z = np.meshgrid(np.arange(-limit_x, limit_x + limit_x/70, limit_x/70),
                                     np.arange(-limit_y, limit_y + limit_y/70, limit_y/70),
                                     np.arange(-limit_z, limit_z + limit_z/70, limit_z/70))
        insertion_positions = np.column_stack([grid_x.ravel(), grid_y.ravel(), grid_z.ravel()])

        boltzmann_factors = []
        
        
        for ts in (u.trajectory[:self.tmax] if self.tmax else u.trajectory):
            beads = u.select_atoms(f"type {int(bead_type)}")
            #com = beads.center_of_mass()
            positions = beads.positions# - com 
            tree = cKDTree(positions)
            boltzmann_factors_timestep = []
                
            distances_to_beads, _ = tree.query(insertion_positions)
            boltzmann_factors_all = (distances_to_beads >= bead_radius + probing_radius)
            #boltzmann_factors_timestep = np.mean(reshaped_boltzmann_factors, axis=1)
            boltzmann_factors_timestep = np.mean(boltzmann_factors_all) #averaging for all slices
            boltzmann_factors.append(boltzmann_factors_timestep)

        boltzmann_factors = np.array(boltzmann_factors)
        return np.mean(boltzmann_factors)
    

