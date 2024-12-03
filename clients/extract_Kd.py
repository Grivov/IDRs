import MDAnalysis as mda
import numpy as np
from scipy.spatial import cKDTree
import argparse
import matplotlib.pyplot as plt
from MDAnalysis.analysis.distances import distance_array

def compute_C(input_file):
    u = mda.Universe(input_file, format='XYZ')

    p_all = []
    total_clients = []
    bound_stickers = []
    total_stickers = []
    
    for ts in u.trajectory:
        beads4 = u.select_atoms("type 4 and prop x > -14 and prop x < 14")
        beads2 = u.select_atoms("type 2 and prop x > -14 and prop x < 14")
        beads3 = u.select_atoms("type 3 and prop x > -14 and prop x < 14")
        #print(len(beads4), len(beads2), len(beads3), ts, "%")
        if len(beads4) == 0: 
            continue
        positions4 = beads4.positions
        positions2 = beads2.positions
        positions3 = beads3.positions
        distances24 = distance_array(positions2, positions4)
        minimal_distances24 = np.min(distances24, axis=1)
        distances34 = distance_array(positions3, positions4)
        minimal_distances34 = np.min(distances34, axis=1)
        bound_stickers.append(sum(minimal_distances24 < 1.8) + sum(minimal_distances34 < 1.8))
        total_stickers.append(len(beads2)+len(beads3))
        total_clients.append(len(beads4))
        #print(bound_stickers[-1], total_stickers[-1], total_clients[-1])
        #print(len(minimal_distances24), len(beads2), len(beads4))
        #distances = distance_array(positions2, positions3)
        #minimal_distances2 = np.min(distances, axis=0)
        #minimal_distances3 = np.min(distances, axis=1)
        #p_all.append(sum(minimal_distances2 > 0.3) + sum(minimal_distances3 > 0.3))
    return np.mean(bound_stickers), np.mean(total_stickers), np.mean(total_clients)

free_energies = []
rho_array = []
Kd_array = []

for epsilon in range(1, 11):
    num_replicas = 11 if (epsilon < 6) else 26
    B_array = []
    T_array = []
    N_array = []

    for replica in range(1, num_replicas):
        path = f"./OutC/FULL_trajectory_clients_eps{epsilon}_rep{replica}.xyz"
        print(path)
        #path = "./OutC/FULL_trajectory_clients_eps9_rep1.xyz"

        B, T, N = compute_C(path)
        B_array.append(B)
        T_array.append(T)
        N_array.append(N)
    
    B, T, N = np.mean(B_array), np.mean(T_array), np.mean(N_array)

    V = 28 * 25 * 25 - (N * 4 * 3.14 * 1.2**3 / 3)
    v = 4 * 3.14 * (1.8**3 - 1.2**3)/ 3

    Kd = (N*T - N*B) / (V*B - N*v*B)
    rho = T/V
    
    Kd_array.append(Kd)
    rho_array.append(rho)
    free_energies.append(-rho / Kd)

    print(f'for epsilon {epsilon} Kd = ', Kd)
    print(f'for epsilon {epsilon} rho = ', rho)
np.savetxt('Kd_array_2.txt', Kd_array)
np.savetxt('rho_array_2.txt', rho_array)
np.savetxt('F_recruit_2.txt', free_energies)
print(free_energies)
print(Kd_array)
print(rho_array)
#print((4 * np.pi / 3 * (1.8**3 - 1.2**3)) * compute_C(path) / 25 / 25/ 20)

