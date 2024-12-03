import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def calculate_vf(filename, BoxSize, NPolymer=50, NParticle=50):

    BeadSize = [2.0, 0.6, 0.6, 3]
    DX = BoxSize[0] / 50
    XSize = np.linspace(-BoxSize[0]/2 + DX/2, BoxSize[0]/2 - DX/2, 50)
    NX = len(XSize)

    # Load the trajectory
    universe = mda.Universe(filename, format='XYZ')

    ConcentrationPol = np.zeros(NX)
    ConcentrationPar = np.zeros(NX)
    VolumeFracPol = np.zeros(NX)
    VolumeFracPar = np.zeros(NX)
    
    NT = 0
    #NT = 3000
    # Process each frame
    for ts in universe.trajectory:
        # Center the positions
        type1 = universe.select_atoms("type 1")
        polymer_positions = type1.positions
        type4 = universe.select_atoms("type 4")
        particle_positions = type4.positions

        mean_x = np.mean(polymer_positions[:, 0])
        polymer_positions[:, 0] -= mean_x
        particle_positions[:, 0] -= mean_x

        # Periodic boundary conditions
        polymer_positions[:, 0] = np.where(polymer_positions[:, 0] >= BoxSize[0]/2, polymer_positions[:, 0] - BoxSize[0],
                                           np.where(polymer_positions[:, 0] <= -BoxSize[0]/2, polymer_positions[:, 0] + BoxSize[0],
                                                    polymer_positions[:, 0]))
        particle_positions[:, 0] = np.where(particle_positions[:, 0] >= BoxSize[0]/2, particle_positions[:, 0] - BoxSize[0],
                                            np.where(particle_positions[:, 0] <= -BoxSize[0]/2, particle_positions[:, 0] + BoxSize[0],
                                                     particle_positions[:, 0]))

        # Check for condensate split
        if np.std(polymer_positions[:, 0]) > 50:
            print('condensate split')

        # Calculate histograms
        CountPol, _ = np.histogram(polymer_positions[:, 0], bins=NX, range=(-BoxSize[0]/2, BoxSize[0]/2))
        CountPar, _ = np.histogram(particle_positions[:, 0], bins=NX, range=(-BoxSize[0]/2, BoxSize[0]/2))

        # Update concentration and volume fraction
        ConcentrationPol += CountPol / (6.02e23 * BoxSize[1] * BoxSize[2] * DX * 1e-27)  # mM
        ConcentrationPar += CountPar / (6.02e23 * BoxSize[1] * BoxSize[2] * DX * 1e-27)  # mM
        VolumeFracPol += CountPol * 4 * np.pi / 3 * (BeadSize[0] / 2)**3 / (BoxSize[1] * BoxSize[2] * DX)
        VolumeFracPar += CountPar * 4 * np.pi / 3 * (BeadSize[3] / 2)**3 / (BoxSize[1] * BoxSize[2] * DX)
        NT += 1
        
    # Normalize by number of frames
    ConcentrationPol /= NT
    ConcentrationPar /= NT
    VolumeFracPol /= NT
    VolumeFracPar /= NT

    return XSize, ConcentrationPol, ConcentrationPar, VolumeFracPol, VolumeFracPar

def plot_vf(VolumeFrac):

    # Plot the volume fraction profile
    plt.figure(figsize=(10, 6))
    plt.plot(VolumeFrac)
    plt.ylabel('Volume Fraction')
    plt.grid(True)
    plt.show()

    print(f"Mean Volume Fraction Polymer: {np.mean(VolumeFrac[int(len(VolumeFrac)//2-5):int(len(VolumeFrac)//2+5)])}")


