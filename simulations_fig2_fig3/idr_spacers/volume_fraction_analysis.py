import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def calculate_vf(filename, BoxSize, Chain=450, Spacing=8, NPolymer=50, NParticle=50):
    Sticker = int(np.ceil(Chain / (2 * Spacing)) * 2)
    NPolymerBeads = NPolymer * (Sticker + Chain)
    NAtom = NParticle + NPolymerBeads

    # Bead sizes and volumes
    R = 0.6 / 2
    H = R - 0.38 / 2
    V1 = np.pi / 3 * (3 * R - H) * H**2
    V2 = 4 * np.pi / 3 * R**3 - 2 * V1

    #BeadSize = [0.6, 0.6, 0.6, 3]
    BeadSize = [0.6, 0.6, 0.6, 4] #FOR CROWDING
    DX = BoxSize[0] / 50
    XSize = np.linspace(-BoxSize[0]/2 + DX/2, BoxSize[0]/2 - DX/2, 50)
    NX = len(XSize)

    # Load the trajectory
    universe = mda.Universe(filename, format='XYZ')

    # Initialize arrays for concentration and volume fraction
    ConcentrationPol = np.zeros(NX)
    ConcentrationPar = np.zeros(NX)
    VolumeFracPol = np.zeros(NX)
    VolumeFracPar = np.zeros(NX)
    
    NT = 0
    # Process each frame
    for ts in universe.trajectory:
        # Center the positions
        positions = universe.atoms.positions
        polymer_positions = positions[:NPolymerBeads]
        particle_positions = positions[NPolymerBeads:]

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
        if np.std(polymer_positions[:, 0]) > 30:
            print('condensate split')

        # Calculate histograms
        CountPol, _ = np.histogram(polymer_positions[:, 0], bins=NX, range=(-BoxSize[0]/2, BoxSize[0]/2))
        CountPar, _ = np.histogram(particle_positions[:, 0], bins=NX, range=(-BoxSize[0]/2, BoxSize[0]/2))

        # Update concentration and volume fraction
        ConcentrationPol += CountPol / (6.02e23 * BoxSize[1] * BoxSize[2] * DX * 1e-27)  # mM
        ConcentrationPar += CountPar / (6.02e23 * BoxSize[1] * BoxSize[2] * DX * 1e-27)  # mM
        VolumeFracPol += V2 * Chain / (Sticker + Chain) * CountPol / (BoxSize[1] * BoxSize[2] * DX)
        VolumeFracPar += CountPar * 4 * np.pi / 3 * (BeadSize[3] / 2)**3 / (BoxSize[1] * BoxSize[2] * DX)
        NT += 1
        
    # Normalize by number of frames
    ConcentrationPol /= NT
    ConcentrationPar /= NT
    VolumeFracPol /= NT
    VolumeFracPar /= NT

    return XSize, ConcentrationPol, ConcentrationPar, VolumeFracPol, VolumeFracPar

def plot_vf(XSize, ConcentrationPol, VolumeFracPol, VolumeFracPar, Sticker, Chain):
    # Plot the concentration profile
    plt.figure(figsize=(10, 6))
    plt.semilogy(XSize, ConcentrationPol / (Sticker + Chain), 's-')
    plt.xlabel('X coordinate (nm)')
    plt.ylabel('Polymer concentration (mM)')
    plt.title('Polymer Concentration Profile')
    plt.grid(True)
    plt.savefig('concentration_profile.png')
    plt.show()

    # Plot the volume fraction profile
    plt.figure(figsize=(10, 6))
    plt.plot(XSize, VolumeFracPol, 's-', label='Polymer')
    plt.plot(XSize, VolumeFracPar, 's-', label='Particle')
    plt.plot(XSize, VolumeFracPol + VolumeFracPar, 's-', label='Total')
    plt.xlabel('X coordinate (nm)')
    plt.ylabel('Volume Fraction')
    plt.title('Volume Fraction Profile')
    plt.legend()
    plt.grid(True)
    plt.savefig('volume_fraction_profile12.png')
    plt.show()

    print(f"Mean Volume Fraction Polymer: {np.mean(VolumeFracPol[(XSize > -20) & (XSize < 20)])}")

    # Save the data
    np.savetxt('concentration_profile.txt', np.column_stack((XSize, ConcentrationPol / (Sticker + Chain))),
               header='X coordinate (nm)\tPolymer concentration (mM)', delimiter='\t')
    np.savetxt('volume_fraction_profile.txt', np.column_stack((XSize, VolumeFracPol, VolumeFracPar, VolumeFracPol + VolumeFracPar)),
               header='X coordinate (nm)\tPolymer Volume Fraction\tParticle Volume Fraction\tTotal Volume Fraction', delimiter='\t')

    print("Calculation complete. Results saved to PNG and TXT files.")

