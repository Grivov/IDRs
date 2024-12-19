import os
import numpy as np
from volume_fraction_analysis import calculate_vf, plot_vf

# Directory containing the XYZ files
directory = '/home/vgrigor2/scratch4-yzhan567/vgrigor/idr/figure4_crowding/StickerSpacer_Chain450/Out8'
# Box dimensions
BoxSize = [100, 25, 25]  # Adjust as necessary

# Initialize arrays for averaging
all_ConcentrationPol = []
all_ConcentrationPar = []
all_VolumeFracPol = []
all_VolumeFracPar = []

# Iterate through all XYZ files in the directory
for filename in os.listdir(directory):
    if filename.endswith('.xyz'):
        filepath = os.path.join(directory, filename)
        XSize, ConcentrationPol, ConcentrationPar, VolumeFracPol, VolumeFracPar = calculate_vf(filepath, BoxSize)
        all_ConcentrationPol.append(ConcentrationPol)
        all_ConcentrationPar.append(ConcentrationPar)
        all_VolumeFracPol.append(VolumeFracPol)
        all_VolumeFracPar.append(VolumeFracPar)

# Convert lists to numpy arrays for averaging
all_ConcentrationPol = np.array(all_ConcentrationPol)
all_ConcentrationPar = np.array(all_ConcentrationPar)
all_VolumeFracPol = np.array(all_VolumeFracPol)
all_VolumeFracPar = np.array(all_VolumeFracPar)

# Calculate the average profiles
avg_ConcentrationPol = np.mean(all_ConcentrationPol, axis=0)
avg_ConcentrationPar = np.mean(all_ConcentrationPar, axis=0)
avg_VolumeFracPol = np.mean(all_VolumeFracPol, axis=0)
avg_VolumeFracPar = np.mean(all_VolumeFracPar, axis=0)

# Plot the average profiles
#plot_vf(XSize, avg_ConcentrationPol, avg_VolumeFracPol, avg_VolumeFracPar, Sticker=58, Chain=450)

# Save the data
np.savetxt('volume_frac_pol_corwding_sp8.txt', avg_VolumeFracPol)
np.savetxt('volume_frac_par_corwding_sp8.txt', avg_VolumeFracPar)
