import os
import glob
import numpy as np

averages = []
# Iterate through all files matching the pattern
for file_path in glob.glob(f'pressure_vf*'):
        # Load the data from the file
        with open(file_path, 'r') as f:
            data = np.loadtxt(f, skiprows=2)
            averages.append(np.mean(data[1:]))
print(np.mean(averages))
