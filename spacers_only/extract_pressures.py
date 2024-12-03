import os
import glob
import numpy as np
import matplotlib.pyplot as plt

def pressure(VOLUME_FRACTION):
    averages = []
    # Iterate through all files matching the pattern
    for file_path in glob.glob(f'pressure_vf{int(VOLUME_FRACTION*1000)}_rep*'):
        # Load the data from the file
        with open(file_path, 'r') as f:
            data = np.loadtxt(f, skiprows=2)
            averages.append(np.mean(data[1:]))
    return np.mean(averages)


pressures = []
for VOLUME_FRACTION in np.exp(np.linspace(np.log(0.03), np.log(0.08),10)):
    path = f"./pressure_vf{int(VOLUME_FRACTION*1000)}"
    os.chdir(path)
    press = pressure(VOLUME_FRACTION)

    print(VOLUME_FRACTION, press)
    pressures.append(press)
    os.chdir("..")
np.savetxt('pressures.txt',pressures)
