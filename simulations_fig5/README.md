# Simulations Figure 5

## Overview
Simulations of IDR polymers (Chain450, 50 stickers) interacting with point client particles at variable interaction strengths (epsilon). Density profiles are measured to characterize polymer partitioning.

## File Locations

Initial configurations: `initials/Sticker58_Chain450_NP50_Particle50_Rep${index}.initial`
Input script: `input_script.in`

## How to Run

Single simulation with specific parameters:
```
mpirun -np 48 lmp -in input_script.in -var index_arg <replica> -var epsilon_value <eps>
```

Arguments:
- `index_arg`: Replica number (1-30)
- `epsilon_value`: Interaction strength (1-10)

Using SLURM batch script:
```
sbatch run_patchy.sh
```

## Output

Density profiles: `Out_profile/density_clients_eps${eps}_rep${replica}.dat`
Trajectories: `OutC/FULL_trajectory_clients_eps${eps}_rep${replica}.xyz`

## Analysis

Plot density profiles: `python plot_profiles.py`
Extract binding affinity: `python extract_Kd.py`
