# IDR Simulations

## Overview
LAMMPS simulations of intrinsically disordered region (IDR) polymers with stickers under various conditions. Includes equilibration, osmotic pressure measurements, and client protein interactions.

## Simulations Figure 2 and 3

LAMMPS equilibration and production runs of IDR polymer systems (Chain450 with 58 stickers) and compact spacer polymers under confinement.

File locations:
- Compact spacers: `simulations_fig2_fig3/compact_spacers/`
- IDR spacers: `simulations_fig2_fig3/idr_spacers/`

Run MATLAB scripts to generate initial configurations, then run LAMMPS input files in each directory.

## Simulations Figure 4

Pressure measurements of IDR polymers at different volume fractions (vf=0.29-0.79). Used to calculate osmotic pressure and phase behavior.

File locations:
- Initial configurations: `simulations_fig4/initial_vf${vf}/polymer_system_vf${vf}_rep1.data`
- Input script: `simulations_fig4/input_script.in`

Generate initial configurations:
```
cd simulations_fig4
module load anaconda
python initial_generator.py
```

Run simulation:
```
mpirun -np 48 lmp -in input_script.in -var replica <rep> -var volume_fraction <vf>
```

Output: `simulations_fig4/pressure_vf${vf}/pressure_vf${vf}_rep${replica}.txt`

## Simulations Figure 5

IDR polymers interacting with point client particles at variable interaction strengths (epsilon=1-10). Density profiles characterize polymer partitioning.

File locations:
- Initial configurations: `simulations_fig5/initials/Sticker58_Chain450_NP50_Particle50_Rep${index}.initial`
- Input script: `simulations_fig5/input_script.in`
- Interaction tables: `simulations_fig5/tables3.dat`

Run simulation:
```
mpirun -np 48 lmp -in simulations_fig5/input_script.in -var index_arg <replica> -var epsilon_value <eps>
```

Arguments:
- `index_arg`: Replica number (1-30)
- `epsilon_value`: Interaction strength (1-10)

Output:
- Density profiles: `simulations_fig5/Out_profile/density_clients_eps${eps}_rep${replica}.dat`
- Trajectories: `simulations_fig5/OutC/FULL_trajectory_clients_eps${eps}_rep${replica}.xyz`

Analysis: `python simulations_fig5/plot_profiles.py` or `python simulations_fig5/extract_Kd.py`
