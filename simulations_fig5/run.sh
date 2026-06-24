#!/bin/bash
#SBATCH --job-name=client
#SBATCH --array=1-2                             # To match with epsilon values
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --cpus-per-task=1
#SBATCH --partition=parallel
#SBATCH --mem-per-cpu=2GB
#SBATCH --account=yzhan567

# Load necessary modules
module load gcc/11.4.0 lammps/20231121

# Receive replica index as first argument
replica_index=$1

# Define epsilon values array
epsilon_values=(4 5)

# Get epsilon value using SLURM array ID
index=$(($SLURM_ARRAY_TASK_ID - 1))
EV=${epsilon_values[$index]}

# Echo for debugging purposes
echo "Running simulation with epsilon_value: $EV and replica_index: $replica_index"

# Run LAMMPS
mpirun -np 48 lmp -in input_script_surf.in -var epsilon_value $EV -var index_arg $replica_index
