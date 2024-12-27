#!/bin/bash
#SBATCH --job-name=volume_fraction_simulations  # Job name
#SBATCH --array=1-50                            # 50 jobs, 10 volume fractions x 5 replicas
#SBATCH --time=3-00:00:00                         # Time limit hh:mm:ss
#SBATCH --nodes=1                               # Each job runs on one node
#SBATCH --ntasks=48                             # Number of tasks per job
#SBATCH --cpus-per-task=1                       # Number of cores per task
#SBATCH --partition=parallel                     # Used partition (-p defq)
#SBATCH --mem-per-cpu=2GB                    # Define memory per core
#SBATCH --account=yzhan567

# Load necessary modules
module load gcc/11.4.0 lammps/20231121

# Define volume fractions (you need to ensure that these are matched by directory or input names)

volume_fractions=(29 33 37 41 46 51 57 64 71 79)


# Calculate index for volume fraction and replica
vol_index=$(( (SLURM_ARRAY_TASK_ID - 1) / 5 ))
replica_index=$(( (SLURM_ARRAY_TASK_ID - 1) % 5 + 1 ))

# Set the volume fraction based on the job array task ID
VF=${volume_fractions[$vol_index]}

# Print job details
echo "Running simulation for volume fraction ${VF} - Replica ${replica_index}"

# Run LAMMPS using the same input file for all simulations
mpirun -np 48 lmp -in input_script.in -var volume_fraction $VF -var replica $replica_index
