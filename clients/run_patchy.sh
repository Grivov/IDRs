#!/bin/bash
#SBATCH --array=8-16
#SBATCH --job-name=clients                       # Job name (-J MyTest)
#SBATCH --time=2-00:00:00                      # Time limit (-t 12:00:00)
#SBATCH --nodes=1                            # Number of nodes (-N 1)
#SBATCH --ntasks=48                          # Number of processors (-n 2)
#SBATCH --cpus-per-task=1                    # Threads per process (-c 6)
#SBATCH --partition=parallel                     # Used partition (-p defq)
#SBATCH --mem-per-cpu=2GB                    # Define memory per core
#SBATCH --account=yzhan567

module load gcc/11.4.0 lammps/20231121

#n_clients=$(( (SLURM_ARRAY_TASK_ID) / 10  + 1))
#replica_index=$(( (SLURM_ARRAY_TASK_ID) % 10 + 1))

n_clients=$(((SLURM_ARRAY_TASK_ID)))
replica_index=1

echo "Running simulation for ${n_clients} clients - Replica ${replica_index}"

mpirun -np 48 lmp -in input_script_patchy.in -var index_arg $replica_index -var n_clients $n_clients
