#!/bin/bash
#SBATCH --array=1-10
#SBATCH --job-name=70                     # Job name (-J MyTest)
#SBATCH --time=3-00:00:00                      # Time limit (-t 12:00:00)
#SBATCH --nodes=1                            # Number of nodes (-N 1)
#SBATCH --ntasks=48                          # Number of processors (-n 2)
#SBATCH --cpus-per-task=1                    # Threads per process (-c 6)
#SBATCH --partition=parallel                     # Used partition (-p defq)
#SBATCH --mem-per-cpu=2GB                    # Define memory per core
#SBATCH --account=yzhan567

#module load gcc/9.3.0 openmpi/3.1.6 lammps/20200721

module load gcc/11.4.0 lammps/20231121

VF=70

mpirun -np 48 lmp -in input_script.in -var volume_fraction $VF -var replica $SLURM_ARRAY_TASK_ID
