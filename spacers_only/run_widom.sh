#!/bin/bash
#SBATCH --array=0
#SBATCH --job-name=widoms_insertion
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --partition=shared
#SBATCH -A yzhan567

module load anaconda
module load mdanalysis

#DIRS=("./traj_vf29" "./traj_vf33" "./traj_vf37" "./traj_vf41" "./traj_vf46" "./traj_vf51" "./traj_vf57" "./traj_vf64" "./traj_vf71" "./traj_vf79") 

# Calculate suffix based on the SLURM_ARRAY_TASK_ID (adjust numbers to match job indices)
#SUFFIXES=(29 33 37 41 46 51 57 64 71 79)
#SUFFIX=${SUFFIXES[$SLURM_ARRAY_TASK_ID]}

# Access the directory using the SLURM_ARRAY_TASK_ID
#DIRECTORY_PATH="${DIRS[$SLURM_ARRAY_TASK_ID]}"

# Execute the Python script with the directory path and output file suffix as arguments
#python extract_prob_distr.py $DIRECTORY_PATH $SUFFIX
python extract_prob_distr.py "./traj_vf70" 70
