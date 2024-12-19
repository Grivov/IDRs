#!/bin/bash
#SBATCH --array=0
#SBATCH --job-name=widoms_insertion
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --partition=shared
#SBATCH -A yzhan567

module load anaconda
module load mdanalysis

#DIRS=(
#      "/home/vgrigor2/scratch4-yzhan567/idr/figure3_matlab_scripts/StickerSpacer_Chain450/Out6"
#      "/home/vgrigor2/scratch4-yzhan567/idr/figure3_matlab_scripts/StickerSpacer_Chain450/Out7"
#      "/home/vgrigor2/scratch4-yzhan567/idr/figure2_matlab_scripts/StickerSpacer_Chain450/Out"
#      "/home/vgrigor2/scratch4-yzhan567/idr/figure3_matlab_scripts/StickerSpacer_Chain450/Out9"
#)

# Calculate suffix based on the SLURM_ARRAY_TASK_ID (adjust numbers to match job indices)
#SUFFIXES=(6 7 8 9)
#SUFFIX=${SUFFIXES[$SLURM_ARRAY_TASK_ID]}

# Access the directory using the SLURM_ARRAY_TASK_ID
#DIRECTORY_PATH="${DIRS[$SLURM_ARRAY_TASK_ID]}"

# Execute the Python script with the directory path and output file suffix as arguments
python extract_prob_distr.py "/home/vgrigor2/scratch4-yzhan567/idr/figure4_crowding/StickerSpacer_Chain450/Out8" 8 
