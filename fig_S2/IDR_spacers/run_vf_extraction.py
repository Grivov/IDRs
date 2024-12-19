#!/bin/bash
#SBATCH --array=1
#SBATCH --job-name=315
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --partition=shared
#SBATCH -A yzhan567

module load anaconda
module load mdanalysis

python average_volume_fraction.py
