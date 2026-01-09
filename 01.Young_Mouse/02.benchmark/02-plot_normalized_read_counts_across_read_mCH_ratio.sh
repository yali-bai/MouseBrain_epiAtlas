#!/usr/bin/bash
#SBATCH -J Integration
#SBATCH --mem=400G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu2
#SBATCH -o ./logs/Integration.log
#SBATCH -e ./logs/Integration.log
#conda activate analysis

python 02-plot_normalized_read_counts_across_read_mCH_ratio.py
