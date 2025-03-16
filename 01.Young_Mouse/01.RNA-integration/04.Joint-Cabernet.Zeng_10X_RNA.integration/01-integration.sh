#!/usr/bin/bash
#SBATCH -J Integration
#SBATCH --mem=400G
#SBATCH --cpus-per-task=6
#SBATCH --partition=cpu2
#SBATCH -o ./logs/Integration.log
#SBATCH -e ./logs/Integration.log
conda activate analysis
Rscript 01-integration.R 01-meta.txt
