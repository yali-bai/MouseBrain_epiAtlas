#!/usr/bin/bash
#SBATCH -J mcool_process
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute_fat
#SBATCH -o ./logs/02-mcool_process.log
#SBATCH -e ./logs/02-mcool_process.log
#SBATCH --time=150:00:00

python 02-mcool_process.py

