#!/usr/bin/bash
#SBATCH -J process
#SBATCH --mem=200G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu2
#SBATCH -o ./logs/RNA_process.log
#SBATCH -e ./logs/RNA_process.log

Rscript 01-scRNAseq.R 01-meta.txt mouse
