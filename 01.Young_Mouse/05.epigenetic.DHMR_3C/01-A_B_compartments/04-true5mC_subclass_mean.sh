#!/usr/bin/bash
#SBATCH -J true5mC
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute_pro
#SBATCH -o ./logs/04-true5mC.log
#SBATCH -e ./logs/04-true5mC.log
#SBATCH --time=150:00:00

if [ "$#" -ne 2 ]; then
    echo "Usage: sbatch 04-true5mC_subclass_mean.sh <group> <element>"
    exit 1
fi

group=$1 # CG
element=$2 # chrom100k

Rscript 04-true5mC_subclass_mean.R "$group" "$element"