#!/usr/bin/bash
#SBATCH -J DNA_methylation_ratio_subclass_mean
#SBATCH --mem=60G
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute_fat
##SBATCH --nodelist=fatnode04
#SBATCH -o ./logs/03-DNA_methylation_ratio_subclass_mean.log
#SBATCH -e ./logs/03-DNA_methylation_ratio_subclass_mean.log
#SBATCH --time=150:00:00

if [ "$#" -ne 3 ]; then
    echo "Usage: sbatch 03-DNA_methylation_ratio_subclass_mean.sh <omics> <element> <group>"
    exit 1
fi

omics=$1 # hmC
element=$2 # chrom100k
group=$3 # CG

set -x
python 03-DNA_methylation_ratio_subclass_mean.py "$omics" "$element" "$group"
set +x

