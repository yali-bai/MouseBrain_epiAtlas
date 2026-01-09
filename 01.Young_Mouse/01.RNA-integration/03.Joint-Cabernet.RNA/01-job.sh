#!/usr/bin/bash
#SBATCH -J process
#SBATCH --mem=200G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu2
#SBATCH -o ./logs/RNA_process.log
#SBATCH -e ./logs/RNA_process.log

Rscript 01-scRNAseq.R 01-meta.txt mouse ../../../03.data/02.metainfo/01.Young_Mouse/TSO-joint.RNA_QC_stat.young.without_filter_nFeature_nCount.csv
