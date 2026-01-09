#!/bin/sh
sbatch << RUN
#!/usr/bin/bash
#SBATCH -J QC
#SBATCH -o add_QC.log
#SBATCH -e add_QC.log
#SBATCH --partition=compute_fat
#SBATCH --cpus-per-task=2
##SBATCH --mem=40g
#SBATCH --mem=500g
set -x

Rscript 02-change_QC_before_integration.R --input TSO-joint.RNA_QC_stat.young.without_filter_nFeature_nCount.csv --in_dir ../../../03.data/02.metainfo/01.Young_Mouse/ --metainfo rds/Joint_Cabernet_seurat.rds

set +x
RUN
