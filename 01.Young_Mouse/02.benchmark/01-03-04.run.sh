#!/usr/bin/bash
#SBATCH -J Integration
#SBATCH --mem=400G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu2
#SBATCH -o ./logs/Integration.log
#SBATCH -e ./logs/Integration.log
#conda activate analysis

Rscript 01-QC_plot.r

Rscript 03-QC_boxplot_and_violin_plot.R

Rscript 04-UMAP_plot_of_read_counts.R