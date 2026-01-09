The files are for overlapping DMRs and DHMRs with genebodies and promoters, and generating the relevant lists and heatmaps.

The "run_plot.py" script submits jobs to the slurm system by calling the "do_plot_DMR_RNA.sh", which calls the "plot_DMR_RNA.py" to actually do the plotting. The exact overlapping thresholds are explicitly shown in the codes.