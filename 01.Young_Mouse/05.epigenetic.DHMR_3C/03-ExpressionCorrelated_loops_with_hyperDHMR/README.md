### "01_filter_intergenic_to_selectedGenes_dots.ipynb", "02_check_dots_PCC_3C_RNA.ipynb", and "03_get_Distal30kbRegion_gene_pair_PCC.ipynb" 
These scripts are used one by one in order to select the loops with one anchor overlapping genes with expression levels over 0.2 in at least 1 subclass and the coefficient of variation (CV) of expression levels among subclasses over 0.5 (the genic anchor), and the other anchor not overlapping with any genes but having at least 1 hyper-DHMR within the 30-kb regions centering at this anchor (the intergenic anchor with hyper-DHMR)

### "Etl4_3C_plots.ipynb" and "Vwc2_3C_plots.ipynb" 
These scripts are used for generating the 3C contact heatmaps for the gene Etl4 and gene Vwc2, respectively.
The files in "plotting_for_supplementary" are for plotting the figures in the related Extended Data Fig. 9. The exact functions for the ipynb files in the folder are self-explanatory by their names.