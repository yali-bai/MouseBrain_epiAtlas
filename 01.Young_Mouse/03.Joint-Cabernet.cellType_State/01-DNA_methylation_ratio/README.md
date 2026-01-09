#### 01-color_global_CG_CH_methylation.R: 
The code uses violin plots to display the global levels of DNA methylation across different cell types, including various data types (5hmC, 5mC, 5mC+5hmC) and different methylation contexts (CG, CH).

1. Firstly, we reads DNA methylation data and quality control information, merges and filters the data to include only qualified samples.
2. Secondly, we associate the methylation information with a Seurat object. 
3. Finally, the script uses the ggplot2 library to plot violin plots, providing a visual representation of methylation levels across different subclasses. 


#### 02-UMAP_plots_show_Cdh18_expr_and_DNA_methylation_ratio.R: 
Primarily displays the RNA expression and DNA methylation ratio (5hmCG, 5mCG, 5mCG+5hmCG, 5mCH) of selected Cdh18 genes on the UMAP plot.
