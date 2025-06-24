### 01-DNA_methylation_ratio_show_in_Joint-Cabernet_RNA_tSNE.ipynb: 
This script is used to visualize global DNA methylation ratios (including mCG, 5hmCG, and true5mC) on RNA t-SNE plots to show the methylation distribution characteristics across different samples in single-cell RNA sequencing data. 
1. The script first reads DNA methylation data and quality control information, merges and filters the data to include only qualified samples.
2. And then we associate the methylation information with a Seurat object. 
3. Finally, we plot t-SNE plots to visually represent the methylation levels across different cell clusters or subtypes.

### 02-DNA_methylation_ratio_show_in_violinplot.ipynb: 
This script is used to visualize global DNA methylation ratios (including mC, 5hmC, and true5mC) across different subclasses using violin plots to compare the differences between aging and young mice. 
1. Firstly, we reads DNA methylation data and quality control information, merges and filters the data to include only qualified samples.
2. Secondly, we associate the methylation information with a Seurat object. 
3. Finally, the script uses the ggplot2 library to plot violin plots, providing a visual representation of methylation levels across different subclasses. 

### 03-UMAP_plots_show_selected_genes_expr_and_DNA_methylation_ratio.ipynb: 
This script is used to visualize the joint distribution of RNA expression levels and DNA methylation ratios (including 5mC, 5hmC, and true5mC) of selected genes on UMAP plots, revealing the heterogeneity of different cell subtypes in both epigenetic and transcriptomic landscapes. 
1. The script first reads single-cell RNA sequencing data and DNA methylation data, combines them with UMAP coordinates, and normalizes the expression ranges of specific genes. 
2. It then uses the ggplot2 library to plot two-dimensional scatter plots, visually demonstrating the spatial relationship between gene expression and methylation status.

#### 04-RNA_var_Cor_modification.ipynb: 
This script is used to plot a scatter plot illustrating the correlation between RNA variance and DNA methylation ratios across different subtypes. 
1. The script first reads single-cell RNA sequencing data and associated metadata, calculates the variance of gene expression for each gene in different subtypes, and associates it with methylation ratios. 
2. And then it uses the ggplot2 library to plot scatter plots, visually demonstrating the relationship between RNA expression variability and methylation status.