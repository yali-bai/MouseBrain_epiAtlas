### 01-compute_true5mC.py:  
We have 5hmC and 5hmC+5mC sequencing data, and we compute 5mC data by '5hmC+5mC' -5hmC.

### 02-all_cell_correlation.py:  
Calculating Pearson correlation coefficients and P-values of RNA expression and DNA methylation of each gene in all cells with 5hmC/5hmC+5mC/5mC, genebody, CG/CH data.

### 03-shuffled_correlation.py:  
Randomly perturbing DNA sample names 5000 times to calculate the correlation between gene RNA expression and DNA methylation in the random state.

### 04-divided_results_combination.r:  
Combining the calculated all cell correlation and 100 times shuffled correlation for each type of data.

### 05-metainfo.r:  
Generating the information of all genes on the genebody, including gene id, gene name, chromosome location, cpg quantity, etc.

### 06-allgene_results_combination.py:  
Integrating the correlation of all genes in all cells into one, which contains the metainfo of all genes.

### 07-all_gene_results_statistics.r:  
Counting the number of genes under different correlation coefficient gradients.

### 08-density_and_upset_plot.r:  
Plotting density plots and upset plots of correlation coefficient of different data.

### 09-one_gene_subclass_dotplot_data_prepare.py:  
RNA expression and DNA methylation dotplot plot data preparation in each subclass of gene AC132685.1.

### 10-one_gene_subclass_RNA-DNA_dotplot.r:  
Plotting RNA expression and DNA methylation dotplot plots of gene AC132685.1 in different subclasses.

### 11-3D_plot.ipynb: 
This script is used to generate a 3D scatter plot illustrating the correlation between CG methylation ratios (5mC, 5hmC, and 5mC+5hmC) in gene body regions. The script first reads metadata files, filters out significant correlations, and then uses the matplotlib library to plot 3D scatter plots, visually representing the spatial distribution of methylation levels across different subtypes. 

