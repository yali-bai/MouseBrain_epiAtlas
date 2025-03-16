01-compute_true5mC.py  :  We have 5hmC and 5hmC+5mC sequencing data, and we compute 5mC data by '5hmC+5mC' -5hmC.

02-all_cell_correlation.py  :  Calculating Pearson correlation coefficients and P-values of RNA expression and DNA methylation of each gene in all cells with 5hmC/5hmC+5mC/5mC, genebody, CG/CH data.

03-shuffled_correlation.py  :  Randomly perturbing DNA sample names 5000 times to calculate the correlation between gene RNA expression and DNA methylation in the random state.

04-divided_results_combination.r  :  Combining the calculated all cell correlation and 100 times shuffled correlation for each type of data.

05-metainfo.r  :  Generating metainfo of all genes on genebody and RNA.

06-allgene_results_combination.py  :  Integrating the correlation of all genes in all cells into one, which contains the metainfo of all genes.

07-all_gene_results_statistics.r  :  Counting the number of genes under different correlation coefficient gradients.

08-density_and_upset_plot.r  :  Plotting density plots and upset plots of correlation coefficient of different data.

09-one_gene_subclass_dotplot_data_prepare.py  :  RNA expression and DNA methylation dotplot plot data preparation in each subclass of gene AC132685.1.

10-one_gene_subclass_RNA-DNA_dotplot.r  :  Plotting RNA expression and DNA methylation dotplot plots of gene AC132685.1 in different subclasses.

11-3D_plot: Describes the 3D plot illustrating the correlation of CG methylation on the genebody.


