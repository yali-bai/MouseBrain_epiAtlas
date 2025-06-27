01-preprocessing.r:  Preprocessing the data. Generating our young and aged RNA data after integration with zeng and filter out low-quality genes.
The filter condition of genes is that the sum of count values is greater than 500 and the number of cells with cpm value greater than 0.5 is greater than 2.

02-dividing_metacells.r:  Dividing single cells into metacells.
Note: If an error occurs when loading the hdWGCNA Package, for example, Package 'Rcpp' version 1.0.11 cannot be unloaded, re-open an R window.

03-screen_aged_DEG.r:  Screening aging related differential genes by DESeq2.

04-Zeng_RNA_logFC_DNA_diff.separate_DNA_into_three_group_by_young_mean.point_plot.NN_select_logFC_by_ncell.r: Plotting the dotplot of Zeng DEG log2FC and Joint Cabernet 5hmCG diff (aged - young) and plotting barplot of the correlation of Zeng DEG RNA log2FC and Joint Cabernet 5hmCG diff.