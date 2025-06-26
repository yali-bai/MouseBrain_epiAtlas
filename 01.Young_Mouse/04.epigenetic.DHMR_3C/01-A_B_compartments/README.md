### 01-mcool_process: 
This script processes downloaded .mcool files for different subclasses, including matrix balancing, PCA analysis, and other related operations.

### 02-addCompartment_subclass: 
This script assigns A and B compartment labels to 100kb-binned intervals for different subclasses.

### 03-plot_AB_lineplot: 
This script is used to plot line graphs of DNA methylation ratios in A and B compartments across different subclasses to compare the differences between aging and young mice. The script first reads a CSV file containing methylation data and compartment information, extracts the mean values for A and B compartments for each subclass, and uses the ggplot2 library to plot line graphs, visually representing the methylation levels in A and B compartments across different subclasses.

### 04-barplot_AB_difference_value: 
This script plots bar charts of the difference in DNA methylation ratios between A and B compartments for different subclasses.

### 05-barplot_AB_logFC: 
This script plots bar charts of the log fold change (logFC) of DNA methylation ratios between A and B compartments for different subclasses.

### 06-PCC_between_DNA_Methylation_Ratio_and_3C_AB_Compartment_Score_across_100kb_bins: 
This script calculates the correlation (PCC) between DNA methylation ratios and 3C A/B compartment scores across 100kb bins, and generates the corresponding plots. The script first reads DNA methylation data and 3C compartment score files, extracts the mean values for each 100kb window, and associates them with methylation ratios. It then uses statistical methods to calculate the correlation and plots the results.

### 07-PCC_between_DNA_Methylation_Ratio_and_3C_AB_Compartment_Score_across_subclasses: 
This script calculates the correlation (PCC) between DNA methylation ratios and 3C A/B compartment scores across subclasses, and generates the corresponding plots.

### 08-generate_DNA_methylation_ratio_100kb: 
This script generates DNA methylation ratios within 100kb bins.

### 09-plot_3C_RNA_UMAP: 
This script generates UMAP plots related to 3C and RNA data, and then first reads relevant data files, calculates the correlation (PCC) between 5hmC and 5mC ratios and 3C A/B compartment scores, and converts the data into a wide format. It then uses the ggplot2 library to plot scatter plots, showing the distribution of correlations across different subtypes, while incorporating UMAP coordinates to visually represent the spatial relationship between gene expression and epigenetic states. 