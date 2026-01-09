### 01-generating_subclass_name.Joint_Cabernet_corresponding_to_3C.ipynb: 
This notebook generates a CSV file that maps the subclass labels in our dataset to the corresponding subclass labels in the Liu dataset.

### 02-mcool_process: 
This script processes downloaded .mcool files for different subclasses, including matrix balancing, PCA analysis, and other related operations.

### 03-DNA_methylation_ratio_subclass_mean: 
This script is used to generate DNA methylation ratio matrices at the subclass level for different combinations, including CG and CH promoter/chrom100k regions [5hmC/total5mC].

### 04-true5mC_subclass_mean: 
This script is used to generate DNA methylation ratio matrices at the subclass level for different combinations, including CG and CH promoter/chrom100k regions [5mC].

### 05-addCompartment_subclass.sh: 
This script assigns A and B compartment labels to 100kb-binned intervals for different subclasses.

### 06-barplot_AB.ipynb: 
This script plots bar charts of the log fold change (logFC) of DNA methylation ratios between A and B compartments for different subclasses.

### 07-Pearson_Correlation_between_DNA_Methylation_Ratio_and_3C_AB_Compartment_Score_across_subclass.ipynb: 
This script calculates the correlation (PCC) between DNA methylation ratios and 3C A/B compartment scores across subclasses, and generates the corresponding plots.

### 08-generate_DNA_methylation_ratio_100kb.ipynb: 
This script generates DNA methylation ratios within 100kb bins.

### 09-plot_UMAP_chr7_903_show_3C_DNA_methylation.ipynb: 
This notebook visualizes DNA methylation patterns of the chr7:903 region on a UMAP embedding.

### 10-Pearson_Correlation_between_DNA_Methylation_Ratio_and_3C_AB_Compartment_Score_across_100K.ipynb: 
This script calculates the correlation (PCC) between DNA methylation ratios and 3C A/B compartment scores across 100kb bins, and generates the corresponding plots. The script first reads DNA methylation data and 3C compartment score files, extracts the mean values for each 100kb window, and associates them with methylation ratios. It then uses statistical methods to calculate the correlation and plots the results.

### 11-plot_RNA_DNA_3C_UMAP.ipynb: 
This script generates UMAP plots related to 3C and RNA data, and then first reads relevant data files, calculates the correlation (PCC) between 5hmC and 5mC ratios and 3C A/B compartment scores, and converts the data into a wide format. It then uses the ggplot2 library to plot scatter plots, showing the distribution of correlations across different subtypes, while incorporating UMAP coordinates to visually represent the spatial relationship between gene expression and epigenetic states. 