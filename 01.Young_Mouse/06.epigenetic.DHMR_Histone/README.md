### 00-DHMR: 
This script is used to generate files related to DHMR.

### 01-process_getMean.R:
Calculate the methylation mean of given bed region

### 02-generate_DHMR_hyper_bed:
This script generates hyper-DHMRs files for different subclasses.

### [03-07]
The shell scripts from 03-07 (including cmd or scale) are used to generate DNA methylation ratio for specified bed intervals [different histone modification bed intervals.].

### 08-plot_heatmap_show_H3K4me1_modification_in_hyperDHMR_merge: 
This script is used to plot a heatmap of H3K4me1 modification signals in hyperactive DNA methylation regions (hyperDHMR) across different cell types, illustrating the relationship between DNA methylation ratios and histone modifications. The script first reads multiple input files, extracts modification data for each cell type, and merges them into a unified matrix format. It then uses the pheatmap library to visualize the data, providing a clear representation of modification patterns in hyperactive DNA methylation regions across different cell types. [08-12]
