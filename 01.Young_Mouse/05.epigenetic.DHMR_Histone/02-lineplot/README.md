### 01_generate_subclass_bw: 
This script is used to generate single-base resolution BW files for different subclasses.

### 02-generate_DHMR_hyper_bed: 
This script generates hyper-DHMRs files for different subclasses.

### [03-07]
The shell scripts from 03-07 (including cmd or scale) are used to generate DNA methylation ratio for specified bed intervals [different histone modification bed intervals.].

### 08-plot_hyper_DHMR_histone_modification: 
This script is used to plot line graphs of DNA methylation ratios in hyperactive DNA methylation regions (hyperDHMR) across different cell types, illustrating the dynamic relationship between histone modifications and DNA methylation. The script first reads multiple histone modification data files, extracts modification information for each cell type, and associates it with hyperactive DNA methylation regions. It then uses the ggplot2 library to plot line graphs, visually representing the methylation levels across different histone modifications and cell types.