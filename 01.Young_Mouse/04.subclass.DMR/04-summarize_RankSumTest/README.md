### "01-call_DMR.ipynb" 
This script is used for DMR-calling with results from the Wilcoxon rank-sum tests and the absolute differences of the mean methylation levels of the group differed from the robust means of all groups. Number of cells during DMR calling, and FDR correction process, are implemented.

### "02-check_DMR_number.ipynb" 
This script is used for count the numbers of the DMRs.

### 03-volcano_plot_neuron_1vsOthers
The files in "03-volcano_plot_neuron_1vsOthers" are for the plotting of volcano plots related to DMR-calling, with the x-axis being the absolute differences, and the y-axis being the adjusted p-values.

### 04-heatmap_for_DMRs
The files in "04-heatmap_for_DMRs" are for the plotting of the heatmaps of the DMRs. The heatmaps are in the format of segment by subclass.

### 05-stacked_barplot_for_DMR_number
The files in "05-stacked_barplot_for_DMR_number" are for the plotting of the stacked barplots of the numbers of the DMRs, utilizing the result from "02-check_DMR_number.ipynb".

### 06-annotation_of_DMRs
The files in "06-annotation_of_DMRs" are for the annotation of the DMRs and the plotting of the relevant pie-charts. Allcools is used for the annotation process before the plotting.

### 07-hyper_DHMR_bw
The files in "07-hyper_DHMR_bw" are for generating BigWig files of the hyper DHMRs, which can be visualized in IGV. The hyper DHMR segments are first converted into a bedgraph format, before converting into the BigWig format.