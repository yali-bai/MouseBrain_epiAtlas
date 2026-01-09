### 00.subset_zhuang_data_by_z_with_threshold_of_7.33_and_major_brain_region.R: 
This script is used to subsetting Zhuang MERFISH data, because the z-coordinates of the Joint Cabernet brain slice dataset are concentrated between 7.3 and 7.4 and the brain regions in the Joint Cabernet data are primarily the Isocortex and hippocampus.

### 01.integrated_by_three_classes.by_twice_transfer_label.R: 
This script separately integrates Joint-Cabernet slice RNA and Zhuang MERFISH RNA according to Exc, Inh, Non Neuron, which is decided by integration of Joint Cabernet dataset and Zeng 10X dataset.

### 02.MERFISH_mapped.R: 
This script mappes Joint-Cabernet slice cells to the nearst Zhuang MERFISH cell, and assign its locus to Joint-Cabernet slice cells

### 03.combine_integration_of_three_classes.R: 
This script integrates all Joint-Cabernet slice cells with Zhuang MERFISH cells, and assigns celltype label of step 01. At the same time, assign the loci calculated in step 02 to all cells for plotting  spatial distribution.

### 04.three_class_UMAP_plot.R：
UMAP plot of three classes (Exc, Inh, Non-neuron)

### 05.RNA_fill.R: 
fill zero of RNA expression matrix by subclass mean

### 06.DNA_fill_na.R: 
fill 5hmCG, 5mCG, 5mCG_5hmCG methlyation NA value by subclass mean

### 07.whole_brain_mapped_plot.R: 
This script plots loci and subclasses information of all Joint Cabernet cells in whole brain region

### 08.Pseudotime_analysis_of_RNA_after_fill_na.R: 
After getting gradually changed genes across IT and CA subclasses, we assign them as pseudotime markers and then visualize RNA expression levels and DNA methylation levels of one pseudotime marker
