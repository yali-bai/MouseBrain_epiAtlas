00.subset_zhuang_data_by_z_with_threshold_of_7.33_and_major_brain_region.R: script for subsetting Zhuang MERFISH data 
01.integrate_by_three_classes_respectively.R: separately integrate Joint-Cabernet slice RNA and Zhuang MERFISH RNA according to Exc, Inh, Non Neuron
02.MERFISH_mapped.R: map Joint-Cabernet slice cells to the nearst Zhuang MERFISH cell, and assign its locus to Joint-Cabernet slice cells
03.combine_integration_of_three_classes.R: integration all Joint-Cabernet slice cells with Zhuang MERFISH cells, assign spatial loci to 
04.fill_RNA_and_DNA_na.R: fill NA by subclass mean
05.Pseudotime.R: pseudotime analysis
06.RNA_Pseudotime_subclass_markers.color_RNA_and_DNA.R: plot RNA expression levels and DNA methylation levels pseudotime marker
07.whole_brain_mapped_plot.R: plot loci of subclasses in whole region