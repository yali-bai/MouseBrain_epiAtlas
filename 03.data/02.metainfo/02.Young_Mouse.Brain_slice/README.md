### subset.z_axis_located_on_7.33.cortex_and_hippo.rds: 
The RNA raw count matrix from Zhuang's dataset after Z-coordinate filtering (7.33 < z < 7.34) and brain region filtering (hippocampus and cortex only).

### RNA_DNA_match_name_QC_class_label_young.brain_slice.add_celltype.csv: 
Slice data QC information, including unique_id, SampleID, and QC metrics for hmC, mC, and RNA, as well as total_QC, and labels for major class, subclass, three_class, Brain_Region, and neuron_type.
The column "unique_id" gives a unified identifier shared by the three omics layers for the same biological cell. 
The 'hmC_SampleID', 'mC_SampleID', and 'RNA_SampleID' columns correspond to the sample names for hmC, mC, and RNA of the same cell, respectively. 
Quality control (QC) information is provided through the RNA_QC, hmC_QC, and mC_QC columns, where 1 indicates the cell passed QC for that modality and 0 indicates it did not. The final column "total_QC" summarizes whether a cell passed QC across all three modalities simultaneously (1 for passed in all, 0 otherwise).
The 9th to 11th columns are class, subclass, three class label after RNA integration and label transfer.
The last four columns are detailed brain regions, brain region categories, whether it is slice data, and neuronal/non-neuronal information.
