### Joint_Cabernet_brain_slice_raw_count.without_QC_filter.rds: 
Joint-Cabernet brain slice RNA raw count matrix without QC filter. The rows are genes and the columns are 11520 samples.

### subset.z_axis_located_on_7.33.cortex_and_hippo.rds: 
The RNA raw count matrix from Zhuang's dataset after Z-coordinate filtering (7.33 < z < 7.34) and brain region filtering (hippocampus and cortex only).

### RNA_DNA_match_name_QC_Joint_Cabernet_brain_slice.csv: 
This CSV file provides metadata for matched single cells profiled across three omics modalities—RNA, 5-hydroxymethylcytosine (5hmC), and 5-methylcytosine (5mC)—from the Joint Cabernet brain slice dataset. Each row represents a single cell that has been uniquely matched across all three modalities, with corresponding cell IDs listed in the RNA, hmC, and mC columns. 
The column "Unique_ID_match" gives a unified identifier shared by the three omics layers for the same biological cell. 
Quality control (QC) information is provided through the RNA_filter4_QC, hmC_QC, and mC_QC columns, where 1 indicates the cell passed QC for that modality and 0 indicates it did not. The final column "total_QC" summarizes whether a cell passed QC across all three modalities simultaneously (1 for passed in all, 0 otherwise).
