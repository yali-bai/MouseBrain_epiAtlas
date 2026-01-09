### TSO-joint.RNA_QC_stat.young.without_filter_nFeature_nCount.csv: 
QC stat of Joint-Cabernet RNA which was filtered by 'Align_rate', 'RNA_reads_ratio', 'Qualimap_Intronic', 'Qualimap_Intergenic', 'chrM_readsN' threshold

### TSO-joint.RNA_QC_stat.young.csv  :  
The final QC stat of young Joint-Cabernet RNA which was filtered by nFeature and nCount.
The first column is the cell ID.
The second column is the plate ID. 
The "species" column is the species of the cell. 
The "QC" column is the quality control information of the cell for RNA, 1 means pass, and 0 means fail. 
The "QC_after_integration" column is the quality control information of the cell for RNA after filtering nFeature and nCount thresholds.

### TSO-joint.DNA_QC_stat.young.csv: 
hmC, mC QC stat of young Joint-Cabernet DNA. 
The first column is the cell ID.
The "QC" column is the quality control information of hmC and mC after applying their respective filtering conditions. 
The 'unique_id' column is the id shared among RNA, hmC, mC samples.
The 'total_QC' column is the quality control information after filtering RNA, hmC, mC threshold.
The Library column indicates which data type the sample belongs to, either hmC or mC.

### RNA_DNA_match_name_QC_class_label_young.csv  :  
RNA, hmC, and mC matched cell names, unique_id, single and paired QC, and cell type labels.
The first 'unique_id' column is the id shared among RNA, hmC, mC samples.
The seond 'hmC_SampleID' column is hmC sample name of all 47712 Joint Cabernet young cells, and the forth and the sixth columns are corresponding mC and RNA cell names.  
The third, fifth, seventh columns are hmC, mC, RNA quality control information. 
The "total_QC" will be 1 (indicating QC pass) only when all three individual QC measures are The 9th to 11th columns are class, subclass, three class label after RNA integration and label transfer.
The last three columns are detailed brain regions, brain region categories, and neuronal/non-neuronal information.

### zeng_v3_metadata_downsample_1000.rds : 
metainfo of downsampled Zeng cells which is the reference of integration.
The "cell_new" column is the cell name of Zeng 10X data.
The 'region' column is the brain region infomation of each cell.
The last two columns 'class_label' and 'subclass_label' are the major class and subclass information of each cell.

### TSO-joint.5hmCH.global_methy.txt: 
global 5hmCH methylation level of Joint-Cabernet.
The first column is cell name.
The last columns are methylated number, coverage and global methylate fraction of each cell.

### TSO-joint.5mCH_5hmCH.global_methy.txt: 
global 5mCH_5hmCH methylation level of Joint-Cabernet.
The column names are the same as TSO-joint.5hmCH.global_methy.txt.

### subclass_name.Joint_Cabernet_corresponding_to_3C.csv: 
The corresponding name of subclass labels after integration of Joint-Cabernet RNA and Zeng 10X RNA with 3C. 
The first column is the subclass label of Joint-Cabernet RNA.
The second column is the subclass label of Zeng 10X RNA.