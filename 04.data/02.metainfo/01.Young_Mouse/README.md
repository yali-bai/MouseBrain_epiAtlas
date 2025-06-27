### RNA_DNA_match_name_QC_class_label.csv: 
The first is RNA cell names of all 47712 Joint Cabernet young cells, and the third and the forth columns are corresponding hmC and mC cell names. 
The second "Unique_ID_match" column is the id shared among RNA, hmC, mC samples. 
The fifth, sixth and the seventh columns are RNA, hmC, mC quality control information. 
The total QC will be 1 (indicating QC pass) only when all three individual QC measures are 1. The last three columns are class, subclass, three class label after RNA integration and label transfer.

### subclass_corresponding_name.csv: 
The corresponding name of subclass labels after integration of Joint-Cabernet RNA and Zeng 10X RNA with 3C. 
The first column is the subclass label of Joint-Cabernet RNA.
The second column is the subclass label of Zeng 10X RNA.

### TSO-joint.hmC_QC_stat.young.csv: 
hmC metainfo of young 47712 Joint-Cabernet cells. 
The first column is the cell ID.
The second column is the plate ID.
The third column is the cell ID without barcode.
The fourth column is the unique ID sharing among RNA, 5hmC and 5mC. 
The "species" column is the aligned reference genome. 
The "datatype" column represents the cell sequencing technology as TSO-joint-RNA. 
The "Library" column seperates RNA, hmC and 5mC. 
The "Region" and "Brain.Region" columns are the brain region information of the cell. 
The "Batch" column is the batch information of the cell. 
The "spatial_type" column represents whether it is sliced data or not. 
The "Age" column is the age of the mouse of the cell. 
The "QC" column is the quality control information of the cell for hmC, 1 means pass, and 0 means fail. 
The "class_label", "subclass_label" and "three_class_label" columns are the cell type labels of the cell. 
The "Neuron_non_neuron" column is the label of whether the cell is a neuron or not. 

### TSO-joint.mC_QC_stat.young.csv: 
mC metainfo of young 47712 Joint-Cabernet cells. The column names are the same as TSO-joint-hmC_QC_stat.young.csv.

### TSO-joint.RNA_QC_stat.young.csv: 
RNA metainfo of young 47712 Joint-Cabernet cells. 
The first column is the cell ID.
The second column is the plate ID. 
The "species" column is the species of the cell. 
The "QC" column is the quality control information of the cell for hmC, 1 means pass, and 0 means fail. 
The "Unique_ID" column is the unique ID sharing among RNA, 5hmC and 5mC. 
The "Region" and "Brain.Region" columns are the brain region information of the cell. 
The "SampleID_without_barcode" column is the cell ID without barcode. 
The "Batch" column is the batch information of the cell. 
The "class_label", "subclass_label" and "three_class_label" columns are the cell type labels of the cell. 
The "Neuron_non_neuron" column is the label of whether the cell is a neuron or not. 


