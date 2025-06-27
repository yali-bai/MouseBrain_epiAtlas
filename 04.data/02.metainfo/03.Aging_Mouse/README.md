### old_mouse.raw_count.without_QC_filter.rds: 
Joint-Cabernet RNA raw count matrix without QC filter of aged mouse. The rows are genes and the columns are samples.

### TSO-joint.RNA_QC_stat.aged.csv: 
RNA metainfo of 86496 Joint-Cabernet cells, containing: 
(1) QC stat information used for quality control filter, QC (1 means pass, and 0 means fail), (2) Unique_ID for pairing with the other samples of the same cell, 
(3) brain region information, 
(4) batch information, 
(5) celltype (major class, subclass, three class, neuron or non-neuron) information after integration with Zeng 10X dataset and label transfer 
(6) age (mice which are before postnatal day 70 (P70) are classified as young, while those obtained after postnatal day 540 (P540) are designated as old).

### TSO-joint.hmC_QC_stat.aged.csv: 
hmC metainfo of young 86496 Joint-Cabernet cells, containing 
(1) QC stat information used for quality control filter, QC (1 means pass, and 0 means fail), 
(2) Unique_ID for pairing with the other samples of the same cell, 
(3) brain region information, 
(4) batch information, 
(5) celltype (major class, subclass, three class, neuron or non-neuron) information derived from RNA,
(6) age (mice which are before postnatal day 70 (P70) are classified as young, while those obtained after postnatal day 540 (P540) are designated as old).

### TSO-joint.mC_QC_stat.aged.csv: 
mC metainfo of young 86496 Joint-Cabernet cells, containing 
(1) QC stat information used for quality control filter, QC (1 means pass, and 0 means fail), 
(2) Unique_ID for pairing with the other samples of the same cell, 
(3) brain region information, 
(4) batch information, 
(5) celltype (major class, subclass, three class, neuron or non-neuron) information derived from RNA and age (mice which are before postnatal day 70 (P70) are classified as young, while those obtained after postnatal day 540 (P540) are designated as old).

### RNA_DNA_match_name_QC.aged.csv: 
The first is RNA cell names of all 86496 Joint Cabernet young cells, and the third and the forth columns are corresponding hmC and mC cell names. 
The second "Unique_ID_match" column is the id shared among RNA, hmC, mC samples. 
The fifth, sixth and the seventh columns are RNA, hmC, mC quality control information. 
The "total_QC" will be 1 (indicating QC pass) only when all three individual QC measures are 1. 
The "class_label", "subclass_label", "three_class_label" are celltype label after RNA integration and label transfer. 
The "age" is postnatal days of mice. Mice which are before postnatal day 70 (P70) are classified as young, while those obtained after postnatal day 540 (P540) are designated as old, which is the "old_young" column.