### all_gene_upset_group_list_new.csv  : 
The upset groups of all genes in which RNA is significantly related to DNA. 
The first column is the gene id. 
The second column is the group name, for example,'[5mCG -]' means that RNA of this gene is significant negative correlation to 5mCG DNA methylation of this gene. The plus or minus sign represents a positive or negative correlation, and colons separate multiple DNA methylation types that are significantly associated with the RNA of this gene simultaneously. 
The third column is the methylation type, including CG and CH. 
The fourth column is the gene region. 
The fifth column represents all the significantly relevant genes used in this table. 
Columns six to eight show the correlation between RNA and 5mC, p-value and p-adjust value of the correlation. 
Columns nine to eleven show the correlation between RNA and 5hmC, p-value and p-adjust value of the correlation. 
Columns twelve to fourteen show the correlation between 5mC and 5hmC, p-value and p-adjust value of the correlation. 
The column fifteen is the gene name. 
The column sixteen is the length of gene. 
The column seventeen is the number of CpG sites in the gene.

### correspondence_of_subclass_and_class_for_integration_with_zeng.csv  :  
Mapping between subclass and class label after integration of Joint-Cabernet RNA and Zeng 10X RNA. The "subclass" column is the subclass and the "class" column is the corresponding class of the subclass.

### gene_metainfo.of_group1_group2.csv: 
The "genename" column is the genename. 
The "length" column is the length of gene. 
The "group" column is the gruop of gene. 
The "Cpg_number" is CpG number of gene. 
The "log_gene_length" is log10(gene length). 
The "log_CpG" column is the log10(gene CpG number). 
The "promoter_CpG_number" column is the CpG number of the promoter. 
The "log_promoter_CpG" is log10(promoter CpG number). 

### Integrated_Joint_Cabernet_zeng_markerGenes.rds: 
subclass marker genes of integration of Joint-Cabernet and Zeng 10X, we use it as features of integration. This is the result of FindAllMarkers of integration of Joint-Cabernet and Zeng 10X.

### zeng_subclass_mean_dat_final.csv  :  
The average of RNA expression in different subclasses of zeng data. The row names of the table are the gene ids, the column names are the subclass labels, and the value is the average RNA expression of the genes in the subclass of 10X Zeng data.




