### HighExprHighCV_genes.txt: 
A list of genes with high expression and high coefficient of variation.

### top1000_hyper_DHMRs_intersected_gene.intersected_subclass_marker.top1_log2FC_in_subclass_markers.csv: 
First, we obtained the genes from the intersection of all hyper DHMRs. After gene filtering (within the list of highly variable and highly expressed genes), we intersected this file with the top 1000 subclass markers based on the "subclass" and "geneid" columns. Finally, we selected the top 1 gene with the highest log2FC (i.e., the "avg_log2FC" column) as the marker for that subclass. 
The "subclass" column represents the subclass whose markers intersect with the hyper DHMRs of that subclass and belong to genes which are both highly variable and highly expressed. 
The "geneid" column is the top1 marker of all subclass. 
The "DHMR_chr", "DHMR_start", "DHMR_end" and "DHMR_region" columns are location information of the hyper DHMR intersected with the top1 marker. 
The "gene_chr", "gene_start", "gene_end" columns are location information of the top1 marker gene. 
The "length" column is the length of the DHMR. 
The "uniq" column is the combination of the "subclass" column and the "DHMR_region" column. 
The "diff" column represents the difference in 5hmCG methylation levels between aged and young samples within the current subclass. 
The "DHMR_region_index" and "sort" column contain the ranked index assigned to each DHMR. 
The "p_val" and "p_val_adj" columns are the result of signicance test. 
The "avg_log2FC" is the log2FoldChange of the gene in the "subclass" versus the others. 
The "pct.1" is the percentage of cell expressed the gene in the "subclass" and the "pct.2" is the percentage of cell expressed the gene in the other subclasses. 
The "genename" is the name of gene.