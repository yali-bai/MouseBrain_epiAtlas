### 01. allcools_merge_subclass_allc.sh: 
merge allc files by subclass and age

### 02. generate_segment
#### 01. generate_shell_for_each_merged_allc.sh: 
generate shell for each allc merged by subclass and age
#### 02. work.sh: 
run wgbstools to find candidate DMRs or DHMRs 
#### 03. get_segment_bed.sh: 
generate 0-base bed for calling mcds
#### merge.strand.py: 
combine methylated and coverage count of both strands into plus strand

### 03. run_segment_mcds.in_subclass_level.sh: 
run mcds by choosing region as segment to calculate 5mCG_5hmCG and 5hmCG methylation level in all segments

### 04. segment_significance_test
#### 01.allcools_to_get_DNA_raw_fraction.py: 
script of converting mcds to anndata 
#### 01.allcools_to_get_DNA_raw_fraction.sh: 
shell script to excute 01.allcools_to_get_DNA_raw_fraction.py
#### 02. DMRs_DHMRs_significance_test.py: 
do significance test (old versus young) for all candidate segments across all subclasses
#### 02. DMRs_DHMRs_significance_test.sh: 
shell script to excute 02. DMRs_DHMRs_significance_test.sh

### 05. DMR_DHMR_filter.by_diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.R: 
select significant DMRs and DHMRs by diff (old minus young), adjusted p-value, the number of significant digits, and length

### 06.DMR_DHMR_chrom_density_plot.R: 
plot the distribution of DMRs/DHMRs using 1-Mb binsize

### 07.aging_DMR_methyl_mean_subclass_diff_age.py: 
calculate the average levels of 5mC, 5hmC, and 5mC_5hmC across subclasses and age groups

### 08.aging_DMR_DNA_methyl_level_boxplot.R: 
boxplot of 5hmCG, 5mCG, 5mCG_5hmCG mean diff value of significant hyper DMRs, hypo DMRs, hyper DHMRs, hypo DHMRs

### 09.aging_DHMR_heatmap.R: 
heatmap of 5hmCG, 5mCG, 5mCG_5hmCG mean diff value of top significant hyper DHMRs

### 10. top_hyper_DHMRs_GO
#### 01. prepare_DHMRs_bed_for_GO.R: 
generate gradient top hyper DHMRs bed
#### 02. extract_DHMRs_intersected_genes.sh: 
generate genes intersected with gradient top hyper DHMRs
#### 03. GO.R: 
GO analysis of top hyper DHMRs covering genes.

### 11.heatmap_of_top1000_hyper_DHMRs_non_uniq.intersected_gene_with_top1_marker_labeled.r: 
plot heatmap of DHMRs which belong to top 1000 hyper DHMRs and are intersected with high expression and high coefficient of variation. After sorting DHMRs by subclass and difference (old minus young), the order of intersecting genes was correspondingly arranged, with subclass top 1 markers annotated.
