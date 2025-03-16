01. allcools_merge_subclass_allc.sh: merge allc files by subclass and age
02. generate_segment/01. generate_shell_for_each_merged_allc.sh: generate shell for each allc merged by subclass and age
02. generate_segment/02. work.sh: run wgbstools to find candidate DMRs or DHMRs 
02. generate_segment/03. get_segment_bed.sh: generate 0-base bed for call mcds
03. generate_segment/merge.strand.py: combine methylated and coverage count of both strands into plus strand
03. run_segment_mcds.in_subclass_level.sh: run mcds by choosing region as segment to calculate 5mCG_5hmCG and 5hmCG methylation level in all segments
04. segment_significance_test/01.allcools_to_get_DNA_raw_fraction.py: script of converting mcds to anndata 
04. segment_significance_test/01.allcools_to_get_DNA_raw_fraction.sh: shell script to excute 01.allcools_to_get_DNA_raw_fraction.py
04. segment_significance_test/02. DMRs_DHMRs_significance_test.py: do significance test to all segments
04. segment_significance_test/02. DMRs_DHMRs_significance_test.sh: shell script to excute 02. DMRs_DHMRs_significance_test.sh
05. DMR_DHMR_filter.by_diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.R: select significant DMRs and DHMRs
06.DMR_DHMR_chrom_density_plot.R: plot DMRs DHMRs distribution per MB
07.aging_DMR_methy_mean_subclass_diff_age.py: generate mean 5hmCG, 5mCG, 5mCG_5hmCG value of all significant segments in all subclasses
08.aging_DMR_DNA_methyl_level_boxplot.R: boxplot of 5hmCG, 5mCG, 5mCG_5hmCG mean diff value of significant hyper DMRs, hypo DMRs, hyper DHMRs, hypo DHMRs
09.aging_DHMR_heatmap.R: heatmap of 5hmCG, 5mCG, 5mCG_5hmCG mean diff value of top significant hyper DHMRs
10.GO analysis of top hyper DHMRs covering genes.

