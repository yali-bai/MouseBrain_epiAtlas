"01-wgbstools_generate_segment" is for generating segments with the merged ALLC files using wgbstools.
"02-generate_segment_MCDS" is then used for calculating the mean CG methylation levels for the segments in individual cells.
"03-stratified_RankSumTest_1vsOthers" conducts the Wilcoxon rank-sum test for DMR-calling.
"04-summarize_RankSumTest" takes the result from "04_stratified_RankSumTest_1vsOthers" to call out the DMRs, and conducts some basic characterization of the DMRs.
"05-DMR_overlapping_genes" is for correlating the methylation status of the DMRs overlapping with genes (genebodies or promoters).

