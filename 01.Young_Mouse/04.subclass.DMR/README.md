### "01-wgbstools_generate_segment" 
This script is for generating segments with the merged ALLC files using wgbstools, which are further used for DMR calling and further analyses.

### "02-generate_segment_MCDS" 
This script is then used for calculating the mean CG methylation levels for the segments in individual cells. The mean CG/CH methylation levels for the segments in individual cells or merged allc files are all conducted similarly.

### "03-stratified_RankSumTest_1vsOthers" 
This script conducts the Wilcoxon rank-sum test for DMR-calling. Requirement for minimum cell number, and FDR correction, are implemented.

### "04-summarize_RankSumTest" 
This script takes the result from "04_stratified_RankSumTest_1vsOthers" to call out the DMRs, and conducts some basic characterization of the DMRs.

### "05-DMR_overlapping_genes" 
This script is for correlating the methylation status of the DMRs overlapping with genes (genebodies or promoters).

