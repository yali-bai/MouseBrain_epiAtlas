01-data_merge.py: Calculating the average value of DNA methylation rate of genes in each subclasses.
01-merge_sbtch.sh: Job submission for 01-data_merge.py.
01-RNA_merge.py: Calculating the average value of RNA expression of genes in each subclasses.
01-RNA_merge_sbtch.sh: Job submission for 01-RNA_merge.py.
02-merge_expr.r: The wide matrix of the average values of genes in subclasses becomes a long matrix.
03-dotplot_corrected_x-log1p.R: Plotting the RNA expression and DNA  methylation rate dotplots of genes in aging and youth respectively.
03-old_and_young_combined_dotplot.R: Plotting dotplot for the combination of RNA expression and DNA  methylation rate of genes in aging and youth.
04-subclass_correlation_barplot.r: Plotting the barplot of the correlation coefficients between RNA expression and DNA  methylation rate in subclasses for youth and aging.
04-subclass_correlation_barplot_Optimized.r: Plotting the barplot of the correlation coefficients between RNA expression and DNA  methylation rate in three classes for youth, and sortting them in the order of mCG+hmCG, mCG, hmCG, mCH+hmCH, mCH and hmCH.