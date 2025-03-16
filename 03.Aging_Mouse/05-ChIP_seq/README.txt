01-get_allc_txt_for_subclass_young: Generate a series of .txt files for merging allc data by subclass for young cortex mice.
02_allcools_merge_subclass_5hmC_young: Use the allcools software to merge individual cell allc files by subclass (5hmC data for young mice).
02_allcools_merge_subclass_5mC_young: Use the allcools software to merge individual cell allc files by subclass (5mC data for young mice).
03_get_allc_txt_for_subclass_aging: Generate a series of .txt files for merging allc data by subclass for aged cortex mice.
04_allcools_merge_subclass_5hmC_aging: Use the allcools software to merge individual cell allc files by subclass (5hmC data for aged mice).
04_allcools_merge_subclass_5mC_aging: Use the allcools software to merge individual cell allc files by subclass (5mC data for aged mice).
05-allc2bw_5hmC_young: Use the allcools software to convert subclass-level allc files into bw files (5hmC data for young mice).
05-allc2bw_5mC_young: Use the allcools software to convert subclass-level allc files into bw files (5mC data for young mice).
06-allc2bw_5hmC_aging: Use the allcools software to convert subclass-level allc files into bw files (5hmC data for aged mice).
06-allc2bw_5mC_aging: Use the allcools software to convert subclass-level allc files into bw files (5mC data for aged mice).
07-generate_true5mC_bw_aging: Generate subclass-level true5mC bw files (data for aged mice).
07-generate_true5mC_bw_young: Generate subclass-level true5mC bw files (data for young mice).
09-plot_TET1_peaks_CG_methylation_ratio: This script plots the fitted curves of DNA methylation ratios on TET1 peaks (sorted by signal).
The shell scripts from 08 (including cmd or scale) are used to generate DNA methylation ratio for specified bed intervals (TET1 ChIP-seq peaks).