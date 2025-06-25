### "01-wgbstools_generate_beta_files.sh" 
This script is first used to generate bash files for converting the merged allc files to beta files.

### "02-wgbstools_generate_segment.sh" 
This script is used to generate the segments based on the beta files,
before the 0/1-base problem is handled by "03-get_segment_bed.sh". 

This whole process is for the generation of DNA segments with relatively homogeneous DNA methylation modification rates,
which are then used for DMR calling and further analyses.