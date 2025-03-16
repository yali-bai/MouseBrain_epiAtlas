"01-wgbstools_generate_beta_files.sh" is first used to generate bash files for converting the merged allc to beta files.
Then "02-wgbstools_generate_segment.sh" is used to generate the segments based on the beta files,
before the 0/1-base problem is handled by "03-get_segment_bed.sh".