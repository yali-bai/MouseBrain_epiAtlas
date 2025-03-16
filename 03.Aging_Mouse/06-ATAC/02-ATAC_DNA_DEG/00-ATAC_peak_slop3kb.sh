#!/bin/bash

# "outdir" is a custom output path.
# outdir=""

input_dir=""
output_dir="${outdir}/00_ATAC_open_peak/"

for bed_file in "$input_dir"/*.bed
do
  output_file="${output_dir}/$(basename "$bed_file")"
  awk '{
    center = int(($2 + $3) / 2);
    start = center - 3000;
    end = center + 3000;
    if (start < 0) start = 0;
    print $1, start, end, $1 ":" start "-" end
  }' OFS='\t' "$bed_file" > "$output_file"
done


