#!/usr/bin/bash

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

input_dir="${indir}/input/"
output_dir="${outdir}/01_bw/"

if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

if [ ! -d "./logs" ]; then
    mkdir -p "./logs"
fi

for bw_file in "$input_dir"/*.bw; do
    filename=$(basename "$bw_file")
    output_file="$output_dir/$filename"

    sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J generate_bw_${filename}
#SBATCH --mem=30G
#SBATCH --cpus-per-task=10
#SBATCH --partition=compute_fat
#SBATCH -o ${outdir}/logs/01-${filename}.log
#SBATCH -e ${outdir}/logs/01-${filename}.log
#SBATCH --time=150:00:00

set -x 

python 01_generate_subclass_bw.py "$bw_file" "$output_file"

set +x
RUN

done

