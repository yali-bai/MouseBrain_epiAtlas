#!/bin/bash

# Usage: bash ./06-H3K27ac_scale_Allsubclass.sh <bed_filename> <bw_filename> <outdir_suffix>

# "indir" is a custom input path, and "your_outdir" is a custom output path.
# indir=""
# your_outdir=""

bed_filename=$1
bw_filename=$2
outdir_suffix=$3

cpus=10

bed_dir="${indir}/02_DHMR_hyper_bed/"
bw_dir="${indir}/01_bw/"
base_outdir="${your_outdir}/06_H3K27ac_Allsubclass/"
bed_file="${bed_dir}${bed_filename}"
bw_file="${bw_dir}${bw_filename}"
outdir="${base_outdir}${outdir_suffix}/"

mkdir -p "${outdir}"

subclass=$(basename "${bw_filename}" .bw)

echo "Processing subclass: ${subclass}"
echo "BED file: ${bed_file}"
echo "BW file: ${bw_file}"
echo "Output directory: ${outdir}"

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J scale_${subclass}_H3K27ac
#SBATCH --mem=30G
#SBATCH --cpus-per-task=10
#SBATCH --partition=compute_fat
#SBATCH -o ${your_outdir}/logs/06-${subclass}_H3K27ac.scale.log
#SBATCH -e ${your_outdir}/logs/06-${subclass}_H3K27ac.scale.log
#SBATCH --time=150:00:00

set -x

# Debugging output to confirm input and output paths
echo "Running computeMatrix for subclass: ${subclass}"
echo "Input file: ${bw_file}"
echo "annotation: ${bed_file}"
echo "Output file: ${outdir}/${subclass}.mm10.gz"

computeMatrix scale-regions \
      -S ${bw_file} \
      -R ${bed_file} \
      --binSize 5 \
      --numberOfProcessors ${cpus} \
      -o ${outdir}/${subclass}.mm10.gz

if [ \$? -eq 0 ]; then
  echo "computeMatrix completed successfully."

  touch ${outdir}/${subclass}.mm10.tmp

  echo "Running R script for subclass: ${subclass}"
  echo "Rscript input: ${outdir}/${subclass}.mm10.gz"
  echo "Rscript output: ${outdir}/${subclass}.mm10.tmp"

  Rscript ./01-process_getMean.R \
  ${outdir}/${subclass}.mm10.gz \
  ${outdir}/${subclass}.mm10.tmp

  if [ \$? -eq 0 ]; then
    echo "R script completed successfully."
  else
    echo "R script encountered an error."
    exit 1  # Exit if R script fails
  fi

else
  echo "computeMatrix encountered an error."
  exit 1  # Exit if computeMatrix fails
fi

set +x

RUN


