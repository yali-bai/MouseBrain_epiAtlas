#!/bin/bash

# Usage: bash ./01-true5mCH_scale_open_peak.sh <bed_filename> <bw_filename> <age>

# "indir" is a custom input path, and "your_outdir" is a custom output path.
# indir=""
# your_outdir=""

bed_filename=$1
bw_filename=$2
age=$3

cpus=10

bed_dir="${indir}/00_ATAC_open_peak/"
bw_file="${indir}/${age}/true5mC/${bw_filename}/CH/merge_bw/subclass.${bw_filename}_${age}_true5mC.Merge.allc.tsv.gz.CH.bw"
outdir="${your_outdir}/01_true5mCH_open_peak/"
bed_file="${bed_dir}${bed_filename}"

mkdir -p "${outdir}"

subclass=${bw_filename}

echo "Processing subclass: ${subclass}"
echo "BED file: ${bed_file}"
echo "BW file: ${bw_file}"
echo "Output directory: ${outdir}"

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J scale_${bed_filename}_{bw_filename}
#SBATCH --mem=30G
#SBATCH --cpus-per-task=10
#SBATCH --partition=compute_pro
#SBATCH -o ${your_outdir}/logs/01-${bed_filename}_${bw_filename}_${age}_true5mCH.scale.log
#SBATCH -e ${your_outdir}/logs/01-${bed_filename}_${bw_filename}_${age}_true5mCH.scale.log
#SBATCH --time=150:00:00

set -x

# Debugging output to confirm input and output paths
echo "Running computeMatrix for subclass: ${subclass}"
echo "Input file: ${bw_file}"
echo "annotation: ${bed_file}"
echo "Output file: ${outdir}/${subclass}_${age}_true5mCH.mm10.gz"

computeMatrix scale-regions \
      -S ${bw_file} \
      -R ${bed_file} \
      --binSize 5 \
      --numberOfProcessors ${cpus} \
      -o ${outdir}/${subclass}_${age}_true5mCH.mm10.gz

if [ \$? -eq 0 ]; then
  echo "computeMatrix completed successfully."

  touch ${outdir}/${subclass}_${age}_true5mCH.mm10.tmp

  echo "Running R script for subclass: ${subclass}"
  echo "Rscript input: ${outdir}/${subclass}_${age}_true5mCH.mm10.gz"
  echo "Rscript output: ${outdir}/${subclass}_${age}_true5mCH.mm10.tmp"

  Rscript ./01-process_getMean.R \
  ${outdir}/${subclass}_${age}_true5mCH.mm10.gz \
  ${outdir}/${subclass}_${age}_true5mCH.mm10.tmp 

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




