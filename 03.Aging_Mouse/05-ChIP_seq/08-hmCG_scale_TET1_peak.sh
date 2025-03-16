#!/bin/bash

# Usage: bash ./08-hmCG_scale_TET1_peak.sh <bw_filename> <age>

# "indir" is a custom input path, and "your_outdir" is a custom output path.
# indir=""
# your_outdir=""


bw_filename=$1
age=$2

cpus=10

bed_file="../../output/03-aging/05-ChIP-seq/01-TET1_5D6_rep1_filter.bed"

outdir="${your_outdir}/08_hmCG_TET1_peak/"

bw_file="${indir}/${age}/02_allc_to_bw/subclass/5hmC/CG/${bw_filename}.CGN-both.frac.bw"

mkdir -p "${outdir}"

subclass=${bw_filename}

echo "Processing subclass: ${subclass}"
echo "BED file: ${bed_file}"
echo "BW file: ${bw_file}"
echo "Output directory: ${outdir}"

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J scale_${bw_filename}
#SBATCH --mem=30G
#SBATCH --cpus-per-task=10
#SBATCH --partition=compute_fat
#SBATCH -o ${your_outdir}/logs/08-${bw_filename}_${age}_5hmCG.scale.log
#SBATCH -e ${your_outdir}/logs/08-${bw_filename}_${age}_5hmCG.scale.log
#SBATCH --time=150:00:00

set -x

# Debugging output to confirm input and output paths
echo "Running computeMatrix for subclass: ${subclass}"
echo "Input file: ${bw_file}"
echo "annotation: ${bed_file}"
echo "Output file: ${outdir}/${subclass}_${age}_5hmCG.mm10.gz"

computeMatrix scale-regions \
      -S ${bw_file} \
      -R ${bed_file} \
      --binSize 5 \
      --numberOfProcessors ${cpus} \
      -o ${outdir}/${subclass}_${age}_5hmCG.mm10.gz

if [ \$? -eq 0 ]; then
  echo "computeMatrix completed successfully."

  touch ${outdir}/${subclass}_${age}_5hmCG.mm10.tmp

  echo "Running R script for subclass: ${subclass}"
  echo "Rscript input: ${outdir}/${subclass}_${age}_5hmCG.mm10.gz"
  echo "Rscript output: ${outdir}/${subclass}_${age}_5hmCG.mm10.tmp"

  Rscript ./08-process_getMean.R \
  ${outdir}/${subclass}_${age}_5hmCG.mm10.gz \
  ${outdir}/${subclass}_${age}_5hmCG.mm10.tmp 

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




