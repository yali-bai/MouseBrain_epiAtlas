#!/usr/bin/bash

# "indir" is a custom input path, and "outdir" is a custom output path.
indir=""
outdir=""
function run_allc2bw {

chrom_size=../../input/reference_genome/mm10.chrom.sizes.nochrM.txt

ls $indir | grep ".Merge.allc.tsv.gz$" | sed 's/5mC_//g' | sed 's/.Merge.allc.tsv.gz$//g' | while read ff
do 
    echo ${1} $ff

    sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J allc2bw_${1}_${ff}
#SBATCH --mem=200G
#SBATCH --partition=compute_pro
#SBATCH -o ${outdir}/logs/06-${ff}.${1}.allc2bw_aging_5mC.log
#SBATCH -e ${outdir}/logs/06-${ff}.${1}.allc2bw_aging_5mC.log
set -x

allcools allc-to-bigwig --allc_path ${indir}/5mC_${ff}.Merge.allc.tsv.gz --output_prefix ${outdir}/${ff} --mc_contexts ${1}N --chrom_size_path ${chrom_size} --bin_size 1

set +x
RUN

done
}

run_allc2bw CG
