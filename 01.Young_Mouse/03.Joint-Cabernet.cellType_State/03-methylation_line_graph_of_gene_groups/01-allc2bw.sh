#!/usr/bin/bash
indir=../02-merge_allc/04.subclass_merged_allc
outdir=./

function run_allc2bw {

chrom_size=../../../03.data/01.ref/mm10.chrom.sizes.nochrM.txt

ls $indir|grep "Merge\.allc\.tsv\.gz$"| sed 's/.Merge.allc.tsv.gz$//g' |while read ff

do 
echo ${1} $ff

sbatch <<RUN 
#!/usr/bin/bash 
#SBATCH -J allc2bw_${1}_${ff} 
#SBATCH --mem=200G 
#SBATCH --partition=compute_fat 
#SBATCH -o logs/${ff}.${1}.allc2bw.log 
#SBATCH -e logs/${ff}.${1}.allc2bw.log 
set -x 

allcools allc-to-bigwig --allc_path ${indir}/${ff}.Merge.allc.tsv.gz --output_prefix ${outdir}/${ff} --mc_contexts ${1}N --chrom_size_path ${chrom_size} --bin_size 1

set +x

RUN

done
}

run_allc2bw CG
run_allc2bw CH
