function run_merge {
indir=../../output/03-aging/02-merge_allc
outdir=""

ls  $indir|grep txt|sed 's/.txt//g' |while read ff  
do
echo $ff

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J hmC_allc_merge
#SBATCH --mem=500G
#SBATCH --cpus-per-task=50
#SBATCH --partition=compute_pro
#SBATCH -o logs/allc_5mc.${ff}.log
#SBATCH -e logs/allc_5mc.${ff}.log
set -x

allcools merge-allc  \
	--allc_paths $indir/${ff}.txt \
	--output_path ${outdir}/merge_allc/${ff}.Merge.allc.tsv.gz \
	--chrom_size_path ../../input/reference_genome/mm10.chrom.sizes.nochrM.txt \
	--cpu 50 \
	--bin_length 10000000 

set +x
RUN
done
}

run_merge 
