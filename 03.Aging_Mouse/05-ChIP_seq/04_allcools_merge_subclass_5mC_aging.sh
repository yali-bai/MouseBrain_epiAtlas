# "indir" is a custom input path, and "outdir" is a custom output path.

function run_merge {
indir=""
outdir=""
ls  $indir | sed 's/_old_5mC.txt$//g' | sed 's/subclass.//g' | while read ff
do

  echo $ff
sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J mC_allc_merge
#SBATCH --mem=500G
#SBATCH --cpus-per-task=60
#SBATCH --partition=compute_pro
#SBATCH -o ${outdir}/logs/allc_5mC_aging.${ff}.log
#SBATCH -e ${outdir}/logs/allc_5mC_aging.${ff}.log
set -x

allcools merge-allc  \
	--allc_paths $indir/subclass.${ff}_old_5mC.txt \
	--output_path ${outdir}/subclass/5mC/5mC_${ff}.Merge.allc.tsv.gz \
	--chrom_size_path ../../04.data/01.ref/mm10.chrom.sizes.nochrM.txt \
	--cpu 60 \
	--bin_length 10000000 

set +x
RUN
done
}

run_merge 
