# "indir" is a custom input path, and "outdir" is a custom output path.
function run_merge_big {
indir=""
outdir=""
wc -l ${indir}/* | grep -v total | sort -gr | awk '{if ( $1 >3300 ) print $2}'| awk -F "/" '{print $NF}' |grep "5hmC_"|sed 's/.txt.*//g' |while read ff
do

  echo $ff
sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J hmC_allc_merge_${ff}
#SBATCH --mem=700G
#SBATCH --cpus-per-task=105
#SBATCH --partition=compute_fat
#SBATCH -o ${outdir}/logs/subclass/allc_5hmC.${ff}.log
#SBATCH -e ${outdir}/logs/subclass/allc_5hmC.${ff}.log
set -x

allcools merge-allc  \
	--allc_paths $indir/${ff}.txt \
	--output_path ${outdir}/${ff}.Merge.allc.tsv.gz \
	--chrom_size_path ../../input/reference_genome/mm10.chrom.sizes.nochrM.txt \
	--cpu 105 \
	--bin_length 10000000 

set +x
RUN
done
}

function run_merge_small {
indir=""
outdir=""
wc -l ${indir}/* | grep -v total | sort -gr | awk '{if ( $1 <=3300 ) print $2}'| awk -F "/" '{print $NF}' |grep "5hmC_"|sed 's/.txt.*//g' |while read ff
#cat ./requeue_hmC_samples.txt | while read ff
do

  echo $ff
sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J hmC_allc_merge_${ff}
#SBATCH --mem=220G
#SBATCH --cpus-per-task=35
#SBATCH --partition=compute_fat
#SBATCH -o ${outdir}/logs/subclass/allc_5hmC.${ff}.log
#SBATCH -e ${outdir}/logs/subclass/allc_5hmC.${ff}.log
set -x

allcools merge-allc  \
        --allc_paths $indir/${ff}.txt \
        --output_path ${outdir}/${ff}.Merge.allc.tsv.gz \
        --chrom_size_path ../../input/reference_genome/mm10.chrom.sizes.nochrM.txt \
        --cpu 35 \
        --bin_length 10000000 
set +x
RUN
done


}

if [ $1 == "merge_big" ]; then
run_merge_big
elif [ $1 == "merge_small" ]; then
run_merge_small
else
echo ${1} not accepted
fi
