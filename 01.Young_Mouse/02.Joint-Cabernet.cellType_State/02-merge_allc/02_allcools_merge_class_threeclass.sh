function run_merge_outdir_indir_modification_level {
outdir=$1
indir=$2
modification=$3 #5mC 5hmC
level=$4 #class threeclass

wc -l ${indir}/* | grep -v total | sort -gr | awk -F "/" '{print $NF}' |grep ${modification}_ | sed 's/.txt.*//g' |while read ff
do

  echo $ff

if [ `wc -l ${indir}/${ff}.txt | awk '{print $1}'` == 1 ]; then
    order="cp `cat ${indir}/${ff}.txt` ${outdir}/${level}/${modification}/${ff}.Merge.allc.tsv.gz; \
           cp `cat ${indir}/${ff}.txt`.tbi ${outdir}/${level}/${modification}/${ff}.Merge.allc.tsv.gz.tbi"
else
    order="allcools merge-allc  \
        --allc_paths ${indir}/${ff}.txt \
        --output_path ${outdir}/${level}/${modification}/${ff}.Merge.allc.tsv.gz \
        --chrom_size_path ../../../04.data/01.ref/mm10.chrom.sizes.nochrM.txt \
        --cpu 40 \
        --bin_length 10000000"
fi


sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J ${level}_${modification}_allc_merge_${ff}
#SBATCH --mem=300G
#SBATCH --cpus-per-task=44
#SBATCH --partition=compute_fat
#SBATCH -o logs/${level}/allc_${modification}.${ff}.240904.log
#SBATCH -e logs/${level}/allc_${modification}.${ff}.240904.log
set -x

        ${order}

##allcools merge-allc  \
##	--allc_paths $indir/${ff}.txt \
##	--output_path ${outdir}/${level}/${modification}/${ff}.Merge.allc.tsv.gz \
##	--chrom_size_path ../../../04.data/01.ref/mm10.chrom.sizes.nochrM.txt \
##	--cpu 40 \
##	--bin_length 10000000 

set +x
RUN
done
}
