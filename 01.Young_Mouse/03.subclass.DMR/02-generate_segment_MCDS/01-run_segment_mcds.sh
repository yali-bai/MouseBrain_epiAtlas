# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

for modification in 5mC 5hmC;
do

ls ${indir} | grep txt | grep ${modification} |sed 's/.txt.*//g'|while read batch;
do
allc_path=${indir}/${batch}.txt

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J MCDS_${batch}
#SBATCH --mem=1500G
#SBATCH --cpus-per-task=120
#SBATCH --partition=compute_fat
#SBATCH -o logs/run.segment.mcds.${batch}.log
#SBATCH -e logs/run.segment.mcds.${batch}.log
set -x


allcools generate-dataset  \
	--allc_table ${allc_path} \
	--output_path ${outdir}/${batch}.mcds \
	--chrom_size_path ../../../04.data/01.ref/mm10.chrom.sizes.nochrM.txt \
	--obs_dim cell  \
	--cpu 70 \
	--chunk_size 50 \
	--region segment ../../../output/01.Young_Mouse/03-wgbstools_generate_segment/TSO_${modification}_segment_3CpG.bed  \
	--quantifiers segment count CGN

set +x
RUN

done

done
