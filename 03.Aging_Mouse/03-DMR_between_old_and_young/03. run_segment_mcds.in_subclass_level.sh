# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

for modification in 5hmC 5mC 
do

for start in $(seq 1 5000 90000)
do
end=`echo "$start+4999"|bc|xargs`

allc_path=${indir}/${start}_${end}.${modification}_allc_table.all_age.txt

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J MCDS_class_${start}_${end}.${modification}
#SBATCH --mem=500G
#SBATCH --cpus-per-task=60
#SBATCH --partition=compute_fat
#SBATCH -o logs/run.segment.mcds.class.${start}_${end}.${modification}.log
#SBATCH -e logs/run.segment.mcds.class.${start}_${end}.${modification}.log
set -x


allcools generate-dataset  \
	--allc_table ${allc_path} \
	--output_path ${outdir}/TSO-joint-RNA_Mouse_${modification}_all_cells.${start}_${end}.mcds \
	--chrom_size_path ../../input/reference_genome/mm10.chrom.sizes.nochrM.txt \
	--obs_dim cell  \
	--cpu 60 \
	--chunk_size 50 \
	--region segment ../../output/03-aging/03-DMRs_DHMRs/segments/tso.${modification}.segment.3cpg.0_base.bed  \
	--quantifiers segment count CGN

set +x
RUN

done

done
