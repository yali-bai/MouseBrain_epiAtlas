for modification in mC hmC 
do
for start in $(seq 45001 5000 47712)
do
max=`echo "$start+4999"|bc|xargs`;
if [$max -gt 47712]; then 
  end=47712; 
else 
  end=`echo "$start+4999"|bc|xargs`; 
fi
end=`echo "$start+4999"|bc|xargs`

allc_path=${modification}_allc.${start}_${end}.txt

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J MCDS_class_${start}_${end}.${modification}
#SBATCH --mem=750G
#SBATCH --cpus-per-task=60
#SBATCH --partition=compute_fat
#SBATCH -o logs/run.segment.mcds.class.${start}_${end}.${modification}.log
#SBATCH -e logs/run.segment.mcds.class.${start}_${end}.${modification}.log
set -x


allcools generate-dataset  \
	--allc_table ${allc_path} \
	--output_path TSO-joint-RNA_Mouse_${modification}_all_cells.${start}_${end}.mcds \
	--chrom_size_path ../../../03.data/01.ref/mm10.chrom.sizes.nochrM.txt \
	--obs_dim cell  \
	--cpu 60 \
	--chunk_size 50 \
	--region segment ../01-wgbstools_generate_segment/tso.5${modification}.segment.3cpg.0_base.bed  \
	--quantifiers segment count CGN,CHN

set +x
RUN

done
done

