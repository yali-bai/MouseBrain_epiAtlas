function run_merge {
indir=./03.subclass_samplelist


for datatype in mC hmC
do
ls  $indir|awk -F "\\." '{print $1}'|while read subclass
do

  echo -e "${subclass}_${datatype}"
sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J hmC_allc_merge
#SBATCH --mem=500G
#SBATCH --cpus-per-task=60
#SBATCH --partition=compute_pro
#SBATCH -o logs/04.merge_allc_subclass.${subclass}_${datatype}.log
#SBATCH -e logs/04.merge_allc_subclass.${subclass}_${datatype}.log
set -x

allcools merge-allc  \
	--allc_paths $indir/${subclass}.${datatype}.samplelist.txt \
	--output_path 04.subclass_merged_allc/${subclass}.${datatype}.Merge.allc.tsv.gz \
	--chrom_size_path ../../../03.data/01.ref/mm10.chrom.sizes.nochrM.txt \
	--cpu 60 \
	--bin_length 10000000 

set +x
RUN
done
done
}

run_merge 
