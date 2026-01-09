function run_merge {
indir=./01.major_class_samplelist


for datatype in hmC mC
do
ls  $indir|awk -F "\\." '{print $1}'|while read major_class
do

  echo -e "${major_class}_${datatype}"
sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J hmC_allc_merge
#SBATCH --mem=500G
#SBATCH --cpus-per-task=60
#SBATCH --partition=compute_fat
#SBATCH -o logs/02.merge_allc_major_class.${major_class}_${datatype}.log
#SBATCH -e logs/02.merge_allc_major_class.${major_class}_${datatype}.log
set -x

allcools merge-allc  \
	--allc_paths $indir/${major_class}.${datatype}.samplelist.txt \
	--output_path 02.major_class_merged_allc/${major_class}.${datatype}.Merge.allc.tsv.gz \
	--chrom_size_path ../../../03.data/01.ref/mm10.chrom.sizes.nochrM.txt \
	--cpu 60 \
	--bin_length 10000000 

set +x
RUN
done
done
}

run_merge 
