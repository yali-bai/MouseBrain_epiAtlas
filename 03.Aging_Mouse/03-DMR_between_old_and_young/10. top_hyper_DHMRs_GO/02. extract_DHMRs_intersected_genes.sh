#!/usr/bin/bash

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

for dir in subclass three_class
do
for subdir in all top500 top1000
do
ls ${indir}/$dir/$subdir/*.bed |while read file
do

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J ss
#SBATCH --mem=150g
#SBATCH --partition=compute_pro
##SBATCH --nodelist=fatnode02
##SBATCH -o logs/bedtools.$file.log
##SBATCH -e ${indir}/logs/subclass.DHMR_${datatype}.chr${col2}_${start}_${end}.log
#SBATCH --cpus-per-task=1
#SBATCH --time=150:00:00

set -x

bedtools intersect -a $file -b ../../../output/03.Aging_Mouse/mm10_gene.bed -wb >${file}.intersect_gene.bed

set +x
RUN

done 
done
done
