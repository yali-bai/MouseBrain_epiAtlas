#!/usr/bin/bash

for dt in mC hmC
do

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J sep_h5ad_by_chr
#SBATCH --mem=1500g
#SBATCH --partition=compute_fat
##SBATCH --nodelist=fatnode02
#SBATCH -o logs/sep_h5ad_by_chr.${dt}.log
#SBATCH -e logs/sep_h5ad_by_chr.${dt}.log
#SBATCH --cpus-per-task=1
#SBATCH --time=150:00:00

set -x

python 02-sep_h5ad_by_chr.py $dt

set +x
RUN

done




