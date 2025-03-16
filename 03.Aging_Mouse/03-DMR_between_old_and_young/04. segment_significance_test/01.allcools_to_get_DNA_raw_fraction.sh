#!/usr/bin/bash

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J ss
#SBATCH --mem=1200g
#SBATCH --partition=compute_fat
##SBATCH --nodelist=fatnode02
#SBATCH -o logs/subclass.mcds_to_adata.log
#SBATCH -e logs/subclass.mcds_to_adata.log
#SBATCH --cpus-per-task=1
#SBATCH --time=150:00:00

set -x

python 01.allcools_to_get_DNA_raw_fraction.py

set +x
RUN





