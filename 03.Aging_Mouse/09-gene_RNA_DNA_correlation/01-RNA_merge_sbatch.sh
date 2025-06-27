#!/usr/bin/bash
# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

ages=("old" "young")
class="subclass"

for age in "${ages[@]}"; do
    # for class in "${classes[@]}"; do
        # 输出组合
        echo "$age,$class"
        sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J RNA_merge
#SBATCH --mem=100G
#SBATCH --partition=compute_pro
#SBATCH --time=150:00:00
#SBATCH -o ${outdir}/logs/RNA_${age}_${class}_corrected.log
#SBATCH -e ${outdir}/logs/RNA_${age}_${class}_corrected.log


set -x
python 01-RNA_merge.py "$age" "$class"
set +x

RUN

done
# done


