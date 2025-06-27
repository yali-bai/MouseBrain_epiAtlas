#!/usr/bin/bash

ages=("old" "young")
# 
datatype=("5mC" "true_5mC")
elements=('genebody' 'promoter')
mc_type=("CG" "CH")
class="subclass"
# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""
for dtype in "${datatype[@]}"; do
    for mtype in "${mc_type[@]}"; do
        for age in "${ages[@]}"; do
            # for class in "${classes[@]}"; do
                for ele in "${elements[@]}"; do
                    echo "$dtype,$ele, $mtype,$age,$class"
                    sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J ${dtype}_merge
#SBATCH --mem=100G
#SBATCH --partition=compute_pro
#SBATCH --time=150:00:00
#SBATCH -o ${outdir}/logs/${dtype}_${ele}_${mtype}_${age}_${class}_corrected.log
#SBATCH -e ${outdir}/logs/${dtype}_${ele}_${mtype}_${age}_${class}_corrected.log


set -x
python 01-data_merge.py "$age" "$class" "$dtype" "$ele" "$mtype" 
set +x

RUN

done
done
done
done
# done

