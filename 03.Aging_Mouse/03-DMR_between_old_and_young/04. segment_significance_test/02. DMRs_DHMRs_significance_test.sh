#!/usr/bin/bash

for datatype in 5mC 5mC
do

while IFS=$' ' read -r col1 col2
do
for start in $(seq 1 500 $col1)
do
max=`echo "$start+499"|bc|xargs`
if [ $max -gt $col1 ]; then
   end=$col1
else
   end=`echo "$start+499"|bc|xargs`
fi


sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J ss
#SBATCH --mem=150g
#SBATCH --partition=compute_pro
##SBATCH --nodelist=fatnode02
#SBATCH -o logs/subclass.DHMR_${datatype}.chr${col2}_${start}_${end}.log
#SBATCH -e logs/subclass.DHMR_${datatype}.chr${col2}_${start}_${end}.log
#SBATCH --cpus-per-task=1
#SBATCH --time=150:00:00

set -x

python DMR_subclass.py ${datatype}G $col2 $start $end

set +x
RUN

done
done < ../../../input/03-aging/subclass_${datatype}_chr_number.txt
done



