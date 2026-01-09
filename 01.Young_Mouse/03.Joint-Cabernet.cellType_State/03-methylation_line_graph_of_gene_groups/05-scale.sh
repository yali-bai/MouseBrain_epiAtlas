#!/usr/bin/bash

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

datatype=("hm" "m")
mc_type=("CG" "CH") #"CH"
for dt in "${datatype[@]}";do
for mc in "${mc_type[@]}";do

bw_indir=./
bed_indir=./gene_list_bed
out_dir=./${dt}C/${mc}
mkdir -p ${out_dir}

ls ${bw_indir}/*.${dt}C.CGN-both.frac.bw|awk -F "/" '{print $NF}'|sed "s/.${dt}C.//"|sed "s/CGN-both.frac.bw//g" |while read class;

do
#ff="${class//-/.}"
ff="${class//[-|_]/.}"
ls ${bed_indir}/ |grep ${ff}|sed 's/_bed.bed//g'|while read group;
do 
echo ${group}
       
sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J ${dt}C_${group}
#SBATCH --mem=200G
#SBATCH --partition=compute_fat
#SBATCH -o logs/${dt}C_${group}.scale.${mc}.log
#SBATCH -e logs/${dt}C_${group}.scale.${mc}.log
#SBATCH --time=150:00:00

set -x

computeMatrix scale-regions \
      -S ${bw_indir}/${class}.${dt}C.${mc}N-both.frac.bw \
      -R ${bed_indir}/${group}_bed.bed \
      --regionBodyLength 5000 \
      --binSize 5 \
      --upstream 2000 \
      --downstream 2000 \
      --numberOfProcessors 50 \
      -o ${out_dir}/${dt}C_${group}_${mc}.genebody.gz

get tmp
Rscript ./04-GetMeans.R \
        ${out_dir}/${dt}C_${group}_${mc}.genebody.gz \
        ${out_dir}/${dt}C_${group}_${mc}.genebody.tmp \
        ${dt}C_${group}_${mc}




set +x
RUN
done
done
done
done


