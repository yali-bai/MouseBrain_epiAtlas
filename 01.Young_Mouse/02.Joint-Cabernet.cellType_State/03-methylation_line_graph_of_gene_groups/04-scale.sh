#!/usr/bin/bash

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

datatype=("5hm" "5m")
mc_type=("CG" "CH")
for dt in "${datatype[@]}";do
for mc in "${mc_type[@]}";do

bw_indir=${indir}/${dt}c/${mc}/merge_bw
bed_indir=${indir}/gene_list_bed
out_dir=${outdir}/${dt}c/${mc}


ls ${bw_indir}/|sed "s/${dt}C_//"|sed "s/.Merge.${mc}.bw//g" |while read class;

do
ff="${class//-/.}"
ls ${bed_indir}/ |grep ${ff}|sed 's/_bed.bed//g'|while read group;
do 
echo ${group}
       
sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J ${dt}c_${group}
#SBATCH --mem=200G
#SBATCH --partition=compute_fat
#SBATCH -o ${outdir}/logs/${dt}c_${group}.scale.${mc}.log
#SBATCH -e ${outdir}/logs/${dt}c_${group}.scale.${mc}.log
#SBATCH --time=150:00:00

set -x

computeMatrix scale-regions \
      -S ${bw_indir}/${dt}C_${class}.Merge.${mc}.bw \
      -R ${bed_indir}/${group}_bed.bed \
      --regionBodyLength 5000 \
      --binSize 5 \
      --upstream 2000 \
      --downstream 2000 \
      --numberOfProcessors 50 \
      -o ${out_dir}/${dt}c_${group}_${mc}.genebody.gz

# get tmp
Rscript ./03-GetMeans.R \
        ${out_dir}/${dt}c_${group}_${mc}.genebody.gz \
        ${out_dir}/${dt}c_${group}_${mc}.genebody.tmp \
        ${dt}c_${group}_${mc}




set +x
RUN
done
done
done
done


