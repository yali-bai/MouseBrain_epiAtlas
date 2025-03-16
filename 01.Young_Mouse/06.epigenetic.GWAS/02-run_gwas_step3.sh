#!/bin/bash

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

ldsc_software=${indir}/LDSC/ldsc
input=${indir}/03_output/liftOver
output=${indir}/03_output
gwas_ref=${indir}/gwas_ref
annot_dir=$output/ld_score
result_dir=$output/results
sumstats=$gwas_ref/sumstats
script_dir=$output/scripts
log_dir=$output/log

mkdir -p $result_dir

#
ls ${input}/*_hyper_DHMR_final_hg19.bed | grep -v "final_merged_subclass" | sed 's/_hyper_DHMR_final_hg19.bed//' | while read i
do 
    j=$(basename $i)
    printf "${j}\t${annot_dir}/${j}.,${annot_dir}/final_merged_subclass.\n"
done > ${result_dir}/hyperDHMR.ldcts


#
mkdir -p $gwas_ref/1000G_EUR/1000G_Phase3_weights_hm3_no_MHC_renamed

for i in {1..22}; do
    cp $gwas_ref/1000G_EUR/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.${i}.l2.ldscore.gz \
    $gwas_ref/1000G_EUR/1000G_Phase3_weights_hm3_no_MHC_renamed/weights.${i}.l2.ldscore.gz
done


#
for file in ${sumstats}/*.gz; do
    if [[ "$(basename "$file")" == *UKB_460K* ]]; then
        continue
    fi

    group=$(basename "$file" .sumstats.gz | sed 's/^LDSCORE_all_sumstats_PASS_//')

    cat <<EOT > ${script_dir}/step3_ldsc_slurm_${group}.sh
#!/bin/bash
#SBATCH --job-name=step3_ldsc_${group}
#SBATCH --partition=compute_fat
#SBATCH --output=${log_dir}/step3_ldsc_${group}.out
#SBATCH --error=${log_dir}/step3_ldsc_${group}.err
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=30GB

echo "Running LDSC for ${group}"

python ${ldsc_software}/ldsc.py \\
    --h2-cts ${file} \\
    --ref-ld-chr ${gwas_ref}/1000G_EUR/baseline_v1.1/baseline. \\
    --out ${result_dir}/${group}.hyperDHMR \\
    --ref-ld-chr-cts ${result_dir}/hyperDHMR.ldcts \\
    --w-ld-chr ${gwas_ref}/1000G_EUR/1000G_Phase3_weights_hm3_no_MHC_renamed/weights.

echo "Finished LDSC for ${group}"
EOT

    sbatch ${script_dir}/step3_ldsc_slurm_${group}.sh
done
