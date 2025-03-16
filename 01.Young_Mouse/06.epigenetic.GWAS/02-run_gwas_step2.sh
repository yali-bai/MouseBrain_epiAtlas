#!/bin/bash

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

ldsc_software=${indir}/LDSC/ldsc
gwas_ref=${indir}/gwas_ref
input=${indir}/03_output/liftOver
output=${indir}/03_output
script_dir=$output/scripts
log_dir=$output/log
plink=$gwas_ref/1000G_EUR/1000G_EUR_Phase3_plink
annot_dir=$output/ld_score

mkdir -p $script_dir $log_dir $annot_dir

for bed_file in ${input}/*_hyper_DHMR_final_hg19.bed; do
    cell_type=$(basename $bed_file _hyper_DHMR_final_hg19.bed)
    cat <<EOT > ${script_dir}/step2_slurm_${cell_type}.sh
#!/bin/bash
#SBATCH --job-name=step2_${cell_type}_ldsc
#SBATCH --partition=compute_fat
#SBATCH --output=${log_dir}/step2_${cell_type}_ldsc.out
#SBATCH --error=${log_dir}/step2_${cell_type}_ldsc.err
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem=30GB

echo "Running LD Score Regression for ${cell_type}"

# annot
for i in {1..22}; do
    python ${ldsc_software}/make_annot.py \\
        --bed-file ${bed_file} \\
        --bimfile ${plink}/1000G.EUR.QC.\${i}.bim \\
        --annot-file ${annot_dir}/${cell_type}.\${i}.annot.gz
done

# LD Score Regression
for i in {1..22}; do
    python ${ldsc_software}/ldsc.py \\
        --l2 \\
        --bfile ${plink}/1000G.EUR.QC.\${i} \\
        --ld-wind-cm 1 \\
        --annot ${annot_dir}/${cell_type}.\${i}.annot.gz \\
        --thin-annot \\
        --out ${annot_dir}/${cell_type}.\${i} \\
#        --print-snps ${gwas_ref}/hapmap3_snps/hm.\${i}.snp
        --print-snps ${indir}/listHM3.txt
done

echo "Finished LD Score Regression for ${cell_type}"
EOT

    sbatch ${script_dir}/step2_slurm_${cell_type}.sh
done

