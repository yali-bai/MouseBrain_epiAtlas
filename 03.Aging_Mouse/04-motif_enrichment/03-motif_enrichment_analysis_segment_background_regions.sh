#!/usr/bin/bash

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

if [ "$#" -ne 4 ]; then
    echo "Usage: ${indir}/03-motif_enrichment_analysis_segment_background_regions.sh <input_csv> <background_file> <Specific_Hyper_Region> <group>"
    exit 1
fi

INPUT_CSV=$1 # /share/analysisdata/Methyl/public/TSO/yangfa/analysis/aging_mouse/motif_enrichment_latest_20250219/output/01-dhmr_hyper_hypo_DHMR_include_region.csv
BACKGROUND_FILE=$2 #  /share/analysisdata/Methyl/workflow/TSO_HT/Datadir/Mouse_Brain/data/RNA/integration/all_age/20241011_integration_by_subclass_marker/DMR/run_mcds.by_3cpg_segment_cell.all_age.20250115/06.chrom_density_plot/20250217_update.filter_before_calculating_p_value_adjust/DHMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.sorted.bed
SPECIFIC_REGION=$3  # Specific_Hyper_Region or Specific_Hypo_Region
GROUP=$4 # DHMR or DMR

MOTIF_DB=${indir}/select_TF_motif_db.txt

OUTPUT_DIR=${outdir}

REGION_TYPE=$(basename "$SPECIFIC_REGION" | grep -o 'Hyper\|Hypo')

mkdir -p ./logs

while IFS=: read -r motif_prefix motif_name; do

    sanitized_motif_name=$(echo "$motif_name" | sed 's/::/_/g' | xargs)
    
    echo "Motif: '${sanitized_motif_name}', Region: '${REGION_TYPE}'"
    
    MOTIF_TSV_GZ="${indir}/${motif_prefix}.tsv.gz"
    OUTPUT_CSV="${OUTPUT_DIR}/03-${sanitized_motif_name}_${GROUP}_${REGION_TYPE}_motif_enrichment_analysis_subclass.csv"

    sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J motif_enrichment_${motif_prefix}_${GROUP}_${SPECIFIC_REGION}
#SBATCH --mem=60G
#SBATCH --cpus-per-task=10
#SBATCH --partition=compute_pro
#SBATCH -o ${outdir}/logs/03-motif_enrichment_${motif_prefix}_${GROUP}_${SPECIFIC_REGION}.log
#SBATCH -e ${outdir}/logs/03-motif_enrichment_${motif_prefix}_${GROUP}_${SPECIFIC_REGION}.log
#SBATCH --time=150:00:00

set -x 

python 03-motif_enrichment_analysis_segment_background_regions.py "$INPUT_CSV" "$MOTIF_TSV_GZ" "$OUTPUT_CSV" "$SPECIFIC_REGION" "$BACKGROUND_FILE"

set +x
RUN

done < "$MOTIF_DB"







