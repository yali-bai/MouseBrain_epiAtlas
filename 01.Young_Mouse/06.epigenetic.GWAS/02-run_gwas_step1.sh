#!/bin/bash

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

gwas_ref="${indir}/gwas_ref"
input="${indir}/01_DHMR_hyper_bed"
output="${outdir}/03_output"
chain_dir=$gwas_ref/chains
liftOver_dir=$output/liftOver

mkdir -p $chain_dir $liftOver_dir

wget -P $chain_dir http://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg19.over.chain.gz
wget -P $chain_dir http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToMm10.over.chain.gz

mm10ToHg19_chain="$chain_dir/mm10ToHg19.over.chain.gz"
hg19ToMm10_chain="$chain_dir/hg19ToMm10.over.chain.gz"


for bed_file in ${input}/*.bed; do
    filename=$(basename "$bed_file" .bed)
    echo "Processing: $bed_file"

    # Step 1: mm10 to hg19
    echo "Step 1: Converting mm10 to hg19..."
    liftOver "$bed_file" "$mm10ToHg19_chain" \
        "${liftOver_dir}/${filename}_hg19.bed" \
        "${liftOver_dir}/${filename}_unmapped_hg19.bed" -minMatch=0.5
    if [ $? -ne 0 ]; then
        echo "Error: liftOver mm10 to hg19 failed for $bed_file" >&2
        continue
    fi

    # Step 2: hg19 back to mm10
    echo "Step 2: Converting hg19 back to mm10..."
    liftOver "${liftOver_dir}/${filename}_hg19.bed" "$hg19ToMm10_chain" \
        "${liftOver_dir}/${filename}_back_to_mm10.bed" \
        "${liftOver_dir}/${filename}_unmapped_mm10.bed" -minMatch=0.5
    if [ $? -ne 0 ]; then
        echo "Error: liftOver hg19 to mm10 failed for $bed_file" >&2
        continue
    fi

    # Step 3: Check consistency (at least 50% bases mapped back to mm10)
    echo "Step 3: Checking consistency..."
    bedtools intersect -a "$bed_file" \
                       -b "${liftOver_dir}/${filename}_back_to_mm10.bed" \
                       -f 1 > "${liftOver_dir}/${filename}_consistent_mm10.bed"
    if [ $? -ne 0 ]; then
        echo "Error: bedtools intersect failed for $bed_file" >&2
        continue
    fi

    # Step 4: Convert consistent mm10 regions to final hg19 coordinates
    echo "Step 4: Generating final hg19 coordinates..."
    liftOver "${liftOver_dir}/${filename}_consistent_mm10.bed" "$mm10ToHg19_chain" \
        "${liftOver_dir}/${filename}_final_hg19_unsorted.bed" \
        "${liftOver_dir}/${filename}_unmapped_final_hg19.bed" -minMatch=0.5
    if [ $? -ne 0 ]; then
        echo "Error: liftOver consistent mm10 to hg19 failed for $bed_file" >&2
        continue
    fi

    # Sort final hg19 regions
    echo "Sorting final hg19 regions..."
    sort -k1,1 -k2,2n "${liftOver_dir}/${filename}_final_hg19_unsorted.bed" \
        > "${liftOver_dir}/${filename}_final_hg19.bed"

    # Remove intermediate files
    echo "Cleaning up intermediate files..."
    rm -f "${liftOver_dir}/${filename}_hg19.bed" \
          "${liftOver_dir}/${filename}_unmapped_hg19.bed" \
          "${liftOver_dir}/${filename}_back_to_mm10.bed" \
          "${liftOver_dir}/${filename}_unmapped_mm10.bed" \
          "${liftOver_dir}/${filename}_consistent_mm10.bed" \
          "${liftOver_dir}/${filename}_final_hg19_unsorted.bed" \
          "${liftOver_dir}/${filename}_unmapped_final_hg19.bed"

    echo "Finished: ${liftOver_dir}/${filename}_final_hg19.bed"
done

echo "All processing completed."



#
input="${outdir}/03_output/liftOver"

temp_merged_bed="${input}/merged_subclass_bed_temp.bed"

echo "Merging subclass BED files..."
for bed_file in ${input}/*hyper_DHMR_final_hg19.bed; do
    cat ${bed_file} >> ${temp_merged_bed}
done

sort -k1,1 -k2,2n ${temp_merged_bed} > ${input}/final_merged_subclass_hyper_DHMR_final_hg19.bed
rm ${temp_merged_bed}
echo "Finished merging and sorting. Output saved to: ${input}/final_merged_subclass_hyper_DHMR_intersect_bins.bed"



