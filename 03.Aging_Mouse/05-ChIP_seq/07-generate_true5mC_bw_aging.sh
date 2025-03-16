#!/usr/bin/bash

# "indir" is a custom input path, and "your_outdir" is a custom output path.
# indir=""
# your_outdir=""

indir_5hmC="${indir}/subclass/5hmC/CG"
indir_5mC="${indir}/subclass/5mC/CG"
outdir="${your_outdir}/subclass/true5mC/CG"
chrom_size="../../input/reference_genome/mm10.chrom.sizes.txt"

mkdir -p ${outdir}/logs

bedGraphToBigWig=${indir}/UCSC_tools/bedGraphToBigWig
bigWigToBedGraph=${indir}/UCSC_tools/bigWigToBedGraph

for file_5mC in ${indir_5mC}/*CGN-both.frac.bw; do

    base_name=$(basename "$file_5mC")
    
    file_5hmC="${indir_5hmC}/${base_name}"

    output_bw=$(basename "$file_5mC")
    
    output_name=$(basename "$file_5mC" | sed 's/.bw//g')

    sample_name=$(basename "$file_5mC" | sed 's/.CGN-both.cov.bw//')

    if [[ -f "$file_5hmC" ]]; then
        echo "Processing ${sample_name}..."

        sbatch <<EOF
#!/usr/bin/bash
#SBATCH -J true5mC_${sample_name}
#SBATCH --mem=60G
#SBATCH --partition=compute_pro
#SBATCH --cpus-per-task=30
#SBATCH -o ${outdir}/logs/07-${sample_name}.log
#SBATCH -e ${outdir}/logs/07-${sample_name}.log

set -x

bigwigCompare -b1 "${file_5mC}" -b2 "${file_5hmC}" \
             --operation subtract \
             -o "${outdir}/${output_name}.bed" \
             --binSize 1 \
             --numberOfProcessors 30 \
             --skipNonCoveredRegions \
             --outFileFormat bedgraph

tmp_bedgraph="${outdir}/${output_name}_sorted.bed"

awk '{ if (\$4 < 0) \$4=0; print }' "${outdir}/${output_name}.bed" | sort -k1,1 -k2,2n > "${tmp_bedgraph}"

${bedGraphToBigWig} "${tmp_bedgraph}" ${chrom_size} "${outdir}/${output_bw}"

set +x
EOF

    else
        echo "Warning: Missing 5hmC file for ${file_5mC}, skipping..."
    fi
done
