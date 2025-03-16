#!/usr/bin/bash
# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""
bedGraphToBigWig=${indir}/UCSC_tools/bedGraphToBigWig
chrom_size=../../../../input/reference_genome/mm10.chrom.sizes.txt

cat ${indir}/hyper_DHMR.bedgraph | sort -k1,1 -k2,2n  > ${indir}/hyper_DHMR_sorted.bedgraph

${bedGraphToBigWig} ${indir}/hyper_DHMR_sorted.bedgraph ${chrom_size} ${outdir}/hyper_DHMR.bw