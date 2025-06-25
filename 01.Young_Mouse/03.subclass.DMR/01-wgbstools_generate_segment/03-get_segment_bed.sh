# start is 0-based, end is 1-based

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

cat ${indir}/tso.5hmc.segment.3cpg.bed | awk -v OFS="\t" '{print $1,$2-1,$3,$1"_"$2-1"_"$3}' > ../../../output/01.Young_Mouse/03-wgbstools_generate_segment/TSO_5hmC_segment_3CpG.bed
cat ${indir}/tso.5hmc.segment.5cpg.bed | awk -v OFS="\t" '{print $1,$2-1,$3,$1"_"$2-1"_"$3}' > ../../../output/01.Young_Mouse/03-wgbstools_generate_segment/TSO_5hmC_segment_5CpG.bed
cat ${indir}/tso.5mc.segment.3cpg.bed | awk -v OFS="\t" '{print $1,$2-1,$3,$1"_"$2-1"_"$3}' > ../../../output/01.Young_Mouse/03-wgbstools_generate_segment/TSO_5mC_segment_3CpG.bed
cat ${indir}/tso.5mc.segment.5cpg.bed | awk -v OFS="\t" '{print $1,$2-1,$3,$1"_"$2-1"_"$3}' > ../../../output/01.Young_Mouse/03-wgbstools_generate_segment/TSO_5mC_segment_5CpG.bed
