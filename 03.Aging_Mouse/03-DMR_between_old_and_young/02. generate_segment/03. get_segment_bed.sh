#!/bin/sh

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

# start is 0-based, end is 1-based

cat ${indir}/tso.5hmc.segment.3cpg.bed | awk -v OFS="\t" '{print $1,$2-1,$3,$1"_"$2-1"_"$3}' > ../../../output/03.Aging_Mouse/03-DMRs_DHMRs/segments/tso.5hmC.segment.3cpg.0_base.bed
cat ${indir}/tso.5mc.segment.3cpg.bed | awk -v OFS="\t" '{print $1,$2-1,$3,$1"_"$2-1"_"$3}' > ../../../output/03.Aging_Mouse/03-DMRs_DHMRs/segments/tso.5mC.segment.3cpg.0_base.bed
