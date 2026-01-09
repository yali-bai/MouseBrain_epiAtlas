#!/bin/sh
# start is 0-based, end is 1-based

cat tso.5hmc.segment.3cpg.bed | awk -v OFS="\t" '{print $1,$2-1,$3,$1"_"$2-1"_"$3}' > ./tso.5hmC.segment.3cpg.0_base.bed
#cat /share/home/wanght/work/TSO_mm10_DMR_0905/alloc_to_beta/5hmc/call_dmr/tso.5hmc.segment.5cpg.bed | awk -v OFS="\t" '{print $1,$2-1,$3,$1"_"$2-1"_"$3}' > ./TSO_5hmC_segment_5CpG.bed
cat tso.5mc.segment.3cpg.bed | awk -v OFS="\t" '{print $1,$2-1,$3,$1"_"$2-1"_"$3}' > ./tso.5mC.segment.3cpg.0_base.bed
#cat /share/home/wanght/work/TSO_mm10_DMR_0905/alloc_to_beta/5mc/call_dmr/tso.5mc.segment.5cpg.bed | awk -v OFS="\t" '{print $1,$2-1,$3,$1"_"$2-1"_"$3}' > ./TSO_5mC_segment_5CpG.bed
