#!/bin/bash
# "indir" is a custom input path.
# indir=""

python ${indir}/wgbs_tools/src/python/wgbs_tools.py segment --betas ../*5mC*beta --min_cpg 3 --max_bp 2000 -o tso.5mc.segment.3cpg.bed 
python ${indir}/wgbs_tools/src/python/wgbs_tools.py segment --betas ../*5hmC*beta --min_cpg 3 --max_bp 2000 -o tso.5hmc.segment.3cpg.bed

