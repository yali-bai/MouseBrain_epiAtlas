#!/bin/bash

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

python ${indir}/wgbs_tools/src/python/wgbs_tools.py segment --betas *beta --min_cpg 3 --max_bp 2000 -o ${outdir}/tso.5hmc.segment.3cpg.bed 

python ${indir}/wgbs_tools/src/python/wgbs_tools.py segment --betas *beta --min_cpg 5 --max_bp 2000 -o ${outdir}/tso.5hmc.segment.5cpg.bed

python ${indir}/wgbs_tools/src/python/wgbs_tools.py beta_to_table ${outdir}/tso.5hmc.segment.3cpg.bed --betas *beta -o ${outdir}/tso.5hmc.segment.3cpg.bed.beta.txt

python ${indir}/wgbs_tools/src/python/wgbs_tools.py beta_to_table ${outdir}/tso.5hmc.segment.5cpg.bed --betas *beta -o ${outdir}/tso.5hmc.segment.5cpg.bed.beta.txt

#python filter.py tso.5hmc.segment.3cpg.bed.beta.txt > tso.5hmc.segment.3cpg.bed.beta.txt.calc.txt
