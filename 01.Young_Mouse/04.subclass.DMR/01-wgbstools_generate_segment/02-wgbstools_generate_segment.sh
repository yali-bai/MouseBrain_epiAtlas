#!/bin/bash
python /share/home/wanght/SOFTWARE/wgbs_tools/src/python/wgbs_tools.py segment --betas ../*.mC*beta --min_cpg 3 --max_bp 2000 --genome mm10 -o tso.5mc.segment.3cpg.bed 
python /share/home/wanght/SOFTWARE/wgbs_tools/src/python/wgbs_tools.py segment --betas ../*.hmC*beta --min_cpg 3 --max_bp 2000 --genome mm10 -o tso.5hmc.segment.3cpg.bed

