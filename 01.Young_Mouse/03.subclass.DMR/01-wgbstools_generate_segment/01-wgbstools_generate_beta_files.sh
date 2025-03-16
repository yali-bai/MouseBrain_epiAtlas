# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

for a in `ls *tsv.gz`;do 
    echo "#!/bin/bash" > run.$a.sh
    echo "less $a | awk '\$4 ~ /^CG/' > $a.cg.txt" >> run.$a.sh
    echo "python ${indir}/merge.strand.py $a.cg.txt $a.cg.txt.bed" >> run.$a.sh 
    echo "python ${indir}/wgbs_tools/src/python/wgbs_tools.py bed2beta $a.cg.txt.bed --outdir ./ --genome mm10 -f" >> run.$a.sh
done
