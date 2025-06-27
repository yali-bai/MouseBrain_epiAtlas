# "indir" is a custom input path.
# indir=""

ls ../../../output/03.Aging_Mouse/02-merge_allc/*tsv.gz|sed 's/merge_allc\///g'|while read a;do 
    echo "#!/bin/bash" > run.$a.sh
    echo "less merge_allc/$a | awk '\$4 ~ /^CG/' > $a.cg.txt" >> run.$a.sh
    echo "python ${indir}/merge.strand.py $a.cg.txt $a.cg.txt.bed" >> run.$a.sh 
    echo "python ${indir}/wgbs_tools/src/python/wgbs_tools.py bed2beta $a.cg.txt.bed --outdir ./ --genome mm10 -f" >> run.$a.sh
done
