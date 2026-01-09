# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

ls 04.subclass_merged_allc/*tsv.gz|sed 's/04.subclass_merged_allc\///g'|while read a;do 
    echo "#!/bin/bash" > run.$a.sh
    echo "less 04.subclass_merged_allc/$a | awk '\$4 ~ /^CG/' > $a.cg.txt" >> run.$a.sh
    echo "python ${indir}/merge.strand.py $a.cg.txt $a.cg.txt.bed" >> run.$a.sh 
    echo "python ${indir}/wgbs_tools.py bed2beta $a.cg.txt.bed --outdir ./ --genome mm10 -f" >> run.$a.sh
done
