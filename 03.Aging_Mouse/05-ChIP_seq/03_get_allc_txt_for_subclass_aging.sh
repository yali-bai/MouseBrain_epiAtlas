# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

mkdir -p ${outdir}/subclass/5hmC

for file in ${indir}/*_old_5hmC.txt
do
    grep "Cortex" "$file" > "${outdir}/subclass/5hmC/$(basename $file)"
done

mkdir -p ${outdir}/subclass/5mC

for file in ${indir}/*_old_5mC.txt
do
    grep "Cortex" "$file" > "${outdir}/subclass/5mC/$(basename $file)"
done
