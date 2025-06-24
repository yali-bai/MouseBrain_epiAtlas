#!/usr/bin/bash
#Convert subclass merge allc files to bed and bw files
datatype=("5hm" "5m")

for dt in "${datatype[@]}";do
# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

bedGraphToBigWig= ${indir}/software/UCSC_tools/bedGraphToBigWig
chrom_size= ../../../04.data/01.ref/mm10.chrom.sizes.txt

ls $indir|grep "\.gz$"| sed 's/.allc.tsv.gz$//g' | while read ff 
do 
echo $ff 
sbatch <<RUN 
#!/usr/bin/bash 
#SBATCH -J ${dt}C_allc_bed_bw 
#SBATCH --mem=200G 
#SBATCH --partition=compute_fat 
#SBATCH -o ${outdir}/logs/${ff}.log 
#SBATCH -e ${outdir}/logs/${ff}.elog 

set -x 
zcat ${indir}/${ff}.allc.tsv.gz | awk '\$4 ~ "^CG" {OFS="\t"; print \$1, \$2, \$2+1, \$5/\$6}' | sort -k1,1 -k2,2n  > ${outdir}/CG/merge_bed/${ff}.CG.bed 
zcat ${indir}/${ff}.allc.tsv.gz | awk '\$4 !~ "^CG" {OFS="\t"; print \$1, \$2, \$2+1, \$5/\$6}' | sort -k1,1 -k2,2n  > ${outdir}/CH/merge_bed/${ff}.CH.bed 

${bedGraphToBigWig} ${outdir}/CG/merge_bed/${ff}.CG.bed ${chrom_size} ${outdir}/CG/merge_bw/${ff}.CG.bw 
${bedGraphToBigWig} ${outdir}/CH/merge_bed/${ff}.CH.bed ${chrom_size} ${outdir}/CH/merge_bw/${ff}.CH.bw

set +x

RUN
done
done



