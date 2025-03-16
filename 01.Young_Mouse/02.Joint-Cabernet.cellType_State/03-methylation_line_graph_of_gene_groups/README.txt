01-allc_to_bed_bw.sh  :  Converting subclass merge allc files to bed and bw files. 
Run method: 
sed 's/\r//' 01-allc_to_bed_bw.sh
bash 01-allc_to_bed_bw.sh

02-gene_list_bed.r  :  Generating bed files for gene groups.

03-GetMeans.R  :  Calculating the average methylation rate of all genes in the gene groups.

04-scale.sh  : Calculating The methylation levels of different gene groups in genebody and upstream and downstream 2kb regions by computeMatrix.
Run method: 
conda activate /share/analysisdata/Methyl/public/rna
sed 's/\r//' 04-scale.sh
bash 04-scale.sh

05-tmp_turn.r  : Flipping the negative strand over.

06-plot.r  :  Drawing line diagrams of methylation level in different subclasses.


