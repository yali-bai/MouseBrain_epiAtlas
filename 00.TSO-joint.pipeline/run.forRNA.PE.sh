#!/usr/bin/bash


src_dir=/share/analysisdata/Methyl/workflow/TSO_HT/src/Pipeline/TSO_joint_RNA_src_240809/scripts
run_script=${src_dir}/TSO_JointRNA_forRNA.20231115.sh


function run_pip {
species=$1
prefix=$2

$run_script --indir $indir \
       --outdir $outdir \
       --prefix $prefix \
       --species $species \
       --seq_type nextera \
       --partition compute \
       --thread 10 

}


########RNA data analysis by RNA pipeline

outdir=Output_s1
indir=Output_s1/sep_cell


ls Output_s1/sep_cell/|grep R1|grep "_RNA_"|grep -v F|sed 's/.R1.fastq.gz//g'|while read prefix;

do
	# select by genome
	#run_pip hg38 $prefix
	run_pip mm10 $prefix

done

