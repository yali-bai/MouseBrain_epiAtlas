#!/usr/bin/bash

## software location ##
seqkit='$PATH/miniconda3/bin/seqkit'
fastq_multx='/share/analysisdata/Methyl/Pipeline/TSO_2.0/scripts/fastq-multx'
python='/share/home/fany/miniconda3/envs/python2/bin/python'
src=../00.TSO-joint.pipeline

function sep_reads {
id=$1
barcode=$2
outdir=$3
indir=$4

R1=$indir/${id}.R1.fastq.gz
R2=$indir/${id}.R2.fastq.gz
cpus=30

sbatch <<RUN
#!/usr/bin/bash
#SBATCH -J $id
#SBATCH -e Output_s1/run_logs/${id}.sep.log
#SBATCH -o Output_s1/run_logs/${id}.sep.log
#SBATCH --cpus-per-task=$cpus
##SBATCH --exclude=node07
#SBATCH --mem=100G
##SBATCH --qos=temp
#SBATCH --partition=compute
#SBATCH --time=150:00:00
set -x
date


$fastq_multx \\
    -B $barcode \\
    $R1 $R2 \\
    -o $outdir/$id.%.R1.fastq -o $outdir/$id.%.R2.fastq \\
    -m 2 \\
    -d 0 \\
    -t $cpus

cat $barcode | cut -f 1 | while read bar
do
  echo $outdir/${id}.\${bar}.R1.fastq
  $seqkit subseq -r 20:-1 -j $cpus $outdir/${id}.\${bar}.R1.fastq >$outdir/${id}.\${bar}.R1.t.fastq
  mv $outdir/${id}.\${bar}.R1.t.fastq $outdir/${id}.\${bar}.R1.fastq
  gzip $outdir/${id}.\${bar}.R1.fastq
  gzip $outdir/${id}.\${bar}.R2.fastq
done

set +x
RUN
}

mkdir -p Output_s1/run_logs
mkdir -p Output_s1/sep_cell

run_script=$src/EM_main.sh

function run_sep {
indir=raw_data

ls  $indir|grep "R1.fastq"|sed 's/.R1.fastq.*//g' |grep -v "Balance"|while read ff # extract plate info
do
  echo $ff # for example: P56_Male_230131Mouse1_DeepCortex_EXP230202_TSO_5hmC_plate9
  sep_reads $ff $src/data/V2_TSO_Barcode_for_fastq_multx.220108.txt Output_s1/sep_cell $indir

done
}

function run_pip {
species=$1
prefix=$2

$run_script --indir $indir \
       --outdir $outdir \
       --prefix $prefix \
       --trim_type basic_trim \
       --read_len 20 \
       --species $species \
       --seq_type nextera \
       --partition compute \
       --thread 10 

}

run_sep # sep plate into cells









