#!/usr/bin/bash

indir=""
outdir="."
prefix="*"
thread=10
read_len=20
species=hg38

ddate=`date -R|awk '{print $2$3}'`

src_dir=/share/analysisdata/Methyl/workflow/TSO_HT/src/Pipeline/TSO_joint_RNA_src_240809


temp=`getopt -o i:o:s:p:m:t: --long indir:,outdir:,sep_type:,prefix:,thread:,read_len:,species:,seq_type:,partition: -- "$@"`
if [ $? != 0 ] ; then echo "terminating..." >&2 ; exit 1 ; fi
eval set -- "$temp"

while true ; do
        case "$1" in
                -i|--indir) echo "Input dir is $2"; indir=$2; shift 2;;
                -o|--outdir) echo "Output dir is $2" ;outdir=$2; shift 2;;
                -p|--prefix) prefix=$2; echo "Reads prefix is $prefix" ; shift 2;;
            		--read_len) read_len=$2; echo "minReadlen is $2" ; shift 2 ;;
            		--species) species=$2; echo "Species is $2"; shift 2 ;;
            		--partition) partition=$2; echo "partition is $2"; shift 2 ;;
            		--seq_type) seq_type=$2; echo "seq_type is $2"; shift 2 ;;
                -t|--thread) thread=$2; echo "multi-thread is $thread"; shift 2 ;;

                --) shift ; break ;;
                *) echo "internal error!" ; exit 1 ;;
        esac
done

if [[ ! -n $indir ]];
then
	echo "Please give input dir"
	exit
fi
if [[ ! -n $seq_type ]];then
        seq_type="nextera"
fi
if [[ ! -n $prefix ]];then prefix="*";fi

#set -x

mkdir -p $outdir
mkdir -p $outdir/run_logs
LOG_DIR=$outdir/run_logs

function pip_run {
pip_sh=$src_dir/scripts/pip.forDNA.20231128.sh

find $indir -maxdepth 1 -name "${prefix}*.f*q.gz" -print |sort --version-sort | xargs -L2 |column -t |while read reads
do
        sampleid=`echo $reads|sed 's/.*\///g'|sed 's/.R[12].*//g'|xargs`
        arr=(${reads// / })
        read1=${arr[0]}
        read2=${arr[1]}
	if [[ ! $read1 =~ "R1" ]];then
                echo $read1 "Read1 not exist"
                exit
        fi
	if [[ ! $read2 =~ "R2" ]];then
                echo $read2 "Read2 not exist"
                exit
        fi

	sample_size=$(stat -c "%s" $read1)
	sample_size=`echo "($sample_size * 2 + 50000000)/100000+40000"|bc|xargs`

	echo "   ---"$sampleid memery: ${sample_size}M

	sbatch -J $sampleid \
        -o $LOG_DIR/${sampleid}.pip.dna.${ddate}.log \
        -e $LOG_DIR/${sampleid}.pip.dna.${ddate}.log \
        --mem=${sample_size}M \
	--partition=${partition} \
        --cpus-per-task=$thread \
        --time=150:00:00 \
	      $pip_sh --read1 $read1 \
	      --read2 $read2 \
	      $temp

done
}

temp=`echo $temp|sed "s/'//g"`
echo $temp


pip_run


set +x
