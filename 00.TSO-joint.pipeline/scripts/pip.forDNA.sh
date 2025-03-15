#!/usr/bin/bash
## langurage
python=$PATH/python
Rscript=$PATH/Rscript
python3=$PATH/python

## soft Path
bismark=$PATH/miniconda3/bin/bismark
bowtie2=$PATH/miniconda3/bin/
picard=$PATH/miniconda3/bin/picard
bismark_methylation_extractor=$PATH/miniconda3/bin/bismark_methylation_extractor
samtools=$PATH/miniconda3/bin/samtools
bedtools=$PATH/miniconda3/bin/bedtools
sambamba=$PATH/miniconda3/bin/sambamba
qualimap=$PATH/miniconda3/bin/qualimap
featureCounts=$PATH/miniconda3/bin/featureCounts
trim_galore=$PATH/miniconda3/bin/trim_galore
cutadapt=$PATH/miniconda3/bin/cutadapt
allcools=$PATH/miniconda3/bin/allcools

## home made code
SRC_DIR=../../00.TSO-joint.pipeline

Remove_Gap_sepRNAreads=$SRC_DIR/scripts/RemoveGap_seperateRNAreadsbymCs.v2.py
Remove_Gap=$SRC_DIR/scripts/Remove_mC_filledGap.py
merge_CpG=$SRC_DIR/scripts/merge_PE_SE_report.py
## reference
#----------------------- Bismark Ref --------------------
hg38_bis=/share/analysisdata/Methyl/database/human/hg38/Bismark_hg38
mm10_bis=/share/analysisdata/Methyl/database/mouse/mm10/Index/Bismark_mm10
full_puc19_bis=/share/analysisdata/Methyl/database/puC19/puc19_bismark
Ecoli_bis=/share/analysisdata/Methyl/database/RNA/ecoli_bismark
lambda_bis=/share/analysisdata/Methyl/database/NEB_lambda/lambda_Bismark
msp1_bis=/share/analysisdata/Methyl/database/RNA/msp1_bwaindex
hmC5_bis=/share/analysisdata/Methyl/database/RNA/hmC5_bismark
clai_bis=/share/analysisdata/Methyl/database/clai/bismark_clai

#----------------------- Bwa Ref ----------------------
hg38_ref_fa=/share/analysisdata/Methyl/database/human/hg38/Index/bwa/hg38.fa
mm10_ref_fa=/share/analysisdata/Methyl/database/mouse/mm10/Index/bwa/mm10.genome.fa
lambda_ref_fa=/share/analysisdata/Methyl/database/NEB_lambda/BWAIndex/genome.fa
hmC5_ref_fa=/share/analysisdata/Methyl/database/RNA/5hmC_zymo/BWAIndex/5hmC_zymo.fa
full_puc19=/share/analysisdata/Methyl/database/puC19/BWAIndex/genome.fa
clai_fa=/share/analysisdata/Methyl/database/clai/clai.fa

#----------------------- Blacklist --------------------
hg38_blacklist=/share/analysisdata/Methyl/database/human/blacklist/hg38.blacklist.bed
mm10_blacklist=/share/analysisdata/Methyl/database/mouse/Blacklist/mm10.blacklist.bed

#----------------------- Gtf Ref --------------------
hg38_gtf=/share/analysisdata/Methyl/database/human/hg38/Annotation/gencode.v35.annotation.gtf
mm10_gtf=/share/analysisdata/Methyl/database/mouse/mm10/Annotation/gencode.vM18.annotation.gtf



## variable
declare -A Ref_bis
Ref_bis["hg38"]=$hg38_bis
Ref_bis["mm10"]=$mm10_bis
Ref_bis["hg38_mm10"]=$hg38_mm10_bis
Ref_bis["hg38_ecoli_puc19"]=$hg38_ecoli_puc19_bis
Ref_bis["mm10_ecoli_puc19"]=$mm10_ecoli_puc19_bis
Ref_bis["fullpuc19"]=$full_puc19_bis
Ref_bis["puc19"]=$song_puc19_bis
Ref_bis["ecoli"]=$Ecoli_bis
Ref_bis["lambda"]=$lambda_bis
Ref_bis["msp1"]=$msp1_bis
Ref_bis["hmC5"]=$hmC5_bis
Ref_bis["clai"]=$clai_bis

declare -A Ref_fa
Ref_fa["hg38"]=$hg38_ref_fa
Ref_fa["mm10"]=$mm10_ref_fa
Ref_fa["puc19"]=$cutpuc19_ref_fa
Ref_fa["lambda"]=$lambda_ref_fa
Ref_fa["ecoli"]=$cutecoli_ref_fa
Ref_fa["hg38_mm10"]=$hg38_mm10_ref
Ref_fa["hg38_mm10_ecoli_puc19"]=$hg38_mm10_ecoli_puc19
Ref_fa["hg38_ecoli_puc19"]=$hg38_ecoli_puc19
Ref_fa["mm10_ecoli_puc19"]=$mm10_ecoli_puc19
Ref_fa["fullpuc19"]=$full_puc19
Ref_fa["hmC5"]=$hmC5_ref_fa
Ref_fa["clai"]=$clai_fa

declare -A Blacklist
Blacklist["hg38"]=$hg38_blacklist
Blacklist["mm10"]=$mm10_blacklist

declare -A Ref_gtf
Ref_gtf["hg38"]=$hg38_gtf
Ref_gtf["mm10"]=$mm10_gtf

declare -A Ref_bed
Ref_bed["hg38"]=$hg38_gene_bed

m_threshold=0.5

module load fastqc
#module load bismark
module load java/1.8.0
module load samtools
module load R
module load pcre2


## Args
temp=`getopt -o i:s:p:m:o:t: --long indir:,prefix:,outdir:,read1:,read2:,read_len:,thread:,species:,seq_type:,partition:  -- "$@"`
if [ $? != 0 ] ; then echo "terminating..." >&2 ; exit 1 ; fi
eval set -- "$temp"

while true ; do
        case "$1" in
		-i|--indir) echo "Input dir is $2"; indir=$2; shift 2;;
                -p|--prefix) prefix=$2; echo "Reads prefix is $prefix" ; shift 2;;
                -o|--outdir) echo "Output dir is $2" ;outdir=$2; shift 2;;
                --read1)
                        case "$2" in
                                "") echo "option read1, no argument"; shift 2 ;;
                                *) Read1=$2; echo "Read1 is $Read1" ; shift 2 ;;
                        esac ;;
                --read2)
                        case "$2" in
                                "") echo "option read2, no argument"; shift 2 ;;
                                *) Read2=$2; echo "Read2 is $2" ; shift 2 ;;
                        esac ;;
                --read_len) read_len=$2; echo "minReadlen is $2" ; shift 2 ;;
            		--species) species=$2; echo "Species is $2"; shift 2 ;;
            		--partition) partition=$2; echo "partition is $2"; shift 2 ;;
            		--seq_type) seq_type=$2; echo "seq_type is $2"; shift 2 ;;
            		-t|--thread) thread=$2; echo "multi-thread is $thread"; shift 2 ;;

                --) shift ; break ;;
                *) echo "internal error!" ; exit 1 ;;
        esac
done


### -------------------------------------------------------------------
###### judge the args
###
if [ ! -s $Read1 ];then
        echo $Read1 "Read1 not exist"
        exit
fi
if [ ! -s $Read2 ];then
        echo $Read2 "Read2 not exist"
        exit
fi

if [[ ! -n $read_len ]];then read_len=20;fi
if [[ ! -n $trim_type ]];then trim_type=basic_trim;fi
if [[ ! -n $outdir ]];then  outdir=".";fi
if [[ ! -n $species ]];then echo "Please choose a species";exit;fi
if [[ ! -n $seq_type ]];then seq_type=nextera;fi
#if [[ ! -n $library ]];then library=hmC;fi



cores=`expr $SLURM_CPUS_PER_TASK / 5`
if [[ $cores == 0 ]];then
	 cores=1
fi

ref_fa=${Ref_fa[$species]}
ref_bis=${Ref_bis[$species]}
Blacklist=${Blacklist[$species]}
Ref_gtf=${Ref_gtf[$species]}
Ref_bed=${Ref_bed[$species]}

####chrom_size
if [[ $species == "hg38" ]];then
 chrom_size=/share/analysisdata/Methyl/workflow/TSO_HT/src/Pipeline/TSO_joint_RNA_src_240117/data/hg38_chrom_size.txt
fi
if [[ $species == "mm10" ]];then
 chrom_size=/share/analysisdata/Methyl/database/mouse/mm10/chrom_size/mm10.chrom.sizes.nochrM.txt
fi

sampleid=`echo $Read1|sed 's/.*\///g'|sed 's/[._]R1.f.*//g'|sed 's/.R1_val_1.fq.gz//g'`
echo $sampleid
id=$sampleid

if [[ ! -f $SRC_DIR/data/genome.length.txt ]];then
        genome_len=`less $ref_fa |grep -v ">"|awk 'BEGIN{sum=0}{sum+=length($1)}END{print sum}'|xargs`
        echo -e "${species}\t$genome_len" >> $SRC_DIR/data/genome.length.txt
else
        genome_len=`grep -w ${species} $SRC_DIR/data/genome.length.txt|awk '{print $2}'|uniq|xargs`

        if [[ $genome_len == "" ]];then
                genome_len=`less $ref_fa |grep -v ">"|awk 'BEGIN{sum=0}{sum+=length($1)}END{print sum}'|xargs`
                echo -e "${species}\t$genome_len" >> $SRC_DIR/data/genome.length.txt

        fi
fi
CpG_num=`grep -w ${species} $SRC_DIR/data/genome.CpG_num.txt|awk '{print $2}'|uniq|xargs`

###
### --------------- functions -----------------
###
function trim_array {
mkdir -p $outdir/${sampleid}/trim
t1=`ls $outdir/${sampleid}/trim/${sampleid}[._]R1_val_1.fq.gz`
t2=`ls $outdir/${sampleid}/trim/${sampleid}[._]R2_val_2.fq.gz`
echo "Reads need to be trimed: $Read1, $Read2"

if [[ -f $t1 ]] && [[ -f $t2 ]];then
        echo "Trim already finished"
        return
else
       rm -fr $outdir/${sampleid}/trim/${sampleid}*
fi

$trim_galore --paired \
      --phred33 \
      --length 20 \
      --retain_unpaired \
      --path_to_cutadapt $cutadapt \
      --output_dir $outdir/${sampleid}/trim \
      $Read1 $Read2

#if [[ ! -f $outdir/${sampleid}/trim/${sampleid}.R2_val_2_fastqc.zip ]];then
fastqc -t $SLURM_CPUS_PER_TASK $outdir/${sampleid}/trim/${sampleid}[._]R1_val_1.fq.gz
fastqc -t $SLURM_CPUS_PER_TASK $outdir/${sampleid}/trim/${sampleid}[._]R2_val_2.fq.gz
#fi

}

function trim_stat {
  mkdir -p $outdir/${sampleid}/stat
  Total_read_pair=`grep "Total reads processed" $outdir/${sampleid}/trim/${sampleid}[._]R1*_trimming_report.txt | sed 's/\s\s*/ /g' | cut -d " " -f 4 | sed 's/,//g'|xargs`
  Total_reads=$(($Total_read_pair * 2 ))
  Total_bases=`echo "$Total_reads * 150"|bc`
  unzip -n $outdir/${sampleid}/trim/${sampleid}[._]R*_val_1_fastqc.zip -d $outdir/${sampleid}/trim
  unzip -n $outdir/${sampleid}/trim/${sampleid}[._]R*_val_2_fastqc.zip -d $outdir/${sampleid}/trim

  Clean_reads=$(($(grep "Total Sequences" $outdir/${sampleid}/trim/${sampleid}.R1_val_1_fastqc/fastqc_data.txt |sed 's/Total Sequences\s//g'|xargs) + $(grep "Total Sequences" $outdir/${sampleid}/trim/${sampleid}.R2_val_2_fastqc/fastqc_data.txt |sed 's/Total Sequences\s//g'|xargs)))
  Clean_read_base=$(($(grep "Total written (filtered)" $outdir/${sampleid}/trim/${sampleid}[._]R1*_trimming_report.txt |sed 's/.*:\s//g'|awk '{print $1}'|sed 's/,//g'|xargs) + $(grep "Total written (filtered)" $outdir/${sampleid}/trim/${sampleid}[._]R2*_trimming_report.txt |sed 's/.*:\s//g'|awk '{print $1}'|sed 's/,//g'|xargs)))

  Filter_reads_r=`echo "scale=2;100*(1-$Clean_reads/$Total_reads)"|bc`
  Filter_base_r=`echo "scale=2;100*(1-$Clean_read_base/$Total_bases)"|bc`

  raw_depth=`echo "scale=6;$Total_bases / $genome_len "|bc|xargs`
  clean_depth=`echo "scale=6;$Clean_read_base / $genome_len "|bc|xargs`
  echo -e "${sampleid}\t" \
    "${Total_reads}\t" \
    "${Clean_reads}\t" \
    "${Filter_reads_r}%\t" \
    "${Filter_base_r}%\t" \
    "${raw_depth}\t" \
    "${clean_depth}"  > $outdir/${sampleid}/stat/${sampleid}.trim.stat.txt
}

function mapping {
trim_Read1=$outdir/${sampleid}/trim/${sampleid}.R1_val_1.fq.gz
trim_Read2=$outdir/${sampleid}/trim/${sampleid}.R2_val_2.fq.gz
mkdir -p $outdir/${sampleid}/align

$bismark --multicore $cores \
        --fastq --non_directional --unmapped \
        --nucleotide_coverage \
        --path_to_bowtie $bowtie2 --bowtie2 \
        --genome_folder $ref_bis \
        --output_dir $outdir/${sampleid}/align \
        --temp_dir $outdir/${sampleid}/align/temp_bismark_${sampleid} \
        -1 $trim_Read1 -2 $trim_Read2

outbam=`ls $outdir/${sampleid}/align/${sampleid}.R1_val_1_bismark_bt2_pe.bam|xargs`
if [[ $species == "hg38" ]] || [[ $species == "mm10" ]]
then
$samtools view -bh -q 1 -F 4 $outbam |$samtools sort -@ $SLURM_CPUS_PER_TASK |$bedtools intersect -a - -b $Blacklist -v > $outdir/${sampleid}/align/${sampleid}.${species}.PE.bam
else
$samtools view -bh -q 1 -F 4 $outbam |$samtools sort -@ $SLURM_CPUS_PER_TASK > $outdir/${sampleid}/align/${sampleid}.${species}.PE.bam
fi
$samtools index $outdir/${sampleid}/align/${sampleid}.${species}.PE.bam
$samtools stats $outdir/${sampleid}/align/${sampleid}.${species}.PE.bam > $outdir/${sampleid}/align/${sampleid}.${species}.PE.bam.stat
$samtools depth $outdir/${sampleid}/align/${sampleid}.${species}.PE.bam > $outdir/${sampleid}/align/${sampleid}.${species}.PE.depth

##########SE
umapped_read1=$outdir/${sampleid}/align/${sampleid}.R1_val_1.fq.gz_unmapped_reads_1.fq.gz
umapped_read2=$outdir/${sampleid}/align/${sampleid}.R2_val_2.fq.gz_unmapped_reads_2.fq.gz

$bismark --multicore $cores \
    --fastq --non_directional --unmapped \
    --nucleotide_coverage \
    --path_to_bowtie $bowtie2 --bowtie2 \
    --genome_folder $ref_bis \
    --output_dir $outdir/${sampleid}/align \
    --temp_dir $outdir/${sampleid}/align/temp_bismark_${sampleid} \
    $umapped_read1

$bismark --multicore $cores \
    --fastq --non_directional --unmapped \
    --nucleotide_coverage \
    --path_to_bowtie $bowtie2 --bowtie2 \
    --genome_folder $ref_bis \
    --output_dir $outdir/${sampleid}/align \
    --temp_dir $outdir/${sampleid}/align/temp_bismark_${sampleid} \
    $umapped_read2
$samtools merge $outdir/${sampleid}/align/${sampleid}.${species}.merged.SE.bam $outdir/${sampleid}/align/${sampleid}*_unmapped_reads_*_bismark_bt2.bam
outbam=$outdir/${sampleid}/align/${sampleid}.${species}.merged.SE.bam
if [[ $species == "hg38" ]] || [[ $species == "mm10" ]]
then
$samtools view -bh -q 1 -F 4 $outbam |$samtools sort -@ $SLURM_CPUS_PER_TASK |$bedtools intersect -a - -b $Blacklist -v > $outdir/${sampleid}/align/${sampleid}.${species}.SE.bam
else
$samtools view -bh -q 1 -F 4 $outbam |$samtools sort -@ $SLURM_CPUS_PER_TASK > $outdir/${sampleid}/align/${sampleid}.${species}.SE.bam
fi

$samtools index $outdir/${sampleid}/align/${sampleid}.${species}.SE.bam
$samtools stats $outdir/${sampleid}/align/${sampleid}.${species}.SE.bam > $outdir/${sampleid}/align/${sampleid}.${species}.SE.bam.stat
$samtools depth $outdir/${sampleid}/align/${sampleid}.${species}.SE.bam > $outdir/${sampleid}/align/${sampleid}.${species}.SE.depth


echo "Calculate the fragment length"
$picard CollectInsertSizeMetrics \
        I=$outdir/${sampleid}/align/${sampleid}.${species}.PE.bam \
        O=$outdir/${sampleid}/align/${sampleid}.${species}.PE.insert_size_metrics.txt \
        H=$outdir/${sampleid}/align/${sampleid}.${species}.PE.insert_size_histogram.pdf \
        M=0.5
date

$samtools merge $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.bam $outdir/${sampleid}/align/${sampleid}.${species}.PE.bam $outdir/${sampleid}/align/${sampleid}.${species}.merged.SE.bam
$samtools sort -@ $SLURM_CPUS_PER_TASK $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.bam > $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam
$samtools index $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam
$samtools stats $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam > $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam.stat
$samtools depth $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam > $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam.depth

}

function remove_Gap_sepReads {
  Read_type=$1  ## PE or SE
  use_species=$2 ## hg38/mm10/
  inbam=$outdir/${sampleid}/align/${sampleid}.${use_species}.${Read_type}.bam
  out_prefix=$outdir/${sampleid}/align/${sampleid}.${use_species}.${Read_type}
  date
  $python $Remove_Gap_sepRNAreads \
          -i $inbam \
          -o ${out_prefix}.bam \
          -r $ref_fa \
          --thread $SLURM_CPUS_PER_TASK \
          --Read_type $Read_type
  date

  ## dna reads
  $samtools sort -@ $SLURM_CPUS_PER_TASK ${out_prefix}.dna.bam > ${out_prefix}.dna.sorted.bam
  $samtools index ${out_prefix}.dna.sorted.bam
  $samtools stats ${out_prefix}.dna.sorted.bam > ${out_prefix}.dna.sorted.bam.stat
  ## rna reads
  $samtools sort -@ $SLURM_CPUS_PER_TASK ${out_prefix}.rna.bam > ${out_prefix}.rna.sorted.bam
  $samtools index ${out_prefix}.rna.sorted.bam
  $samtools stats ${out_prefix}.rna.sorted.bam > ${out_prefix}.rna.sorted.bam.stat

echo "check error"
  date

}

function Remove_Duplicates {
  Read_type=$1  ## PE or SE
  use_species=$2 ## hg38/mm10/
  data_type=$3 ## dna/rna
  out_prefix=$outdir/${sampleid}/align/${sampleid}.${use_species}.${Read_type}.${data_type}

  $picard MarkDuplicates \
      I=${out_prefix}.sorted.bam \
      O=${out_prefix}.sorted.rmdup.bam \
      M=${out_prefix}.sorted.rmdup.txt \
      READ_NAME_REGEX=null \
      REMOVE_DUPLICATES=true
  date
  $samtools index ${out_prefix}.sorted.rmdup.bam
  $samtools sort -n -@ $SLURM_CPUS_PER_TASK ${out_prefix}.sorted.rmdup.bam > ${out_prefix}.sorted.rmdup.sortn.bam
  date
}

function call_methyl {
  mkdir -p $outdir/${sampleid}/methyl
  species=$1 ## hg38/mm10/
  data_type=$2 ## dna/rna
  echo "Starting call the CpG methylation..."
  date
#:<<MULTILINECOMMENT
  $bismark_methylation_extractor -p \
   --multicore $cores \
   --comprehensive \
   --no_overlap \
   --bedGraph \
   --counts \
   --buffer_size 20G \
   --report \
   --cytosine_report \
   --CX_context \
   --genome_folder $ref_bis \
   -o $outdir/${sampleid}/methyl \
   $outdir/${sampleid}/align/${sampleid}.${species}.PE.${data_type}.sorted.rmdup.sortn.bam

  date

  $bismark_methylation_extractor -s \
     --multicore $cores \
     --comprehensive \
     --no_overlap \
     --bedGraph \
     --counts \
     --buffer_size 20G \
     --report \
     --cytosine_report \
     --CX_context \
     --genome_folder $ref_bis \
     -o $outdir/${sampleid}/methyl \
     $outdir/${sampleid}/align/${sampleid}.${species}.SE.${data_type}.sorted.rmdup.sortn.bam

  date

#MULTILINECOMMENT
  awk -F "\t" '$6=="CG"' $outdir/${sampleid}/methyl/${sampleid}.${species}.PE.${data_type}.sorted.rmdup.sortn.CX_report.txt > $outdir/${sampleid}/methyl/${sampleid}.${species}.PE.${data_type}.sorted.rmdup.sortn.CpG_report.txt
  awk -F "\t" '$6=="CG"' $outdir/${sampleid}/methyl/${sampleid}.${species}.SE.${data_type}.sorted.rmdup.sortn.CX_report.txt > $outdir/${sampleid}/methyl/${sampleid}.${species}.SE.${data_type}.sorted.rmdup.sortn.CpG_report.txt

  sort --parallel=8 -k1,1 -k2n,2 -k3 -k6 -k7 $outdir/${sampleid}/methyl/${sampleid}.${species}.PE.${data_type}.sorted.rmdup.sortn.CpG_report.txt > $outdir/${sampleid}/methyl/${sampleid}.${species}.PE.${data_type}.sorted.rmdup.sortn.CpG_report.sort.txt
  sort --parallel=8 -k1,1 -k2n,2 -k3 -k6 -k7 $outdir/${sampleid}/methyl/${sampleid}.${species}.SE.${data_type}.sorted.rmdup.sortn.CpG_report.txt > $outdir/${sampleid}/methyl/${sampleid}.${species}.SE.${data_type}.sorted.rmdup.sortn.CpG_report.sort.txt

  $python3 $merge_CpG $outdir/${sampleid}/methyl/${sampleid}.${species}.SE.${data_type}.sorted.rmdup.sortn.CpG_report.sort.txt $outdir/${sampleid}/methyl/${sampleid}.${species}.PE.${data_type}.sorted.rmdup.sortn.CpG_report.sort.txt  $outdir/${sampleid}/methyl/${sampleid}.${species}.${data_type}.merged.sortn.CpG_report.txt
  CpG_f=$outdir/${sampleid}/methyl/${sampleid}.${species}.${data_type}.merged.sortn.CpG_report.txt

  awk 'BEGIN{chr="Chr";pos="Pos";methyl="Methyl";cov="UMethyl"}{OFS="\t";if($3=="-"){$2=$2-1};if(pos==$2){methyl+=$4;cov+=$5}else{sum=methyl+cov;if(sum==0){rr=""}else{rr=methyl/sum};print chr,pos,methyl,cov,rr;chr=$1;pos=$2;methyl=$4;cov=$5}}' $CpG_f  > $outdir/${sampleid}/methyl/${sampleid}.${species}.${data_type}.CpG_methy.txt


  #bis_stat=$outdir/${sampleid}/methyl/${sampleid}.${species}.correctCH.sorted.rmdup.sortn_splitting_report.txt
  #CGmethy_stat=`grep "CpG" $bis_stat |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{CpG=$1+$2; cpg_r=100*$1/CpG;print CpG,cpg_r"%"}'|sed 's/ /\t/g'`
  #CHmethy_stat=`grep "CH[HG]" $bis_stat |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{CHG=$1+$3;CHH=$2+$4;chg_r=100*$1/CHG;chh_r=100*$2/CHH;print chg_r"%",chh_r"%"}'|sed 's/ /\t/g'`
  bis_stat_SE=$outdir/${sampleid}/methyl/${sampleid}.${species}.PE.${data_type}.sorted.rmdup.sortn_splitting_report.txt
  bis_stat_PE=$outdir/${sampleid}/methyl/${sampleid}.${species}.SE.${data_type}.sorted.rmdup.sortn_splitting_report.txt

  CGmethy_SE_num=`grep "CpG" $bis_stat_SE |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{print $1}'|sed 's/ /\t/g'`
  CG_total_SE_num=`grep "CpG" $bis_stat_SE |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{print $1+$2}'|sed 's/ /\t/g'`

  CGmethy_PE_num=`grep "CpG" $bis_stat_PE |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{print $1}'|sed 's/ /\t/g'`
  CG_total_PE_num=`grep "CpG" $bis_stat_PE |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{print $1+$2}'|sed 's/ /\t/g'`

  CGmethy_num=`echo "$CGmethy_SE_num+$CGmethy_PE_num"|bc`
  CG_total_num=`echo "$CG_total_SE_num+$CG_total_PE_num"|bc`
  CGmethy_stat=`echo "scale=4;100*$CGmethy_num/$CG_total_num"|bc`

  CHGmethy_SE_num=`grep "CHG" $bis_stat_SE |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{print $1}'|sed 's/ /\t/g'`
  CHG_total_SE_num=`grep "CHG" $bis_stat_SE |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{print $1+$2}'|sed 's/ /\t/g'`

  CHGmethy_PE_num=`grep "CHG" $bis_stat_PE |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{print $1}'|sed 's/ /\t/g'`
  CHG_total_PE_num=`grep "CHG" $bis_stat_PE |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{print $1+$2}'|sed 's/ /\t/g'`

  CHGmethy_num=`echo "$CHGmethy_SE_num+$CHGmethy_PE_num"|bc`
  CHG_total_num=`echo "$CHG_total_SE_num+$CHG_total_PE_num"|bc`
  CHGmethy_stat=`echo "scale=4;100*$CHGmethy_num/$CHG_total_num"|bc`

  CHHmethy_SE_num=`grep "CHH" $bis_stat_SE |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{print $1}'|sed 's/ /\t/g'`
  CHH_total_SE_num=`grep "CHH" $bis_stat_SE |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{print $1+$2}'|sed 's/ /\t/g'`

  CHHmethy_PE_num=`grep "CHH" $bis_stat_PE |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{print $1}'|sed 's/ /\t/g'`
  CHH_total_PE_num=`grep "CHH" $bis_stat_PE |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{print $1+$2}'|sed 's/ /\t/g'`

  CHHmethy_num=`echo "$CHHmethy_SE_num+$CHHmethy_PE_num"|bc`
  CHH_total_num=`echo "$CHH_total_SE_num+$CHH_total_PE_num"|bc`
  CHHmethy_stat=`echo "scale=4;100*$CHHmethy_num/$CHH_total_num"|bc`

  CpG_n=`awk '($3+$4)>0' $outdir/${sampleid}/methyl/${sampleid}.${species}.${data_type}.CpG_methy.txt|wc -l|xargs`
  CpG_Cov=`echo "scale=4;100*$CpG_n/$CpG_num"|bc`
  echo -e "${CpG_Cov}%\t${CG_total_num}\t${CGmethy_stat}%\t${CHGmethy_stat}%\t${CHHmethy_stat}%"  > $outdir/${sampleid}/stat/${sampleid}.${species}.${data_type}.Methyl.stat.txt

}

function FeatureCount_Qualimap {
species=$1 ## hg38/mm10/
data_type=$2 ## dna/rna
mkdir -p $outdir/${sampleid}/feature_count

$featureCounts -p --countReadPairs \
-T $SLURM_CPUS_PER_TASK \
-a $Ref_gtf \
-t exon \
-g gene_id \
-R BAM \
-o $outdir/${sampleid}/feature_count/${sampleid}.${species}.PE.${data_type}.exon_counts.txt \
$outdir/${sampleid}/align/${sampleid}.${species}.PE.${data_type}.sorted.rmdup.bam

$featureCounts -p --countReadPairs \
-T $SLURM_CPUS_PER_TASK \
-a $Ref_gtf \
-t transcript \
-g gene_id \
-R BAM \
-o $outdir/${sampleid}/feature_count/${sampleid}.${species}.PE.${data_type}.transcript_counts.txt \
$outdir/${sampleid}/align/${sampleid}.${species}.PE.${data_type}.sorted.rmdup.bam

$featureCounts -p --countReadPairs \
-T $SLURM_CPUS_PER_TASK \
-a $Ref_gtf \
-t gene \
-g gene_id \
-R BAM \
-o $outdir/${sampleid}/feature_count/${sampleid}.${species}.PE.${data_type}.gene_counts.txt \
$outdir/${sampleid}/align/${sampleid}.${species}.PE.${data_type}.sorted.rmdup.bam

grep -v "Geneid" $outdir/${sampleid}/feature_count/${sampleid}.${species}.PE.${data_type}.gene_counts.txt |awk '$7>0'|awk '{print $6}' > $outdir/${sampleid}/feature_count/${sampleid}.${species}.PE.${data_type}.mapped_gene_length.txt


mkdir -p $outdir/${sampleid}/Qualimap_for_QC

read1=$outdir/sep_cell/${sampleid}.R1.fastq.gz
sample_size=$(stat -c "%s" $read1)
#sample_size=`echo "($sample_size * 2 + 1000000000)/100000"|bc|xargs`
sample_size=`echo "($sample_size * 2 + 50000000)/100000+40000"|bc|xargs`


${qualimap} rnaseq -bam $outdir/${sampleid}/align/${sampleid}.${species}.PE.${data_type}.sorted.rmdup.bam \
        --java-mem-size=${sample_size}"M" \
        -gtf $Ref_gtf \
        -pe \
        -outdir $outdir/${sampleid}/Qualimap_for_QC/ \
        -outfile ${sampleid}.${species}.PE.${data_type}.rnaseq_report.pdf \
        -oc ${sampleid}.${species}.PE.${data_type}.count.txt
#
date
FeatureCounts_Exon=$(grep Assigned $outdir/${sampleid}/feature_count/${sampleid}.${species}.PE.${data_type}.exon_counts.txt.summary|cut -f2|awk '{print $1*2}')
FeatureCounts_Gene=$(grep Assigned $outdir/${sampleid}/feature_count/${sampleid}.${species}.PE.${data_type}.gene_counts.txt.summary|cut -f2|awk '{print $1*2}')
FeatureCounts_Intron=`echo "$FeatureCounts_Gene - $FeatureCounts_Exon"|bc`
FeatureCounts_IntergenicRegion=$(grep "Unassigned_NoFeatures" $outdir/${sampleid}/feature_count/${sampleid}.${species}.PE.${data_type}.gene_counts.txt.summary|cut -f2|awk '{print $1*2}')

FeatureCounts_Total=`echo "$FeatureCounts_Gene + $FeatureCounts_IntergenicRegion"|bc`
FeatureCounts_Exonic=`echo "scale=2;100*$FeatureCounts_Exon/$FeatureCounts_Total"|bc`
FeatureCounts_Intronic=`echo "scale=2;100*$FeatureCounts_Intron/$FeatureCounts_Total"|bc`
FeatureCounts_Intergenic=`echo "scale=2;100*$FeatureCounts_IntergenicRegion/$FeatureCounts_Total"|bc`

Exon_gene_number=`grep -v Geneid $outdir/${sampleid}/feature_count/${sampleid}.${species}.PE.${data_type}.exon_counts.txt|awk '$7 > 0'|wc -l`
Gene_gene_number=`grep -v Geneid $outdir/${sampleid}/feature_count/${sampleid}.${species}.PE.${data_type}.gene_counts.txt |awk '$7 > 0'|wc -l`

Qualimap_Exon=`grep exonic $outdir/${sampleid}/Qualimap_for_QC/rnaseq_qc_results.txt|cut -d '(' -f1|sed 's/.*=[[:space:]]*//g'|sed 's/,//g'`
Qualimap_Exonic=`grep exonic $outdir/${sampleid}/Qualimap_for_QC/rnaseq_qc_results.txt|cut -d '(' -f2|sed 's/%.*//g'`
Qualimap_Intron=`grep intronic $outdir/${sampleid}/Qualimap_for_QC/rnaseq_qc_results.txt|cut -d '(' -f1|sed 's/.*=[[:space:]]*//g'|sed 's/,//g'`
Qualimap_Intronic=`grep intronic $outdir/${sampleid}/Qualimap_for_QC/rnaseq_qc_results.txt|cut -d '(' -f2|sed 's/%.*//g'`
Qualimap_IntergenicRegion=`grep intergenic $outdir/${sampleid}/Qualimap_for_QC/rnaseq_qc_results.txt|cut -d '(' -f1|sed 's/.*=[[:space:]]*//g'|sed 's/,//g'`
Qualimap_Intergenic=`grep intergenic $outdir/${sampleid}/Qualimap_for_QC/rnaseq_qc_results.txt|cut -d '(' -f2|sed 's/%.*//g'`
Qualimap_overlap_Exon=`grep overlapping $outdir/${sampleid}/Qualimap_for_QC/rnaseq_qc_results.txt|cut -d '(' -f1|sed 's/.*=[[:space:]]*//g'|sed 's/,//g'`
Qualimap_overlap_Exonic=`grep overlapping $outdir/${sampleid}/Qualimap_for_QC/rnaseq_qc_results.txt|cut -d '(' -f2|sed 's/%.*//g'`

echo -e "$sampleid\t${species}\t" \
"${FeatureCounts_Exon}\t" \
"${FeatureCounts_Exonic}%\t" \
"${FeatureCounts_Intron}\t" \
"${FeatureCounts_Intronic}%\t" \
"${FeatureCounts_IntergenicRegion}\t" \
"${FeatureCounts_Intergenic}%\t" \
"${Exon_gene_number}\t" \
"${Gene_gene_number}\t" \
"${Qualimap_Exon}\t" \
"${Qualimap_Exonic}%\t" \
"${Qualimap_Intron}\t" \
"${Qualimap_Intronic}%\t" \
"${Qualimap_IntergenicRegion}\t" \
"${Qualimap_Intergenic}%\t" \
"${Qualimap_overlap_Exon}\t" \
"${Qualimap_overlap_Exonic}%\t" \
"${bias}\t" \
> $outdir/${sampleid}/stat/FeatureCount_Qualimap.${species}.${data_type}.txt

}

function spike_in {
  s_species=$1
  ref_bis=${Ref_bis[$s_species]}
  ref_fa=${Ref_fa[$s_species]}
  echo $ref_bis
  echo $ref_fa
  mkdir -p $outdir/${sampleid}/align/$s_species
  Read1=$outdir/${sampleid}/align/${sampleid}.R1_val_1.fq.gz_unmapped_reads_1.fq.gz
  Read2=$outdir/${sampleid}/align/${sampleid}.R2_val_2.fq.gz_unmapped_reads_2.fq.gz

  $bismark --multicore $cores \
      --fastq --non_directional \
      --nucleotide_coverage \
      --path_to_bowtie $bowtie2 --bowtie2 \
      --genome_folder $ref_bis \
      --output_dir $outdir/${sampleid}/align/$s_species \
      --temp_dir $outdir/${sampleid}/align/${s_species}/temp_bismark_${sampleid} \
      -1 $Read1 -2 $Read2

  outbam=`ls $outdir/${sampleid}/align/${s_species}/${sampleid}.R1_val_1*_bismark_bt2_pe.bam|xargs`
  $samtools view -bh -q 1 -F 4 $outbam |$samtools sort -@ $SLURM_CPUS_PER_TASK -O BAM -o $outdir/${sampleid}/align/${sampleid}.${s_species}.PE.bam
  $samtools index $outdir/${sampleid}/align/${sampleid}.${s_species}.PE.bam
  $samtools stats $outdir/${sampleid}/align/${sampleid}.${s_species}.PE.bam > $outdir/${sampleid}/align/${sampleid}.${s_species}.bam.stat


  ## Only Remove filledGap
    $python $Remove_Gap \
            -i $outdir/${sampleid}/align/${sampleid}.${s_species}.PE.bam \
            -o $outdir/${sampleid}/align/${sampleid}.${s_species}.PE.removeGap.bam \
            -r $ref_fa \
            --thread $SLURM_CPUS_PER_TASK \
            --Read_type PE

    $samtools sort -@ $SLURM_CPUS_PER_TASK $outdir/${sampleid}/align/${sampleid}.${s_species}.PE.removeGap.bam > $outdir/${sampleid}/align/${sampleid}.${s_species}.PE.removeGap.sorted.bam
    $samtools index $outdir/${sampleid}/align/${sampleid}.${s_species}.PE.removeGap.sorted.bam
    $samtools stats $outdir/${sampleid}/align/${sampleid}.${s_species}.PE.removeGap.sorted.bam > $outdir/${sampleid}/align/${sampleid}.${s_species}.PE.removeGap.sorted.bam.stat

    out_prefix=$outdir/${sampleid}/align/${sampleid}.${s_species}.PE.removeGap
    $picard MarkDuplicates \
        I=${out_prefix}.sorted.bam \
        O=${out_prefix}.sorted.rmdup.bam \
        M=${out_prefix}.sorted.rmdup.txt \
        READ_NAME_REGEX=null \
        REMOVE_DUPLICATES=true

    date
    $samtools index ${out_prefix}.sorted.rmdup.bam
    $samtools sort -n -@ $SLURM_CPUS_PER_TASK ${out_prefix}.sorted.rmdup.bam > ${out_prefix}.sorted.rmdup.sortn.bam

    $bismark_methylation_extractor -p \
       --multicore $cores \
       --comprehensive \
       --no_overlap \
       --bedGraph \
       --counts \
       --buffer_size 20G \
       --report \
       --cytosine_report \
       --genome_folder $ref_bis \
       -o $outdir/${sampleid}/methyl \
       $outdir/${sampleid}/align/${sampleid}.${s_species}.PE.removeGap.sorted.rmdup.sortn.bam
    date

    if [[ $s_species == "clai" ]];then
      cat $outdir/${sampleid}/methyl/C*context_${sampleid}.${s_species}.PE.removeGap.sorted.rmdup.sortn.txt|cut -f 3-|awk '$2>18 && $2<68'|awk 'BEGIN{sum=0;umethy=0}{sum+=1;if($3~"X|H|Z"){methy+=1}}END{print sum,100*methy/sum"%","*","*"}'|sed 's/ /\t/g' > $outdir/${sampleid}/stat/${sampleid}.${s_species}.Methyl.stat.txt
    else
      ## Remove Gap and filter DNA reads
      remove_Gap_sepReads PE ${s_species}

      Remove_Duplicates PE ${s_species} dna

      $bismark_methylation_extractor -p \
         --multicore $cores \
         --comprehensive \
         --no_overlap \
         --bedGraph \
         --counts \
         --buffer_size 20G \
         --report \
         --cytosine_report \
         --genome_folder $ref_bis \
         -o $outdir/${sampleid}/methyl \
         $outdir/${sampleid}/align/${sampleid}.${s_species}.PE.dna.sorted.rmdup.sortn.bam

      date

      bis_stat_pre=$outdir/${sampleid}/methyl/${sampleid}.${s_species}.PE.removeGap.sorted.rmdup.sortn_splitting_report.txt
      bis_stat=$outdir/${sampleid}/methyl/${sampleid}.${s_species}.PE.dna.sorted.rmdup.sortn_splitting_report.txt
      bis_stat_pre_stat=`grep "CpG" $bis_stat_pre |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{CpG=$1+$2; cpg_r=100*$1/CpG;print CpG,cpg_r"%"}'|sed 's/ /\t/g'`
      bis_stat_after_stat=`grep "CpG" $bis_stat |sed 's/.*:\t//g'|awk -v ORS="\t" '{print}'|awk '{CpG=$1+$2; cpg_r=100*$1/CpG;print CpG,cpg_r"%"}'|sed 's/ /\t/g'`
      echo -e "${bis_stat_pre_stat}\t${bis_stat_after_stat}" > $outdir/${sampleid}/stat/${sampleid}.${s_species}.Methyl.stat.txt
  fi

}


function merge_CH_report {
species=$1
mkdir -p $outdir/${sampleid}/methyl/sort_tmp

export TMPDIR=$outdir/${sampleid}/methyl/sort_tmp

awk -F "\t" '$6~"CH"' $outdir/${sampleid}/methyl/${sampleid}.${species}.SE.dna.sorted.rmdup.sortn.CX_report.txt > $outdir/${sampleid}/methyl/${sampleid}.${species}.SE.dna.sorted.rmdup.sortn.CH_report.txt
awk -F "\t" '$6~"CH"' $outdir/${sampleid}/methyl/${sampleid}.${species}.PE.dna.sorted.rmdup.sortn.CX_report.txt > $outdir/${sampleid}/methyl/${sampleid}.${species}.PE.dna.sorted.rmdup.sortn.CH_report.txt

module load methylpy
mkdir -p $outdir/bismark_allc

awk -v OFS='\t' '{if (($4+$5)>0) print $1,$2,$3,$7,$4,$4+$5,1}' $outdir/${sampleid}/methyl/${sampleid}.${species}.SE.dna.sorted.rmdup.sortn.CH_report.txt > $outdir/bismark_allc/allc_${sampleid}.${species}.dna.SE.tsv
awk -v OFS='\t' '{if (($4+$5)>0) print $1,$2,$3,$7,$4,$4+$5,1}' $outdir/${sampleid}/methyl/${sampleid}.${species}.SE.dna.sorted.rmdup.sortn.CpG_report.txt >> $outdir/bismark_allc/allc_${sampleid}.${species}.dna.SE.tsv

awk -v OFS='\t' '{if (($4+$5)>0) print $1,$2,$3,$7,$4,$4+$5,1}' $outdir/${sampleid}/methyl/${sampleid}.${species}.PE.dna.sorted.rmdup.sortn.CH_report.txt > $outdir/bismark_allc/allc_${sampleid}.${species}.dna.PE.tsv
awk -v OFS='\t' '{if (($4+$5)>0) print $1,$2,$3,$7,$4,$4+$5,1}' $outdir/${sampleid}/methyl/${sampleid}.${species}.PE.dna.sorted.rmdup.sortn.CpG_report.txt >> $outdir/bismark_allc/allc_${sampleid}.${species}.dna.PE.tsv


sort --parallel=8 -k1,1 -k2n,2 $outdir/bismark_allc/allc_${sampleid}.${species}.dna.SE.tsv|gzip -c  > $outdir/bismark_allc/allc_${sampleid}.${species}.dna.SE.tsv.gz
sort --parallel=8 -k1,1 -k2n,2 $outdir/bismark_allc/allc_${sampleid}.${species}.dna.PE.tsv|gzip -c  > $outdir/bismark_allc/allc_${sampleid}.${species}.dna.PE.tsv.gz

methylpy merge-allc --allc-files $outdir/bismark_allc/allc_${sampleid}.${species}.dna.SE.tsv.gz $outdir/bismark_allc/allc_${sampleid}.${species}.dna.PE.tsv.gz --output-file $outdir/bismark_allc/allc_${sampleid}.${species}.dna.tsv.gz --num-procs 8 --compress-output True
$allcools standardize-allc \
    --chrom_size_path $chrom_size \
    --remove_additional_chrom \
    --allc_path  $outdir/bismark_allc/allc_${sampleid}.${species}.dna.tsv.gz

rm $outdir/bismark_allc/allc_${sampleid}.${species}.dna.PE.* $outdir/bismark_allc/allc_${sampleid}.${species}.dna.SE.*
}


function Final_stat {
#  mkdir -p $outdir/stat
  species=$1
  Total_read_pair=`grep "Total reads processed" $outdir/${sampleid}/trim/${sampleid}[._]R1*_trimming_report.txt | sed 's/\s\s*/ /g' | cut -d " " -f 4 | sed 's/,//g'|xargs`
  if [[ $Total_read_pair == "" ]];then
  Total_read_pair=`grep "Total reads processed" $outdir/${sampleid}/trim/${sampleid}[._]R2*_trimming_report.txt | sed 's/\s\s*/ /g' | cut -d " " -f 4 | sed 's/,//g'|xargs`
  fi

  Total_reads=$(($Total_read_pair * 2 ))
  Total_bases=`echo "$Total_reads * 150"|bc`

  if [[ ! -s $outdir/${sampleid}/trim/${sampleid}.R1_val_1_fastqc/fastqc_data.txt ]];then
    unzip -n $outdir/${sampleid}/trim/${sampleid}[._]R*_val_1_fastqc.zip -d $outdir/trim
    unzip -n $outdir/${sampleid}/trim/${sampleid}[._]R*_val_2_fastqc.zip -d $outdir/trim
  fi

  Clean_reads=$(($(grep "Total Sequences" $outdir/${sampleid}/trim/${sampleid}.R1_val_1_fastqc/fastqc_data.txt |sed 's/Total Sequences\s//g'|xargs) + $(grep "Total Sequences" $outdir/${sampleid}/trim/${sampleid}.R2_val_2_fastqc/fastqc_data.txt |sed 's/Total Sequences\s//g'|xargs)))
  Clean_read_base=$(($(grep "Total written (filtered)" $outdir/${sampleid}/trim/${sampleid}[._]R1*_trimming_report.txt |sed 's/.*:\s//g'|awk '{print $1}'|sed 's/,//g'|xargs) + $(grep "Total written (filtered)" $outdir/${sampleid}/trim/${sampleid}[._]R2*_trimming_report.txt |sed 's/.*:\s//g'|awk '{print $1}'|sed 's/,//g'|xargs)))

  Filter_reads_r=`echo "scale=2;100*(1-$Clean_reads/$Total_reads)"|bc`
  Filter_base_r=`echo "scale=2;100*(1-$Clean_read_base/$Total_bases)"|bc`



  reads_aligned=`grep "reads mapped:" $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam.stat|cut -f 3|xargs`
  base_aligned=`grep "bases mapped (cigar):" $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam.stat|cut -f 3|xargs`
  Align_rate=`echo "scale=2;100*$reads_aligned/$Clean_reads"|bc|xargs`
  if [[ ! -s $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam.depth ]];then
    $samtools depth $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam > $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam.depth
  fi

  COVERED_REGION=$(wc -l $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam.depth|awk '{print $1}')
  COVERAGE=`echo "scale=4;100*$COVERED_REGION/$genome_len"|bc|xargs`
  #chrM_readsN=`$samtools index $outdir/${sampleid}/align/${sampleid}.${species}.dna.merged.PE_SE.bam | $samtools idxstats $outdir/${sampleid}/align/${sampleid}.${species}.dna.merged.PE_SE.bam|awk '{print $3,$1}'|sort|grep -v "_"|grep chr|awk 'BEGIN{sum=0;chrM=0}{if($2=="chrM"){chrM=$1};sum+=$1}END{print chrM/sum*100}'`
  chrM_readsN=`$samtools idxstats $outdir/${sampleid}/align/${sampleid}.${species}.merged.PE_SE.sorted.bam |awk '{print $3,$1}'|sort|grep -v "_"|grep chr|awk 'BEGIN{sum=0;chrM=0}{if($2=="chrM"){chrM=$1};sum+=$1}END{print chrM/sum*100}'`


  echo "chrM" $chrM_readsN
  Dup_READS=`echo "$(awk -F "\t" '{print $6}' $outdir/${sampleid}/align/${sampleid}.${species}.SE.dna.sorted.rmdup.txt | grep -A 1 "UNPAIRED_READ_DUPLICATES" | sed -n 2p) + $(awk -F "\t" '{print $7}' $outdir/${sampleid}/align/${sampleid}.${species}.PE.dna.sorted.rmdup.txt | grep -A 1 "READ_PAIR_DUPLICATES" | sed -n 2p)*2"|bc|xargs`
  Dup_rate=`echo "scale=4;100*$Dup_READS/$reads_aligned"|bc|xargs`

  Mean_Insertion_size=`grep -A 1 MEAN_INSERT_SIZE $outdir/${sampleid}/align/${sampleid}.${species}.PE.insert_size_metrics.txt|cut -f 6|tail -n 1|xargs`

  Depth_after_align=`echo "scale=6;$base_aligned / $genome_len "|bc|xargs`
  lambda_aligned=`grep "reads mapped:" $outdir/${sampleid}/align/${sampleid}.lambda.bam.stat|cut -f 3|xargs`
  puc19_aligned=`grep "reads mapped:" $outdir/${sampleid}/align/${sampleid}.fullpuc19.bam.stat|cut -f 3|xargs`
  clai_aligned=`grep "reads mapped:" $outdir/${sampleid}/align/${sampleid}.clai.bam.stat|cut -f 3|xargs`

  dna_reads_num_PE=`grep "reads mapped:" $outdir/${sampleid}/align/${sampleid}.${species}.PE.dna.sorted.bam.stat|cut -f 3|xargs`
  dna_reads_num_SE=`grep "reads mapped:" $outdir/${sampleid}/align/${sampleid}.${species}.SE.dna.sorted.bam.stat|cut -f 3|xargs`
  rna_reads_num_PE=`grep "reads mapped:" $outdir/${sampleid}/align/${sampleid}.${species}.PE.rna.sorted.bam.stat|cut -f 3|xargs`
  rna_reads_num_SE=`grep "reads mapped:" $outdir/${sampleid}/align/${sampleid}.${species}.SE.rna.sorted.bam.stat|cut -f 3|xargs`
  dna_reads_num=`echo "$dna_reads_num_PE + $dna_reads_num_SE"|bc|xargs`
  rna_reads_num=`echo "$rna_reads_num_PE + $rna_reads_num_SE"|bc|xargs`
  dna_reads_ratio=`echo "scale=4; $dna_reads_num/$reads_aligned*100"|bc|xargs`
  rna_reads_ratio=`echo "scale=4; $rna_reads_num/$reads_aligned*100"|bc|xargs`

  $samtools merge $outdir/${sampleid}/align/${sampleid}.${species}.dna.merged.PE_SE.bam $outdir/${sampleid}/align/${sampleid}.${species}.PE.dna.sorted.bam $outdir/${sampleid}/align/${sampleid}.${species}.SE.dna.sorted.bam
  $samtools index $outdir/${sampleid}/align/${sampleid}.${species}.dna.merged.PE_SE.bam
  $samtools stats $outdir/${sampleid}/align/${sampleid}.${species}.dna.merged.PE_SE.bam > $outdir/${sampleid}/align/${sampleid}.${species}.dna.merged.PE_SE.bam.stat
  $samtools depth $outdir/${sampleid}/align/${sampleid}.${species}.dna.merged.PE_SE.bam > $outdir/${sampleid}/align/${sampleid}.${species}.dna.merged.PE_SE.bam.depth

  dna_COVERED_REGION=$(wc -l $outdir/${sampleid}/align/${sampleid}.${species}.dna.merged.PE_SE.bam.depth|awk '{print $1}')
  dna_coverage=`echo "scale=4;100*$dna_COVERED_REGION/$genome_len"|bc|xargs`

  echo -e "${reads_aligned}\t" \
  "${Align_rate}%\t" \
  "${COVERAGE}%\t" \
  "${dna_reads_num}\t" \
  "${rna_reads_num}\t" \
  "${dna_reads_ratio}%\t" \
  "${rna_reads_ratio}%\t" \
  "${dna_coverage}%\t" \
  "${chrM_readsN}%\t" \
  "${Dup_rate}%\t" \
  "${Mean_Insertion_size}\t" \
  "${Depth_after_align}\t" \
  "${lambda_aligned}\t" \
  "${puc19_aligned}\t" \
  "${clai_aligned}" \
  > $outdir/${sampleid}/stat/${sampleid}.Align.stat.txt

echo -e \
  "SampleID\tTotal_reads\tClean_reads\tFilter_reads_r\tFilter_base_r\traw_depth\tclean_depth\t" \
  "Aligned_Reads\tAlign_rate\tCOVERAGE\t" \
  "dna_reads_num\trna_reads_num\tdna_reads_ratio\trna_reads_ratio\tdna_reads_coverage\t" \
  "chrM_readsN\tDup_rate\tMean_Insertion_size\tDepth_after_align\tlambda_aligned\tfullpuc19_aligned\tclai_aligned\t" \
  "dna_CpG_Cov\tdna_CGn\tdna_mCG_R\tdna_mCHG_R\tdna_mCHH_R\t" \
  "rna_CpG_Cov\trna_CGn\trna_mCG_R\trna_mCHG_R\trna_mCHH_R\t" \
  "fullpuc19_CGn_pre\tfullpuc19_mCG_R_pre\tfullpuc19_CGn_after\ttfullpuc19_mCG_R_after\t" \
  "lambda_CGn_pre\tlambda_mCG_R_pre\tlambda_CGn_after\tlambda_mCG_R_after\t" \
  "cla1_CGn_pre\tcla1_mCG_R_pre\tcla1_CGn_after\tcla1_mCG_R_after\t" \
  > $outdir/${sampleid}/stat/${sampleid}.total.stat.txt

paste $outdir/${sampleid}/stat/${sampleid}.trim.stat.txt $outdir/${sampleid}/stat/${sampleid}.Align.stat.txt $outdir/${sampleid}/stat/${sampleid}.${species}.dna.Methyl.stat.txt $outdir/${sampleid}/stat/${sampleid}.${species}.rna.Methyl.stat.txt $outdir/${sampleid}/stat/${sampleid}.fullpuc19.Methyl.stat.txt $outdir/${sampleid}/stat/${sampleid}.lambda.Methyl.stat.txt $outdir/${sampleid}/stat/${sampleid}.clai.Methyl.stat.txt >> $outdir/${sampleid}/stat/${sampleid}.total.stat.txt

#rm $outdir/${sampleid}/stat/${sampleid}*


}

function del_intermediate_file {
####Delete intermediate files
echo "Delete intermediate files"

rm $outdir/${sampleid}/align/*.depth
rm $outdir/${sampleid}/methyl/*SE.*report.*txt
rm $outdir/${sampleid}/methyl/*PE.*report.*txt
rm $outdir/${sampleid}/methyl/C*context*${species}*txt
rm $outdir/${sampleid}/methyl/C*context*lambda*txt
rm $outdir/${sampleid}/methyl/C*context*puc19*txt

}

set -x

if [[ $species == "hg38" ]] || [[ $species == "mm10" ]];then
        echo "hg38/mm10 DNA pipeline"
#:<<MULTILINECOMMENT
trim_array

trim_stat

mapping

remove_Gap_sepReads PE $species

remove_Gap_sepReads SE $species

Remove_Duplicates PE $species dna

Remove_Duplicates PE $species rna

Remove_Duplicates SE $species dna

Remove_Duplicates SE $species rna

call_methyl $species dna

call_methyl $species rna

FeatureCount_Qualimap $species dna

FeatureCount_Qualimap $species rna

spike_in fullpuc19

spike_in lambda

spike_in clai

Final_stat $species

merge_CH_report $species
#MULTILINECOMMENT
del_intermediate_file

echo "All pipeline Finished"
date

set +x
fi
