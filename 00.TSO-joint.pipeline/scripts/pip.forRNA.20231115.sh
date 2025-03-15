#!/usr/bin/bash
## language
python=/share/home/fany/miniconda3/envs/python2/bin/python
Rscript=/share/home/fany/software/conda_software/Rscript
python3=/share/analysisdata/Methyl/workflow/software/miniconda/envs/bioinfo/bin/python

## soft Path
picard=/share/home/fany/miniconda3/bin/picard
samtools=/share/home/fany/miniconda3/bin/samtools
bedtools=/share/home/fany/miniconda3/bin/bedtools
sambamba=/share/home/fany/miniconda3/bin/sambamba
STAR=/share/analysisdata/Methyl/workflow/software/miniconda/envs/bioinfo/bin/STAR
qualimap=/share/analysisdata/Methyl/workflow/software/miniconda/envs/bioinfo/bin/qualimap
featureCounts=/share/analysisdata/Methyl/workflow/software/miniconda/envs/bioinfo/bin/featureCounts
trim_galore=/share/home/fany/miniconda3/bin/trim_galore
cutadapt=/share/home/fany/miniconda3/bin/cutadapt
## home made code
SRC_DIR=`pwd`/../
SRC_DIR=/share/analysisdata/Methyl/workflow/TSO_HT/src/Pipeline/TSO_joint_RNA_src_240809


Remove_mC_filledGap=$SRC_DIR/scripts/Remove_mC_filledGap.py
select_rna=${SRC_DIR}/scripts/mct_star_bam_filter.split_bys.py

## reference

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

#----------------------- STAR Index --------------------
hg38_STAR=/share/analysisdata/Methyl/database/human/hg38/Index/STAR_hg38_Gencode_release_35/STAR_index
mm10_STAR=/share/analysisdata/Methyl/database/mouse/mm10/STAR_index
lambda_STAR=/share/analysisdata/Methyl/database/NEB_lambda/STAR_index/  
puc19_STAR=/share/analysisdata/Methyl/database/puC19/STAR_index
clai_STAR=/share/analysisdata/Methyl/database/clai/STAR_index

#----------------------- Gtf Ref --------------------
hg38_gtf=/share/analysisdata/Methyl/database/human/hg38/Annotation/gencode.v35.annotation.gtf
mm10_gtf=/share/analysisdata/Methyl/database/mouse/mm10/Annotation/gencode.vM18.annotation.gtf



## variable

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

declare -A STAR_Ref
STAR_Ref["hg38"]=$hg38_STAR
STAR_Ref["mm10"]=$mm10_STAR
STAR_Ref["lambda"]=$lambda_STAR
STAR_Ref["fullpuc19"]=$puc19_STAR
STAR_Ref["clai"]=$clai_STAR


declare -A Ref_gtf
Ref_gtf["hg38"]=$hg38_gtf
Ref_gtf["mm10"]=$mm10_gtf

declare -A Ref_bed
Ref_bed["hg38"]=$hg38_gene_bed



m_threshold=0.5


module load java/1.8.0
## Args
temp=`getopt -o i:s:p:m:o:t: --long indir:,prefix:,outdir:,read1:,read2:,read_len:,thread:,species:,seq_type:,partition:  -- "$@"`
if [ $? != 0 ] ; then echo "terminating..." >&2 ; exit 1 ; fi
eval set -- "$temp"

while true ; do
        case "$1" in
          -i|--indir) echo "Input dir is $2"; indir=$2; shift 2;;
          -s|--sep_type) sep_type=$2; echo "Separate the reads with barcodes. Default: meta_mark" ; shift 2;;
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


cores=`expr $SLURM_CPUS_PER_TASK / 5`
if [[ $cores == 0 ]];then
	 cores=1
fi

ref_fa=${Ref_fa[$species]}
Blacklist=${Blacklist[$species]}
STAR_Ref=${STAR_Ref[$species]}
Ref_gtf=${Ref_gtf[$species]}
Ref_bed=${Ref_bed[$species]}


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
mkdir -p $outdir/trim
t1=`ls $outdir/trim/${sampleid}[._]R1_val_1.fq.gz`
t2=`ls $outdir/trim/${sampleid}[._]R2_val_2.fq.gz`
echo "Reads need to be trimed: $Read1, $Read2"

if [[ -f $t1 ]] && [[ -f $t2 ]];then
        echo "Trim already finished"
        return
else
       rm -fr $outdir/trim/${sampleid}*
fi

$trim_galore --fastqc --fastqc_args "-t $SLURM_CPUS_PER_TASK" --paired \
      --phred33 \
      --length 20 \
      --retain_unpaired \
      --path_to_cutadapt $cutadapt \
      --output_dir $outdir/trim \
      $Read1 $Read2

if [[ ! -f $outdir/trim/${sampleid}.R2_val_2_fastqc.zip ]];then
fastqc -t $SLURM_CPUS_PER_TASK $outdir/trim/${sampleid}[._]R1_val_1.fq.gz
fastqc -t $SLURM_CPUS_PER_TASK $outdir/trim/${sampleid}[._]R2_val_2.fq.gz
fi

}

function trim_stat {
  mkdir -p $outdir/trim_stat
  Total_read_pair=`grep "Total reads processed" $outdir/trim/${sampleid}[._]R1*_trimming_report.txt | sed 's/\s\s*/ /g' | cut -d " " -f 4 | sed 's/,//g'|xargs`
  Total_reads=$(($Total_read_pair * 2 ))
  Total_bases=`echo "$Total_reads * 150"|bc`
  unzip -n $outdir/trim/${sampleid}[._]R*_val_1_fastqc.zip -d $outdir/trim
  unzip -n $outdir/trim/${sampleid}[._]R*_val_2_fastqc.zip -d $outdir/trim

  Clean_reads=$(($(grep "Total Sequences" $outdir/trim/${sampleid}.R1_val_1_fastqc/fastqc_data.txt |sed 's/Total Sequences\s//g'|xargs) + $(grep "Total Sequences" $outdir/trim/${sampleid}.R2_val_2_fastqc/fastqc_data.txt |sed 's/Total Sequences\s//g'|xargs)))
  Clean_read_base=$(($(grep "Total written (filtered)" $outdir/trim/${sampleid}[._]R1*_trimming_report.txt |sed 's/.*:\s//g'|awk '{print $1}'|sed 's/,//g'|xargs) + $(grep "Total written (filtered)" $outdir/trim/${sampleid}[._]R2*_trimming_report.txt |sed 's/.*:\s//g'|awk '{print $1}'|sed 's/,//g'|xargs)))

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
    "${clean_depth}"  > $outdir/trim_stat/${sampleid}.trim.stat.txt
}

function mapping {
trim_Read1=$outdir/trim/${sampleid}.R1_val_1.fq.gz
trim_Read2=$outdir/trim/${sampleid}.R2_val_2.fq.gz
mkdir -p $outdir/align


$STAR --runThreadN $SLURM_CPUS_PER_TASK \
  --genomeDir $STAR_Ref \
  --readFilesIn $trim_Read1 $trim_Read2 \
  --readFilesCommand gunzip -c \
  --outSAMtype BAM Unsorted \
  --outFileNamePrefix $outdir/align/${sampleid}.${species}. \
  --alignEndsType Local \
  --genomeLoad NoSharedMemory \
  --outSAMunmapped None \
  --outSAMattributes NH HI AS NM MD \
  --sjdbOverhang 100 \
  --outFilterType BySJout \
  --outFilterMultimapNmax 1 \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverLmax 0.04 \
  --alignIntronMin 20 \
  --alignIntronMax 1000000 \
  --alignMatesGapMax 1000000


outbam=`ls $outdir/align/${sampleid}.${species}.Aligned.out.bam|xargs`
###
if [[ $species == "hg38" ]] || [[ $species == "mm10" ]]
then
  $samtools view -bh -q 1 -F 4 $outbam |$samtools sort -@ $SLURM_CPUS_PER_TASK |$bedtools intersect -a - -b $Blacklist -v > $outdir/align/${sampleid}.${species}.bam
else
  $samtools view -bh -q 1 -F 4 $outbam |$samtools sort -@ $SLURM_CPUS_PER_TASK > $outdir/align/${sampleid}.${species}.bam
fi

$samtools index $outdir/align/${sampleid}.${species}.bam
$samtools stats $outdir/align/${sampleid}.${species}.bam > $outdir/align/${sampleid}.${species}.bam.stat
$samtools depth $outdir/align/${sampleid}.${species}.bam > $outdir/align/${sampleid}.${species}.depth

}

function select_rna {

  $python3 $select_rna --sampleid ${sampleid} --species ${species} --sequence PE --input_bam $outdir/align/${sampleid}.${species}.bam --output_bam $outdir/align/${sampleid}.${species}.rna.bam

  $samtools sort -@ $SLURM_CPUS_PER_TASK $outdir/align/${sampleid}.${species}.rna.bam > $outdir/align/${sampleid}.${species}.rna.sorted.bam

  $samtools index $outdir/align/${sampleid}.${species}.rna.sorted.bam
  $samtools stats $outdir/align/${sampleid}.${species}.rna.sorted.bam > $outdir/align/${sampleid}.${species}.rna.sorted.bam.stat
  $samtools depth $outdir/align/${sampleid}.${species}.rna.sorted.bam > $outdir/align/${sampleid}.${species}.rna.sorted.depth

  echo "Calculate the fragment length"
  $picard CollectInsertSizeMetrics \
          I=$outdir/align/${sampleid}.${species}.rna.sorted.bam \
          O=$outdir/align/${sampleid}.${species}.rna.insert_size_metrics.txt \
          H=$outdir/align/${sampleid}.${species}.rna.insert_size_histogram.pdf \
          M=0.5

}

function feature_count.stranded {

mkdir -p $outdir/feature_count.stranded/

$featureCounts -p --countReadPairs \
-T $SLURM_CPUS_PER_TASK \
-a $Ref_gtf \
-t exon \
-g gene_id \
-R BAM \
-s 1 \
-o $outdir/feature_count.stranded/${sampleid}.${species}.exon_counts.txt \
$outdir/align/${sampleid}.${species}.rna.sorted.bam

$featureCounts -p --countReadPairs \
-T $SLURM_CPUS_PER_TASK \
-a $Ref_gtf \
-t transcript \
-g gene_id \
-R BAM \
-s 1 \
-o $outdir/feature_count.stranded/${sampleid}.${species}.transcript_counts.txt \
$outdir/align/${sampleid}.${species}.rna.sorted.bam

$featureCounts -p --countReadPairs \
-T $SLURM_CPUS_PER_TASK \
-a $Ref_gtf \
-t gene \
-g gene_id \
-R BAM \
-s 1 \
-o $outdir/feature_count.stranded/${sampleid}.${species}.gene_counts.txt \
$outdir/align/${sampleid}.${species}.rna.sorted.bam


grep -v "Geneid" $outdir/feature_count.stranded/${sampleid}.${species}.gene_counts.txt |awk '$7>0'|awk '{print $6}' > $outdir/feature_count.stranded/${sampleid}.${species}.mapped_gene_length.stranded.txt

}

function Qualimap_for_QC_stranded {
  mkdir -p $outdir/Qualimap_for_QC_stranded/${sampleid}/
  read1=Output_s1/sep_cell/${sampleid}.R1.fastq.gz
  sample_size=$(stat -c "%s" $read1)
  sample_size=`echo "($sample_size * 2 + 50000000)/100000+40000"|bc|xargs`

  ${qualimap} rnaseq -bam $outdir/align/${sampleid}.${species}.rna.sorted.bam \
        --java-mem-size=${sample_size}"M" \
        -gtf $Ref_gtf \
        -pe \
        -p strand-specific-forward \
        -outdir $outdir/Qualimap_for_QC_stranded/${sampleid}/ \
        -outfile ${sampleid}.${species}.rnaseq_report.pdf \
        -oc ${sampleid}.${species}.count.txt
}

function Final_stat {
  mkdir -p $outdir/stat

  Total_read_pair=`grep "Total reads processed" $outdir/trim/${sampleid}[._]R1*_trimming_report.txt | sed 's/\s\s*/ /g' | cut -d " " -f 4 | sed 's/,//g'|xargs`
  if [[ $Total_read_pair == "" ]];then
    Total_read_pair=`grep "Total reads processed" $outdir/trim/${sampleid}[._]R2*_trimming_report.txt | sed 's/\s\s*/ /g' | cut -d " " -f 4 | sed 's/,//g'|xargs`
  fi

  Total_reads=$(($Total_read_pair * 2 ))
  Total_bases=`echo "$Total_reads * 150"|bc`

  if [[ ! -s $outdir/trim/${sampleid}.R1_val_1_fastqc/fastqc_data.txt ]];then
    unzip -n $outdir/trim/${sampleid}[._]R*_val_1_fastqc.zip -d $outdir/trim
    unzip -n $outdir/trim/${sampleid}[._]R*_val_2_fastqc.zip -d $outdir/trim
  fi

  Clean_reads=$(($(grep "Total Sequences" $outdir/trim/${sampleid}.R1_val_1_fastqc/fastqc_data.txt |sed 's/Total Sequences\s//g'|xargs) + $(grep "Total Sequences" $outdir/trim/${sampleid}.R2_val_2_fastqc/fastqc_data.txt |sed 's/Total Sequences\s//g'|xargs)))
  Clean_read_base=$(($(grep "Total written (filtered)" $outdir/trim/${sampleid}[._]R1*_trimming_report.txt |sed 's/.*:\s//g'|awk '{print $1}'|sed 's/,//g'|xargs) + $(grep "Total written (filtered)" $outdir/trim/${sampleid}[._]R2*_trimming_report.txt |sed 's/.*:\s//g'|awk '{print $1}'|sed 's/,//g'|xargs)))

  Filter_reads_r=`echo "scale=2;100*(1-$Clean_reads/$Total_reads)"|bc`
  Filter_base_r=`echo "scale=2;100*(1-$Clean_read_base/$Total_bases)"|bc`

  reads_aligned=`grep "reads mapped:" $outdir/${sampleid}/align/${sampleid}.${species}.bam.stat|cut -f 3|xargs`
  base_aligned=`grep "bases mapped (cigar):" $outdir/${sampleid}/align/${sampleid}.${species}.bam.stat|cut -f 3|xargs`
  #Align_rate=`echo "scale=2;100*$reads_aligned/$Clean_reads"|bc|xargs`

  ##reads_aligned=`grep "reads mapped:" $outdir/align/${sampleid}.${species}.bam.stat|cut -f 3|xargs`
  ##base_aligned=`grep "bases mapped (cigar):" $outdir/align/${sampleid}.${species}.bam.stat|cut -f 3|xargs`
  ##Align_rate=`echo "scale=2;100*$reads_aligned/$Clean_reads"|bc|xargs`

  reads_aligned=`grep -E "Number.*loci|Uniquely mapped reads number" $outdir/align/${sampleid}.${species}.Log.final.out|cut -f2|awk '{align += $1} END {print 2*align}'`
  Align_rate=`grep -E "%.*loci|Uniquely mapped reads %" $outdir/align/${sampleid}.${species}.Log.final.out|cut -f2|sed 's/%//g'|awk '{align += $1} END {print align}'`
  unique_reads_aligned=`grep "Uniquely mapped reads number" $outdir/align/${sampleid}.${species}.Log.final.out|cut -f2|awk '{print 2*$1}'`
  unique_Align_rate=`grep "Uniquely mapped reads %" $outdir/align/${sampleid}.${species}.Log.final.out|cut -f2 |sed 's/%//g' `

  ##after remove blacklist and filter mCH
  RNA_reads=`$samtools view -@ $SLURM_CPUS_PER_TASK $outdir/align/${sampleid}.${species}.rna.sorted.bam|wc -l`
  RNA_reads_ratio=`echo "scale=2;100*${RNA_reads}/${reads_aligned}"|bc`
  rna_COVERED_REGION=$(wc -l $outdir/align/${sampleid}.${species}.rna.sorted.depth|awk '{print $1}')
  rna_coverage=`echo "scale=4;100*$rna_COVERED_REGION/$genome_len"|bc|xargs`
  Mean_Insertion_size=`grep -A 1 MEAN_INSERT_SIZE $outdir/align/${sampleid}.${species}.rna.insert_size_metrics.txt|cut -f 6|tail -n 1|xargs`
  #Mean_Insertion_size=`grep -A 1 MEAN_INSERT_SIZE $outdir/align/${sampleid}.${species}.rna.insert_size_metrics.txt|cut -f 5|tail -n 1|xargs`

  FeatureCounts_Exon=$(grep Assigned $outdir/feature_count.stranded/${sampleid}.${species}.exon_counts.txt.summary|cut -f2|awk '{print $1*2}')
  FeatureCounts_Gene=$(grep Assigned $outdir/feature_count.stranded/${sampleid}.${species}.gene_counts.txt.summary|cut -f2|awk '{print $1*2}')
  FeatureCounts_Intron=`echo "$FeatureCounts_Gene - $FeatureCounts_Exon"|bc`
  FeatureCounts_IntergenicRegion=$(grep "Unassigned_NoFeatures" $outdir/feature_count.stranded/${sampleid}.${species}.gene_counts.txt.summary|cut -f2|awk '{print $1*2}')

  FeatureCounts_Total=`echo "$FeatureCounts_Gene + $FeatureCounts_IntergenicRegion"|bc`
  FeatureCounts_Exonic=`echo "scale=2;100*$FeatureCounts_Exon/$FeatureCounts_Total"|bc`
  FeatureCounts_Intronic=`echo "scale=2;100*$FeatureCounts_Intron/$FeatureCounts_Total"|bc`
  FeatureCounts_Intergenic=`echo "scale=2;100*$FeatureCounts_IntergenicRegion/$FeatureCounts_Total"|bc`

  Exon_gene_number=`grep -v Geneid $outdir/feature_count.stranded/${sampleid}.${species}.exon_counts.txt|awk '$7 > 0'|wc -l`
  Gene_gene_number=`grep -v Geneid $outdir/feature_count.stranded/${sampleid}.${species}.gene_counts.txt |awk '$7 > 0'|wc -l`

  Qualimap_Exon=`grep exonic $outdir/Qualimap_for_QC_stranded/${sampleid}/rnaseq_qc_results.txt|cut -d '(' -f1|sed 's/.*=[[:space:]]*//g'|sed 's/,//g'`
  Qualimap_Exonic=`grep exonic $outdir/Qualimap_for_QC_stranded/${sampleid}/rnaseq_qc_results.txt|cut -d '(' -f2|sed 's/%.*//g'`
  Qualimap_Intron=`grep intronic $outdir/Qualimap_for_QC_stranded/${sampleid}/rnaseq_qc_results.txt|cut -d '(' -f1|sed 's/.*=[[:space:]]*//g'|sed 's/,//g'`
  Qualimap_Intronic=`grep intronic $outdir/Qualimap_for_QC_stranded/${sampleid}/rnaseq_qc_results.txt|cut -d '(' -f2|sed 's/%.*//g'`
  Qualimap_IntergenicRegion=`grep intergenic $outdir/Qualimap_for_QC_stranded/${sampleid}/rnaseq_qc_results.txt|cut -d '(' -f1|sed 's/.*=[[:space:]]*//g'|sed 's/,//g'`
  Qualimap_Intergenic=`grep intergenic $outdir/Qualimap_for_QC_stranded/${sampleid}/rnaseq_qc_results.txt|cut -d '(' -f2|sed 's/%.*//g'`
  Qualimap_overlap_Exon=`grep overlapping $outdir/Qualimap_for_QC_stranded/${sampleid}/rnaseq_qc_results.txt|cut -d '(' -f1|sed 's/.*=[[:space:]]*//g'|sed 's/,//g'`
  Qualimap_overlap_Exonic=`grep overlapping $outdir/Qualimap_for_QC_stranded/${sampleid}/rnaseq_qc_results.txt|cut -d '(' -f2|sed 's/%.*//g'`

  ###5_bias=`grep "5' bias =" $outdir/Qualimap_for_QC_stranded/${sampleid}/rnaseq_qc_results.txt |sed 's/.*=[[:space:]]*//g'`
  ###3_bias=`grep " 3' bias =" $outdir/Qualimap_for_QC_stranded/${sampleid}/rnaseq_qc_results.txt |sed 's/.*=[[:space:]]*//g'`
  ###5_3_bias=`grep "5'-3' bias =" $outdir/Qualimap_for_QC_stranded/${sampleid}/rnaseq_qc_results.txt |sed 's/.*=[[:space:]]*//g'`

  bias=`grep ".*bias =" $outdir/Qualimap_for_QC_stranded/${sampleid}/rnaseq_qc_results.txt|sed 's/.*=[[:space:]]*//g'|awk BEGIN{RS=EOF}'{gsub(/\n/,"\t");print $0}'`

  if [[ ! -f ${outdir}/stat/STAR_align_stat.${species}.forRNA.txt ]];then
          echo -e "SampleID\tSpecies\tTotal_reads\tClean_Reads\tReads_Clean_rate\tBase_Clean_rate\treads_aligned\tAlign_rate\t" \
            "unique_reads_aligned\tunique_Align_rate\t" \
            "RNA_reads\tRNA_reads_ratio\tCOVERAGE\tMean_Insertion_size\t" \
            "FeatureCounts_Exon\tFeatureCounts_Exonic\tFeatureCounts_Intron\tFeatureCounts_Intronic\tFeatureCounts_IntergenicRegion\tFeatureCounts_Intergenic\t" \
            "Exon_gene_number\tGene_gene_number\t" \
            "Qualimap_Exon\tQualimap_Exonic\tQualimap_Intron\tQualimap_Intronic\tQualimap_IntergenicRegion\tQualimap_Intergenic\t" \
            "Qualimap_overlap_Exon\tQualimap_overlap_Exonic\t" \
            "5_bias\t3_bias\t5_3_bias\t" \
                  > ${outdir}/stat/STAR_align_stat.${species}.forRNA.txt
  fi

  echo -e "$sampleid\t${species}\t" \
    "${Total_reads}\t" \
    "${Clean_reads}\t" \
    "${Filter_reads_r}%\t" \
    "${Filter_base_r}%\t" \
    "${reads_aligned}\t" \
    "${Align_rate}%\t" \
    "${unique_reads_aligned}\t" \
    "${unique_Align_rate}%\t" \
    "${RNA_reads}\t" \
    "${RNA_reads_ratio}%\t" \
    "${rna_coverage}%\t" \
    "${Mean_Insertion_size}\t" \
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
    >> $outdir/stat/STAR_align_stat.${species}.forRNA.txt

}



set -x
if [[ $species == "hg38" ]] || [[ $species == "mm10" ]];then
        echo "hg38/mm10 DNA pipeline"
trim_array

trim_stat

mapping

select_rna

feature_count.stranded

Qualimap_for_QC_stranded

Final_stat

fi
echo "All pipeline Finished"
date

set +x
