### 01-generate_DHMR_hyper_bed.ipynb: 
This script generates hyper-DHMRs for different subclasses.

### 02-run_gwas_step1.sh: 
This script is used to perform the first stage of GWAS analysis, including converting BED files from mm10 chromosome coordinates to hg19 chromosome coordinates and validating the consistency of the conversion results. 
The script first downloads the conversion chain files from UCSC, then processes each BED file through four steps: 
(1) converting mm10 coordinates to hg19 coordinates, 
(2) converting hg19 coordinates back to mm10, 
(3) checking the consistency of the converted regions (at least 50% of the regions can be mapped back to mm10), 
(4) and finally generating the final hg19 coordinates and sorting them. 

### 02-run_gwas_step2-3.sh: 
This script is used to perform the second and third stages of GWAS analysis, focusing on the annotation, heritability estimation, and statistical analysis of hyperactive DNA methylation regions (hyperDHMR) using LD Score Regression (LDSC). 
The script first annotates each BED file for different cell types, generates annotation files, and then uses the LDSC tool to estimate heritability and potential epigenetic regulatory loci.