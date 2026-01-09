args <- commandArgs(trailingOnly = TRUE)
options(stringsAsFactors=FALSE)
print(args)

group <- args[1] # CG
element <- args[2] # chrom100k

library(Seurat)
library(readr)

true5mC <- function(group,element){ ## "CG","chrom100k"
mC <- read.csv(sprintf("./01_3C/output/03-5mC_%s_%s_subclass_mean_dat_final.csv",group,element),row.names=1,check.names = FALSE) 
hmC <- read.csv(sprintf("./01_3C/output/03-5hmC_%s_%s_subclass_mean_dat_final.csv",group,element),row.names=1,check.names = FALSE)
true5mC <- mC-hmC
true5mC[true5mC < 0] <- 0
write.csv(true5mC,sprintf("./01_3C/output/04-true5mC_%s_%s_subclass_mean_dat_final.csv",group,element),quote=F)
}

true5mC(group,element)




