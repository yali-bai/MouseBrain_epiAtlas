library(getopt)
library(stringr)
library(Seurat)

arguments = matrix(c(
  'help', 'h', 0, "logical", "",
  'input', 'i', 1, "character", "",
  'in_dir', 'd', 1, "character", "",
  'metainfo','m', 1, "character", ""
), byrow=TRUE, ncol=5)
args = getopt(arguments)

# if help was asked for print a friendly message
# and exit with a non-zero error code
#if ( !is.null(args$help) ) {
#  cat(getopt(arguments, usage=TRUE))
#  q(status=1)
#}

if (!is.null(args$help) || is.null(args$input) || is.null(args$in_dir) || is.null(args$metainfo)) {
  cat(paste(getopt(arguments, usage = T), "\n"))
  q()
}


stat.df = read.csv(paste0(args$in_dir,"/",args$input),header=T)

metainfo = readRDS(args$metainfo)
metainfo@meta.data$sampleid = unlist(lapply(rownames(metainfo@meta.data), function(x) strsplit(x,"@@_")[[1]][2]))

stat.df$QC_after_integration = 0
stat.df$QC_after_integration[match(intersect(metainfo@meta.data$sampleid,stat.df$SampleID),stat.df$SampleID)] = 1

write.csv(stat.df,file = paste0(args$in_dir,"/","TSO-joint.RNA_QC_stat.young.csv"),quote=F,row.names=F,col.names=T,sep="\t")


