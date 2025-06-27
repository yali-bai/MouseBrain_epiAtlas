##### 01.import packages #####
library(data.table)
library(stringr)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

##### 02.set working path #####
# setwd("./")

##### 03.filter DMRs #####
result = fread(paste0(indir,"/DMR.statistic.csv"),data.table=F,header=T)
dim(result)

## calculate diff old_mean - young_mean ##
result$diff = result$old_mean - result$young_mean

## calculate segment length ##
result$start = unlist(lapply(as.character(result$chrom), function(x) strsplit(x,'_')[[1]][2]))
result$end = unlist(lapply(as.character(result$chrom), function(x) strsplit(x,'_')[[1]][3]))
result$length = as.numeric(result$end) - as.numeric(result$start)

## filter before adjusting p-value  ##
result.subset = subset(result, young_number >= 10 & old_number >= 10 & abs(diff) > 0.05 & length >= 200)
dim(result.subset)

## adjusting p-value in each subclass ##
result.subset$mannwhitneyu_p_adj = NA
for(cl in unique(result.subset$cluster)){
    result.subset$mannwhitneyu_p_adj[which(result.subset$cluster == cl)] = p.adjust(result.subset$mannwhitneyu_p[which(result.subset$cluster == cl)], method = "BH")
}
length(which(is.na(result.subset$mannwhitneyu_p_adj)))

## DMRs filter ##
result.sig = subset(result.subset, mannwhitneyu_p_adj < 0.05)
dim(result.sig)[1]
table(result.sig$cluster)

## save result ##
saveRDS(result.sig,file="../../output/03.Aging_Mouse/03-DMRs_DHMRs/DMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.rds")

## hyper DMRs ##
up_DMR = subset(result.sig, diff > 0.05 & young_number >= 10 & old_number >= 10 & mannwhitneyu_p_adj < 0.05 & length >= 200)
up_DMR_bed = data.frame(unlist(lapply(as.character(up_DMR$chrom), function(x) strsplit(x,'_')[[1]][1])),unlist(lapply(as.character(up_DMR$chrom), function(x) strsplit(x,'_')[[1]][2])),unlist(lapply(as.character(up_DMR$chrom), function(x) strsplit(x,'_')[[1]][3])))
#write.table(up_DMR_bed,file=paste0("../../output/02-aging/02-DMRs_DHMRs/DMR_upregulated_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.bed"),sep="\t",row.names=F,col.names=F,quote=F)

## hypo DMRs ##
down_DMR = subset(result.sig, diff < -0.05 & young_number >= 10 & old_number >= 10 & mannwhitneyu_p_adj < 0.05 & length >= 200)
down_DMR_bed = data.frame(unlist(lapply(as.character(down_DMR$chrom), function(x) strsplit(x,'_')[[1]][1])),unlist(lapply(as.character(down_DMR$chrom), function(x) strsplit(x,'_')[[1]][2])),unlist(lapply(as.character(down_DMR$chrom), function(x) strsplit(x,'_')[[1]][3])))
#write.table(down_DMR_bed,file=paste0("../../output/02-aging/02-DMRs_DHMRs/DMR_downregulated_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.bed"),sep="\t",row.names=F,col.names=F,quote=F)

## output all DMRs: 1 means hyper, and 0 means hypo ##
DMR_subset = result.sig[,c("cluster","chrom","diff")]
DMR_subset$hyper_hypo = 0
DMR_subset$hyper_hypo[which(DMR_subset$diff > 0)] = 1
head(DMR_subset)
head(DMR_subset[which(DMR_subset$diff < 0),])
#write.csv(DMR_subset,file="../../output/02-aging/02-DMRs_DHMRs/DMR.subclass_diff_hyper_hypo_info.csv",quote=F,row.names=F,col.names=T)

##### 04.filter DHMRs #####
result = fread(paste0(indir,"/DHMR.statistic.csv"),data.table=F,header=T)

## calculate diff old_mean - young_mean ##
result$diff = result$old_mean - result$young_mean

## calculate length ##
result$start = unlist(lapply(as.character(result$chrom), function(x) strsplit(x,'_')[[1]][2]))
result$end = unlist(lapply(as.character(result$chrom), function(x) strsplit(x,'_')[[1]][3]))
result$length = as.numeric(result$end) - as.numeric(result$start)

## filter before adjusting p-value ##
result.subset = subset(result, young_number >= 10 & old_number >= 10 & abs(diff) > 0.05 & length >= 200)

## adjusting p-value in each subclass ##
result.subset$mannwhitneyu_p_adj = NA
for(cl in unique(result.subset$cluster)){
    result.subset$mannwhitneyu_p_adj[which(result.subset$cluster == cl)] = p.adjust(result.subset$mannwhitneyu_p[which(result.subset$cluster == cl)], method = "BH")
}
length(which(is.na(result.subset$mannwhitneyu_p_adj)))

## DHMRs filter ##
result.sig = subset(result.subset, mannwhitneyu_p_adj < 0.05)
dim(result.sig)[1]
table(result.sig$cluster)

## save result ##
saveRDS(result.sig,file="../../output/03.Aging_Mouse/03-DMRs_DHMRs/DHMR_significant_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.rds")

## hyper DHMRs ##
up_DHMR = subset(result.sig, diff > 0.05 & young_number >= 10 & old_number >= 10 & mannwhitneyu_p_adj < 0.05 & length >= 200)
up_DHMR_bed = data.frame(unlist(lapply(as.character(up_DHMR$chrom), function(x) strsplit(x,'_')[[1]][1])),unlist(lapply(as.character(up_DHMR$chrom), function(x) strsplit(x,'_')[[1]][2])),unlist(lapply(as.character(up_DHMR$chrom), function(x) strsplit(x,'_')[[1]][3])))
write.table(up_DHMR_bed,file=paste0("../../output/03.Aging_Mouse/03-DMRs_DHMRs/DHMR_upregulated_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.bed"),sep="\t",row.names=F,col.names=F,quote=F)

## hypo DHMRs ##
down_DHMR = subset(result.sig, diff < -0.05 & young_number >= 10 & old_number >= 10 & mannwhitneyu_p_adj < 0.05 & length >= 200)
down_DHMR_bed = data.frame(unlist(lapply(as.character(down_DHMR$chrom), function(x) strsplit(x,'_')[[1]][1])),unlist(lapply(as.character(down_DHMR$chrom), function(x) strsplit(x,'_')[[1]][2])),unlist(lapply(as.character(down_DHMR$chrom), function(x) strsplit(x,'_')[[1]][3])))
write.table(down_DHMR_bed,file=paste0("../../output/03.Aging_Mouse/03-DMRs_DHMRs/DHMR_downregulated_in_old.diff_0.05_p_adjusted_0.05.non_na_number_10.length_200_2000.bed"),sep="\t",row.names=F,col.names=F,quote=F)

## output all DHMRs: 1 means hyper, and 0 means hypo ##
DHMR_subset = result.sig[,c("cluster","chrom","diff")]
DHMR_subset$hyper_hypo = 0
DHMR_subset$hyper_hypo[which(DHMR_subset$diff > 0)] = 1
write.csv(DHMR_subset,file="../../output/03.Aging_Mouse/03-DMRs_DHMRs/DHMR.subclass_diff_hyper_hypo_info.csv",quote=F,row.names=F,col.names=T)


