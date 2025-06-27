#########    All "our" in the following code refers to Joint Cabernet.
library(DESeq2)
library(dplyr)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""


metacell_obj<-readRDS(paste0(indir,"/metacells_obj_k15_mincell40_maxsh3_corrected.rds"))
count<-metacell_obj@assays$RNA$count
filtered_genes<-read.csv("../../output/03.Aging_Mouse/07-aged_DEG/count_cpm_filtered_gene.csv")
count<-count[filtered_genes$x,]
info<-data.frame(row.names=colnames(metacell_obj),class=metacell_obj$lt_twice_subclass,age=metacell_obj$age)
subclass<-unique(info$class)
total_result<-data.frame()
statistics<-data.frame()
filter_fc<-as.numeric(c(0.5,1,2))
for(cl in subclass){
    cat(cl,"\n")
    countData<-count[,rownames(info)[info$class==cl]]
    condition <- factor(info$age[info$class==cl])
    colData <- data.frame(row.names=colnames(countData), condition)
    dds <- DESeqDataSetFromMatrix(countData = round(countData), colData = colData, design = ~ condition)
    dds1 <- DESeq(dds) 
    res <- results(dds1, contrast=c("condition", "old", "young"))%>%as.data.frame()
    res<-res[!is.na(res$log2FoldChange),]
    res$class<-cl
    total_result<-rbind(total_result,res)

    #statistics
    for(fc in filter_fc){
        sig_res<-res[(abs(as.numeric(res$log2FoldChange))>as.numeric(fc)&as.numeric(res$padj)<0.05),]
        sig_res<-sig_res[complete.cases(sig_res),]
        up_num<-sum(sig_res$log2FoldChange>0)
        down_num<-sum(sig_res$log2FoldChange<0)
        statistics_dt<-data.frame(level="subclass",class=cl,filter_condition=paste0("abs(log2FC)>",fc,"&padj<0.05"),
            sig_change=nrow(sig_res),sig_up=up_num,sig_down=down_num)
        statistics<-rbind(statistics,statistics_dt)
    }   

}

write.csv(total_result,"../../output/03.Aging_Mouse/07-aged_DEG/Deseq_total_result_k15_corrected.csv")
write.csv(statistics,"../../output/03.Aging_Mouse/07-aged_DEG/Deseq_statistics_k15_corrected.csv",row.names=F)

filter_result<-total_result[abs(as.numeric(total_result$log2FoldChange))>1&as.numeric(total_result$padj)<0.05,]
filter_result<-filter_result[complete.cases(filter_result),]
filter_result$class<-factor(filter_result$class,level=subclass_order)
filter_result<-filter_result%>%arrange(class)
write.csv(filter_result,paste0(outdir,"/Deseq_total_result_k15_filter_genes_corrected.csv"))



