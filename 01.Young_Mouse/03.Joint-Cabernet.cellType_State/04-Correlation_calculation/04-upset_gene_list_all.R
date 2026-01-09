library(dplyr)
library(UpSetR)
library(RColorBrewer)
library(ggplot2)
library(data.table)

setwd("./all_cell_correlation")
result<-data.frame()
for(dt in c("5hmC","5mC","true_5mC")){
  df<-fread(paste0("all_cells_",dt,"_CG_genebody_gene_correlation_results.csv"),data.table = F)
  result<-rbind(result,df)
}
colnames(result)[1]<-"gene_id"
data<-result[result$`Adjusted P-value`<0.05,]
data$corr_direction<-ifelse(data$Correlation>0,"Positive correlation","Negative correlation")
mc='CG'
vardim='genebody'
genetypes<-c("all_gene","abs(correlation)>0.05")
select_gene<-data.frame()
for(gt in genetypes){
                print(gt)
            if (gt=="all_gene"){
                aa <- subset(data,mc_type==mc&var_dim==vardim) 
            }else{
                aa <- subset(data,mc_type==mc&var_dim==vardim&abs(Correlation)>0.05)  
            }
            set1 <- subset(aa,datatype == "true_5mC"&corr_direction=="Positive correlation")$gene_id  #mc
            set2 <- subset(aa,datatype == "true_5mC"&corr_direction=="Negative correlation")$gene_id  
            set3 <- subset(aa,datatype == "5hmC"&corr_direction=="Positive correlation")$gene_id  #hmc
            set4<- subset(aa,datatype == "5hmC"&corr_direction=="Negative correlation")$gene_id
            set5 <- subset(aa,datatype == "5mC"&corr_direction=="Positive correlation")$gene_id  #mc+hmc
            set6 <- subset(aa,datatype == "5mC"&corr_direction=="Negative correlation")$gene_id
            name<-c(paste0("[5m",mc," +]"),paste0("[5m",mc," -]"),paste0("[5hm",mc," +]"),paste0("[5hm",mc," -]"),paste0("[(5m",mc,"+5hm",mc,") +]"),paste0("[(5m",mc,"+5hm",mc,") -]"))
            no_filter_aa<-subset(result,mc_type==mc&var_dim==vardim)
            #独有基因
            uniqueGeneset1 <- setdiff(set1, Reduce(union, list(set2, set3, set4, set5, set6)))
            uniqueGeneset2 <- setdiff(set2, Reduce(union, list(set1, set3, set4, set5, set6)))
            uniqueGeneset3 <- setdiff(set3, Reduce(union, list(set2, set1, set4, set5, set6)))
            uniqueGeneset4 <- setdiff(set4, Reduce(union, list(set2, set3, set1, set5, set6)))
            uniqueGeneset5 <- setdiff(set5, Reduce(union, list(set2, set3, set4, set1, set6)))
            uniqueGeneset6 <- setdiff(set6, Reduce(union, list(set2, set3, set4, set5, set1)))
            uniqueGenesets<-c(uniqueGeneset1,uniqueGeneset2,uniqueGeneset3,uniqueGeneset4,uniqueGeneset5,uniqueGeneset6)
            df<-data.frame(gene_id=uniqueGenesets,
                            group=c(rep(name[1],length(uniqueGeneset1)),rep(name[2],length(uniqueGeneset2)),rep(name[3],length(uniqueGeneset3)),rep(name[4],length(uniqueGeneset4)),rep(name[5],length(uniqueGeneset5)),rep(name[6],length(uniqueGeneset6))),
                            mc_type=rep(mc,length(uniqueGenesets)),
                            var_dim=rep(vardim,length(uniqueGenesets)),
                            gene_type=rep(gt,length(uniqueGenesets)),
                            Correlation_5mC=no_filter_aa$Correlation[no_filter_aa$datatype=="true_5mC"][match(uniqueGenesets,no_filter_aa$gene_id[no_filter_aa$datatype=="true_5mC"])],
                            new.P.value_5mC=no_filter_aa$`P-value`[no_filter_aa$datatype=="true_5mC"][match(uniqueGenesets,no_filter_aa$gene_id[no_filter_aa$datatype=="true_5mC"])],
                            new.P.adjust_5mC=no_filter_aa$`Adjusted P-value`[no_filter_aa$datatype=="true_5mC"][match(uniqueGenesets,no_filter_aa$gene_id[no_filter_aa$datatype=="true_5mC"])],
                            Correlation_5hmC=no_filter_aa$Correlation[no_filter_aa$datatype=="5hmC"][match(uniqueGenesets,no_filter_aa$gene_id[no_filter_aa$datatype=="5hmC"])],
                            new.P.value_5hmC=no_filter_aa$`P-value`[no_filter_aa$datatype=="5hmC"][match(uniqueGenesets,no_filter_aa$gene_id[no_filter_aa$datatype=="5hmC"])],
                            new.P.adjust_5hmC=no_filter_aa$`Adjusted P-value`[no_filter_aa$datatype=="5hmC"][match(uniqueGenesets,no_filter_aa$gene_id[no_filter_aa$datatype=="5hmC"])],
                            Correlation_mC_hmC=no_filter_aa$Correlation[no_filter_aa$datatype=="5mC"][match(uniqueGenesets,no_filter_aa$gene_id[no_filter_aa$datatype=="5mC"])],
                            new.P.value_mC_hmC=no_filter_aa$`P-value`[no_filter_aa$datatype=="5mC"][match(uniqueGenesets,no_filter_aa$gene_id[no_filter_aa$datatype=="5mC"])],
                            new.P.adjust_mC_hmC=no_filter_aa$`Adjusted P-value`[no_filter_aa$datatype=="5mC"][match(uniqueGenesets,no_filter_aa$gene_id[no_filter_aa$datatype=="5mC"])], 
                            stringsAsFactors = FALSE  
                            ) 
            df<-df[!duplicated(df),]
            select_gene<-rbind(select_gene,df)

            sets<-list(set1 = set1,set2 = set2,set3 = set3,set4 = set4,set5 = set5,set6 = set6)
            names(sets)<-name
 
            combos <- combn(names(sets), 3, simplify = FALSE)    
            overlapGenes <- lapply(combos, function(combo) {Reduce(intersect, sets[combo])}) 
            names(overlapGenes) <- sapply(combos, paste, collapse = ":")   
            result_list <- list()    
            for (na in names(overlapGenes)) {  
            if (!is.null(overlapGenes[[na]]) && length(overlapGenes[[na]]) > 0) {  
                temp_df <- data.frame(  
                    gene_id = overlapGenes[[na]],
                group = rep(na, length(overlapGenes[[na]])),  
                mc_type=mc,
                var_dim=vardim,
                gene_type=gt,
                Correlation_5mC=no_filter_aa$Correlation[no_filter_aa$datatype=="true_5mC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="true_5mC"])],
                new.P.value_5mC=no_filter_aa$`P-value`[no_filter_aa$datatype=="true_5mC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="true_5mC"])],
                new.P.adjust_5mC=no_filter_aa$`Adjusted P-value`[no_filter_aa$datatype=="true_5mC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="true_5mC"])],
                Correlation_5hmC=no_filter_aa$Correlation[no_filter_aa$datatype=="5hmC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="5hmC"])],
                new.P.value_5hmC=no_filter_aa$`P-value`[no_filter_aa$datatype=="5hmC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="5hmC"])],
                new.P.adjust_5hmC=no_filter_aa$`Adjusted P-value`[no_filter_aa$datatype=="5hmC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="5hmC"])],
                Correlation_mC_hmC=no_filter_aa$Correlation[no_filter_aa$datatype=="5mC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="5mC"])],
                new.P.value_mC_hmC=no_filter_aa$`P-value`[no_filter_aa$datatype=="5mC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="5mC"])],
                new.P.adjust_mC_hmC=no_filter_aa$`Adjusted P-value`[no_filter_aa$datatype=="5mC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="5mC"])],               
                stringsAsFactors = FALSE  
                )  
                result_list[[na]] <- temp_df    
            }  }
            result_df <- do.call(rbind, result_list)%>%as.data.frame() 
            if(nrow(result_df)!=0){
                result_df<-result_df[!duplicated(result_df),] 
                select_gene<-rbind(select_gene,result_df)
            }
            # 
            bb<-aa[!aa$gene_id%in%result_df$gene_id,]
            set1 <- subset(bb,datatype == "true_5mC"&corr_direction=="Positive correlation")$gene_id  #mc
            set2 <- subset(bb,datatype == "true_5mC"&corr_direction=="Negative correlation")$gene_id  
            set3 <- subset(bb,datatype == "5hmC"&corr_direction=="Positive correlation")$gene_id  #hmc
            set4<- subset(bb,datatype == "5hmC"&corr_direction=="Negative correlation")$gene_id
            set5 <- subset(bb,datatype == "5mC"&corr_direction=="Positive correlation")$gene_id  #mc+hmc
            set6 <- subset(bb,datatype == "5mC"&corr_direction=="Negative correlation")$gene_id
            sets<-list(set1 = set1,set2 = set2,set3 = set3,set4 = set4,set5 = set5,set6 = set6)
            names(sets)<-name
            combos <- combn(names(sets), 2, simplify = FALSE)  #  
            #  
            overlapGenes <- lapply(combos, function(combo) {  
            Reduce(intersect, sets[combo])  
            })  
            names(overlapGenes) <- sapply(combos, paste, collapse = ":")  # 
            result_list <- list()   # 
            # 
            for (na in names(overlapGenes)) {  
            # 
            if (!is.null(overlapGenes[[na]]) && length(overlapGenes[[na]]) > 0) {  
                # 
                    temp_df <- data.frame(  
                gene_id = overlapGenes[[na]],
                group = rep(na, length(overlapGenes[[na]])),  
                mc_type=mc,
                var_dim=vardim,
                gene_type=gt,
                Correlation_5mC=no_filter_aa$Correlation[no_filter_aa$datatype=="true_5mC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="true_5mC"])],
                new.P.value_5mC=no_filter_aa$`P-value`[no_filter_aa$datatype=="true_5mC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="true_5mC"])],
                new.P.adjust_5mC=no_filter_aa$`Adjusted P-value`[no_filter_aa$datatype=="true_5mC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="true_5mC"])],
                Correlation_5hmC=no_filter_aa$Correlation[no_filter_aa$datatype=="5hmC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="5hmC"])],
                new.P.value_5hmC=no_filter_aa$`P-value`[no_filter_aa$datatype=="5hmC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="5hmC"])],
                new.P.adjust_5hmC=no_filter_aa$`Adjusted P-value`[no_filter_aa$datatype=="5hmC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="5hmC"])],
                Correlation_mC_hmC=no_filter_aa$Correlation[no_filter_aa$datatype=="5mC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="5mC"])],
                new.P.value_mC_hmC=no_filter_aa$`P-value`[no_filter_aa$datatype=="5mC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="5mC"])],
                new.P.adjust_mC_hmC=no_filter_aa$`Adjusted P-value`[no_filter_aa$datatype=="5mC"][match(overlapGenes[[na]],no_filter_aa$gene_id[no_filter_aa$datatype=="5mC"])],
                stringsAsFactors = FALSE  # 
                    )
                result_list[[na]] <- temp_df  # 
            }  }
            result_df <- do.call(rbind, result_list)%>%as.data.frame()  # 
            if(nrow(result_df)!=0){  
            result_df<-result_df[!duplicated(result_df),] 
            select_gene<-rbind(select_gene,result_df)}
}         

gene_metainfo<-read.csv("../../../03.data/01.ref/gene_CpG_number_metainfo.csv")
select_gene<-merge(select_gene,gene_metainfo[,c("gene_id","gene_name","Cpg_number","Gene_length")],by="gene_id",all.x=TRUE)
all_gene<-select_gene[select_gene$gene_type=="all_gene",]
filter_gene<-select_gene[select_gene$gene_type!="all_gene",]
print(paste0("dim(all_gene)=",dim(all_gene)))
print(paste0("dim(filter_gene)=",dim(filter_gene)))

write.csv(all_gene,"all_gene_upset_group_list.csv",row.names=F)
write.csv(filter_gene,"filter_gene_upset_group_list.csv",row.names=F)



######  extract group1/2 gene list
# all_gene<-read.csv("all_gene_upset_group_list.csv")
group1<-all_gene$gene_name[all_gene$group=="[5mCG -]:[5hmCG +]:[(5mCG+5hmCG) -]"]
group2<-all_gene$gene_name[all_gene$group=="[5mCG -]:[5hmCG +]"]
write.table(group1,"group1_gene_list.txt",row.names=F,col.names=F,quote=F)
write.table(group2,"group2_gene_list.txt",row.names=F,col.names=F,quote=F)
