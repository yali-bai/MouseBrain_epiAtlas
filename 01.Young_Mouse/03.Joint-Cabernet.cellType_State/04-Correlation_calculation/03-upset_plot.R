library(UpSetR)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(dplyr)

setwd("./all_cell_correlation")
result<-data.frame()
for(dt in c("5hmC","5mC","true_5mC")){
  df<-fread(paste0("all_cells_",dt,"_CG_genebody_gene_correlation_results.csv"),data.table = F)
  result<-rbind(result,df)
}
data<-result[result$`Adjusted P-value`<0.05,]
data$corr_direction<-ifelse(data$Correlation>0,"Positive correlation","Negative correlation")
colnames(data)[1]<-"gene_id"
mc='CG'
vardim='genebody'
aa <- subset(data,mc_type==mc&var_dim==vardim)  
              set1 <- subset(aa,datatype == "true_5mC"&corr_direction=="Positive correlation")$gene_id  #mc
              set2 <- subset(aa,datatype == "true_5mC"&corr_direction=="Negative correlation")$gene_id  
              set3 <- subset(aa,datatype == "5hmC"&corr_direction=="Positive correlation")$gene_id  #hmc
              set4<- subset(aa,datatype == "5hmC"&corr_direction=="Negative correlation")$gene_id
              set5 <- subset(aa,datatype == "5mC"&corr_direction=="Negative correlation")$gene_id  #mc+hmc
              set6 <- subset(aa,datatype == "5mC"&corr_direction=="Positive correlation")$gene_id

              name<-c(paste0("[5m",mc," +]"),paste0("[5m",mc," -]"),paste0("[5hm",mc," +]"),paste0("[5hm",mc," -]"),paste0("[unmotified ",mc," +]"),paste0("[unmotified ",mc," -]"))

              color<-c(rep('#F2CD5C', 2), rep("#F8766D", 2),rep("#aa96da",2))

              length1<-c(length(set1),length(set2),length(set3),length(set4),length(set5),length(set6))
              list1 <- list(set1 = set1,set2 = set2,set3 = set3,set4 = set4,set5 = set5,set6 = set6)
              names(list1) <- name
              non_empty_dfs <- lapply(list1, function(x) { 
              if (length(x) > 0) {  
                  return(x)  
              }  }) 
              upset_list1 <- non_empty_dfs[!sapply(non_empty_dfs, is.null)]      
              color_vector1 <- rev(color[length1!=0])
p1 <- upset(fromList(upset_list1),
                    nsets = length(upset_list1), 
                    nintersects = 80, 
                    sets = rev(c(name[length1!=0])), 
                    keep.order = TRUE, 
                    number.angles = 0, 
                    point.size = 4, 
                    line.size = 1, 
                    mainbar.y.label = "Intersection size", 
                    main.bar.color = 'black', 
                    matrix.color = "black", 
                    sets.x.label = "Set size", 
                    sets.bar.color = color_vector1, 
                    mb.ratio = c(0.7, 0.3), 
                    order.by = "freq", 
                    text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.6), 
                    shade.color = "#12507B", 
              )
pdf("all_cells_correlation_upset.pdf",width =8, height =5)
print(p1)
dev.off()


