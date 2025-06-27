library(Seurat)
library(dplyr)
library(future)
library(presto)
library(stringr)
library(getopt)
library(data.table)
library(reshape2)
library(ggpubr)
library(ggplot2)
library(ggunchained)
library(cowplot)
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(scales)
library(RColorBrewer)
library(ggpointdensity) 
library(ggtext)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

metainfo = readRDS("../../04data/02.metainfo/03.Aging_Mouse/RNA_DNA_match_name_QC.aged.csv")
head(metainfo)
colnames(metainfo)

subclass_order = readRDS("../../04.data/04.config_files/order.subclass.rds")

####### process data  #####
#subclass mean
promoter_genebody_subclass_mean<-function(age,region,datatype){
    data<-read.csv(paste0(age,"_",region,"_",datatype,".mean_methyl.csv"),row.names=1,header=T)
    data$age = age
    data$sampleid = str_replace(data$sampleid ,"allc_","")
    data$sampleid = str_replace(data$sampleid,".mm10.dna.tsv.gz","")
    if(datatype=="hmCG"){
        data$subclass = metainfo[match(data$sampleid,metainfo$hmC_ID),"subclass_label"]
    }
    if(datatype=="mCG"){
        data$subclass = metainfo[match(data$sampleid,metainfo$mC_ID),"subclass_label"]
    }
    subclass_mean<-data %>%
        group_by(subclass) %>%
        summarise(mean_level = mean(gene_mean_level, na.rm = TRUE))   %>%
        as.data.frame() 
    return(subclass_mean)
}
other_subclass_mean<-function(region,datatype){
    data<-read.csv(paste0(region,"_",datatype,".mean_methyl.csv"),row.names=1,header=T)
    data$sampleid = str_replace(str_replace(data$sampleid,".mm10.dna.tsv.gz",""),"allc_","")
    if(datatype=="hmCG"){
        data$age = metainfo[match(data$sampleid ,metainfo$hmC_ID),"old_young"]
        data$subclass = metainfo[match(data$sampleid,metainfo$hmC_ID),"subclass_label"]
    }
    if(datatype=="mCG"){
        data$age = metainfo[match(data$sampleid ,metainfo$mC_ID),"old_young"]
        data$subclass = metainfo[match(data$sampleid,metainfo$mC_ID),"subclass_label"]
    }
    subclass_mean<-data %>%
        group_by(subclass,age) %>%
        summarise(mean_level = mean(gene_mean_level, na.rm = TRUE))   %>%
        as.data.frame() 
    return(subclass_mean)
}
true_subclass_mean<-function(age,region){
    if(region%in%c("promoter","genebody")){
        hmCG<-promoter_genebody_subclass_mean(age,region,"hmCG")
        mCG<-promoter_genebody_subclass_mean(age,region,"mCG")
        true_5hmCG<-data.frame(subclass=hmCG$subclass,mean_level=mCG$mean_level-hmCG$mean_level)
    }else{
        hmCG<-other_subclass_mean(region,"hmCG")
        hmCG<-hmCG[hmCG$age==age,]
        mCG<-other_subclass_mean(region,"mCG")
        mCG<-mCG[mCG$age==age,]
        true_5hmCG<-data.frame(subclass=hmCG$subclass,mean_level=mCG$mean_level-hmCG$mean_level)
    }
    return(true_5hmCG)
}
#old-young
old_young_diff<-function(region,datatype){
    if(datatype=="true_5mCG"){
        old_region<-true_subclass_mean("old",region)
        young_region<-true_subclass_mean("young",region)
    }else{
        if(region%in%c("promoter","genebody")){
            old_region<-promoter_genebody_subclass_mean("old",region,datatype)
            young_region<-promoter_genebody_subclass_mean("young",region,datatype)
        }else{
            region<-other_subclass_mean(region,datatype)
            old_region<-region[region$age=="old",]
            young_region<-region[region$age=="young",]
    }}
    diff<-data.frame(subclass=unique(metainfo$subclass_label),diff=old_region$mean_level-young_region$mean_level)
    return(diff)
}


create_diff_dataframe<-function(regions,datatype){
    total_diff_matrix<-data.frame(subclass=unique(metainfo$subclass_label))
    for(re in regions){
        diff<-old_young_diff(re,datatype)
        total_diff_matrix<-merge(total_diff_matrix,diff,by="subclass")
        if(re=="Intergenetic"){
            colnames(total_diff_matrix)[ncol(total_diff_matrix)]<-"Intergenic"
        }else{
            colnames(total_diff_matrix)[ncol(total_diff_matrix)]<-re
        }
    }
    total_diff_matrix$subclass<-factor(total_diff_matrix$subclass,levels=subclass_order)
    total_diff_matrix<-total_diff_matrix[order(total_diff_matrix$subclass),]
    rownames(total_diff_matrix)<-total_diff_matrix$subclass
    total_diff_matrix<-total_diff_matrix[,-1]
    return(total_diff_matrix)
}
    






#################  plot
unnormalized_boxplot_adjust_with_significant<-function(dataframe,datatype){
    inte.col = readRDS("../../04.data/04.config_files/subclass_new.col_latest.rds")
    inte.col = c(inte.col,c("CA2-FC-IG Glut"="#89C75F","L6b CTX Glut"="#0C727C","Lamp5 Lhx6 Gaba"="#90D5E4","PAL-STR Gaba-Chol" ="#00ae9d","DG-PIR Ex IMN"="#1d953f",
        "Vip Gaba"="#009ad6","L2/3 IT PIR-ENTl Glut"="#6E4B9E","L6 IT CTX Glut"="#AA0DFE","HPF CR Glut"='#e74c3c',"Pvalb chandelier Gaba"="#A6BDD7","STR D1 Sema5a Gaba" ="#B32851",
        "OB-mi Frmd7 Gaba"='#5AC2F1FF',"OB Trdn Gaba"="#e4c6d0","OB Meis2 Thsd7b Gaba"="#f9906f","VLMC NN"="#ffc773","Sst Chodl Gaba"="#88c4e8","STR Prox1 Lhx6 Gaba"="#eb7f54",
        "OT D3 Folh1 Gaba"="#815463","ABC NN" ="#253494","BAM NN" ="#FFFF00","Endo NN"="#d6ecf0","Peri NN"="#DEA0FD","Lymphoid NN"="#808080","SMC NN"="#bce672",'zeng'= 'lightgrey',
        "our"='lightgrey',"old"="lightgrey","young"="lightgrey","IT AON-TT-DP Glut"="#F6768E","LA-BLA-BMA-PA Glut"="#ff3300","COAa-PAA-MEA Barhl2 Glut"="#801dae"))
    inte.col<-inte.col[rownames(dataframe)]

    if(datatype=="hmCG"){dtype="5hmCG"
        limits=c(0,0.065)}
    if(datatype=="mCG"){dtype="5mCG+5hmCG"
        limits=c(-0.01,0.06)}
    if(datatype=="true_5mCG"){dtype="5mCG"
        limits=c(-0.04,0.04)}

    # 转换数据为长格式
    df_long <- reshape2::melt(as.matrix(dataframe))
    colnames(df_long) <- c("subclass", "gene_element", "value")
    elements<-colnames(dataframe)

    # 绘制图形
    p<-ggplot(df_long, aes(x = gene_element, y = value)) +
    geom_boxplot(color = "black", width = 0.8, outlier.shape = NA,linewidth=0.7) + 
     geom_line(
            aes( group = subclass), color = "lightgrey",
            alpha = 0.5, size = 0.8
        ) +
        geom_point(
            aes( color = subclass),
            alpha = 0.8, size = 1.5
        ) +
    stat_compare_means(comparisons = combn(elements, 2, simplify = FALSE), 
                         method = "t.test",
                         paired = T,     
                         symnum.args = list(cutpoints=c(0,0.001,0.01,0.05,1),
                                                 symbols=c("***","**","*","ns")),
                         label = "p.signif",
                         size=3.5,step.increase = 0.12)+ 
    scale_y_continuous(limits=limits)+
    scale_color_manual(values = inte.col) +
    theme_minimal() +
    theme(
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.5),
        axis.ticks = element_line(color = "black", linewidth = 0.5),
        legend.position = "none",
        # legend.key.size = unit(0.5, "cm"),
        axis.text.x = element_text(hjust =0.5, color = "black"),#angle=45,hjust =1
        axis.text.y = element_text(color = "black"),
        axis.title = element_text(color = "black", face = "bold",size=10)
    ) +
    labs(x = "Gene Elements", y = paste0(dtype," diff"))

    pdf(paste0("../../output/03.Aging_Mouse/08-gene_elements_diff_boxplot/promoter_genebody_intergenic_enhancer_",dtype,"_diff_boxplot_with_significant.pdf"),width=4.5,height=4)
    plot(p)
    dev.off()
}

regions<-c("promoter","genebody","Intergenetic","Enhancer")
datatypes<-c("mCG","hmCG")
for(datatype in datatypes){
    data<-create_diff_dataframe(regions,datatype)
    unnormalized_boxplot_adjust_with_significant(data,datatype)
}



