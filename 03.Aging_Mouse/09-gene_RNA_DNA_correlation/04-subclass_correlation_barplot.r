library(ggplot2)
library(Seurat)
library(dplyr)
library(stringr)
library(reshape2)
library(ggpubr)
library(ggpointdensity) 
library(RColorBrewer)
library(cowplot)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

metainfo<-read.csv("../../04data/02.metainfo/03.Aging_Mouse/RNA_DNA_match_name_QC.aged.csv")
old_subclass<-unique(metainfo$subclass_label)
subclass_order = readRDS("../../04.data/04.config_files/order.subclass.rds")
subclass_order<-subclass_order[subclass_order%in%old_subclass]
subclass_order<-gsub(" ",".",subclass_order)
subclass_order<-gsub("-",".",subclass_order)
subclass_order<-gsub("/",".",subclass_order)
setwd(paste0(indir,"/Joint_Cabernet_old+Joint_Cabernet_young"))

subclass_color<-readRDS("../../04.data/04.config_files/color.subclass.rds")

##############################   our
correlation_subclass_RNA_DNA_mean<-function(age,varim,datatype,mc_type){
    indir=paste0("./",age,"_corrected/02data_expr/subclass/")
    #读取RNA文件
    class_expr <- read.csv(paste0(indir,"RNA_subclass_expr.csv"))
    colnames(class_expr) <- c("gene","class","RNA")
    #读取DNA文件
    tmp<- read.csv(paste0(indir,datatype,"_",mc_type,"_",varim,"_subclass_expr.csv"))
    colnames(tmp) <- c("gene","class","DNA")
    tmp$gene<-gsub("\\..*","",tmp$gene)
    all_data <- merge(class_expr,tmp,by = c("gene","class"))
    all_data <- subset(all_data ,RNA != 0)
    subclass<-unique(all_data$class)
    result<-data.frame()
    for(cl in subclass){
        print(paste0(datatype,"-",mc_type,"-",cl))
        cor_data<-all_data[all_data$class==cl,]
        cor_data<-cor_data[complete.cases(cor_data),]
        cor<-cor.test(as.numeric(cor_data$RNA),as.numeric(cor_data$DNA),alternative = "two.sided",method = "pearson",conf.level = 0.95)
        res<-data.frame(age=age,datatype=datatype,var_dim=varim,mc_type=mc_type,subclass=cl,correlation=cor[["estimate"]][["cor"]],p_value=cor$p.value)
        result<-rbind(result,res)
    }  
    return(result)
}

# data<-correlation_subclass_RNA_DNA_mean("old","genebody","5hmC","CG")

#####old and young combination
draw_barplot<-function(correlation_dataframe){
    outdir=outdir
    correlation_dataframe<-merge(correlation_dataframe,subclasses,by.x="subclass",by.y="subclass_order_change")
    correlation_dataframe$age<-factor(correlation_dataframe$age,levels=c("young","old"))
    correlation_dataframe$subclass_order<-factor(correlation_dataframe$subclass_order,levels=subclass_order)
    correlation_dataframe$sig<-ifelse(correlation_dataframe$p_value < 0.001, "***",
                                    ifelse(correlation_dataframe$p_value  < 0.01, "**",
                                        ifelse(correlation_dataframe$p_value  < 0.05, "*","ns")))
    if(max(correlation_dataframe$correlation)<0){
        y_max <-0.01
    }else{
        y_max <- max(correlation_dataframe$correlation) + max(correlation_dataframe$correlation) * 0.1  
    }
    datatype<-unique(correlation_dataframe$datatype)
    if(datatype=="true_5mC"){dataty="5mC"
        dtype<-paste0(dataty,substr(unique(correlation_dataframe$mc_type),2,2))}
    if(datatype=="5hmC"){dataty="5hmC"
        dtype<-paste0(dataty,substr(unique(correlation_dataframe$mc_type),2,2))}
    if(datatype=="5mC"){dtype<-paste0("5mC",substr(unique(correlation_dataframe$mc_type),2,2),"+5hmC",substr(unique(correlation_dataframe$mc_type),2,2))}

    p<-ggplot(correlation_dataframe, aes(x = subclass_order, y = correlation, fill = age,group=age)) +
        geom_bar(stat = 'identity', position = 'dodge') + 
        scale_fill_manual(values = c('young'='#3498db',"old"='#e74c3c')) +     
        # facet_wrap(~ mc_type, ncol = 2) +                 
        labs(x="Subclasses",y = 'Correlation', fill = 'age') +
        ggtitle(paste0(dtype," genebody subclass correlation barplot"))+
        theme(legend.position="none",
        # strip.position="bottom",
        strip.background = element_blank(), 
        strip.text = element_text(face = "bold", colour = "black"), 
        plot.title = element_text(hjust = 0.5,family="ArialMT"),
        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
        axis.title.y = element_text(color="black", size=10,family="ArialMT"),
        axis.text.x = element_text(angle=30,hjust=1,color="black", size=8,family="ArialMT"),
        axis.text.y = element_text(color="black", size=8,family="ArialMT"),
        axis.line=element_line(size=0.6),
        axis.ticks.y=element_line(size=0.6),
        axis.ticks.x = element_blank() ,
        panel.background = element_blank(),  
        panel.border=element_blank(),
        panel.grid=element_blank(),
        element_line(linetype = "dashed"), 
        panel.grid.minor = element_blank())
    ggsave(paste0(outdir,dtype,"_Joint_Cabernet_RNA_genebody_subclass_correlation_barplot.pdf"),plot=p,width=8,height=3,dpi=300) 
}


    


datatypes = c("5mC","true_5mC","5hmC")
varim = c("genebody")
mc_types=c("CG","CH")
for(datatype in datatypes){
    for(mc_type in mc_types){
        old_data<-correlation_subclass_RNA_DNA_mean("old","genebody",datatype,mc_type)
        young_data<-correlation_subclass_RNA_DNA_mean("young","genebody",datatype,mc_type)
        all_data<-rbind(old_data,young_data)
        draw_barplot(all_data)
    }
}

