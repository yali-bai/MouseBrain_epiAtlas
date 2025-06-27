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

subclass_color<-readRDS("../../../04.data/04.config_files/color.subclass.rds")

subclass_order = read.table("../../../04.data/04.config_files/subclass_order_for_integration_with_zeng.txt",sep="\t",header=F,stringsAsFactors=F)
subclass_order<-subclass_order$V1
subclass_order_change<-gsub(" ",".",subclass_order)
subclass_order_change<-gsub("-",".",subclass_order_change)
subclass_order_change<-gsub("/",".",subclass_order_change)
subclasses<-data.frame(subclass_order,subclass_order_change)


setwd(indir)


##############################   our
correlation_subclass_RNA_DNA_mean<-function(varim,datatype,mc_type){
    indirs=paste0(indir,"/merge_expr/subclass/")
    #读取RNA文件
    class_expr <- read.csv(paste0(indirs,"Joint_Cabernet_RNA_subclass_expr.csv"))
    colnames(class_expr) <- c("gene","class","RNA")
    #读取DNA文件
    tmp<- read.csv(paste0(indirs,datatype,"_",mc_type,"_",varim,"_subclass_expr.csv"))
    colnames(tmp) <- c("gene","class","DNA")
    # tmp$gene<-gsub("\\..*","",tmp$gene)
    all_data <- merge(class_expr,tmp,by = c("gene","class"))
    all_data <- subset(all_data ,RNA != 0)
    subclass<-unique(all_data$class)
    result<-data.frame()
    for(cl in subclass){
        print(paste0(datatype,"-",mc_type,"-",cl))
        cor_data<-all_data[all_data$class==cl,]
        cor_data<-cor_data[complete.cases(cor_data),]
        cor<-cor.test(as.numeric(cor_data$RNA),as.numeric(cor_data$DNA),alternative = "two.sided",method = "pearson",conf.level = 0.95)
        res<-data.frame(datatype=datatype,var_dim=varim,mc_type=mc_type,subclass=cl,correlation=cor[["estimate"]][["cor"]],p_value=cor$p.value)
        res$type<-ifelse(datatype=="5hmC",paste0("hmC",substr(unique(mc_type),2,2)),
                                    ifelse(datatype=="5mC",paste0("mC",substr(unique(mc_type),2,2),"+hmC",substr(unique(mc_type),2,2)),
                                        paste0("mC",substr(unique(mc_type),2,2))))
        result<-rbind(result,res)
    }  
    return(result)
}

# data<-correlation_subclass_RNA_DNA_mean("genebody","5hmC","CG")

#####
draw_subclass_barplot<-function(correlation_dataframe){
    outdir=outdir
    correlation_dataframe<-merge(correlation_dataframe,subclasses,by.x="subclass",by.y="subclass_order_change")
    correlation_dataframe$subclass_order<-factor(correlation_dataframe$subclass_order,levels=subclass_order)
    correlation_dataframe$type<-factor(correlation_dataframe$type,levels=c("mCG+hmCG","mCG","hmCG","mCH+hmCH","mCH","hmCH"))
    correlation_dataframe$sig<-ifelse(correlation_dataframe$p_value < 0.001, "***",
                                    ifelse(correlation_dataframe$p_value  < 0.01, "**",
                                        ifelse(correlation_dataframe$p_value  < 0.05, "*","ns")))
    if(max(correlation_dataframe$correlation)<0){
        y_max <-0.01
    }else{
        y_max <- max(correlation_dataframe$correlation) + max(correlation_dataframe$correlation) * 0.1  
    }
    p<-ggplot(correlation_dataframe, aes(x = subclass_order, y = correlation, fill = type,group=type)) +
        geom_bar(stat = 'identity', position = 'dodge') + 
        labs(x="Subclasses",y = 'Correlation') +
        ggtitle("genebody subclass correlation barplot")+
        theme(legend.position="right",
        strip.background = element_blank(), 
        strip.text = element_text(face = "bold", colour = "black"), 
        plot.title = element_text(hjust = 0.5,family="ArialMT"),
        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
        axis.title.y = element_text(color="black", size=10,family="ArialMT"),
        axis.text.x = element_text(angle=90,hjust=1,color="black", size=8,family="ArialMT"),
        axis.text.y = element_text(color="black", size=8,family="ArialMT"),
        axis.line=element_line(size=0.6),
        axis.ticks.y=element_line(size=0.6),
        axis.ticks.x = element_blank() ,
        panel.background = element_blank(),  
        panel.border=element_blank(),
        panel.grid=element_blank(),
        element_line(linetype = "dashed"), 
        panel.grid.minor = element_blank())

    # ggsave(paste0(outdir,dtype,"_Joint_Cabernet_RNA_genebody_subclass_correlation_barplot.png"),plot=p,width=8,height=3,dpi=300)     
    ggsave(paste0(outdir,"Joint_Cabernet_RNA_genebody_subclass_correlation_barplot_color_by_type.pdf"),plot=p,width=10,height=3,dpi=300) 
# correlation_dataframe$subclass_color<-subclass_color[match(correlation_dataframe$subclass_order,names(subclass_color))]
p2<-ggplot(correlation_dataframe, aes(x = subclass_order, y = correlation, fill = subclass_order,group=type)) +
        geom_bar(stat = 'identity', position = 'dodge') + 
        scale_fill_manual(values = subclass_color) +  
        labs(x="Subclasses",y = 'Correlation') +
        ggtitle("genebody subclass correlation barplot")+
        theme(legend.position="none",
        strip.background = element_blank(), 
        strip.text = element_text(face = "bold", colour = "black"), 
        plot.title = element_text(hjust = 0.5,family="ArialMT"),
        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
        axis.title.y = element_text(color="black", size=10,family="ArialMT"),
        axis.text.x = element_text(angle=90,hjust=1,color="black", size=8,family="ArialMT"),
        axis.text.y = element_text(color="black", size=8,family="ArialMT"),
        axis.line=element_line(size=0.6),
        axis.ticks.y=element_line(size=0.6),
        axis.ticks.x = element_blank() ,
        panel.background = element_blank(),  
        panel.border=element_blank(),
        panel.grid=element_blank(),
        element_line(linetype = "dashed"), 
        panel.grid.minor = element_blank())
    ggsave(paste0(outdir,"Joint_Cabernet_RNA_genebody_subclass_correlation_barplot_color_by_subclass.pdf"),plot=p2,width=10,height=3,dpi=300) 
}

draw_threeclass_barplot<-function(correlation_dataframe){
    outdir=outdir
    correlation_dataframe<-merge(correlation_dataframe,subclasses,by.x="subclass",by.y="subclass_order_change")
    correlation_dataframe$subclass_order<-factor(correlation_dataframe$subclass_order,levels=subclass_order)
    correlation_dataframe$three_class<-as.vector(correlation_dataframe$subclass_order)
    correlation_dataframe$three_class[grep("Glut",correlation_dataframe$three_class)]<-"Exc"
    correlation_dataframe$three_class[grep("Gaba",correlation_dataframe$three_class)]<-"Inh"
    correlation_dataframe$three_class[grep("IMN",correlation_dataframe$three_class)]<-"Inh"
    correlation_dataframe$three_class[grep("NN",correlation_dataframe$three_class)]<-"NN"
    summary_data <- correlation_dataframe %>%
        group_by(three_class,type) %>%
        summarise(
            mean_corr = mean(correlation, na.rm = TRUE),
            min_corr = min(correlation, na.rm = TRUE),    
            max_corr = max(correlation, na.rm = TRUE)    
        )%>%as.data.frame()
    summary_data$type<-factor(summary_data$type,levels=c("mCG+hmCG","mCG","hmCG","mCH+hmCH","mCH","hmCH"))
    summary_data$three_class<-factor(summary_data$three_class,levels=c("Exc","Inh","NN"))

    p<-ggplot(summary_data, aes(x = three_class, y = mean_corr, fill = type,group=type)) +
        geom_bar(stat = 'identity', position = position_dodge(width = 0.8)) +  
        geom_errorbar(aes(x = three_class,ymin = min_corr, ymax = max_corr,group=type),  
            width = 0.2,color = "black",linewidth = 0.5,position = position_dodge(width = 0.8))+                      
        labs(x="three classes",y = 'Correlation') +
        ggtitle("genebody three classes correlation barplot")+
        theme(legend.position="right",
        strip.background = element_blank(), 
        strip.text = element_text(face = "bold", colour = "black"),  
        plot.title = element_text(hjust = 0.5,family="ArialMT"),
        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
        axis.title.y = element_text(color="black", size=10,family="ArialMT"),
        axis.text.x = element_text(hjust=0.5,color="black", size=8,family="ArialMT"),
        axis.text.y = element_text(color="black", size=8,family="ArialMT"),
        axis.line=element_line(size=0.6),
        axis.ticks.y=element_line(size=0.6),
        axis.ticks.x = element_blank() ,
        panel.background = element_blank(), 
        panel.border=element_blank(),
        panel.grid=element_blank(),
        element_line(linetype = "dashed"), 
        panel.grid.minor = element_blank())
    ggsave(paste0(outdir,"Joint_Cabernet_RNA_genebody_three_class_correlation_barplot_color_by_type.pdf"),plot=p,width=5,height=3,dpi=300) 
p2<-ggplot(summary_data, aes(x = three_class, y = mean_corr, fill = three_class,group=type)) +
        geom_bar(stat = 'identity', position = position_dodge(width = 0.8)) +  
        geom_errorbar(aes(ymin = min_corr, ymax = max_corr), 
            width = 0.2,color = "black",linewidth = 0.5,position = position_dodge(width = 0.8) ) +
        labs(x="three classes",y = 'Correlation') +
        ggtitle("genebody three classes correlation barplot")+
        theme(legend.position="none",
        strip.background = element_blank(), 
        strip.text = element_text(face = "bold", colour = "black"), 
        plot.title = element_text(hjust = 0.5,family="ArialMT"),
        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
        axis.title.y = element_text(color="black", size=10,family="ArialMT"),
        axis.text.x = element_text(hjust=0.5,color="black", size=8,family="ArialMT"),
        axis.text.y = element_text(color="black", size=8,family="ArialMT"),
        axis.line=element_line(size=0.6),
        axis.ticks.y=element_line(size=0.6),
        axis.ticks.x = element_blank() ,
        panel.background = element_blank(), 
        panel.border=element_blank(),
        panel.grid=element_blank(),
        element_line(linetype = "dashed"), 
        panel.grid.minor = element_blank())
    ggsave(paste0(outdir,"Joint_Cabernet_RNA_genebody_three_class_correlation_barplot_color_by_class.pdf"),plot=p2,width=5,height=3,dpi=300) 
}
    


datatypes = c("5mC","true_5mC","5hmC")
varim = c("genebody")
mc_types=c("CG","CH")
all_data<-data.frame()
for(datatype in datatypes){
    for(mc_type in mc_types){
        data<-correlation_subclass_RNA_DNA_mean("genebody",datatype,mc_type)
        all_data<-rbind(all_data,data)
    }
}
draw_subclass_barplot(all_data)
draw_threeclass_barplot(all_data)
