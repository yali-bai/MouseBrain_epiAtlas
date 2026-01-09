#Draw line diagrams of methylation level in different subclasses.
library(dplyr)
library(readxl)
library(Seurat)
library(cowplot)
library(ggplot2)

# "indir" is a custom input path, and "outdir" is a custom output path.
indir="./"
# outdir=""

###########  hmc_CG  ##################

file_path <- paste0(indir,"/hmC/CG")
#Extract all file names in the folder  
all_files <- list.files(file_path, full.names = FALSE)  
  
#Filter out file names ending in.txt  
tmp_files <- all_files[grepl("\\.tmp$", all_files, ignore.case = TRUE)] 
tmp<-strsplit(tmp_files,"_")
matrix_result <- do.call(rbind, lapply(tmp, function(x) as.character(x)))%>%as.data.frame()

class=matrix_result[!duplicated(matrix_result[,2]),2]
group<-c("=0","<0.7960165",">0.7960165")  # threshold to sep gene as three classes: No Expr., Low Expr., and High Expr.


plots<-list()
for(n in c(1:3)){
    k = c(8,18,24)[n]
    if(k == 8){
        other_theme = theme(legend.position="none",
            plot.title = element_text(hjust = 0.5,family="ArialMT"),
            title = element_text(size=25,hjust = 0.5,family="ArialMT"),
            axis.title.x = element_text(color="black", size=25,family="ArialMT"),
            #axis.title.y = element_text(color="black", size=13,family="ArialMT"),
            axis.text.x = element_blank(),
            axis.text.y = element_text(color="black", size=30,family="ArialMT"),
            axis.line=element_line(size=0.8),
            axis.ticks=element_line(size=0.6),
            panel.border=element_blank(),
            panel.grid=element_blank(),
            plot.margin = margin(5.5, 12.5, 5.5, 5.5, "pt"))
    }else if(k == 18){
        other_theme = theme(legend.position="none",
            plot.title = element_blank(),
            title = element_text(size=25,hjust = 0.5,family="ArialMT"),
            axis.title.x = element_blank(),
            #axis.title.y = element_text(color="black", size=13,family="ArialMT"),
            axis.text.x = element_blank(),
            axis.text.y = element_text(color="black", size=30,family="ArialMT"),
            axis.line=element_line(size=0.8),
            axis.ticks=element_line(size=0.6),
            panel.border=element_blank(),
            panel.grid=element_blank(),
            plot.margin = margin(5.5, 12.5, 5.5, 5.5, "pt"))
    }else{
        other_theme = theme(legend.position="none",
            plot.title = element_blank(),
            title = element_text(size=25,hjust = 0.5,family="ArialMT"),
            axis.title.x = element_blank(),
            #axis.title.y = element_text(color="black", size=13,family="ArialMT"),
            axis.text.x = element_text(color="black", size=25,family="ArialMT"),
            axis.text.y = element_text(color="black", size=30,family="ArialMT"),
            axis.line=element_line(size=0.8),
            axis.ticks=element_line(size=0.6),
            panel.border=element_blank(),
            panel.grid=element_blank(),
            plot.margin = margin(5.5, 12.5, 5.5, 5.5, "pt"))
    }

###########  mc_CG  ##################

file_path <- paste0(indir,"/mC/CG")

#for(k in 1:length(class)){
    cat(k,'\n')
    tmp_df<-data.frame()
    for(i in c(1:3)){
        data<-read.csv(paste0(file_path,"/mC_",class[k],"_group",i,"_CG.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df<-rbind(tmp_df,data1)
    }
    tmp_df$value<-tmp_df$value*100
    p <- ggplot() +
    geom_vline(xintercept = c(0,400,1400,1800), color = "grey", linetype = "solid", size = 0.2) +
    geom_hline(yintercept = c(0,20,40,60,80), color = "grey", linetype = "solid", size = 0.2) + 
    geom_line(data = tmp_df, mapping = aes(x=pos, y=value, color =type,group=type),alpha=1,linewidth=0.3) +
    labs(title = "mCG+hmCG") +
     xlab("") +
     ylab(paste0(class[k]))+ 
    scale_color_manual(values = c("=0"="#0278ae","<0.7960165"="#ffd480",">0.7960165"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,20,40,60,80),limits=c(0,90))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(axis.title.y = element_text(color="black", size=25,family="ArialMT"))+
    other_theme
   plots[[(n-1)*6+1]]<-p
#}

#####################   true5mc_CG    ################################


#for(k in 1:length(class)){
    cat(k,'\n')
    tmp_df_hmcg<-data.frame()
    for(i in 1:3){
        data<-read.csv(paste0(indir,"/hmC/CG/hmC_",class[k],"_group",i,"_CG.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df_hmcg<-rbind(tmp_df_hmcg,data1)
    }
    tmp_df_mcg<-data.frame()
    for(i in 1:3){
        data<-read.csv(paste0(indir,"/mC/CG/mC_",class[k],"_group",i,"_CG.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df_mcg<-rbind(tmp_df_mcg,data1)
    }    
    true_mcg<-data.frame(value=rep(0,nrow(tmp_df_hmcg)),SeqType=tmp_df_hmcg$SeqType,pos=tmp_df_hmcg$pos,type=tmp_df_hmcg$type)
    true_mcg$value<-tmp_df_mcg$value-tmp_df_hmcg$value
    true_mcg$value<-true_mcg$value*100
    true_mcg$value[true_mcg$value<0]<-0

    p <- ggplot() +
    geom_vline(xintercept = c(0,400,1400,1800), color = "grey", linetype = "solid", size = 0.2) + 
    geom_hline(yintercept = c(20,40,60,80), color = "grey", linetype = "solid", size = 0.2) +
    geom_line(data = true_mcg, mapping = aes(x=pos, y=value, color =type,group=type),alpha=1,linewidth=0.3) +
    labs(title = "mCG") +
     xlab("") +
     ylab(paste0(class[k]))+   
    scale_color_manual(values = c("=0"="#0278ae","<0.7960165"="#ffd480",">0.7960165"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,20,40,60,80),limits=c(0,80))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(axis.title.y = element_blank())+
    other_theme
   plots[[(n-1)*6+2]]<-p
#}




#for(k in 1:length(class)){
    file_path <- paste0(indir,"/hmC/CG")
    cat(k,'\n')
    tmp_df<-data.frame()
    for(i in c(1:3)){
        data<-read.csv(paste0(file_path,"/hmC_",class[k],"_group",i,"_CG.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df<-rbind(tmp_df,data1)
    }
    tmp_df$value<-tmp_df$value*100
    p <- ggplot() +
    geom_vline(xintercept = c(0,400,1400,1800), color = "grey", linetype = "solid", size = 0.2) + #Add vertical guides
    geom_hline(yintercept = c(0,15,30,45), color = "grey", linetype = "solid", size = 0.2) + # Add horizontal guides 
    geom_line(data = tmp_df, mapping = aes(x=pos, y=value, color =type,group=type),alpha=1,linewidth=0.3) +
    labs(title = "hmCG") +
     xlab("") +
     ylab(paste0(class[k]))+
    scale_color_manual(values = c("=0"="#0278ae","<0.7960165"="#ffd480",">0.7960165"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,15,30,45),limits=c(0,45))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(axis.title.y = element_blank())+
    other_theme     
   plots[[(n-1)*6+3]]<-p  
#}

############### hmcg/(mcg+hmcg)   ###############################


#for(k in 1:length(class)){
    cat(k,'\n')
    tmp_df_hmcg<-data.frame()
    for(i in 1:3){
        data<-read.csv(paste0(indir,"/hmC/CG/hmC_",class[k],"_group",i,"_CG.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df_hmcg<-rbind(tmp_df_hmcg,data1)
    }
    tmp_df_mcg<-data.frame()
    for(i in 1:3){
        data<-read.csv(paste0(indir,"/mC/CG/mC_",class[k],"_group",i,"_CG.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df_mcg<-rbind(tmp_df_mcg,data1)
    }
    hmcg_mh<-data.frame(value=rep(0,nrow(tmp_df_hmcg)),SeqType=tmp_df_hmcg$SeqType,pos=tmp_df_hmcg$pos,type=tmp_df_hmcg$type)
    hmcg_mh$value<-tmp_df_hmcg$value/tmp_df_mcg$value
    hmcg_mh$value<-hmcg_mh$value*100

    p <- ggplot() +
    geom_vline(xintercept = c(0,400,1400,1800), color = "grey", linetype = "solid", size = 0.2) + 
    geom_hline(yintercept = c(0,20,40), color = "grey", linetype = "solid", size = 0.2) +
    geom_line(data = hmcg_mh, mapping = aes(x=pos, y=value, color =type,group=type),alpha=1,linewidth=0.3) +

    labs(title = "hmCG/(mCG+hmCG)") +
     xlab("") +
     ylab(paste0(class[k]))+ 
    scale_color_manual(values = c("=0"="#0278ae","<0.7960165"="#ffd480",">0.7960165"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,20,40),limits=c(0,56))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(axis.title.y = element_blank())+
    other_theme
   plots[[(n-1)*6+4]]<-p 
#}

####################  true-5mCH   ##############################

#for(k in 1:length(class)){
    cat(k,'\n')
    tmp_df_hmch<-data.frame()
    for(i in 1:3){
        data<-read.csv(paste0(indir,"/hmC/CH/hmC_",class[k],"_group",i,"_CH.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df_hmch<-rbind(tmp_df_hmch,data1)
    }
   tmp_df_mch<-data.frame()
    for(i in 1:3){
        data<-read.csv(paste0(indir,"/mC/CH/mC_",class[k],"_group",i,"_CH.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df_mch<-rbind(tmp_df_mch,data1)
    }   
    true_mch<-data.frame(value=rep(0,nrow(tmp_df_hmch)),SeqType=tmp_df_hmch$SeqType,pos=tmp_df_hmch$pos,type=tmp_df_hmch$type)
    true_mch$value<-tmp_df_mch$value-tmp_df_hmch$value
    true_mch$value<-true_mch$value*100
    true_mch$value[true_mch$value<0]<-0

    p <- ggplot() +
    geom_vline(xintercept = c(0,400,1400,1800), color = "grey", linetype = "solid", size = 0.2) +
    geom_hline(yintercept = c(0,2,4), color = "grey", linetype = "solid", size = 0.2) + 
    geom_line(data = true_mch, mapping = aes(x=pos, y=value, color =type,group=type),alpha=1,linewidth=0.3) +
    labs(title = "mCH") +
     xlab("") +
     ylab(paste0(class[k]))+ 
    scale_color_manual(values = c("=0"="#0278ae","<0.7960165"="#ffd480",">0.7960165"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,2,4),limits=c(0,4))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(axis.title.y = element_blank())+
    other_theme
   plots[[(n-1)*6+5]]<-p

#}

###########  hmc_CH  ##################

file_path <- paste0(indir,"/hmC/CH")

#for(k in 1:length(class)){
    cat(k,'\n')
    tmp_df<-data.frame()
    for(i in c(1:3)){
        data<-read.csv(paste0(file_path,"/hmC_",class[k],"_group",i,"_CH.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df<-rbind(tmp_df,data1)
    }
    tmp_df$value<-tmp_df$value*100
    p <- ggplot() +
    geom_vline(xintercept = c(0,400,1400,1800), color = "grey", linetype = "solid", size = 0.2) + 
    geom_hline(yintercept = c(0,0.25,0.5,0.75,1), color = "grey", linetype = "solid", size = 0.2) +
    geom_line(data = tmp_df, mapping = aes(x=pos, y=value, color =type,group=type),alpha=1,linewidth=0.3) +
    labs(title = "hmCH") +
     xlab("") +
     ylab(paste0(class[k]))+   
    scale_color_manual(values = c("=0"="#0278ae","<0.7960165"="#ffd480",">0.7960165"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,0.5,1),limits=c(0,1))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(axis.title.y = element_blank())+
    other_theme
    plots[[(n-1)*6+6]]<-p
#}
}


pdf("line_diagrams_of_methylation_level_in_different_subclasses.pdf",height = 9,width = 30)
#print(final_plot)
plot_grid(plotlist = plots[1:18], ncol = 6,rel_heights = c(4.5,4,4.5),rel_widths =c(2.2,rep(2,5)))
dev.off()

