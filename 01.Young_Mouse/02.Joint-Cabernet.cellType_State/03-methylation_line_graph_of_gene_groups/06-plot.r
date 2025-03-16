#Draw line diagrams of methylation level in different subclasses.
library(dplyr)
library(readxl)
library(Seurat)
library(cowplot)
library(ggplot2)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

###########  hmc_CG  ##################

file_path <- paste0(indir,"/5hmc/CG")
#Extract all file names in the folder  
all_files <- list.files(file_path, full.names = FALSE)  
  
#Filter out file names ending in.txt  
tmp_files <- all_files[grepl("\\.tmp$", all_files, ignore.case = TRUE)] 
tmp<-strsplit(tmp_files,"_")
matrix_result <- do.call(rbind, lapply(tmp, function(x) as.character(x)))%>%as.data.frame()

class=matrix_result[!duplicated(matrix_result[,2]),2]
group<-c("=0","<0.44",">0.44")


plots<-list()
for(k in 1:length(class)){
    cat(k,'\n')
    tmp_df<-data.frame()
    for(i in c(1:3)){
        data<-read.csv(paste0(file_path,"/5hmc_",class[k],"_group",i,"_CG.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df<-rbind(tmp_df,data1)
    }
    tmp_df$value<-tmp_df$value*100
    p <- ggplot() +
    geom_vline(xintercept = c(0,400,1400,1800), color = "grey", linetype = "solid", size = 0.2) + #Add vertical guides
    geom_hline(yintercept = c(0,15,30,45), color = "grey", linetype = "solid", size = 0.2) + # Add horizontal guides 
    geom_line(data = tmp_df, mapping = aes(x=pos, y=value, color =type,group=type),alpha=1,linewidth=0.3) +
    labs(title = paste0(class[k])) +
     xlab("") +
     ylab("hmCG")+
    scale_color_manual(values = c("=0"="#0278ae","<0.44"="#ffd480",">0.44"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,15,30,45),limits=c(0,45))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,family="ArialMT"),
        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
        axis.title.x = element_text(color="black", size=13,family="ArialMT"),
        axis.title.y = element_text(color="black", size=13,family="ArialMT"),
        axis.text.x = element_text(color="black", size=10,family="ArialMT"),
        axis.text.y = element_text(color="black", size=10,family="ArialMT"),
        axis.line=element_line(size=0.6),
        axis.ticks=element_line(size=0.6),
        panel.border=element_blank(),
        panel.grid=element_blank())     
   plots[[k]]<-p  
}


###########  hmc_CH  ##################

file_path <- paste0(indir,"/5hmc/CH")

for(k in 1:length(class)){
    cat(k,'\n')
    tmp_df<-data.frame()
    for(i in c(1:3)){
        data<-read.csv(paste0(file_path,"/5hmc_",class[k],"_group",i,"_CH.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df<-rbind(tmp_df,data1)
    }
    tmp_df$value<-tmp_df$value*100
    p <- ggplot() +
    geom_vline(xintercept = c(0,400,1400,1800), color = "grey", linetype = "solid", size = 0.2) + 
    geom_hline(yintercept = c(0,0.25,0.5,0.75,1), color = "grey", linetype = "solid", size = 0.2) +
    geom_line(data = tmp_df, mapping = aes(x=pos, y=value, color =type,group=type),alpha=1,linewidth=0.3) +
    labs(title = paste0(class[k])) +
     xlab("") +
     ylab("hmCH")+   
    scale_color_manual(values = c("=0"="#0278ae","<0.44"="#ffd480",">0.44"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1),limits=c(0,1))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,family="ArialMT"),
        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
        axis.title.x = element_text(color="black", size=13,family="ArialMT"),
        axis.title.y = element_text(color="black", size=13,family="ArialMT"),
        axis.text.x = element_text(color="black", size=10,family="ArialMT"),
        axis.text.y = element_text(color="black", size=10,family="ArialMT"),
        axis.line=element_line(size=0.6),
        axis.ticks=element_line(size=0.6),
        panel.border=element_blank(),
        panel.grid=element_blank())
   plots[[k+30]]<-p
}




###########  mc_CG  ##################

file_path <- paste0(indir,"/5mc/CG")

for(k in 1:length(class)){
    cat(k,'\n')
    tmp_df<-data.frame()
    for(i in c(1:3)){
        data<-read.csv(paste0(file_path,"/5mc_",class[k],"_group",i,"_CG.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df<-rbind(tmp_df,data1)
    }
    tmp_df$value<-tmp_df$value*100
    p <- ggplot() +
    geom_vline(xintercept = c(0,400,1400,1800), color = "grey", linetype = "solid", size = 0.2) +
    geom_hline(yintercept = c(0,20,40,60,80), color = "grey", linetype = "solid", size = 0.2) + 
    geom_line(data = tmp_df, mapping = aes(x=pos, y=value, color =type,group=type),alpha=1,linewidth=0.3) +
    labs(title = paste0(class[k])) +
     xlab("") +
     ylab("mCG+hmCG")+ 
    scale_color_manual(values = c("=0"="#0278ae","<0.44"="#ffd480",">0.44"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,20,40,60,80),limits=c(0,90))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,family="ArialMT"),
        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
        axis.title.x = element_text(color="black", size=13,family="ArialMT"),
        axis.title.y = element_text(color="black", size=13,family="ArialMT"),
        axis.text.x = element_text(color="black", size=10,family="ArialMT"),
        axis.text.y = element_text(color="black", size=10,family="ArialMT"),
        axis.line=element_line(size=0.6),
        axis.ticks=element_line(size=0.6),
        panel.border=element_blank(),
        panel.grid=element_blank())
   plots[[k+60]]<-p
}



###########  mc_CH  ##################

file_path <- paste0(indir,"/5mc/CH")

for(k in 1:length(class)){
    cat(k,'\n')
    tmp_df<-data.frame()
    for(i in c(1:3)){
        data<-read.csv(paste0(file_path,"/5mc_",class[k],"_group",i,"_CH.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df<-rbind(tmp_df,data1)
    }
    tmp_df$value<-tmp_df$value*100
    p <- ggplot() +
    geom_vline(xintercept = c(0,400,1400,1800), color = "grey", linetype = "solid", size = 0.2) +
    geom_hline(yintercept = c(0,2,4), color = "grey", linetype = "solid", size = 0.2) +
    geom_line(data = tmp_df, mapping = aes(x=pos, y=value, color =type,group=type),alpha=1,linewidth=0.3) +
    labs(title = paste0(class[k])) +
     xlab("") +
     ylab("mCH+hmCH")+ 
    scale_color_manual(values = c("=0"="#0278ae","<0.44"="#ffd480",">0.44"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,2,4),limits=c(0,4.1))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,family="ArialMT"),
        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
        axis.title.x = element_text(color="black", size=13,family="ArialMT"),
        axis.title.y = element_text(color="black", size=13,family="ArialMT"),
        axis.text.x = element_text(color="black", size=10,family="ArialMT"),
        axis.text.y = element_text(color="black", size=10,family="ArialMT"),
        axis.line=element_line(size=0.6),
        axis.ticks=element_line(size=0.6),
        panel.border=element_blank(),
        panel.grid=element_blank())
   plots[[k+90]]<-p    
}


############### hmcg/(mcg+hmcg)   ###############################


for(k in 1:length(class)){
    cat(k,'\n')
    tmp_df_hmcg<-data.frame()
    for(i in 1:3){
        data<-read.csv(paste0(indir,"/5hmc/CG/5hmc_",class[k],"_group",i,"_CG.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df_hmcg<-rbind(tmp_df_hmcg,data1)
    }
    tmp_df_mcg<-data.frame()
    for(i in 1:3){
        data<-read.csv(paste0(indir,"/5mc/CG/5mc_",class[k],"_group",i,"_CG.genebody.tmp"), sep = "\t")%>%as.data.frame()
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

    labs(title = paste0(class[k])) +
     xlab("") +
     ylab("hmCG/(mCG+hmCG)")+ 
    scale_color_manual(values = c("=0"="#0278ae","<0.44"="#ffd480",">0.44"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,20,40),limits=c(0,56))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,family="ArialMT"),
        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
        axis.title.x = element_text(color="black", size=10,family="ArialMT"),
        axis.title.y = element_text(color="black", size=10,family="ArialMT"),
        axis.text.x = element_text(color="black", size=10,family="ArialMT"),
        axis.text.y = element_text(color="black", size=10,family="ArialMT"),
        axis.line=element_line(size=0.6),
        axis.ticks=element_line(size=0.6),
        panel.border=element_blank(),
        panel.grid=element_blank())
   plots[[k+120]]<-p 
}



#####################   true5mc_CG    ################################


for(k in 1:length(class)){
    cat(k,'\n')
    tmp_df_hmcg<-data.frame()
    for(i in 1:3){
        data<-read.csv(paste0(indir,"/5hmc/CG/5hmc_",class[k],"_group",i,"_CG.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df_hmcg<-rbind(tmp_df_hmcg,data1)
    }
    tmp_df_mcg<-data.frame()
    for(i in 1:3){
        data<-read.csv(paste0(indir,"/5mc/CG/5mc_",class[k],"_group",i,"_CG.genebody.tmp"), sep = "\t")%>%as.data.frame()
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
    labs(title = paste0(class[k])) +
     xlab("") +
     ylab("mCG")+   
    scale_color_manual(values = c("=0"="#0278ae","<0.44"="#ffd480",">0.44"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,20,40,60,80),limits=c(0,80))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,family="ArialMT"),
        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
        axis.title.x = element_text(color="black", size=13,family="ArialMT"),
        axis.title.y = element_text(color="black", size=13,family="ArialMT"),
        axis.text.x = element_text(color="black", size=10,family="ArialMT"),
        axis.text.y = element_text(color="black", size=10,family="ArialMT"),
        axis.line=element_line(size=0.6),
        axis.ticks=element_line(size=0.6),
        panel.border=element_blank(),
        panel.grid=element_blank())
   plots[[k+150]]<-p
}



####################  true-5mCH   ##############################

for(k in 1:length(class)){
    cat(k,'\n')
    tmp_df_hmch<-data.frame()
    for(i in 1:3){
        data<-read.csv(paste0(indir,"/5hmc/CH/5hmc_",class[k],"_group",i,"_CH.genebody.tmp"), sep = "\t")%>%as.data.frame()
        data1 <- data %>% `colnames<-`(c("value"))%>% mutate(SeqType=paste0("group",i))  %>% mutate(pos=as.numeric(rownames(.)))%>%mutate(type=group[i])     
        tmp_df_hmch<-rbind(tmp_df_hmch,data1)
    }
   tmp_df_mch<-data.frame()
    for(i in 1:3){
        data<-read.csv(paste0(indir,"/5mc/CH/5mc_",class[k],"_group",i,"_CH.genebody.tmp"), sep = "\t")%>%as.data.frame()
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
    labs(title = paste0(class[k])) +
     xlab("") +
     ylab("mCH")+ 
    scale_color_manual(values = c("=0"="#0278ae","<0.44"="#ffd480",">0.44"="#e84a5f"))+    
    scale_x_continuous(breaks=c(0,400,1400,1800),
                     labels=c("-2kb","TSS","TES","+2kb"))+
    scale_y_continuous(breaks=c(0,2,4),limits=c(0,4))+
    guides(fill = FALSE,color=FALSE) + 
    theme_bw()+
    theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,family="ArialMT"),
        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
        axis.title.x = element_text(color="black", size=13,family="ArialMT"),
        axis.title.y = element_text(color="black", size=13,family="ArialMT"),
        axis.text.x = element_text(color="black", size=10,family="ArialMT"),
        axis.text.y = element_text(color="black", size=10,family="ArialMT"),
        axis.line=element_line(size=0.6),
        axis.ticks=element_line(size=0.6),
        panel.border=element_blank(),
        panel.grid=element_blank())
   plots[[k+180]]<-p

}



pdf(paste0("../../../output/01-youth/01-methylation_line_graph/hmCG_group3_turn.pdf"),height = 10,width = 16)
plot_grid(plotlist = plots[1:30], ncol = 6)
dev.off()
pdf(paste0("../../../output/01-youth/01-methylation_line_graph/hmCH_group3_turn.pdf"),height = 10,width = 16)
plot_grid(plotlist = plots[31:60], ncol = 6)
dev.off()
pdf(paste0("../../../output/01-youth/01-methylation_line_graph/mCG+hmCG_group3_turn.pdf"),height = 10,width = 16)
plot_grid(plotlist = plots[61:90], ncol = 6)
dev.off()
pdf(paste0("../../../output/01-youth/01-methylation_line_graph/mCH+hmCH_group3_turn.pdf"),height = 10,width = 16)
plot_grid(plotlist = plots[91:120], ncol = 6)
dev.off()
pdf(paste0("../../../output/01-youth/01-methylation_line_graph/hmCG_mCG+hmCG_group3_turn.pdf"),height = 10,width = 16)
plot_grid(plotlist = plots[121:150], ncol = 6)
dev.off()
pdf(paste0("../../../output/01-youth/01-methylation_line_graph/mCG_group3_turn.pdf"),height = 10,width = 16)
plot_grid(plotlist = plots[151:180], ncol = 6)
dev.off()
pdf(paste0("../../../output/01-youth/01-methylation_line_graph/mCH_group3_turn.pdf"),height = 10,width = 16)
plot_grid(plotlist = plots[181:210], ncol = 6)
dev.off()