########################  1.density plot
library(ggplot2)
library(dplyr)
library(data.table)
library(cowplot)

data<-fread("../../../output/01-youth/02-correlation_calculation/shuffled_RNA_DNA_correlation_result_all.csv")%>%as.data.frame()
data_filter<-data[data$new.P.adjust<0.05,]

datatype<-c("5mC","true_5mC","5hmC")
mc_type<-c("CG","CH")
var_dim<-c("genebody")
plot<-list()
g=1
for(mc in mc_type){
    for(varim in var_dim){
        for(dt in datatype){
            cat(g,"\n")
            if(mc=="CG"){
                if(dt=="true_5mC"){
                    xlab="TSS to TES 5mCG"
                }
                if(dt=="5mC"){
                    xlab="TSS to TES 5mCG+5hmCG"

                }
                if(dt=="5hmC"){
                    xlab="TSS to TES 5hmCG"
                }
            }else{
                if(dt=="true_5mC"){
                    xlab="TSS to TES 5mCH"
                }
                if(dt=="5mC"){
                    xlab="TSS to TES 5mCH+5hmCH"
                }
                if(dt=="5hmC"){
                    xlab="TSS to TES 5hmCH"
                }
            }
            #all
            df<-data[data$datatype==dt&data$mc_type==mc&data$var_dim==varim,]
            df_long <- reshape2::melt(df[,c(12,17:5016)])
            df_long$group<-ifelse(df_long$variable == "Correlation", "Observed", "Shuffled")  
            quantiles <- quantile(df_long$value[df_long$group=="Shuffled"], probs = c(0.025, 0.975))
            p<-ggplot(df_long, aes(x = value, group = group, color = group,alpha=group)) +  
                geom_density(aes(alpha=group))+
                geom_area(stat = "density", aes(y = ..density..,alpha=group),fill="#8dc6ff" ,position = "identity") +  # Use geom_area to simulate the shadow
                scale_alpha_manual(values = c("Observed" = 0.4, "Shuffled" = 0.0000001)) +  # Set the transparency of each group separately 
                scale_color_manual(values = c("Observed" = "#3498db", "Shuffled" = "#97a0a6"))+
                scale_x_continuous(breaks=c(-0.1,as.numeric(format(quantiles[1],digits=1)),0,as.numeric(format(quantiles[2],digits=1)),0.1),limit=c(-0.1,0.1))+
                scale_y_continuous(breaks = seq(from = 0, to = 85, by = 20),  
                     limits = c(0, 85))+ 
                geom_vline(xintercept = quantiles, linetype = "dashed", color = "#c9d6df",size=0.5) +                
                labs(title=paste0("RNA expression and ",xlab),
                    x ="Pearson correlation" ,y="Density")+
                theme_bw()+
                theme(legend.position="top",
                plot.title = element_text(hjust = 0.5,family="ArialMT"),
                title = element_text(size=7,hjust = 0.5,family="ArialMT"),
                axis.title.x = element_text(color="black", size=9,family="ArialMT"),
                axis.title.y = element_text(color="black", size=9,family="ArialMT"),
                axis.text.x = element_text(color="black", size=9,family="ArialMT", hjust = 1,),
                axis.text.y = element_text(color="black", size=9,family="ArialMT"),
                axis.line=element_line(size=0.6),
                axis.ticks=element_line(size=0.6),
                panel.border=element_blank(),
                panel.grid=element_blank())
            plot[[g]]<-p      
            g<-g+1
        }
    }
}

pdf("../../../output/01-youth/02-correlation_calculation/plot/Correlation_density_plot_limit_0.1.pdf",width=9,height=12)
plot_grid(plotlist = plot, ncol = 3)
dev.off()



#################   2.upset plot
library(UpSetR)
library(RColorBrewer)
library(ggplot2)
library(data.table)
library(dplyr)


data<-fread("../../../output/01-youth/02-correlation_calculation/shuffled_RNA_DNA_correlation_result_all.csv")%>%as.data.frame()
all<-data[data$new.P.adjust<0.05,]
all_data<-all[,1:16]
all_data$corr_direction<-ifelse(all_data$Correlation>0,"Positive correlation","Negative correlation")
mc_types<-c("CG","CH")
vardims <- c("genebody")
for (mc in mc_types){  
        for (vardim in vardims){
              #所有基因
              aa <- subset(all_data,mc_type==mc&var_dim==vardim)  
              set1 <- subset(aa,datatype == "true_5mC"&corr_direction=="Positive correlation")$gene_id  #mc
              set2 <- subset(aa,datatype == "true_5mC"&corr_direction=="Negative correlation")$gene_id  
              set3 <- subset(aa,datatype == "5hmC"&corr_direction=="Positive correlation")$gene_id  #hmc
              set4<- subset(aa,datatype == "5hmC"&corr_direction=="Negative correlation")$gene_id
              set5 <- subset(aa,datatype == "5mC"&corr_direction=="Positive correlation")$gene_id  #mc+hmc
              set6 <- subset(aa,datatype == "5mC"&corr_direction=="Negative correlation")$gene_id

              name<-c(paste0("[5m",mc," +]"),paste0("[5m",mc," -]"),paste0("[5hm",mc," +]"),paste0("[5hm",mc," -]"),paste0("[(5m",mc,"+5hm",mc,") +]"),paste0("[(5m",mc,"+5hm",mc,") -]"))
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
              
              # Genes with absolute value greater than 0.05
              bb<-subset(aa,abs(Correlation)>0.05)  
              set7 <- subset(bb,datatype == "true_5mC"&corr_direction=="Positive correlation")$gene_id  #mc
              set8 <- subset(bb,datatype == "true_5mC"&corr_direction=="Negative correlation")$gene_id  
              set9 <- subset(bb,datatype == "5hmC"&corr_direction=="Positive correlation")$gene_id  #hmc
              set10<- subset(bb,datatype == "5hmC"&corr_direction=="Negative correlation")$gene_id
              set11 <- subset(bb,datatype == "5mC"&corr_direction=="Positive correlation")$gene_id  #mc+hmc
              set12 <- subset(bb,datatype == "5mC"&corr_direction=="Negative correlation")$gene_id
              length<-c(length(set7),length(set8),length(set9),length(set10),length(set11),length(set12))
              list2 <- list(set1 = set7,set2 = set8,set3 = set9,set4 = set10,set5 = set11,set6 = set12)
              names(list2) <- name
              non_empty_dfs <- lapply(list2, function(x) { 
              if (length(x) > 0) {  
                  return(x)  
              }  })       
              # Remove NULL elements (empty data boxes are converted to NULL by the if statement in lapply)
              upset_list2 <- non_empty_dfs[!sapply(non_empty_dfs, is.null)]
              
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
                    #sets.y.label = "Counts",
                    sets.bar.color = color_vector1, 
                    mb.ratio = c(0.7, 0.3), 
                    order.by = "degree", 
                    text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.6), # 6 parameters intersection size title (y title size),intersection size tick labels (y scale label size), set size title (set title size), set size tick labels (set scale label size), set names (set classification label size), and set numbers above bars
                    shade.color = "#12507B" 
              )
              
              pdf(paste0("../../../output/01-youth/02-correlation_calculation/plot/Upset_",vardim,"_",mc,"_all_gene.pdf"),width =8, height =5)
              print(p1)
              dev.off()
              png(paste0("../../../outpu/01-youtht/02-correlation_calculation/plot/Upset_",vardim,"_",mc,"_all_gene.png"),width =8, height =5,units="in", res = 300)
              print(p1)
              dev.off()
            color_vector2 <- rev(color[length!=0])
              p2 <- upset(fromList(upset_list2),
                    nsets = length(upset_list2), 
                    nintersects = 80, 
                    sets = rev(c(name[length!=0])), 
                    keep.order = TRUE, 
                    number.angles = 0, 
                    point.size = 4, 
                    line.size = 1, #
                    mainbar.y.label = "Intersection size", 
                    main.bar.color = 'black',
                    matrix.color = "black", 
                    sets.x.label = "Set size", 
                    #sets.y.label = "Counts",
                    sets.bar.color = color_vector2, 
                    mb.ratio = c(0.7, 0.3),
                    order.by = "degree",
                    text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.6), 
                    shade.color = "#12507B"
              )
              pdf(paste0("../../../output/01-youth/02-correlation_calculation/plot/Upset_",vardim,"_",mc,"_filterd_gene.pdf"),width =8, height =5)
              print(p2)
              dev.off()
              png(paste0("../../../output/01-youth/02-correlation_calculation/plot/Upset_",vardim,"_",mc,"_filterd_gene.png"),width =8, height =5,units="in", res = 300)
              print(p2)
              dev.off()
          }}




