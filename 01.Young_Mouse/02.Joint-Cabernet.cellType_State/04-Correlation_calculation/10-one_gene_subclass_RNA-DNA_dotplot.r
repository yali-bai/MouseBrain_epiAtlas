#########    All "our" in the following code refers to Joint Cabernet.
# RNA expression and DNA methylation dotplot plots of gene AC132685.1 in different subclasses were mapped.
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)
library(ggpubr)
library(ggpointdensity) 
library(RColorBrewer)
library(cowplot)

# gene_name<-'AC132685.1'
# gene_list<-'ENSMUSG00000068151.7'

datatype<-c("5hmC","true_5mC","5mC")
mc_type<-c("CG")
var_dim<-c("genebody")
palette <- rev(brewer.pal(9,"Greys")[3:9]) #Choose the seven darkest colors
subclass<-read.table("../../../input/01-youth/subclass_order_for_integration_with_zeng.txt",sep='\n')
subclass<-subclass[,1]


plots<-list()
g=1
for(mc in mc_type){
    for(varim in var_dim){
        data<-read.csv("../../../output/01-youth/02-correlation_calculation/genebody_CG_ENSMUSG00000068151.7.csv")
        for(dt in datatype){
            if(dt=="true_5mC"){
                xlab="TSS to TES 5mCG"
            }
            if(dt=="5mC"){
                xlab="TSS to TES 5mCG+5hmCG"
            }
            if(dt=="5hmC"){
                xlab="TSS to TES 5hmCG"
            }
            aa<-data[data$datatype==dt,]
            for(sc in subclass){
                if(sc%in%aa$subclass){
                    bb<-aa[aa$subclass==sc,]
                    bb$log_DNA<-log2(bb$DNA+1)
                    p <- ggplot(bb, aes(x=log_DNA, y=RNA)) +
                        geom_pointdensity(adjust=0.1,show.legend = FALSE,size = 0.8)+
                        geom_smooth(method="lm", formula = y ~ x, linetype=2,color = "#e84545") + 
                        theme_bw() +
                        stat_cor( method = 'pearson',color = "#e84545",size=3) +
                        scale_colour_gradientn(colours = rev(palette))+
                        scale_x_continuous(breaks=c(0,0.25,0.5,0.75,1),limits=c(0,1))+
                        scale_y_continuous(limits=c(0,5))+
                        theme(legend.position="none",
                        plot.title = element_text(hjust = 0.5,family="ArialMT"),
                        title = element_text(size=10,hjust = 0.5,family="ArialMT"),
                        axis.title.x = element_text(color="black", size=10,family="ArialMT"),
                        axis.title.y = element_text(color="black", size=10,family="ArialMT"),
                        axis.text.x = element_text(color="black", size=8,family="ArialMT"),
                        axis.text.y = element_text(color="black", size=8,family="ArialMT"),
                        axis.line=element_line(size=0.6),
                        axis.ticks=element_line(size=0.6),
                        panel.border=element_blank(),
                        panel.grid=element_blank(),
                        element_line(linetype = "dashed"), 
                        panel.grid.minor = element_blank())+   # 去除次栅格线
                        labs(title=sc)+
                        xlab(xlab) +
                        ylab("RNA expression")
                        plots[[g]]<-p
                        g=g+1
            }}
        }
    }
}


pdf("../../../output/01-youth/02-correlation_calculation/plot/AC132685.1_subclass_plots.pdf",width=63,height=8)
plot_grid(plotlist=plots,ncol=25)
dev.off()
