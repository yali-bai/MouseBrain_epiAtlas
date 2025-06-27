library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
library(patchwork)
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



#######
plot_RNA_DNA_correlation <- function(datatype,mc_type , region_type,class) {
  old_RNA <- read.csv(paste0("./old_corrected/02data_expr/",class,"/RNA_",class,"_expr.csv"))
  colnames(old_RNA) <- c("gene","class","RNA")
  old_DNA <- read.csv(paste0("./old_corrected/02data_expr/",class,"/",datatype,"_",mc_type,"_",region_type,"_",class,"_expr.csv"))
  colnames(old_DNA) <- c("gene","class","DNA")
  old_DNA$gene<-gsub("\\..*","",old_DNA$gene)
  old_RNA <- subset(old_RNA, RNA != 0)
  old_RNA_DNA_df <- merge(old_RNA, old_DNA, by = c("gene", "class"))
  old_RNA_DNA_df <- old_RNA_DNA_df[is.finite(old_RNA_DNA_df$RNA) & is.finite(old_RNA_DNA_df$DNA), ]
  old_RNA_DNA_df$group<-"old"

  young_RNA <- read.csv(paste0("./young_corrected/02data_expr/",class,"/RNA_",class,"_expr.csv"))
  colnames(young_RNA) <- c("gene","class","RNA")
  young_DNA <- read.csv(paste0("./young_corrected/02data_expr/",class,"/",datatype,"_",mc_type,"_",region_type,"_",class,"_expr.csv"))
  colnames(young_DNA) <- c("gene","class","DNA")
  young_DNA$gene<-gsub("\\..*","",young_DNA$gene)
  young_RNA <- subset(young_RNA, RNA != 0)
  young_RNA_DNA_df <- merge(young_RNA, young_DNA, by = c("gene", "class"))
  young_RNA_DNA_df <- young_RNA_DNA_df[is.finite(young_RNA_DNA_df$RNA) & is.finite(young_RNA_DNA_df$DNA), ]
  young_RNA_DNA_df$group<-"young"

  combined<-rbind(old_RNA_DNA_df,young_RNA_DNA_df)
  if(datatype=="true_5mC"){dataty="5mC"
    dtype<-paste0(dataty,substr(mc_type,2,2))}
  if(datatype=="5hmC"){dataty="5hmC"
    dtype<-paste0(dataty,substr(mc_type,2,2))}
  if(datatype=="5mC"){dtype<-paste0("5mC",substr(mc_type,2,2),"+5hmC",substr(mc_type,2,2))}

  if(mc_type=="CG"){
    limity=1
  }else{limity=0.1}

  plots_list <- list()
  for (i in seq_along(subclass_order)) {
    subclass_i <- subclass_order[i]
    df <- combined[combined$class == subclass_i, ]
    set.seed(123)
    index<-sample(nrow(df))
    df<-df[index,]
    
    p <- ggplot(df, aes(x = RNA, y = DNA,color=group,fill=group)) +
      geom_point( alpha = 0.2, size = 1.2) +
      geom_smooth(method = "lm", formula = y ~ x, linetype = 2, alpha = 0.3,se = TRUE ) +    
      scale_color_manual(values = c('young'='#3498db',"old"='#e74c3c'))+
      scale_fill_manual(values = c('young' = '#3498db', 'old' = '#e74c3c'),guide = "none") +  
      scale_alpha_continuous(range = c(0.1, 0.4), guide = "none") + 
      scale_x_continuous(trans = "log1p", limits = c(0, 5.1)) +
      scale_y_continuous(limits = c(0, limity)) +
      theme_bw() +
      theme(
          text = element_text(size = 12),
          legend.position = "top",
          plot.title = element_text(hjust = 0.5, family = "ArialMT", size = 14),
          title = element_text(size = 12, hjust = 0.5, family = "ArialMT"),
          axis.title.x = element_text(color = "black", size = 12, family = "ArialMT"),
          axis.title.y = element_text(color = "black", size = 12, family = "ArialMT"),
          axis.text.x = element_text(color = "black", size = 10, family = "ArialMT"),
          axis.text.y = element_text(color = "black", size = 10, family = "ArialMT"),
          axis.line = element_line(size = 0.6),
          axis.ticks = element_line(size = 0.6),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.grid.minor = element_blank()
      ) +
      labs(title = subclass_i, x = "log1p(RNA)", y = paste0(dtype,"\n" ," (", region_type, ")"))
    plots_list[[i]] <- p
  }
  pic <- wrap_plots(plots_list, ncol = 6) + plot_layout(guides = "collect") & theme(plot.margin = margin(5, 5, 5, 5))& theme(legend.position = "bottom")#top图例会显示在第一行title下面 
  # ggsave(paste0("./old_young_dotplot_corrected/combined_adjust/random/",dtype,"_",region_type,"_RNA_",class,"_old_and_young_correlation_pointdensity.pdf"), plot = pic, width = 24, height = 20, units = "in", dpi = 300)
  # message("✅ 图已保存：", paste0("./old_young_dotplot_corrected/combined_adjust/random/",dtype,"_",region_type,"_RNA_",class,"_old_and_young_correlation_pointdensity.pdf"))
  ggsave(paste0("./old_young_dotplot_corrected/combined_adjust/random/",dtype,"_",region_type,"_RNA_",class,"_old_and_young_correlation_pointdensity.png"), plot = pic, width = 24, height = 20, units = "in", dpi = 300)
  message("✅ 图已保存：", paste0("./old_young_dotplot_corrected/combined_adjust/random/",dtype,"_",region_type,"_RNA_",class,"_old_and_young_correlation_pointdensity.png"))
}

datatypes<-c("5hmC","5mC","true_5mC")#
mc_types<-c("CG","CH")#
regions<-c("genebody")#,"promoter"

for(datatype in datatypes){
        for(mc_type in mc_types){
            for(region in regions){
                plot_RNA_DNA_correlation(datatype,mc_type,region,"subclass")
            }
        }
}