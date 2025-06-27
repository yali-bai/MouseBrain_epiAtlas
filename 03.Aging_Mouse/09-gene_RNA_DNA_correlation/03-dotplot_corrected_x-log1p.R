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
setwd(paste0(indir,"/our_old+our_young"))

plot_RNA_DNA_correlation <- function(age,datatype,mc_type , region_type,class) {
  our_RNA <- read.csv(paste0("./",age,"_corrected/02data_expr/",class,"/RNA_",class,"_expr.csv"))
  colnames(our_RNA) <- c("gene","class","RNA")
  our_DNA <- read.csv(paste0("./",age,"_corrected/02data_expr/",class,"/",datatype,"_",mc_type,"_",region_type,"_",class,"_expr.csv"))
  colnames(our_DNA) <- c("gene","class","DNA")
  our_DNA$gene<-gsub("\\..*","",our_DNA$gene)
  our_RNA <- subset(our_RNA, RNA != 0)
  RNA_DNA_df <- merge(our_RNA, our_DNA, by = c("gene", "class"))
  RNA_DNA_df <- RNA_DNA_df[is.finite(RNA_DNA_df$RNA) & is.finite(RNA_DNA_df$DNA), ]

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
    df <- RNA_DNA_df[RNA_DNA_df$class == subclass_i, ]
    
    p <- ggplot(df, aes(x = RNA, y = DNA)) +
      geom_point(color = "lightgrey", alpha = 0.2, size = 1.2) +
      stat_density_2d(aes(fill = ..level..), geom = "polygon", alpha = 0.3, color = "grey") +
      geom_smooth(method = "lm", formula = y ~ x, linetype = 2, color = "#DC3F4E") + 
      scale_fill_gradientn(colors = c("grey"), limits = c(0, 200), oob = scales::squish) +
      stat_cor(method = 'pearson', color = "#DC3F4E", size = 4) +
      scale_x_continuous(trans = "log1p", limits = c(0, 5.1)) +
      scale_y_continuous(limits = c(0, limity)) +
      theme_bw() +
      theme(
          text = element_text(size = 12),
          legend.position = "none",
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

  pic <- wrap_plots(plots_list, ncol = 6) + plot_layout(guides = "collect") & theme(plot.margin = margin(5, 5, 5, 5))

  ggsave(paste0("./",age,"_corrected/03dotplot/",class,"/",dtype,"_",region_type,"_RNA_",class,"_correlation_pointdensity.pdf"), plot = pic, width = 24, height = 20, units = "in", dpi = 300)
  message("✅ 图已保存：", paste0("./",age,"_corrected/dotplot/",class,"/",dtype,"_",region_type,"_RNA_",class,"_correlation_pointdensity.pdf"))
  ggsave(paste0("./",age,"_corrected/03dotplot/",class,"/",dtype,"_",region_type,"_RNA_",class,"_correlation_pointdensity.png"), plot = pic, width = 24, height = 20, units = "in", dpi = 300)
  message("✅ 图已保存：", paste0("./",age,"_corrected/dotplot/",class,"/",dtype,"_",region_type,"_RNA_",class,"_correlation_pointdensity.png"))
}

ages<-c("old","young")
datatypes<-c("5hmC","5mC","true_5mC")
mc_types<-c("CG","CH")
regions<-c("genebody","promoter")
for(age in ages){
    for(datatype in datatypes){
        for(mc_type in mc_types){
            for(region in regions){
                plot_RNA_DNA_correlation(age,datatype,mc_type,region,"subclass")
            }
        }
    }
}









