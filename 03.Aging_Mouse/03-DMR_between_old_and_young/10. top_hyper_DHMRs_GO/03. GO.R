library(biomaRt)
library(cowplot)
library(reshape2)
library(stringr)
library(clusterProfiler)
library(org.Mm.eg.db)

for(dir in c("three_class")){  #"subclass",
  for(subdir in c("top1000","top500")){ #,"top500","top1000","all","top100"
    files = list.files(path = paste0(dir,"/",subdir),pattern = "hyper_DHMRs_per_subclass\\.bed\\.intersect_gene\\.bed$")
    for(file in files){
      rt=read.table(paste0(dir,"/",subdir,"/",file),sep="\t",header=F)       #读取id.txt文件
      genes=as.vector(unlist(lapply(rt[,7], function(x) strsplit(x,"\\.")[[1]][1])))
      mm_symbols <- biomaRt::select(
        biomaRt::useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl"),
        keys = genes,
        columns = c("ensembl_gene_id", "entrezgene_id","go_id"),
        keytype = "ensembl_gene_id"
      )
      if(all(is.na(mm_symbols$entrezgene_id))){
        next;
      }
      # GO - Biological Process
      gene.GO.BP <- enrichGO(gene = unique(mm_symbols$entrezgene_id),
                          OrgDb = org.Mm.eg.db,
                          keyType = "ENTREZID",
                          ont = "BP",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 1, #0.01,
                          qvalueCutoff = 1, #0.05,
                          readable = T)
      
      # GO - Molecular Function
      gene.GO.MF <- enrichGO(gene = unique(mm_symbols$entrezgene_id),
                             OrgDb = org.Mm.eg.db,
                             keyType = "ENTREZID",
                             ont = "MF",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 1, #0.01,
                             qvalueCutoff = 1, #0.05,
                             readable = T)
      # GO - Cell Components
      gene.GO.CC <- enrichGO(gene = unique(mm_symbols$entrezgene_id),
                             OrgDb = org.Mm.eg.db,
                             keyType = "ENTREZID",
                             ont = "CC",
                             pAdjustMethod = "BH",
                             pvalueCutoff = 1, #0.01,
                             qvalueCutoff = 1, #0.05,
                             readable = T)
      
      ## plot
      if(is.null(gene.GO.BP) || is.null(gene.GO.MF) || is.null(gene.GO.CC)){
        next;
      }
      p1 <- barplot(gene.GO.BP, showCategory = 10)
      p2 <- dotplot(gene.GO.BP, showCategory = 10)
      p3 <- barplot(gene.GO.MF, showCategory = 10)
      p4 <- dotplot(gene.GO.MF, showCategory = 10)
      p5 <- barplot(gene.GO.CC, showCategory = 10)
      p6 <- dotplot(gene.GO.CC, showCategory = 10)
      pdf(paste0(dir,"/",subdir,"/",file,".GO.pdf"),height = 18,width= 15)
      p_grid <- plot_grid(p1,p2,p3,p4,p5,p6,ncol = 2,labels = c("GO_BP_barplot", "GO_BP_dotplot", "GO_MF_barplot", "GO_MF_dotplot","GO_CC_barplot", "GO_CC_dotplot"))
      print(p_grid)
      dev.off()
      
    }
  }
}



