library(reticulate)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

sc <- import('scanpy')
scvi_h5ad <- sc$read_h5ad(sprintf('%s/WMB-10Xv3-CB-raw.h5ad',indir)) 
meta = scvi_h5ad$obs
head(meta);dim(meta)
meta$cell_new <- paste(meta$cell_barcode,meta$library_label,sep = "-")

h5ad_cell <- function(sample,name){
  h5ad <- sc$read_h5ad(sprintf('%s/%s.h5ad',indir,sample)) 
  meta = h5ad$obs
  meta$cell_new <- paste(meta$cell_barcode,meta$library_label,sep = "-")
  write.table(meta,sprintf("%s/%s_meta.txt",outdir,name),row.names = T,col.names = T,quote = F,sep = "\t") 
}

h5ad_cell("WMB-10Xv3-CB-raw","CB")
h5ad_cell("WMB-10Xv3-CTXsp-raw","CTXsp")
h5ad_cell("WMB-10Xv3-HPF-raw","HPF")
h5ad_cell("WMB-10Xv3-HY-raw","HY")
h5ad_cell("WMB-10Xv3-Isocortex-1-raw","Isocortex-1")
h5ad_cell("WMB-10Xv3-Isocortex-2-raw","Isocortex-2")
h5ad_cell("WMB-10Xv3-MB-raw","MB")
h5ad_cell("WMB-10Xv3-MY-raw","MY")
h5ad_cell("WMB-10Xv3-OLF-raw","OLF")
h5ad_cell("WMB-10Xv3-PAL-raw","PAL")
h5ad_cell("WMB-10Xv3-P-raw","P")
h5ad_cell("WMB-10Xv3-STR-raw","STR")
h5ad_cell("WMB-10Xv3-TH-raw","TH")

