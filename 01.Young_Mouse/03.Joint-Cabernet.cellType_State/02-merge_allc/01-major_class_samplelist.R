library(stringr)

total_QC = read.csv("../../../03.data/02.metainfo/01.Young_Mouse/RNA_DNA_match_name_QC_class_label_young.csv",header=T)
for(major_class in unique(total_QC$class_label)){
    df = data.frame(allc_path = paste0("allc_",total_QC[total_QC$RNA_QC ==1 & total_QC$hmC_QC == 1 & total_QC$class_label == major_class,"hmC_SampleID"],".mm10.dna.tsv.gz"))
    write.table(df,file = paste0("01.major_class_samplelist/",str_replace_all(str_replace_all(major_class," ","_"),"/","_"),".hmC.samplelist.txt"),quote=F,row.names=F,col.names=F,sep="\t")
    df = data.frame(allc_path = paste0("allc_",total_QC[total_QC$RNA_QC ==1 & total_QC$mC_QC == 1 & total_QC$class_label == major_class,"mC_SampleID"],".mm10.dna.tsv.gz"))
    write.table(df,file = paste0("01.major_class_samplelist/",str_replace_all(str_replace_all(major_class," ","_"),"/","_"),".mC.samplelist.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}

for(subclass in unique(total_QC$subclass_label)){
    df = data.frame(allc_path = paste0("allc_",total_QC[total_QC$RNA_QC ==1 & total_QC$hmC_QC == 1 & total_QC$subclass_label == subclass,"hmC_SampleID"],".mm10.dna.tsv.gz"))
    write.table(df,file = paste0("03.subclass_samplelist/",str_replace_all(str_replace_all(subclass," ","_"),"/","_"),".hmC.samplelist.txt"),quote=F,row.names=F,col.names=F,sep="\t")
    df = data.frame(allc_path = paste0("allc_",total_QC[total_QC$RNA_QC ==1 & total_QC$mC_QC == 1 & total_QC$subclass_label == subclass,"mC_SampleID"],".mm10.dna.tsv.gz"))
    write.table(df,file = paste0("03.subclass_samplelist/",str_replace_all(str_replace_all(subclass," ","_"),"/","_"),".mC.samplelist.txt"),quote=F,row.names=F,col.names=F,sep="\t")
}


