library(data.table)

args = commandArgs(T)
data_f = args[1]
means_o = args[2]
sampleID = args[3]

data_head = read.table(data_f,header=F,sep='\t',skip=0,nrow = 1)
head_tmp = gsub("],sample_labels.*","",gsub(".*group_labels:","",data_head$V1))
group_labels = strsplit(gsub("],group_boundaries.*","",head_tmp),"[[,]")[[1]]
group_labels = group_labels[2:length(group_labels)]

sample_labels = strsplit(gsub(".*group_boundaries:","",head_tmp),"[[,]")[[1]]
sample_labels = sample_labels[2:length(sample_labels)]
sample_labels = as.numeric(sample_labels)

data_data = read.table(data_f,header=F,sep='\t',skip=1)

data_matrix = data_data[,7:dim(data_data)[2]]
data_means = colMeans(data_matrix,na.rm=T)
data_means = as.data.frame(data_means)
colnames(data_means) = sampleID


if(length(group_labels)!=1){
  for(ss in seq(1,length(group_labels))){
    print(ss)
    sub_data = data_data[sample_labels[ss]:sample_labels[ss+1],7:dim(data_data)[2]]
    sub_means = colMeans(sub_data,na.rm=T)
    sub_means = as.data.frame(sub_means)
    colnames(sub_means) = sampleID

    write.table(sub_means, paste0(gsub("csv","",means_o),group_labels[ss],".csv"), col.names = T, row.names=F,sep='\t',append = F,quote = F)
  }
}


write.table(data_means, means_o, col.names = T, row.names=F,sep='\t',append = F,quote = F)

