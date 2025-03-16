# Flip the negative strand over
library(data.table)

# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""
class<-read.table("../../../input/01-youth/subclass_order_for_integration_with_zeng.txt",sep='\n')
class<-gsub("/",".",class[,1])
class<-gsub(" ",".",class)
class<-gsub("-",".",class)
datatype<-c("5hmc","5mc")
mc_type<-c("CG","CH")
for(dt in datatype){
    for(mc in mc_type){
        file_path <- paste0(indir,"/",dt,"/",mc)
        for(k in 1:length(class)){
            cat(k,'\n')
            for(i in c(1:3)){
                data<-read.table(paste0(file_path,"/",dt,"_",class[k],"_group",i,"_",mc,".genebody.gz"),header=F,sep='\t',skip=1)%>%as.data.frame()
                data$V4<-gsub("_.*","",data$V4)
                data[data$V4=="-",7:ncol(data)]<-rev(data[data$V4=="-",7:ncol(data)])
                data_matrix = data[,7:dim(data)[2]]
                data_means = colMeans(data_matrix,na.rm=T)
                data_means = as.data.frame(data_means)
                colnames(data_means) = paste0(dt,"_group",i,"_",mc)
                write.table(data_means, paste0(outdir,"/",dt,"_",class[k],"_group",i,"_",mc,".genebody.tmp"), col.names = T, row.names=F,sep='\t',append = F,quote = F)
            }
        }
    }
}



