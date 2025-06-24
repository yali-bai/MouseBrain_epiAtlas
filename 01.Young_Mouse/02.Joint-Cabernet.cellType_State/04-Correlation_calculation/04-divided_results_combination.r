# "indir" is a custom input path, and "outdir" is a custom output path.
# indir=""
# outdir=""

###########  combination of correlation results 
#Combine the calculated all cell correlation and 100 times shuffled correlation for each type of data
datatype<-c("5hmC","5mC","true_5mC")
mc_type<-c("CH","CG")
var="genebody"
for(mc in mc_type){
    for(dt in datatype){
        Joint_Cabernet_correlation<-read.csv(paste0("../../../output/01.Young_Mouse/02-correlation_calculation/all_cell_correlation/all_cells_",dt,"_",mc,"_",var,"_gene_correlation_results.csv"))
        colnames(Joint_Cabernet_correlation)[1]<-"gene_id"
        data<-Joint_Cabernet_correlation
        for(i in 1:100){
            cat(paste0(dt,"-",mc,"-",i),'\n')
            shuffle_data<-read.csv(paste0(indir,"shuffled_",dt,"_",mc,"_",var,"_gene_correlation_results(",i,").csv"))  
            colnames(shuffle_data)[1]<-"gene_id"     
            data<-merge(data,shuffle_data,by="gene_id") 
        }
        data$new.P.value<-apply(data,1,function(x){sum(abs(as.numeric(x[9:5008]))>=abs(as.numeric(x[6])))/5000})
        data$new.P.adjust<-p.adjust(data$new.P.value, method = "BH")
        write.csv(data,paste0("../../../output/01.Young_Mouse/02-correlation_calculation/shuffled_correlation/shuffled_",dt,"_",mc,"_",var,"_gene_correlation_total_results.csv"),row.names=F)
    }
}












