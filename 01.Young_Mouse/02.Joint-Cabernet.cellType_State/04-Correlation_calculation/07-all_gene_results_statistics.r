#########    All "our" in the following code refers to Joint Cabernet.
####################   8.statistics
library(data.table)
library(dplyr)

datatype<-c("5mC","5hmC","true_5mC")
var_dim<-c("genebody")
mctype<-c("CH","CG")
data<-fread("../../../output/01-youth/02-correlation_calculation/shuffled_RNA_DNA_correlation_result_all.csv")%>%as.data.frame()
data_filter<-data[data$new.P.adjust<0.05,]

result_df <-data.frame()
for (i in 1:3){
  for(k in 1){
    for (j in 1:2) {
        type_data<-data_filter[which(data_filter$datatype==datatype[i]&data_filter$mc_type==mctype[j]&data_filter$var_dim==var_dim[k]),]
        my_vector <- type_data$Correlation
        tmp_df <- data.frame(
           idents=paste0("our_RNA_",datatype[i],"_",mctype[j],"_",var_dim[k]),
           filled_NA="NO",
           Datatype=datatype[i],
           mc_type=mctype[j],
           var_dim=var_dim[k],
           '<-0.5' = sum(my_vector < -0.5),
           '<-0.4' = sum(my_vector < -0.4),
           '<-0.3' = sum(my_vector < -0.3),
           '<-0.2' = sum(my_vector < -0.2),
           '<-0.1' = sum(my_vector < -0.1),
           '<-0.05' = sum(my_vector < -0.05),
           '<-0.01' = sum(my_vector < -0.01),
           '<0' = sum(my_vector < 0),
           '>0' = sum(my_vector > 0),
           '>0.01'  = sum(my_vector > 0.01),
           '>0.05' = sum(my_vector > 0.05),
           '>0.1' = sum(my_vector > 0.1),
           '>0.2' = sum(my_vector > 0.2),
           '>0.3' = sum(my_vector > 0.3),
           '>0.4' = sum(my_vector > 0.4),
           '>0.5' = sum(my_vector > 0.5),
           '25%' = quantile(my_vector, probs = 0.25),
           '75%' = quantile(my_vector, probs = 0.75),
           Median = quantile(my_vector, probs = 0.5),
           Mean = mean(my_vector)
          )
        result_df <- rbind(result_df, tmp_df)
    }
  }
 }

colnames(result_df)[6:23]<-c('<-0.5','<-0.4','<-0.3','<-0.2','<-0.1','<-0.05',
                             '<-0.01','<0','>0','>0.01','>0.05','>0.1','>0.2','>0.3',
                             '>0.4','>0.5','%25','%75')

write.csv(result_df,"../../../output/01-youth/02-correlation_calculation/all_combined_gene_correlation_result.csv",row.names=F)

