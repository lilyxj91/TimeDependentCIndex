####################################################################
##  Name:           DataReformat.R  
##
##  Description:    This code is to reformat the public available 
##                  data, "readmission" in R package "frailtypack", 
##                  for evaluating the time-dependent C-index using
##                  the developed R package for illustration. 
##
## Author:          Jian Wang
## Completion date: May 2023                                          
## Modified:         
## Date             Initials                  Change                           
####################################################################

data(readmission, package = "frailtypack")
readmission_hold<-readmission

## Create data format to be used in the R package
readmission$t.obs<-readmission$Num<-rep(NA,dim(readmission)[1])
ID<-unique(readmission$id)
index.rm<-c();k<-1
for (i in ID){
  tmp<-readmission[readmission$id %in% i,]
  readmission[readmission$id %in% i,'t.obs']<-max(tmp$t.stop)
  readmission[readmission$id %in% i,'Num']<-sum(tmp$event)
  if (dim(tmp)[1]>1){
    index.rm[k]<-rownames(tmp)[tmp$t.stop==max(tmp$t.stop)]
    k<-k+1
  }
}
readmission_hold1<-readmission
readmission<-readmission_hold1[!(rownames(readmission_hold) %in% index.rm),]
readmission$t<-readmission$t.stop

## Re-code the variables
## sex
readmission$sex_num<-rep(NA,dim(readmission)[1])
readmission$sex_num[readmission$sex %in% 'Male']<-1
readmission$sex_num[readmission$sex %in% 'Female']<-2
## chemo
readmission$chemo_num<-rep(NA,dim(readmission)[1])
readmission$chemo_num[readmission$chemo %in% 'NonTreated']<-1
readmission$chemo_num[readmission$chemo %in% 'Treated']<-2
## dukes
readmission$dukes_num<-rep(NA,dim(readmission)[1])
readmission$dukes_num[readmission$dukes %in% 'A-B']<-1
readmission$dukes_num[readmission$dukes %in% 'C']<-2
readmission$dukes_num[readmission$dukes %in% 'D']<-3
## change time from days to years
readmission<-readmission[,c('id','t','t.obs','Num','sex_num','chemo_num','dukes_num')]
readmission$indi<-rep(0,dim(readmission)[1])
readmission$t.obs<-readmission$t.obs/365.25
readmission$t<-readmission$t/365.25

## Split the data into two data sets: training and validating
## For the purpose of illustration, we fix the random number seed so the results can be reproduced
seeds<-c(659)
set.seed(seeds)
ID_train<-sample(ID,200)
ID_valid<-setdiff(ID,ID_train)
data_train<-readmission[readmission$id %in% ID_train,]
data_valid<-readmission[readmission$id %in% ID_valid,]

## Save the training and validating data sets
data_train<-data_train[,c("id","t","t.obs","indi","Num","sex_num","chemo_num","dukes_num")]
names(data_train)<-c("id","t","t.obs","indi","Num","X1","X2","X3")
data_valid<-data_valid[,c("id","t","t.obs","indi","Num","sex_num","chemo_num","dukes_num")]
names(data_valid)<-c("id","t","t.obs","indi","Num","X1","X2","X3")
write.csv(data_train,"TrainData.csv",row.names = F)
write.csv(data_valid,"ValidateData.csv",row.names = F)