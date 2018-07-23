## 2018-06-07
## script to process and look at BWA data sets for paper

folderList<-c("../mergedData", "../mergedData_trim_slWin", "../mergedData_trim_unpaired")

#read in SNP calls aligned to N2 genome
n2list<-list()
for (f in folderList) {
  n2list[f]<-list(list.files(paste0(f,"/finalData/PRJNA13758.WS250"), pattern="vars_.*", full.names =T))
}
n2list<-unlist(n2list)

#read in SNP calls aligned to CB genome
cblist<-list()
for (f in folderList) {
  cblist[f]<-list(list.files(paste0(f,"/finalData/PRJNA275000.WS250"), pattern="vars_.*", full.names =T))
}
cblist<-unlist(cblist)


for (i in 1:length(n2list)) {
  i=1
  n2data<-read.delim(n2list[i])
  cbdata<-read.delim(cblist[i])
}


data<-read.delim(n2list[1])

plot(1:dim(data)[1],data$X1Fmerged_4999_N2counts,type="l",col="red")
lines(1:dim(data)[1],data$X1Fmerged_4999_CBcounts,col="blue")

quantile(data$X1Fmerged_4999_CBcounts,probs=c(0.5,0.9,0.95,0.99,0.999,0.9999,1))

quantile(data$X1Fmerged_4999_N2counts,probs=c(0.5,0.9,0.95,0.99,0.999,0.9999,1))

cbdata<-read.delim(cblist[1])
plot(1:dim(data)[1],data$X1Fmerged_4999_N2counts,type="l",col="red")
lines(1:dim(cbdata)[1],cbdata$X1Fmerged_4999_N2counts,col="blue")

quantile(cbdata$X1Fmerged_4999_CBcounts,probs=c(0.5,0.9,0.95,0.99,0.999,0.9999,1))

quantile(cbdata$X1Fmerged_4999_N2counts,probs=c(0.5,0.9,0.95,0.99,0.999,0.9999,1))

data<-read.delim("../mergedData/finalData/PRJNA13758.WS250/processedData_20170417_RNAi_mel26_MM10_gen5.txt")
head(data)
