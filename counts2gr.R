# 2018-07-10
# convert processedData____noStats.txt to genomic Ranges

setwd("~/Documents/LehnerLab/BSApaper")

fileList<-list.files(path="..",pattern="processedData_.*_noStats.txt",recursive=T,full.names=T)
#remove files where trimmomatic was used to process reads, as it does not chnage the CBfreq difference
i<-grep("_trim_slWin",fileList)
fileList<-fileList[-i]
i<-grep("_trim_unpaired",fileList)
fileList<-fileList[-i]


samples<-read.table("./sampleList.txt",header=T,stringsAsFactors=F)

if (!dir.exists("rds")) {
  dir.create("./rds")
}


for (ff in fileList) {
  processing<-strsplit(ff,"/")[[1]][2]
  data<-read.delim(ff,sep=" ",stringsAsFactors=F)
  CBfreqCol<-grep("_CBfreq",names(data))[1]
  nameFields<-unlist(strsplit(names(data)[CBfreqCol],"_"))
  sampleCode<-gsub("mrg","",nameFields[grep("mrg",nameFields)])
  if (length(sampleCode)==0) {
    sampleCode<-nameFields[grep("[0-9A-Z]{4}",nameFields)]
  }
  gr<-GRanges(seqnames=data$Chr,ranges=IRanges(start=data$Pos, width=1),strand="*",DataFrame(score=data[,CBfreqCol]))
  saveRDS(gr,paste0("./rds/",sampleCode,"_",processing,".RDS"))
}


