# script to process counts from megamate sequencing. Requires functions in MegaMatePlot.R
# modified 20160129 to make it generalizable for non RNAi phenotypes
# added modifications of different genomes
# modified 20160629 to add latest smoothing functions
#modified 20160630 to use separate G and smG functions
#modified 20160818 to remove uncalled SNPs and rename Chromosomes
#modified 20160906 to plot and store smoothed data in processed data.file
#modified 20170101 to deal with reads mapped to N2 and CB genomes
# modified 20180709 to combine files without doing stats calculation, and convert to genomic Ranges

library(signal)
library(RColorBrewer)
library(zoo)
library(dplyr)
options(scipen=10)
library(signal)

#setwd("/media/jenny/670FC52648FA85C4/Dropbox (CRG)/CRG_files/MegaMate/seqResults/2013-14merged/scripts")
setwd("~/Documents/LehnerLab/mergedData/scripts")
setwd("~/Documents/LehnerLab/mergedData_trim_slWin/scripts")
setwd("~/Documents/LehnerLab/mergedData_trim_unpaired/scripts")
setwd("~/Documents/LehnerLab/m201512and201602/scripts")
setwd("~/Documents/LehnerLab/m201512and201602_trim_slWin/scripts")
setwd("~/Documents/LehnerLab/m201512and201602_trim_unpaired/scripts")
setwd("~/Documents/LehnerLab/m201603and06/scripts")
setwd("~/Documents/LehnerLab/m201603and06_trim_slWin/scripts")
setwd("~/Documents/LehnerLab/m201603and06_trim_unpaired/scripts")

source("~/Documents/LehnerLab/BSApaper/MegaMatePlot.R")

genomeVer="PRJNA13758.WS250"
genomeVerCB="PRJNA275000.WS250"
subdir=""
fileNames<-list.files(paste0("../finalData/",genomeVer,subdir),pattern="vars.*txt")
#samples<-read.table(paste0("../finalData/",genomeVer,"/samples.txt"),header=TRUE,stringsAsFactors = FALSE)
samples<-read.table(paste0("./sampleList.txt"),header=TRUE,stringsAsFactors = FALSE)
codeField<-grep("mrg",as.vector(t(data.frame(strsplit(fileNames,"_")))))
codes<-as.vector(t(data.frame(strsplit(fileNames,"_"))))[codeField]
codes<-gsub("mrg","",codes)
codes<-gsub(".txt","",codes)
#sampleName<-as.vector(t(data.frame(strsplit(fileNames,"_"))))[codeField-1]
fileNames<-data.frame(fileName=fileNames,code=codes,stringsAsFactors = FALSE)
fileList<-left_join(fileNames,samples,by="code")
fileList<-arrange(fileList,phenotype,population,treatment)
i<-which(names(fileList)!="code")
fileList<-fileList[,c(names(fileList)[i],"code")]
write.table(fileList,paste0("../finalData/",genomeVer,"/fileList.txt"),quote=FALSE)

# read in list of experiments and reference genome index
fileList<-read.table(paste0("../finalData/",genomeVer,"/fileList.txt"),header=TRUE, colClasses="character")
genomeIdx<-read.table(file=paste0("../../GenomeBuilds/",genomeVer,"/c_elegans.",genomeVer,".genomic.fa.fai"))
#genomeIdxCB<-read.table(file=paste0("../../GenomeBuilds/",genomeVerCB,"/c_elegans.",genomeVerCB,".genomic.fa.fai"))
IDconversion<-read.csv(file=paste0("../../GenomeBuilds/",genomeVer,"/N2-CB_IDs.txt"),
                         stringsAsFactors=FALSE)
names(IDconversion)[c(3,4)]<-paste0("IDtbl_",names(IDconversion)[c(3,4)])

#load granges object with genomic position of genes in this genome build (object is called gff_a)
load(file=paste0("../../GenomeBuilds/",genomeVer,"/Granges_",genomeVer,".RData"))

phenotypes<-list(RNAi_mel26=c("zeel-1","ppw-1","mel-26","fem-1"),
                 RNAi_his55=c("zeel-1","ppw-1","fem-1","his-55"),
                 Longevity=c("zeel-1","npr-1","glb-5","exp-1","nath-10","rec-8","daf-2","daf-16","age-1","akt-1","daf-18","cki-1","unc-31"),
                 L1_survival=c("zeel-1","npr-1","glb-5","exp-1","nath-10","lin-4","daf-2","daf-16","age-1","akt-1","daf-18","cki-1","unc-31"))

# bundle repeatedly used data into list
#chrs=c("CHROMOSOME_I","CHROMOSOME_II","CHROMOSOME_III","CHROMOSOME_IV","CHROMOSOME_V","CHROMOSOME_X") #for chr names in data table
chrs=as.vector(genomeIdx$V1[1:6])
chrNum<-c("ChrI","ChrII","ChrIII","ChrIV","ChrV","ChrX") #for labels on plots
#locNames<-c()#locNames<-as.character(gff_a[gff_a$public %in% phenotypes$Longevity,]$public)
#locs=c()#locs<-c(2342232,4185076,9132545,50551562,56356262)
labelAdj<-c(0)#rep(0,length(locs))#c(0.9,0.2,0.1,NA,NA)
chrEnds<-c(0,cumsum(genomeIdx$V2))[1:(length(genomeIdx$V2))]
midpoints<-chrEnds[1:(length(chrEnds)-1)]+diff(chrEnds)/2
bundle<-list(chrs=chrs, chrNum=chrNum, locNames=c(), locs=c(), labelAdj=labelAdj, chrEnds=chrEnds,
             midpoints=midpoints, mainTitle="",legend="")


d<-format(Sys.time(), "%Y%m%d")

#create a column in which to write name of file with processed data
fileList["processedData"]<-NA

for (ff in seq(1,dim(fileList)[1],by=2)) {
  #ff=1
  # read in control and selected data for single experiment
  contFile=fileList[fileList$treatment=="control",1][(ff+1)/2]
  selectFile=fileList[fileList$treatment=="selected",1][(ff+1)/2]
  cont<-read.table(paste0("../finalData/",genomeVer,subdir,"/",contFile),header=TRUE,stringsAsFactors=F)
  select<-read.table(paste0("../finalData/",genomeVer,subdir,"/",selectFile),header=TRUE,stringsAsFactors=F)
  contCB<-read.table(paste0("../finalData/",genomeVerCB,subdir,"/",contFile),header=TRUE,stringsAsFactors=F)
  selectCB<-read.table(paste0("../finalData/",genomeVerCB,subdir,"/",selectFile),header=TRUE,stringsAsFactors=F)

  names(cont)[1]<-"N2id.WS250"
  #correct N2 and CB labelling of columns
  names(contCB)<-gsub("CB","HAWAII",names(contCB))
  names(contCB)<-gsub("N2","CB",names(contCB))
  names(contCB)<-gsub("HAWAII","N2",names(contCB))
  names(contCB)[1]<-"CBid"
  #merge tables
  newCont<-merge(IDconversion,cont,by=c("N2id.WS250"),all=TRUE, sort=FALSE)
  newCont<-merge(newCont,contCB,by=c("CBid"),all=TRUE,sort=FALSE)
  contSampleName<-gsub("_N2fwd.x","",names(newCont)[11])

  #do same for select
  names(select)[1]<-"N2id.WS250"
  #correct N2 and CB labelling of columns
  names(selectCB)<-gsub("CB","HAWAII",names(selectCB))
  names(selectCB)<-gsub("N2","CB",names(selectCB))
  names(selectCB)<-gsub("HAWAII","N2",names(selectCB))
  names(selectCB)[1]<-"CBid"
  #merge tables
  newSelect<-merge(IDconversion,select,by=c("N2id.WS250"),all=TRUE, sort=FALSE)
  newSelect<-merge(newSelect,selectCB,by=c("CBid"),all=TRUE,sort=FALSE)
  selectSampleName<-gsub("_N2fwd.x","",names(newSelect)[11])

  #compare readDepths to summed N2 and CB counts from their own alignments

  sum(newCont$N2var.x==newCont$IDtbl_N2var,na.rm=TRUE)
  # 326323
  sum(newCont$CBvar.x==newCont$IDtbl_CBvar,na.rm=TRUE)
  # 207375
  sum(newCont$CBvar.y==newCont$IDtbl_CBvar,na.rm=TRUE)
  # 324992
  sum(newCont$N2var.y==newCont$IDtbl_N2var,na.rm=TRUE)
  # 217194

  dim(cont)
  # 326905     21
  dim(contCB)
  # 324992     21


  newCont1<-newCont[,c(1,2,4,5,6,7,8,28,15,35,17,18)]
  # find count columns
  N2countsCol<-grep("_N2counts.x",names(newCont1))
  CBcountsCol<-grep("_CBcounts.y",names(newCont1))
  newCont1[,paste0(contSampleName,"_readDepth.x")]<-newCont1[,N2countsCol]+newCont1[,CBcountsCol]
  newCont1[,paste0(contSampleName,"_CBfreq")]<-newCont1[,CBcountsCol]/newCont1[,paste0(contSampleName,"_readDepth.x")]

  newSelect1<-newSelect[,c(1,2,4,5,6,7,8,28,15,35,17,18)]
  N2countsCol<-grep("_N2counts.x",names(newSelect1))
  CBcountsCol<-grep("_CBcounts.y",names(newSelect1))
  newSelect1[,paste0(selectSampleName,"_readDepth.x")]<-newSelect1[,N2countsCol]+newSelect1[,CBcountsCol]
  newSelect1[,paste0(selectSampleName,"_CBfreq")]<-newSelect1[,CBcountsCol]/newSelect1[,paste0(selectSampleName,"_readDepth.x")]


  names(newCont1)<-gsub(".x$","",names(newCont1))
  names(newCont1)<-gsub(".y$","",names(newCont1))
  cont<-newCont1
  # > names(newCont1)
  # [1] "CBid"                     "N2id.WS250"               "IDtbl_N2var"              "IDtbl_CBvar"
  # [5] "Chr"                      "Pos"                      "N2var"                    "CBvar"
  # [9] "X1Fmerged_4999_N2counts"  "X1Fmerged_4999_CBcounts"  "X1Fmerged_4999_readDepth" "X1Fmerged_4999_CBfreq"
  #

  names(newSelect1)<-gsub(".x$","",names(newSelect1))
  names(newSelect1)<-gsub(".y$","",names(newSelect1))
  select<-newSelect1
  # > names(newSelect1)
  # [1] "CBid"                    "N2id.WS250"              "IDtbl_N2var"             "IDtbl_CBvar"
  # [5] "Chr"                     "Pos"                     "N2var"                   "CBvar"
  # [9] "MEmerged_4998_N2counts"  "MEmerged_4998_CBcounts"  "MEmerged_4998_readDepth" "MEmerged_4998_CBfreq"

  # add cumulative genome coordinates to data table
  cont<-cont[complete.cases(cont),]
  #327028 goes to 324889
  cont<-cumPosition(cont,genomeIdx)
  select<-select[complete.cases(select),]
  #327028 goes to 324530
  select<-cumPosition(select,genomeIdx)

  #change name of chromosomes to "CHROMOSOME_I" format
  cont<-renameChr(cont,bundle$chrs,dataName=paste0(fileList[ff,"phenotype"],"_",fileList[ff,"population"]),genomeVer=genomeVer)
  select<-renameChr(select,bundle$chrs,dataName=paste0(fileList[ff,"phenotype"],"_",fileList[ff,"population"]),genomeVer=genomeVer)

  # clean data by removing data with less than minRead, and whose read numbers are extreme outliers
  cont<-cleanData(cont, paste0(fileList[ff,"phenotype"],"_",fileList[ff,"population"]), genomeVer=genomeVer,
                  MADs=2, minReads=10)
  #down to 297417
  select<-cleanData(select, paste0(fileList[ff,"phenotype"],"_",fileList[ff,"population"]), genomeVer=genomeVer,
                    MADs=2, minReads=10)
  #down to 291030
  # merge the two tables on common columns (chr, position, variant)
  names(cont)[2]<-"ID"
  names(select)[2]<-"ID"
  myExp<-merge(cont[,c(2,4,6:13)],select[,c(2,4,6:13)],by=c("ID","Chr","Pos","Position","N2var","CBvar"),sort=FALSE)
  # remove loci that are not polymorphic in both control and selected
  myExp<-removeNonPolymorphic(myExp, paste0(fileList[ff,"phenotype"],"_",fileList[ff,"population"]),genomeVer=genomeVer)
  # down to 279746

  # order table according to genome coordinates
  myExp<-myExp[order(myExp$Position),]


  ##########

  #output processed data as a table
  processedFileName<-paste0("../finalData/",genomeVer,"/processedData_",d,"_",fileList[ff,"phenotype"],
                           "_",fileList[ff,"population"],"_noStats.txt")
  write.table(myExp,file=processedFileName)

  #add filename of processedData to fileList table
  fileList[ff,"processedData"]<-processedFileName
  fileList[ff+1,"processedData"]<-processedFileName
}

#write new fileList_date.txt file to store processedFileName
write.table(fileList,file=paste0("../finalData/",genomeVer,"/fileList_",d,"_noStats.txt"))

####################

setwd("~/Documents/LehnerLab/BSApaper")
readRDS("./plots/bundledInfo4Plotting.RDS")
fileList<-list.files(path="..",pattern="processedData_.*_noStats.txt",recursive=T,full.names=T)
i<-c(rep(2:5,times=3),rep(6:9,times=3),rep(1,3))
fileList<-fileList[order(i)]

samples<-read.table("./sampleList.txt",header=T,stringsAsFactors=F)

if(!dir.exists("./plots/freqDiff")) {
  dir.create("./plots/freqDiff",recursive=T)
}

pdf(file="./plots/freqDiff/compareProcessing.pdf", paper="a4", height=11, width=8)
par(mfrow=c(3,1))
for (ff in fileList) {
  processing<-strsplit(ff,"/")[[1]][2]
  data<-read.delim(ff,sep=" ",stringsAsFactors=F)
  CBfreqCol<-grep("_CBfreq",names(data))[1]
  nameFields<-unlist(strsplit(names(data)[CBfreqCol],"_"))
  sampleCode<-gsub("mrg","",nameFields[grep("mrg",nameFields)])
  if (length(sampleCode)==0) {
      sampleCode<-nameFields[grep("[0-9A-Z]{4}",nameFields)]
  }
  bundle$mainTitle<-paste(processing,sampleCode)
  plotSmoothedDiff1(data,bundle,sm=101)
}
dev.off()



setwd("~/Documents/LehnerLab/BSApaper")

fileList<-list.files(path="..",pattern="processedData_.*_noStats.txt",recursive=T,full.names=T)
i<-c(rep(2:5,times=3),rep(6:9,times=3),rep(1,3))
fileList<-fileList[order(i)]

samples<-read.table("./sampleList.txt",header=T,stringsAsFactors=F)

if(!dir.exists("./plots/freq2sample")) {
  dir.create("./plots/freq2sample",recursive=T)
}

pdf(file="./plots/freq2sample/compareProcessing_freq.pdf", paper="a4", height=11, width=8)
par(mfrow=c(3,1))
for (ff in fileList) {
  processing<-strsplit(ff,"/")[[1]][2]
  data<-read.delim(ff,sep=" ",stringsAsFactors=F)
  CBfreqCol<-grep("_CBfreq",names(data))[1]
  nameFields<-unlist(strsplit(names(data)[CBfreqCol],"_"))
  sampleCode<-gsub("mrg","",nameFields[grep("mrg",nameFields)])
  if (length(sampleCode)==0) {
    sampleCode<-nameFields[grep("[0-9A-Z]{4}",nameFields)]
  }
  bundle$mainTitle<-paste(processing,sampleCode)
  plotFreq(data,bundle)
}
dev.off()

bundle$mainTitle<-""
saveRDS(bundle,"./plots/bundledInfo4Plotting.RDS")
