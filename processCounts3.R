# script to process counts from megamate sequencing. Requires functions in MegaMatePlot.R
# modified 20160129 to make it generalizable for non RNAi phenotypes
# added modifications of different genomes
# modified 20160629 to add latest smoothing functions
#modified 20160630 to use separate G and smG functions
#modified 20160818 to remove uncalled SNPs and rename Chromosomes
#modified 20160906 to plot and store smoothed data in processed data.file
#modified 20170101 to deal with reads mapped to N2 and CB genomes

library(signal)
library(RColorBrewer)
library(zoo)
library(dplyr)
options(scipen=10)

#setwd("/media/jenny/670FC52648FA85C4/Dropbox (CRG)/CRG_files/MegaMate/seqResults/2013-14merged/scripts")
setwd("/home/jenny/Documents/LehnerLab/mergedData/scripts")

source("./MegaMatePlot.R")
genomeVer="PRJNA13758.WS250"
genomeVerCB="PRJNA275000.WS250"
subdir=""
fileNames<-list.files(paste0("../finalData/",genomeVer,subdir),pattern="vars.*txt")
#samples<-read.table(paste0("../finalData/",genomeVer,"/samples.txt"),header=TRUE,stringsAsFactors = FALSE)
samples<-read.table(paste0("./sampleList.txt"),header=TRUE,stringsAsFactors = FALSE)
codes<-as.vector(t(data.frame(strsplit(fileNames,"_"))[3,]))
codes<-gsub("mrg","",codes)
codes<-gsub(".txt","",codes)
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

### to find distance between two gr objects (my peaks and the possible loci)
# locOfInterest<-gff_a[gff_a$public %in% phenotypes$RNAi_mel26,]
# peak=GRanges(seqnames=c("CHROMOSOME_I"),
#            ranges=IRanges(start=c(5000000),end=c(5001000)),
#            strand=c("+"))
# locdistance<-distanceToNearest(peak,locOfInterest)
# mcols(locdistance)$distance
# locOfInterest[subjectHits(locdistance),]
#########

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

  # cont<-cbind(cont,IDconversion[which(IDconversion$N2id.WS250 %in% cont$ID),2:5])
  # if (sum(cont$ID!=cont$N2id.WS250)>0) {
  #   print("IDs do not match!!!")
  #   break
  # }
  # newCont<-merge(cont,contCB,by=c("CBid"),all=TRUE,sort=FALSE)
  names(cont)[1]<-"N2id.WS250"
  #correct N2 and CB labelling of columns
  names(contCB)<-gsub("CB","HAWAII",names(contCB))
  names(contCB)<-gsub("N2","CB",names(contCB))
  names(contCB)<-gsub("HAWAII","N2",names(contCB))
  names(contCB)[1]<-"CBid"
  #merge tables
  newCont<-merge(IDconversion,cont,by=c("N2id.WS250"),all=TRUE, sort=FALSE)
  newCont<-merge(newCont,contCB,by=c("CBid"),all=TRUE,sort=FALSE)

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

  #compare readDepths to summed N2 and CB counts from their own alignments
  if (!dir.exists("../plots")) {dir.create("../plots")}

  pdf(file="../plots/countsVreadDepths.pdf",paper="a4",height=11,width=8)
  par(mfrow=c(2,1))
  plot(newCont$X1Fmerged_4999_N2counts.x+newCont$X1Fmerged_4999_CBcounts.y,
       newCont$X1Fmerged_4999_readDepth.x+newCont$X1Fmerged_4999_readDepth.y,
       xlab="N2+CB counts",ylab="readDepthN2+readDepthCB", pch=16,cex=0.5,
       col="#00000066", main="Validity of alternative allele counts (control)")

  abline(0,1,col="red")
  abline(0,2,col="blue")
  legend("topleft",legend=c("y=x (alt allele always missed)",
                            "y=2x (alt allele count always right)"),col=c("red","blue"),lty=1)
  plot(newCont$X1Fmerged_4999_N2counts.x+newCont$X1Fmerged_4999_CBcounts.y,
       newCont$X1Fmerged_4999_readDepth.x+newCont$X1Fmerged_4999_readDepth.y,
       xlab="N2+CB counts",ylab="readDepthN2+readDepthCB", pch=16,cex=0.5,
       col="#00000033", main="Validity of alternative allele counts(control)",xlim=c(1,1000),
       ylim=c(1,2000))

  abline(0,1,col="red")
  abline(0,2,col="blue")
  legend("topleft",legend=c("y=x (alt allele always missed)",
                            "y=2x (alt allele count always right)"),col=c("red","blue"),lty=1)


  plot(newSelect$MEmerged_4998_N2counts.x+newSelect$MEmerged_4998_CBcounts.y,
       newSelect$MEmerged_4998_readDepth.x+newSelect$MEmerged_4998_readDepth.y,
       xlab="N2+CB counts",ylab="readDepthN2+readDepthCB", pch=16,cex=0.5,
       col="#00000066", main="Validity of alternative allele counts (mel-26)")

  abline(0,1,col="red")
  abline(0,2,col="blue")
  legend("topleft",legend=c("y=x (alt allele always missed)",
                            "y=2x (alt allele count always right)"),col=c("red","blue"),lty=1)
  plot(newSelect$MEmerged_4998_N2counts.x+newSelect$MEmerged_4998_CBcounts.y,
       newSelect$MEmerged_4998_readDepth.x+newSelect$MEmerged_4998_readDepth.y,
       xlab="N2+CB counts",ylab="readDepthN2+readDepthCB", pch=16,cex=0.5,
       col="#00000033", main="Validity of alternative allele counts(mel-26)",xlim=c(1,1000),
       ylim=c(1,2000))

  abline(0,1,col="red")
  abline(0,2,col="blue")
  legend("topleft",legend=c("y=x (alt allele always missed)",
                            "y=2x (alt allele count always right)"),col=c("red","blue"),lty=1)

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

  par(mfrow=c(2,1))
  plot(log10(newCont$X1Fmerged_4999_readDepth.x),
       log10(newCont$X1Fmerged_4999_readDepth.y),pch=16, col="#00000033",cex=0.5)
  abline(0,1)

  plot(log10(newSelect$MEmerged_4998_readDepth.x),
       log10(newSelect$MEmerged_4998_readDepth.y),pch=16, col="#00000033",cex=0.5)
  abline(0,1)

  dev.off()


  pdf(file="../plots/CBfreqWithDiffGenomes.pdf",paper="a4",height=11,width=8)
  par(mfrow=c(3,1))
  newCont1<-newCont[,c(1,2,4,5,6,7,8,28,15,35,17,18)]
  newCont1$X1Fmerged_4999_readDepth.x<-newCont1$X1Fmerged_4999_N2counts.x+newCont1$X1Fmerged_4999_CBcounts.y
  newCont1$X1Fmerged_4999_CBfreq<-newCont1$X1Fmerged_4999_CBcounts.y/newCont1$X1Fmerged_4999_readDepth.x
  plot(1:dim(newCont)[1],newCont$X1Fmerged_4999_CBfreq,pch=16,cex=0.5,
       col="#00000033",main="Mapping to N2 genome (control)",xlab="SNP",ylab="CB allele freq")
  plot(1:dim(newCont)[1],1-newCont$X1Fmerged_4999_N2freq,pch=16,cex=0.5,
       col="#00000033", main="Mapping to CB genome (control)",xlab="SNP",ylab="CB allele freq")
  plot(1:dim(newCont1)[1],newCont1$X1Fmerged_4999_CBfreq,pch=16,cex=0.5,
       col="#00000033", main="Mapping to both genomes (control)",xlab="SNP",ylab="CB allele freq")


  newSelect1<-newSelect[,c(1,2,4,5,6,7,8,28,15,35,17,18)]
  newSelect1$MEmerged_4998_readDepth.x<-newSelect1$MEmerged_4998_N2counts.x+newSelect1$MEmerged_4998_CBcounts.y
  newSelect1$MEmerged_4998_CBfreq<-newSelect1$MEmerged_4998_CBcounts.y/newSelect1$MEmerged_4998_readDepth.x
  plot(1:dim(newSelect)[1],newSelect$MEmerged_4998_CBfreq,pch=16,cex=0.5,
       col="#00000033",main="Mapping to N2 genome (mel-26)",xlab="SNP",ylab="CB allele freq")
  plot(1:dim(newSelect)[1],1-newSelect$MEmerged_4998_N2freq,pch=16,cex=0.5,
       col="#00000033", main="Mapping to CB genome (mel-26)",xlab="SNP",ylab="CB allele freq")
  plot(1:dim(newSelect1)[1],newSelect1$MEmerged_4998_CBfreq,pch=16,cex=0.5,
       col="#00000033", main="Mapping to both genomes (mel-26)",xlab="SNP",ylab="CB allele freq")

  dev.off()

  names(newCont1)<-gsub(".x$","",names(newCont1))
  names(newCont1)<-gsub(".y$","",names(newCont1))
  cont<-newCont1

  names(newSelect1)<-gsub(".x$","",names(newSelect1))
  names(newSelect1)<-gsub(".y$","",names(newSelect1))
  select<-newSelect1

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
  # remove SNPs where the alternative allele was not properly called
  #myExp<-removeUncalledSNPs(myExp, paste0(fileList[ff,"phenotype"],"_",fileList[ff,"population"]),genomeVer=genomeVer)
  # remove superfluous columns (separate counts of forward and reverse reads)
  #excl<-grep("fwd|rev",names(myExp))
  #myExp<-myExp[,-excl]
  # remove superfluous columns (mapping and base quality scores)
  #excl<-grep("QUAL|RPB|MQB|BQB|MQSB|MQ0F|MQ.|AN.",names(myExp))
  #myExp<-myExp[,-excl]
  # order table according to genome coordinates
  myExp<-myExp[order(myExp$Position),]

  # do fisher.test for each locus
  countCols<-grep("counts",names(myExp))
  myExp["fisherExact.pvals"]<-sapply(1:dim(myExp)[1],function(i){
    fisher.test(matrix(as.matrix(myExp[i,countCols]),ncol=2))$p.value})

  # adjust p values for multiple testing
  myExp["padj.fdr"]<-p.adjust(myExp[,"fisherExact.pvals"],method="fdr")
  # smoothe -log10(pvalues) using tricube kernel
  myExp["sm.logFisherPvals.tck1Mb"]<-smoothelogPbyChr(myExp,rep(1000000,dim(myExp)[1]),myExp$Position,"padj.fdr")


  # do GTest
  myExp["Gval.GTest"]<-getGbyChr(myExp,myExp$Position)
  #smoothe G values with tricube kernel
  myExp["sm.Gvals.tck1Mb"]<-smootheGbyChr(myExp,rep(1000000,dim(myExp)[1]),myExp$Position,"Gval.GTest")

  #do GTest with python function
  #First write a file with counts
  #write.table

  ########
  # open pdf for plots
  pdf(file=paste0("../finalData/",genomeVer,"/",d,"_",fileList[ff,"phenotype"],"_",
                  fileList[ff,"population"],".pdf"), paper="a4", height=11, width=8)
  par(mfrow=c(2,1))
  bundle$mainTitle<-substitute(paste(x," ",y), list(x=fileList[ff,"population"],y=fileList[ff,"phenotype"]))
  bundle$locNames<-getLocNames(gff_a,phenotypes,fileList[ff,"phenotype"],genomeVer)
  bundle$locs<-getLocPositions(gff_a,phenotypes,fileList[ff,"phenotype"],genomeIdx,genomeVer)
  #useLoci<-1:length(bundle$locs)
  useLoci<-which(bundle$locNames=="ppw-1")

  # plot raw Hawaii allele frequencies of control and selected along genome
  plotFreq(myExp,extraData=bundle,useLoci=useLoci)

  # plot the difference between Hawaii allele frequencies in control and selected populations
  # smoothe with sgolay polynomial smoothing (window size=sm) before plotting
  plotSmoothedDiff1(myExp,extraData=bundle,sm=1001,useLoci=useLoci)

  # plot -log10 of p values after polynomial smoothing (sgolay filter)
  plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
  title(sub="sm.logFisherPvals.sg1001")
  myExp["sm.logFisherPvals.sg1001"]<-smootheByChr1(-log10(myExp$padj.fdr), myExp$Chr, sm=1001)
  findDistppw2(myExp,ppw1_loc=c(4185062,4189932),dataCol="sm.logFisherPvals.sg1001")

  # plot -log10 of p values after tricube kernel smoothing
  plotVals(myExp,myExp$sm.logFisherPvals.tck1Mb,"tricube kernel smoothing 1Mb window",bundle,useLoci=useLoci,
           fdr="log10",ylab="-log10(Fisher test p values)")
  findDistppw2(myExp,ppw1_loc=c(4185062,4189932),dataCol="sm.logFisherPvals.tck1Mb")

  #plot sgolay smoothed G
  plotSmGvals_sGolay1(myExp, myExp$Gval.GTest, extraData=bundle, sm=1001, useLoci=useLoci)
  title(sub="sm.Gvals.sg1001")
  myExp["sm.Gvals.sg1001"]<-smootheByChr1(myExp$Gval.GTest, myExp$Chr, sm=1001)
  findDistppw2(myExp,ppw1_loc=c(4185062,4189932),dataCol="sm.Gvals.sg1001")

  #plot G vals smoothed by tricube kernel
  plotVals(myExp,myExp$sm.Gvals.tck1Mb,"tricube kernel smoothing 1Mb window",bundle,useLoci=useLoci,
           fdr="None",ylab="smoothed G values")
  findDistppw2(myExp,ppw1_loc=c(4185062,4189932),dataCol="sm.Gvals.tck1Mb")
  # #plot smG
  # chrList<-levels(myExp$Chr)
  # plot(myExp$Position,smG,type='n',lwd=2,ylab="Smoothed G value",
  #      xlab="Position (bp)",main=bundle$mainTitle,ylim=c(0,max(smG)*1.1))
  # plotByChr(myExp, smG, chrList, chrColumn="Chr",lwd=2)
  #
  # text(bundle$midpoint,max(smG)*1.1,bundle$chrNum,cex=1, col="black")
  # abline(v=bundle$chrEnds,col="dark gray",lty=5)
  # abline(v=bundle$locs[useLoci], col="red",lwd=0.8)
  # mtext(bundle$locNames[useLoci], at=bundle$locs[useLoci], cex=0.7, col="red",
  #       adj=bundle$labelAdj[useLoci])
  # title(sub="2Mb window smoothed GTest statistic (1000 permutations)")
  dev.off()
  ##########

  #output processed data as a table
  processedFileName<-paste0("../finalData/",genomeVer,"/processedData_",d,"_",fileList[ff,"phenotype"],
                           "_",fileList[ff,"population"],".txt")
  write.table(myExp,file=processedFileName)

  #add filename of processedData to fileList table
  fileList[ff,"processedData"]<-processedFileName
  fileList[ff+1,"processedData"]<-processedFileName
}

#write new fileList_date.txt file to store processedFileName
write.table(fileList,file=paste0("../finalData/",genomeVer,"/fileList_",d,".txt"))

#myExp<-read.table(paste0("../finalData/",genomeVer,"/processedData_20170417_RNAi_mel26_MM10_gen5.txt"),header=TRUE,stringsAsFactors=FALSE)
#fileList<-read.table(paste0("../finalData/",genomeVer,"/fileList_20170417.txt"),header=TRUE,stringsAsFactors=FALSE)
#ff=1

#plot histograms of readDepths
pdf(file=paste0("../finalData/",genomeVer,"/",d,"_readDepths.pdf"), paper="a4", height=11, width=8)
par(mfrow=c(4,2))
for (ff in seq(1,dim(fileList)[1],by=2)) {
  # read in control and selected data for single experiment
  myExp<-read.table(as.character(fileList[ff,"processedData"]),header=TRUE)

  # plot histograms of readDepths
  rD<-grep("readDepth",names(myExp))
  hist(myExp[,rD[1]],breaks=100, main=paste("Control",fileList[ff,2],fileList[ff,3]),xlab="Read depth")
  abline(v=mean(myExp[,rD[1]]),col="red",lwd=2)
  hist(myExp[,rD[2]],breaks=100, main=paste("Selected",fileList[ff,2],fileList[ff,3]),xlab="Read depth")
  abline(v=mean(myExp[,rD[2]]),col="red",lwd=2)

  #plot(as.matrix(myExp[rD[1]])~myExp$Position,main=paste("Control",fileList[ff,2],fileList[ff,3]),xlab="Position",type='l')

 }



dev.off()

###
# all this was with myExp that has 279244 SNPs in it.

### check quality metrics
IDmrg<-as.vector(myExp$ID)
i<-which(newCont$N2id.WS250 %in% IDmrg)
param<-names(newCont)[c(10,19:25,30,39:45)]

pdf(file="../plots/qualityMetrics.pdf",width=8,height=11,paper="a4")
par(mfrow=c(3,2))
for (j in param) {
  boxplot(as.numeric(newCont[i,j]),as.numeric(newCont[-i,j]),
          varwidth=TRUE,notch=TRUE, names=c("picked","rejected"),
          col=c("red","blue"),main=j,outline=FALSE)
}
dev.off()

###############################################################################3
#try filtering data with quality scores
removeBQ<-filter(newCont,RPB.x<0.1,RPB.y<0.1,BQB.x<0.1,BQB.y<0.1)
plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeBQ$N2id.WS250)
points(myExp[i,"Position"],myExp[i,"X1Fmerged_4999_CBfreq"],pch=16,col="darkblue")
# <0.6 do not look particularly bad
# <0.3 looks ok
# even <0.1 looks ok (644 points)

# these plot the points at the bottom of the plot
# pdf(file="../plots/mappingMetrics_sm.pdf",width=8,height=11, paper="a4")
# par(mfrow=c(2,1))
# Threshold=0
# removeMQ<-filter(newCont,QUAL.x==Threshold | QUAL.y==Threshold)
# plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
# #plotFreq(myExp,extraData=bundle,useLoci=useLoci)
# i<-which(myExp$ID %in% removeMQ$N2id.WS250)
# length(i)
# points(myExp[i,"Position"],myExp[i,"X1Fmerged_4999_CBfreq"],pch=16,
#        col="#0000ff55")
# title(sub=paste0("QUAL=",Threshold))
#
# Threshold=0
# removeMQ<-filter(newCont,MQB.x==Threshold | MQB.y==Threshold)
# plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
# #plotFreq(myExp,extraData=bundle,useLoci=useLoci)
# i<-which(myExp$ID %in% removeMQ$N2id.WS250)
# length(i)
# points(myExp[i,"Position"],myExp[i,"X1Fmerged_4999_CBfreq"],pch=16,
#        col="#0000ff55")
# title(sub=paste0("MQB=",Threshold))
#
# Threshold=0
# removeMQ<-filter(newCont,MQSB.x==Threshold | MQSB.y==Threshold)
# plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
# #plotFreq(myExp,extraData=bundle,useLoci=useLoci)
# i<-which(myExp$ID %in% removeMQ$N2id.WS250)
# length(i)
# points(myExp[i,"Position"],myExp[i,"X1Fmerged_4999_CBfreq"],pch=16,
#        col="#0000ff55")
# title(sub=paste0("MQSB=",Threshold))
#
# Threshold=10
# removeMQ<-filter(newCont,MQ.x<Threshold | MQ.y<Threshold)
# #plotFreq(myExp,extraData=bundle,useLoci=useLoci)
# plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
# i<-which(myExp$ID %in% removeMQ$N2id.WS250)
# length(i)
# points(myExp[i,"Position"],myExp[i,"X1Fmerged_4999_CBfreq"],pch=16,
#        col="#0000ff55")
# title(sub=paste0("MQ<",Threshold))
#
# dev.off()


#these overlay the plots on the smoothed line
pdf(file="../plots/mappingMetrics_sm.pdf",width=8,height=11, paper="a4")
par(mfrow=c(2,1))
Threshold=0
removeMQ<-filter(newCont,QUAL.x==Threshold | QUAL.y==Threshold)
plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
length(i)
points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=10,
       col="#0000ff55")
title(sub=paste0("QUAL=",Threshold,"  ",length(i)," SNPs"))

Threshold=0
removeMQ<-filter(newCont,MQB.x==Threshold | MQB.y==Threshold)
plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
length(i)
points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=10,
       col="#0000ff55")
title(sub=paste0("MQB=",Threshold,"  ",length(i)," SNPs"))

Threshold=0
removeMQ<-filter(newCont,MQSB.x==Threshold | MQSB.y==Threshold)
plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
length(i)
points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=10,
       col="#0000ff55")
title(sub=paste0("MQSB=",Threshold,"  ",length(i)," SNPs"))

Threshold=10
removeMQ<-filter(newCont,MQ.x<Threshold | MQ.y<Threshold)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
length(i)
points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=10,
       col="#0000ff55")
title(sub=paste0("MQ<",Threshold,"  ",length(i)," SNPs"))

dev.off()

#these remove low quality from smoothed data
pdf(file="../plots/mappingMetrics_sm_sub.pdf",width=8,height=11, paper="a4")
par(mfrow=c(2,1))
Threshold=0
removeMQ<-filter(newCont,QUAL.x==Threshold | QUAL.y==Threshold)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
plotSmPvals1(myExp[-i,], -log10(myExp$padj.fdr[-i]), extraData=bundle, sm=1001, useLoci=useLoci)
length(i)
#points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=16,col="#0000ff55")
title(sub=paste0("QUAL=",Threshold,"  ",length(i)," SNPs"))

Threshold=0
removeMQ<-filter(newCont,MQB.x==Threshold | MQB.y==Threshold)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
plotSmPvals1(myExp[-i,], -log10(myExp$padj.fdr[-i]), extraData=bundle, sm=1001, useLoci=useLoci)
length(i)
#points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=16,col="#0000ff55")
title(sub=paste0("MQB=",Threshold,"  ",length(i)," SNPs"))

Threshold=0
removeMQ<-filter(newCont,MQSB.x==Threshold | MQSB.y==Threshold)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
plotSmPvals1(myExp[-i,], -log10(myExp$padj.fdr[-i]), extraData=bundle, sm=1001, useLoci=useLoci)
length(i)
#points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=16,col="#0000ff55")
title(sub=paste0("MQSB=",Threshold,"  ",length(i)," SNPs"))

Threshold=10
removeMQ<-filter(newCont,MQ.x<Threshold | MQ.y<Threshold)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
plotSmPvals1(myExp[-i,], -log10(myExp$padj.fdr[-i]), extraData=bundle, sm=1001, useLoci=useLoci)
length(i)
#points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=16,col="#0000ff55")
title(sub=paste0("MQ<",Threshold,"  ",length(i)," SNPs"))


removeMQ<-filter(newCont,MQB.x==0 | MQB.y==0 |
                 MQSB.x==0 | MQSB.y==0 | MQ.x<10 | MQ.y<10)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
plotSmPvals1(myExp[-i,], -log10(myExp$padj.fdr[-i]), extraData=bundle, sm=1001, useLoci=useLoci)
length(i)
#points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=16,col="#0000ff55")
title(sub=paste0("combined","  ",length(i)," SNPs"))

removeMQ<-filter(newCont,MQB.x==0 | MQB.y==0 |
                    MQSB.x==0 | MQSB.y==0 | MQ.x<10 | MQ.y<10 |
                    QUAL.x==0 | QUAL.y==0)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
plotSmPvals1(myExp[-i,], -log10(myExp$padj.fdr[-i]), extraData=bundle, sm=1001, useLoci=useLoci)
length(i)
#points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=16,col="#0000ff55")
title(sub=paste0("combined+QUAL0","  ",length(i)," SNPs"))

dev.off()

################################################
### what about when both are bad?? (now that mapping to both genomes this is sure sign that SNP is unreliable)

pdf(file="../plots/mappingMetrics_sm_both.pdf",width=8,height=11, paper="a4")
par(mfrow=c(2,1))
Threshold=0
removeMQ<-filter(newCont,QUAL.x==Threshold & QUAL.y==Threshold)
plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
length(i)
points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=10,
       col="#0000ff55")
title(sub=paste0("QUAL=",Threshold,"  ",length(i)," SNPs"))

Threshold=50
removeMQ<-filter(newCont,QUAL.x<Threshold & QUAL.y<Threshold)
plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
length(i)
points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=10,
       col="#0000ff55")
title(sub=paste0("QUAL<",Threshold,"  ",length(i)," SNPs"))

Threshold=100
removeMQ<-filter(newCont,QUAL.x<Threshold & QUAL.y<Threshold)
plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
length(i)
points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=10,
       col="#0000ff55")
title(sub=paste0("QUAL<",Threshold,"  ",length(i)," SNPs"))

Threshold=150
removeMQ<-filter(newCont,QUAL.x<Threshold & QUAL.y<Threshold)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
plotSmPvals1(myExp, -log10(myExp$padj.fdr), extraData=bundle, sm=1001, useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
length(i)
points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=10,
       col="#0000ff55")
title(sub=paste0("QUAL<",Threshold,"  ",length(i)," SNPs"))

dev.off()


pdf(file="../plots/mappingMetrics_sm_sub_both.pdf",width=8,height=11, paper="a4")
par(mfrow=c(2,1))
Threshold=0
removeMQ<-filter(newCont,QUAL.x==Threshold & QUAL.y==Threshold)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
plotSmPvals1(myExp[-i,], -log10(myExp$padj.fdr[-i]), extraData=bundle, sm=1001, useLoci=useLoci)
length(i)
#points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=16,col="#0000ff55")
title(sub=paste0("QUAL=",Threshold,"  ",length(i)," SNPs"))

Threshold=50
removeMQ<-filter(newCont,QUAL.x<Threshold & QUAL.y<Threshold)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
plotSmPvals1(myExp[-i,], -log10(myExp$padj.fdr[-i]), extraData=bundle, sm=1001, useLoci=useLoci)
length(i)
#points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=16,col="#0000ff55")
title(sub=paste0("QUAL<",Threshold,"  ",length(i)," SNPs"))

Threshold=100
removeMQ<-filter(newCont,QUAL.x<Threshold & QUAL.y<Threshold)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
plotSmPvals1(myExp[-i,], -log10(myExp$padj.fdr[-i]), extraData=bundle, sm=1001, useLoci=useLoci)
length(i)
#points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=16,col="#0000ff55")
title(sub=paste0("QUAL<",Threshold,"  ",length(i)," SNPs"))

Threshold=150
removeMQ<-filter(newCont,QUAL.x<Threshold & QUAL.y<Threshold)
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
plotSmPvals1(myExp[-i,], -log10(myExp$padj.fdr[-i]), extraData=bundle, sm=1001, useLoci=useLoci)
length(i)
#points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=16,col="#0000ff55")
title(sub=paste0("QUAL<",Threshold,"  ",length(i)," SNPs"))


removeMQ<-filter(newCont,(MQB.x==0 | MQB.y==0) |
                    (MQSB.x==0 | MQSB.y==0) | (MQ.x<10 | MQ.y<10))
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
plotSmPvals1(myExp[-i,], -log10(myExp$padj.fdr[-i]), extraData=bundle, sm=1001, useLoci=useLoci)
length(i)
#points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=16,col="#0000ff55")
title(sub=paste0("combined","  ",length(i)," SNPs"))


removeMQ<-filter(newCont,(MQB.x==0 | MQB.y==0) |
                    (MQSB.x==0 | MQSB.y==0) | (MQ.x<10 | MQ.y<10) |
                    (QUAL.x==0 & QUAL.y==0))
#plotFreq(myExp,extraData=bundle,useLoci=useLoci)
i<-which(myExp$ID %in% removeMQ$N2id.WS250)
plotSmPvals1(myExp[-i,], -log10(myExp$padj.fdr[-i]), extraData=bundle, sm=1001, useLoci=useLoci)
length(i)
#points(myExp[i,"Position"],-log10(myExp[i,"padj.fdr"]),pch=16,col="#0000ff55")
title(sub=paste0("combined+QUAL0","  ",length(i)," SNPs"))

dev.off()

