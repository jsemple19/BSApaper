# script to convert genomicRanges of data to smoothed data (mean) in sliding windows

setwd("~/Documents/LehnerLab/BSApaper")
samples<-read.table("./sampleList.txt",header=T,stringsAsFactors=F)

#fileList<-list.files(path="./rds",pattern="[0-9A-Z]{4}_m.*RDS")
#write.table(fileList,"./rdsFileList.txt",row.names=F,col.names=F)

fileList<-read.table("./rdsFileList.txt",header=F,stringsAsFactors=F)[,1]
library("BSgenome.Celegans.UCSC.ce11")
genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6],ranges=IRanges(start=1,
                                                                  end=seqlengths(Celegans)[1:6]), strand="*")
seqlevels(genomeGR)<-gsub("chr","",seqlevels(genomeGR))



#######################################################################
########### sliding windows to create smoothed genomicRanges ##########
#######################################################################

## function to create list of GRanges with sliding windows of different sizes
makeWinList<-function(GRobj,widths,step=0.1){
  winList<-list()
  for (winSize in widths) {
    # get formatted text for window size
    units<-formatWinSize(winSize)
    # create GRanges of windows of that size and add to list of windows
    winList[paste0("win",units)]<-unlist(slidingWindows(x=GRobj,width=winSize,step=winSize*step))
  }
  return(winList)
}

# function to take a numerical window size and convert to nice text for labels etc
formatWinSize<-function(winSize) {
  if (winSize/1000 < 1) {
    units<-paste0(winSize,"bp")
  } else if (winSize/1e6 >= 1) {
    units<-paste0(winSize/1e6,"Mb")
  } else {
    units<-paste0(winSize/1000,"kb")
  }
  return(units)
}

#sizes of windows you want to create
winWidths<-c(1e4,1e5,5e5)
#create create list of GRanges with sliding windows of different sizes
winList<-makeWinList(genomeGR,winWidths)


# read in the bedgraph files in pairs to calculate enrichment counts at different scales
for (f in fileList) {
  #get normalised counts
  data<-readRDS(paste0("./rds/",f))
  sampleCode=strsplit(f,split="_",fixed=T)[[1]][1]

  for (i in 1:length(winList)) {
    #make a GRangesList from GRanges
    windows<-split(winList[[i]],seqnames(winList[[i]]))

    # convert GRanges back to coverage
    rle<-coverage(data,weight="score")

    #create views corresponding to the windows on the coverage rle and average counts
    summedCov<-viewMeans(RleViewsList(rangesList=as(windows,"RangesList"),rleList = rle))

    #now save the values into windows
    windows<-unlist(windows)
    mcols(windows)$score<-unlist(summedCov)
    #windows<-resize(windows,width=1,fix="center")
    saveRDS(windows,paste0("./rds/",names(winList)[i],"__",
                           paste(samples[samples$code==sampleCode,c("code","population","phenotype","days")],collapse="__"),".rds"))
  }
}



###################################################################
######### Plotting data by window size to compare conditions ######
###################################################################

#sizes of windows you want to create
winWidths<-c(1e4,1e5,5e5)

#create create list of sliding windows of different sizes
winSizeList<-unlist(lapply(winWidths,formatWinSize))

for (w in winSizeList) {
  fileList<-list.files(path="./rds",pattern=w,full.names=T)
  data<-lapply(fileList,readRDS)
  names(data)<-basename(fileList)
  dataWide<-data[[1]]
  mcols(dataWide)<-as.data.frame(lapply(data,function(x) {mcols(x)$score}))

}


# function to make genomeGR which is used for plotting and making windows
makeGenomeGR<-function(genome=Celegans) {
  genomeGR<-GRanges(seqnames=seqnames(Celegans)[1:6],ranges=IRanges(start=1,
                                                                    end=seqlengths(Celegans)[1:6]), strand="*")
  seqlengths(genomeGR)<-end(genomeGR)
  cumStart<-1+cumsum(c(0,end(genomeGR)[1:5]))
  cumEnd<-cumsum(end(genomeGR))
  midPoints<-cumStart+floor(end(genomeGR)/2)
  seqlevels(genomeGR)<-gsub("chr","",seqlevels(genomeGR))
  mcols(genomeGR)<-DataFrame(cumStart,cumEnd,midPoints)
  return(genomeGR)
}

genomeGR<-makeGenomeGR()

plotGenomeCanvas4freq<-function(genomeGR){
  options(scipen=8)
  plot(c(1,genomeGR$cumEnd),rep(1,length(genomeGR$cumEnd)+1),type="n",ylim=c(0,1),
       xlab="Position in genome (bp)",ylab="Frequency")
  abline(v=c(genomeGR$cumStart,genomeGR$cumEnd), col="light grey",lty=2)
  text(genomeGR$midPoints,0,labels=paste0("chr",seqnames(genomeGR)))
}

plotGenomeCanvas4freq(genomeGR)

plotGR<-function(gr,genomeGR, samples) {
  resize(gr,width=1,fix="center")
  for (chr in seqlevels(gr)) {
    subGR<-gr[seqnames(gr)==chr,]
    cumPos<-start(subGR)+mcols(genomeGR[seqlevels(genomeGR)==chr,])$cumStart-1
    for (c in names(mcols(subGR))) {
      lines(cumPos,mcols(subGR)[,c])
    }
  }

  lines(start(gr),gr$score)
}
