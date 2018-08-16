###make sure that .seq files aren't duplicated within themselves

###Change working directory

setwd("C:\\Users\\ablanchette\\Documents\\CM Peak Processing\\Cardiotox_Repeat Experiment\\1083")

source("ChiuPeakProcessingAlgorithm.R")

source("ChiuPeakProcessingUtilities.R")

source("tcpl_fit_functions.R")

library(ggplot2)

seq1.file.list<- get.seq1.file.list(".")

n.files<-length(seq1.file.list)

dat.list<-list()

for (n in 1:n.files) {
  
  fname<-seq1.file.list[n]
  
  dat.list[[n]]<-get.seq1.dat(fname)

}

## While file?

filenum<-1

## "dat" is a data frame with Well (e.g., E3) in the first column, and the data in the other columns

dat<-dat.list[[filenum]]

t<-as.numeric(sub("X","",colnames(dat)[-1]))

##Set up loop to export each trace as a .png into new directory folder

setwd("C:\\Users\\ablanchette\\Documents\\CM Peak Processing\\Cardiotox_Repeat Experiment\\1083\\1083_Traces")

for (i in 1:length(dat$Well)) {

  plotname<-paste("1083_Trace_", dat$Well[i], ".png", sep="")
  
  rownum<-i

  y<-dat[rownum,-1]

  well<-dat[rownum,1]

  png(plotname)
  
  plot(t,y,cex=0.5,type="l");

  mtext(paste(seq1.file.list[filenum]," Row ",rownum," (Well ",well,")",sep=""),side=3,cex=0.75)
  
  dev.off()
}
