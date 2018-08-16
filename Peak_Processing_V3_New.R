###make sure that .seq files aren't duplicated within themselves
###Change working directory
setwd("C:\\Users\\ablanchette\\Documents\\CM Peak Processing\\Cardiotox_Repeat Experiment\\1083")
source("ChiuPeakProcessingAlgorithm.R")
source("ChiuPeakProcessingUtilities.R")
source("tcpl_fit_functions.R")
library(tcpl)
seq1.file.list<- get.seq1.file.list(".")
n.files<-length(seq1.file.list)
dat.list<-list()
for (n in 1:n.files) {
  fname<-seq1.file.list[n]
  dat.list[[n]]<-get.seq1.dat(fname)
}
## While file?
filenum<-1
dat<-dat.list[[filenum]]
t<-as.numeric(sub("X","",colnames(dat)[-1]))

## "dat" is a data frame with Well (e.g., E3) in the first column, and the data in the other columns

## Which row?
rownum<-100
y<-dat[rownum,-1]
well<-dat[rownum,1]
plot(t,y,cex=0.5,type="l");
mtext(paste(seq1.file.list[filenum]," Row ",rownum," (Well ",well,")",sep=""),side=3,cex=0.75)
peak.parms.vec.list<-list()
peak.parms.sum.list<-list()
for (n in 1:n.files) {
  peak.parms.vec.list[[n]]<-get.peak.parms(dat.list[[n]],trimsec=10)
  peak.parms.sum.list[[n]]<-summarize.peak.parms(peak.parms.vec.list[[n]])
}
print(t(head(peak.parms.sum.list[[n]],3)))
nplates<-1
plate.map.frame<-get.plate.design(nplates=nplates)
chem.map<-get.chem.map(fname="chemicals.dat")
plate.map.frame$Chemical.name<-chem.map[paste(plate.map.frame$Chemical.Number),]
plate.peak.parms.sum.list<-list()
timepoint<-c("15min","30min","60min","90min","24hr")
for (n in 1:n.files) {
  plate.peak.parms.sum.list[[n]]<-merge(plate.map.frame,
                                        peak.parms.sum.list[[n]],
                                        by.x="rowcol",by.y="Well",sort=FALSE)
  plate.peak.parms.sum.list[[n]]<-plate.peak.parms.sum.list[[n]][
    order(plate.peak.parms.sum.list[[n]]$Chemical.Number,
          plate.peak.parms.sum.list[[n]]$Solvent,
          plate.peak.parms.sum.list[[n]]$Chemical.Concentration),]
  ## Add timepoint to data frame
  plate.peak.parms.sum.list[[n]]$TimePoint<-timepoint[n]
  ncol<-dim(plate.peak.parms.sum.list[[n]])[2]
  plate.peak.parms.sum.list[[n]]<-plate.peak.parms.sum.list[[n]][,c(ncol,1:(ncol-1))]
}
print(t(head(plate.peak.parms.sum.list[[n]],3)))
peak.parms.sum.all<-plate.peak.parms.sum.list[[1]]
if (n.files>1) {
  for (n in 2:n.files) {
    peak.parms.sum.all<-rbind(peak.parms.sum.all,plate.peak.parms.sum.list[[n]])
  }
}
write.csv(peak.parms.sum.all,file="Peak.Parms.Summary.csv",row.names=FALSE)
