###make sure that .seq files aren't duplicated within themselves

###Change working directory

setwd("C:\\Users\\Alex Blanchette\\Documents\\Chiu Lab\\Cardiotox_Repeat Experiment\\1434")

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

nplates=1
plate.map.frame<-get.plate.design(nplates=nplates)
chem.map<-get.chem.map(fname="chemicals.dat")
plate.map.frame$Chemical.name<-chem.map[paste(plate.map.frame$Chemical.Number),]

## While file?

filenum<-1

## "dat" is a data frame with Well (e.g., E3) in the first column, and the data in the other columns

dat<-dat.list[[filenum]]

sec<-matrix(as.numeric(sub("X","",colnames(dat)[-1])))

controls<- subset(plate.map.frame, Solvent=="DMSO" & Chemical.Concentration==0)

##Set up loop to export plots in pdf file

plotname<-paste("1434_Trace_Controls",".pdf", sep="")

pdf(plotname)

par(mfrow=c(3,2))

for (rowcol in controls$rowcol) {

  onedat<- subset(dat, Well==rowcol)

  rfu<-as.numeric(onedat[,-1]) 

  p1= plot(sec,rfu,cex=0.5,type="l", col="red", xlim=c(0,30));

  mtext(paste("Controls"," Row ",rownum," (Well ",rowcol,")",sep=""),side=3,cex=0.75)
  
}

dev.off()
