library(tidyr)

nonblind<- read.table("Nonblinded.peak.parms.sum.dat", header = TRUE, as.is = TRUE)

blind<- read.table("Blinded.peak.parms.stdnum.sum.dat", header = TRUE, as.is = TRUE)

blind<-drop_na(blind)

alldat<-rbind(nonblind, blind)

write.table(alldat, "Nonblinded.blinded.peak.parms.sum.dat", sep="\t",row.names=FALSE)
