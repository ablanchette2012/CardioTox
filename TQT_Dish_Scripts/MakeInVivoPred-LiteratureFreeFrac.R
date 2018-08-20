

invivoparms <- read.csv("QTposnegparameters-Literature.csv",as.is=TRUE)

linearfuncfoldchange <- function(xfree,parms) {
  xtot <- xfree/parms$FracFreePlasma
  y <- 1 + parms$SlopemsecuMtot/parms$QTc0Linmsec*xtot
  y
}

hillfuncfoldchange <- function(xfree,parms) {
  xtot <- xfree/parms$FracFreePlasma
  y <- 1 + parms$Emaxmsec*xfree^parms$nHill/
    (parms$QTc0Hillmsec*(parms$KduMtot^parms$nHill + xfree^parms$nHill))
  y
}

linearfuncfracchange <- function(xfree,parms) {
  xtot <- xfree/parms$FracFreePlasma
  y <- parms$SlopemsecuMtot/parms$QTc0Linmsec*xtot
  y
}

hillfuncfracchange <- function(xfree,parms) {
  xtot <- xfree/parms$FracFreePlasma
  y <- parms$Emaxmsec*xfree^parms$nHill/
    (parms$QTc0Hillmsec*(parms$KduMtot^parms$nHill + xfree^parms$nHill))
  y
}

invivodf <- data.frame()
for (i in 1:nrow(invivoparms)) {
  parmnow<-invivoparms[i,]
  if (!is.na(parmnow$CmaxLinuMtot)) {
    xfree <- 10^seq(-3,log10(parmnow$CmaxLinuMtot*parmnow$FracFreePlasma),0.05)
    nx <- length(xfree)
    dflin <- data.frame(Chemical.name=rep(parmnow$Chemical,nx),
                        Model=rep("in vivo Linear",nx),
                        xfree=xfree,
                        PredFoldChange=linearfuncfoldchange(xfree,parmnow),
                        PredFracChange=linearfuncfracchange(xfree,parmnow))
    invivodf <- rbind(invivodf,dflin)
  }
  if (!is.na(parmnow$CmaxHilluMtot)) {
    xfree <- 10^seq(-3,log10(parmnow$CmaxHilluMtot*parmnow$FracFreePlasma),0.05)
    nx <- length(xfree)
    dfhill <- data.frame(Chemical.name=rep(parmnow$Chemical,nx),
                         Model=rep("in vivo Hill",nx),
                         xfree=xfree,
                         PredFoldChange=hillfuncfoldchange(xfree,parmnow),
                         PredFracChange=hillfuncfracchange(xfree,parmnow))
    invivodf <- rbind(invivodf,dfhill)
  }
}

invivodf.pos <- subset(invivodf,PredFracChange>0)
invivodf.neg <- subset(invivodf,PredFracChange==0)

write.csv(invivodf.pos,file="InVivoPred-LiteratureFreeFrac.pos.csv")
write.csv(invivodf.neg,file="InVivoPred-LiteratureFreeFrac.neg.csv")
