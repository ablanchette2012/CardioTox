library(beepr)
# Reparameterized Hill function
DR_up<-function(x=1, y0=1, x0=1, Emax=1, n=1){
  
  y0*(1+((x/x0)^n)/(1+((x/x0)^n)*(1/Emax))) 
  
}

# Fractional change from control
# without y0, and without "1 +"
DR_up_Frac<- function(x=1, x0=1, Emax=1, n=1){
  (x/x0)^n/(1+((x/x0)^n)*(1/Emax))
}

## Load chemical map - mapping between names and numbers

chemmap<-read.csv("IVIVEchems.csv",as.is=TRUE)

## Load in vivo results
InVivoPos<- read.csv("InVivoPred.pos.v2.csv", header=TRUE)
InVivoNeg<- read.csv("InVivoPred.neg.v2.csv", header=TRUE)
InVivoParm<- read.csv("QTposnegparameters-Literature.csv", header=TRUE)
colnames(InVivoParm)[1]<- "Chemical"
## Concentration range plotted
concframe<-data.frame(x=c(1e-6,10^(seq(-3,2,0.05))))

load("InVivoInVitrodat-ExperimentalFreeFrac_43.Rdata")
for (j in 1:nrow(chemmap)) {
  cat(j,"...\n")
  tmplist<-data.list[[j]]
  ##
  chemical.num <- chemmap$Chemical.num[j]
  chemical.name <- chemmap$Chemical.name[j]
  tmplist$chemical.num<-chemical.num
  tmplist$chemical.name<-chemical.name
  
  ### In vivo data
  InVivo<-subset(InVivoPos, Chemical.name== chemical.name )
  if (nrow(InVivo)==0) {
    InVivo<-subset(InVivoNeg, Chemical.name== chemical.name )
  }
  InVivo$Model <- factor(InVivo$Model,levels=c("in vivo Linear","in vivo Hill"))
  InVivoFree<- subset(InVivoParm, Chemical== chemical.name )
  tmplist$InVivo <- InVivo
  tmplist$InVivoFree <- InVivoFree
  
  data.list[[j]]<-tmplist
  lapply(data.list[[1]],head)
}

save(data.list,chemmap,file="InVivoInVitrodat-LiteratureFreeFrac_43.Rdata")
beep(sound=3)