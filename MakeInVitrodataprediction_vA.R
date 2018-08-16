library(reshape2)
library(beepr)
# Reparameterized Hill function
DR_up<-function(x=1, y0=1, x0=1, Emax=1, n=1){
  
  y0*(1+((x/x0)^n)/(1+((x/x0)^n)*(1/Emax))) 
  
}
## Load chemical map - mapping between names and numbers

chemmap<-read.csv("IVIVEchems.csv",as.is=TRUE)

## Concentration range plotted
concframe<-data.frame(x=c(1e-6,10^(seq(-3,2,0.05))))
data.list <- list()
DR_Frame_All <- data.frame()
for (j in 1:nrow(chemmap)) {
  cat(j,"...\n")
  tmplist<-list()
  ##
  chemical.num <- chemmap$Chemical.num[j]
  chemical.name <- chemmap$Chemical.name[j]
  tmplist$chemical.num<-chemical.num
  tmplist$chemical.name<-chemical.name
  #### in vitro data***
  file.sources=list.files(path="All",pattern=paste("*_",chemical.num,"_dat.R",sep=""))
  sapply(paste("All/",file.sources,sep=""),source,.GlobalEnv)
  invitrodat.list <- list(cell=cell,Ni=Ni,Nj=Nj,Nquants=Nquants,
                          quants=quants,scale_factor=scale_factor,x=x,ys=ys)
  tmplist$invitrodat.list<-invitrodat.list
  DR_Frame<- data.frame(Individual=invitrodat.list$cell, 
                        Dose=invitrodat.list$x, 
                        ScaledResponse=invitrodat.list$ys)
  DR_Frame$Response <- DR_Frame$ScaledResponse * invitrodat.list$scale_factor
  DR_Frame$Chemical.name <- chemical.name
  DR_Frame$Chemical.num <- chemical.num
  ### load in CSV files, R dat file,  set up concentration frame
  files = list.files(path="All",pattern=paste("*_",chemical.num,"_samples_*",sep=""))
  for (l in 1:length(files)) assign(files[l], read.csv(paste("All/",files[l],sep="")))
  #### load parameters***
  y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,125:151], (eval(parse(text=files[2])))[2005:4004,125:151], (eval(parse(text=files[3])))[2005:4004,125:151], (eval(parse(text=files[4])))[2005:4004,125:151]))
  colnames(y0)<- (1:27)
  rownames(y0)<- c(1:8000)
  parmframe<- melt(as.matrix(y0))
  colnames(parmframe)<- c("Iter", "ID", "y0")
  x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, 152:178], (eval(parse(text=files[2])))[2005:4004, 152:178], (eval(parse(text=files[3])))[2005:4004, 152:178], (eval(parse(text=files[4])))[2005:4004, 152:178]))
  colnames(x0)<- paste("ID",1:27, sep="")
  rownames(x0)<- c(1:8000)
  parmframe$x0<- melt(as.matrix(x0))$value
  Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, 179:205], (eval(parse(text=files[2])))[2005:4004, 179:205], (eval(parse(text=files[3])))[2005:4004, 179:205], (eval(parse(text=files[4])))[2005:4004, 179:205]))
  colnames(Emax)<- paste("ID",1:27, sep="")
  rownames(Emax)<- c(1:8000)
  parmframe$Emax<- melt(as.matrix(Emax))$value
  n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, 206:232], (eval(parse(text=files[2])))[2005:4004, 206:232], (eval(parse(text=files[3])))[2005:4004, 206:232], (eval(parse(text=files[4])))[2005:4004, 206:232]))
  colnames(n)<- paste("ID",1:27, sep="")
  rownames(n)<- c(1:8000)
  parmframe$n<- melt(as.matrix(n))$value
  tmplist$parmframe<-parmframe
  ### Concentration-response
  parmframeconc<- data.frame()
  for (i in 1:nrow(concframe)){
    parmframe$concentration<-concframe$x[i]
    parmframeconc<- rbind(parmframeconc, parmframe)
  }
  parmframeconc$Pred<- DR_up(x=parmframeconc$concentration, 
                              y0=parmframeconc$y0, 
                              x0=parmframeconc$x0, 
                              Emax=parmframeconc$Emax, 
                              n=parmframeconc$n)*
    invitrodat.list$scale_factor
  tmplist$parmframeconc<-parmframeconc
  predframe_median<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.5)
  predframe_975<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.975)
  predframe_25<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.025)
  Prediction_Frame<- predframe_median
  Prediction_Frame$p97.5<- predframe_975$Pred
  Prediction_Frame$p2.5<- predframe_25$Pred
  names(Prediction_Frame)[3]<- "p50"
  names(Prediction_Frame)[2]<- "Individual"
  tmplist$Prediction_Frame<-Prediction_Frame
  
  ### Random prediction 
  DR_Frame$Predicted<- t(DR_up(x=DR_Frame$Dose, 
                               x0=x0[8000,DR_Frame$Individual], 
                               y0=y0[8000,DR_Frame$Individual], 
                               n=n[8000,DR_Frame$Individual], 
                               Emax=Emax[8000,DR_Frame$Individual]))*
    invitrodat.list$scale_factor
  tmplist$DR_Frame<-DR_Frame
  DR_Frame_All <- rbind(DR_Frame_All,DR_Frame)
  data.list[[j]]<-tmplist
}
DR_Frame_Controls <- subset(DR_Frame_All,Dose==0)
DR_Frame_Controls$Chemical.name<-"Controls"
DR_Frame_Treated <- subset(DR_Frame_All,Dose>0)

save(data.list,chemmap,DR_Frame_All,
     DR_Frame_Controls,
     DR_Frame_Treated,
     file="InVitrodatandprediction.Rdata")
beep(sound=3)
  
  