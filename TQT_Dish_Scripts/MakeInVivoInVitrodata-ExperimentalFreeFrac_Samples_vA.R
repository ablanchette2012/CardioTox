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
InVivoParm<- read.csv("QTposnegparameters-Experimental.csv", header=TRUE)
colnames(InVivoParm)[1]<- "Chemical"
data.samples.list <- list()
for (j in 1:nrow(chemmap)) {
  cat(j,"...\n")
  tmplist<-list()
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
  
  ### load in CSV files, R dat file,  set up concentration frame
  files = list.files(path="All",pattern=paste("*_",chemical.num,"_samples_*",sep=""))
  for (l in 1:length(files)) assign(files[l], read.csv(paste("All/",files[l],sep="")))
  file.sources=list.files(path="All",pattern=paste("*_",chemical.num,"_dat.R",sep=""))
  sapply(paste("All/",file.sources,sep=""),source,.GlobalEnv)
  #### in vitro data***
  invitrodat.list <- list(cell=cell,Ni=Ni,Nj=Nj,Nquants=Nquants,
                          quants=quants,scale_factor=scale_factor,x=x,ys=ys)
  tmplist$invitrodat.list<-invitrodat.list
  
  ### extract median parameter data from csv files, create new frame with calculated parameter data, do the same when z scores and SD =0 
  #### simframe with random individual ***
  ##### y0
  m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))
  m_y0$z_score_sim<- rnorm(8000,0,1)
  simframe<- as.data.frame(exp(m_y0$m_y0+(m_y0$sd_y0*m_y0$z_score_sim)))
  rownames(simframe)<- c(1:8000)
  colnames(simframe)<- "sim_y0"
  ##### x0
  m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))
  m_x0$z_score_sim<- rnorm(8000,0,1)
  simframe$sim_x0<-exp(m_x0$m_x0+(m_x0$sd_x0*m_x0$z_score_sim))
  ##### Emax
  m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))
  m_Emax$z_score_sim<- rnorm(8000,0,1)
  simframe$sim_Emax<-exp(m_Emax$m_Emax+(m_Emax$sd_Emax*m_Emax$z_score_sim))
  ##### n
  m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))
  m_n$z_score_sim<- rnorm(8000,0,1)
  simframe$sim_n<-exp(m_n$m_n+(m_n$sd_n*m_n$z_score_sim))
  #### Population median
  simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), Emax=(exp(m_Emax$m_Emax)), n=(exp(m_n$m_n)))
  #### Standard donor (cell line 27)
  simframe_27<-data.frame(y0=c((eval(parse(text=files[1])))[2005:4004, "y0.27"], (eval(parse(text=files[2])))[2005:4004, "y0.27"], (eval(parse(text=files[3])))[2005:4004, "y0.27"], (eval(parse(text=files[4])))[2005:4004, "y0.27"]),
                          x0=c((eval(parse(text=files[1])))[2005:4004, "x0.27"], (eval(parse(text=files[2])))[2005:4004, "x0.27"], (eval(parse(text=files[3])))[2005:4004, "x0.27"], (eval(parse(text=files[4])))[2005:4004, "x0.27"]), 
                          Emax=c((eval(parse(text=files[1])))[2005:4004, "Emax.27"], (eval(parse(text=files[2])))[2005:4004, "Emax.27"], (eval(parse(text=files[3])))[2005:4004, "Emax.27"], (eval(parse(text=files[4])))[2005:4004, "Emax.27"]), 
                          n=c((eval(parse(text=files[1])))[2005:4004, "n.27"], (eval(parse(text=files[2])))[2005:4004, "n.27"], (eval(parse(text=files[3])))[2005:4004, "n.27"], (eval(parse(text=files[4])))[2005:4004, "n.27"]))
  
  ## Concentration range plotted (nominal from 1e-3 to 100)
  concframe<-data.frame(x=10^(seq(-3,2,0.05)))*InVivoFree$FracFreePlasma
  
  ### set up working frames with concentrations specified in concframe
  simframeconc<- data.frame()
  simframe_zero_conc<- data.frame()
  simframe_27_conc<- data.frame()
  for (i in 1:nrow(concframe)) {
    simframe$concentration<-concframe$x[i]
    simframe$freeconcentration<-concframe$x[i]*InVivoFree$FracFreeMedia
    simframeconc<- rbind(simframeconc, simframe)
  }
  for (i in 1:nrow(concframe)) {
    simframe_zero$concentration<-concframe$x[i]
    simframe_zero$freeconcentration<-concframe$x[i]*InVivoFree$FracFreeMedia
    simframe_zero_conc<- rbind(simframe_zero_conc, simframe_zero)
  }
  for (i in 1:nrow(concframe)) {
    simframe_27$concentration<-concframe$x[i]
    simframe_27$freeconcentration<-concframe$x[i]*InVivoFree$FracFreeMedia
    simframe_27_conc<- rbind(simframe_27_conc, simframe_27)
  }
  #### Frac change - random indiv
  simframeconc$PredFrac<- DR_up_Frac(x=simframeconc$concentration,  x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n) 
  tmplist$simframeconc<-simframeconc
  #### Frac change - Population median
  simframe_zero_conc$PredFrac<- DR_up_Frac(x=simframe_zero_conc$concentration,  x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n) 
  tmplist$simframe_zero_conc<-simframe_zero_conc
  #### Frac change - Standard donor
  simframe_27_conc$PredFrac<- DR_up_Frac(x=simframe_27_conc$concentration,  x0=simframe_27_conc$x0, Emax=simframe_27_conc$Emax, n=simframe_27_conc$n) 
  tmplist$simframe_27_conc<-simframe_27_conc
  data.samples.list[[j]]<-tmplist
  lapply(data.samples.list[[j]],head)
}

save(data.samples.list,chemmap,file="InVivoInVitrodat-ExperimentalFreeFracSamples.Rdata")
beep(sound=3)
  