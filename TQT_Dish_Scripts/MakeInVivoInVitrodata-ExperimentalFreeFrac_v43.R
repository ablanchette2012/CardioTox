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
## Concentration range plotted
concframe<-data.frame(x=c(1e-6,10^(seq(-3,2,0.05))))

data.list <- list()
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
  file.sources=list.files(path="All",pattern=paste("C",chemical.num, "_standat_v5_do_stan_fit_Peak_Decay_Rise_Ratio_up_dat.R",sep=""))
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
  files = list.files(path="All",pattern=paste("C",chemical.num,"_standat_v5_do_stan_fit_Peak_Decay_Rise_Ratio_up_samples*",sep=""))
  for (l in 1:length(files)) assign(files[l], read.csv(paste("All/",files[l],sep="")))
  
  ### extract median parameter data from csv files, create new frame with calculated parameter data, do the same when z scores and SD =0 
  #### simframe with random individual ***
  ##### y0
  m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001,c(8,12)],(eval(parse(text=files[2])))[2:2001, c(8,12)], (eval(parse(text=files[3])))[2:2001,c(8,12)], (eval(parse(text=files[4])))[2:2001,c(8,12)]))
  m_y0$z_score_sim<- rnorm(8000,0,1)
  simframe<- as.data.frame(exp(m_y0$m_y0+(m_y0$sd_y0*m_y0$z_score_sim)))
  rownames(simframe)<- c(1:8000)
  colnames(simframe)<- "sim_y0"
  ##### x0
  m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, c(9,13)], (eval(parse(text=files[2])))[2:2001, c(9,13)], (eval(parse(text=files[3])))[2:2001, c(9,13)], (eval(parse(text=files[4])))[2:2001, c(9,13)]))
  m_x0$z_score_sim<- rnorm(8000,0,1)
  simframe$sim_x0<-exp(m_x0$m_x0+(m_x0$sd_x0*m_x0$z_score_sim))
  ##### Emax
  m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, c(10,14)], (eval(parse(text=files[2])))[2:2001,c(10,14)], (eval(parse(text=files[3])))[2:2001, c(10,14)], (eval(parse(text=files[4])))[2:2001, c(10,14)]))
  m_Emax$z_score_sim<- rnorm(8000,0,1)
  simframe$sim_Emax<-exp(m_Emax$m_Emax+(m_Emax$sd_Emax*m_Emax$z_score_sim))
  ##### n
  m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2:2001, c(11,15)], (eval(parse(text=files[2])))[2:2001, c(11,15)], (eval(parse(text=files[3])))[2:2001, c(11,15)], (eval(parse(text=files[4])))[2:2001, c(11,15)]))
  m_n$z_score_sim<- rnorm(8000,0,1)
  simframe$sim_n<-exp(m_n$m_n+(m_n$sd_n*m_n$z_score_sim))
  #### Population median
  simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), Emax=(exp(m_Emax$m_Emax)), n=(exp(m_n$m_n)))
  #### Standard donor (cell line 27)
  simframe_27<-data.frame(y0=c((eval(parse(text=files[1])))[2:2001, "y0.28"], (eval(parse(text=files[2])))[2:2001, "y0.28"], (eval(parse(text=files[3])))[2:2001, "y0.28"], (eval(parse(text=files[4])))[2:2001, "y0.28"]),
                          x0=c((eval(parse(text=files[1])))[2:2001, "x0.28"], (eval(parse(text=files[2])))[2:2001, "x0.28"], (eval(parse(text=files[3])))[2:2001, "x0.28"], (eval(parse(text=files[4])))[2:2001, "x0.28"]), 
                          Emax=c((eval(parse(text=files[1])))[2:2001, "Emax.28"], (eval(parse(text=files[2])))[2:2001, "Emax.28"], (eval(parse(text=files[3])))[2:2001, "Emax.28"], (eval(parse(text=files[4])))[2:2001, "Emax.28"]), 
                          n=c((eval(parse(text=files[1])))[2:2001, "n.28"], (eval(parse(text=files[2])))[2:2001, "n.28"], (eval(parse(text=files[3])))[2:2001, "n.28"], (eval(parse(text=files[4])))[2:2001, "n.28"]))
  
  #### Standard donor (cell line 22)
  simframe_22<-data.frame(y0=c((eval(parse(text=files[1])))[2:2001, "y0.22"], (eval(parse(text=files[2])))[2:2001, "y0.22"], (eval(parse(text=files[3])))[2:2001, "y0.22"], (eval(parse(text=files[4])))[2:2001, "y0.22"]),
                          x0=c((eval(parse(text=files[1])))[2:2001, "x0.22"], (eval(parse(text=files[2])))[2:2001, "x0.22"], (eval(parse(text=files[3])))[2:2001, "x0.22"], (eval(parse(text=files[4])))[2:2001, "x0.22"]), 
                          Emax=c((eval(parse(text=files[1])))[2:2001, "Emax.22"], (eval(parse(text=files[2])))[2:2001, "Emax.22"], (eval(parse(text=files[3])))[2:2001, "Emax.22"], (eval(parse(text=files[4])))[2:2001, "Emax.22"]), 
                          n=c((eval(parse(text=files[1])))[2:2001, "n.22"], (eval(parse(text=files[2])))[2:2001, "n.22"], (eval(parse(text=files[3])))[2:2001, "n.22"], (eval(parse(text=files[4])))[2:2001, "n.22"]))
  
  #### Standard donor (cell line 11)
  simframe_11<-data.frame(y0=c((eval(parse(text=files[1])))[2:2001, "y0.11"], (eval(parse(text=files[2])))[2:2001, "y0.11"], (eval(parse(text=files[3])))[2:2001, "y0.11"], (eval(parse(text=files[4])))[2:2001, "y0.11"]),
                          x0=c((eval(parse(text=files[1])))[2:2001, "x0.11"], (eval(parse(text=files[2])))[2:2001, "x0.11"], (eval(parse(text=files[3])))[2:2001, "x0.11"], (eval(parse(text=files[4])))[2:2001, "x0.11"]), 
                          Emax=c((eval(parse(text=files[1])))[2:2001, "Emax.11"], (eval(parse(text=files[2])))[2:2001, "Emax.11"], (eval(parse(text=files[3])))[2:2001, "Emax.11"], (eval(parse(text=files[4])))[2:2001, "Emax.11"]), 
                          n=c((eval(parse(text=files[1])))[2:2001, "n.11"], (eval(parse(text=files[2])))[2:2001, "n.11"], (eval(parse(text=files[3])))[2:2001, "n.11"], (eval(parse(text=files[4])))[2:2001, "n.11"]))
  
  #### Standard donor (cell line 6)
  simframe_6<-data.frame(y0=c((eval(parse(text=files[1])))[2:2001, "y0.6"], (eval(parse(text=files[2])))[2:2001, "y0.6"], (eval(parse(text=files[3])))[2:2001, "y0.6"], (eval(parse(text=files[4])))[2:2001, "y0.6"]),
                         x0=c((eval(parse(text=files[1])))[2:2001, "x0.6"], (eval(parse(text=files[2])))[2:2001, "x0.6"], (eval(parse(text=files[3])))[2:2001, "x0.6"], (eval(parse(text=files[4])))[2:2001, "x0.6"]), 
                         Emax=c((eval(parse(text=files[1])))[2:2001, "Emax.6"], (eval(parse(text=files[2])))[2:2001, "Emax.6"], (eval(parse(text=files[3])))[2:2001, "Emax.6"], (eval(parse(text=files[4])))[2:2001, "Emax.6"]), 
                         n=c((eval(parse(text=files[1])))[2:2001, "n.6"], (eval(parse(text=files[2])))[2:2001, "n.6"], (eval(parse(text=files[3])))[2:2001, "n.6"], (eval(parse(text=files[4])))[2:2001, "n.6"]))
  
  #### Standard donor (cell line 3)
  simframe_3<-data.frame(y0=c((eval(parse(text=files[1])))[2:2001, "y0.3"], (eval(parse(text=files[2])))[2:2001, "y0.3"], (eval(parse(text=files[3])))[2:2001, "y0.3"], (eval(parse(text=files[4])))[2:2001, "y0.3"]),
                         x0=c((eval(parse(text=files[1])))[2:2001, "x0.3"], (eval(parse(text=files[2])))[2:2001, "x0.3"], (eval(parse(text=files[3])))[2:2001, "x0.3"], (eval(parse(text=files[4])))[2:2001, "x0.3"]), 
                         Emax=c((eval(parse(text=files[1])))[2:2001, "Emax.3"], (eval(parse(text=files[2])))[2:2001, "Emax.3"], (eval(parse(text=files[3])))[2:2001, "Emax.3"], (eval(parse(text=files[4])))[2:2001, "Emax.3"]), 
                         n=c((eval(parse(text=files[1])))[2:2001, "n.3"], (eval(parse(text=files[2])))[2:2001, "n.3"], (eval(parse(text=files[3])))[2:2001, "n.3"], (eval(parse(text=files[4])))[2:2001, "n.3"]))
  
  
  ### set up working frames with concentrations specified in concframe
  simframeconc<- data.frame()
  simframe_zero_conc<- data.frame()
  simframe_27_conc<- data.frame()
  simframe_22_conc<- data.frame()
  simframe_11_conc<- data.frame()
  simframe_6_conc<- data.frame()
  simframe_3_conc<- data.frame()
  
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
  for (i in 1:nrow(concframe)) {
    simframe_22$concentration<-concframe$x[i]
    simframe_22$freeconcentration<-concframe$x[i]*InVivoFree$FracFreeMedia
    simframe_22_conc<- rbind(simframe_22_conc, simframe_22)
  }
  for (i in 1:nrow(concframe)) {
    simframe_11$concentration<-concframe$x[i]
    simframe_11$freeconcentration<-concframe$x[i]*InVivoFree$FracFreeMedia
    simframe_11_conc<- rbind(simframe_11_conc, simframe_11)
  }
  for (i in 1:nrow(concframe)) {
    simframe_6$concentration<-concframe$x[i]
    simframe_6$freeconcentration<-concframe$x[i]*InVivoFree$FracFreeMedia
    simframe_6_conc<- rbind(simframe_6_conc, simframe_6)
  }
  for (i in 1:nrow(concframe)) {
    simframe_3$concentration<-concframe$x[i]
    simframe_3$freeconcentration<-concframe$x[i]*InVivoFree$FracFreeMedia
    simframe_3_conc<- rbind(simframe_3_conc, simframe_3)
  }
  
  #### Fold change - random indiv
  simframeconc$PredFold<- DR_up(x=simframeconc$concentration, y0=1,x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n) 
  ##### Quantiles
  simframe_median_fold<- aggregate(PredFold~concentration, data=simframeconc, quantile, prob=0.5)
  simframe_975_fold<- aggregate(PredFold~concentration, data=simframeconc, quantile, prob=0.975)
  simframe_25_fold<- aggregate(PredFold~concentration, data=simframeconc, quantile, prob=0.025)
  Simulation_Frame_Fold<- simframe_median_fold
  Simulation_Frame_Fold$freeconcentration <- Simulation_Frame_Fold$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_Fold$p97.5<- simframe_975_fold$PredFold
  Simulation_Frame_Fold$p2.5<- simframe_25_fold$PredFold
  names(Simulation_Frame_Fold)[2]<- "p50"
  #### Frac change - random indiv
  simframeconc$PredFrac<- DR_up_Frac(x=simframeconc$concentration,  x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n) 
  ##### Quantiles
  simframe_median_frac<- aggregate(PredFrac~concentration, data=simframeconc, quantile, prob=0.5)
  simframe_975_frac<- aggregate(PredFrac~concentration, data=simframeconc, quantile, prob=0.975)
  simframe_25_frac<- aggregate(PredFrac~concentration, data=simframeconc, quantile, prob=0.025)
  Simulation_Frame_Frac<- simframe_median_frac
  Simulation_Frame_Frac$freeconcentration <- Simulation_Frame_Frac$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_Frac$p97.5<- simframe_975_frac$PredFrac
  Simulation_Frame_Frac$p2.5<- simframe_25_frac$PredFrac
  names(Simulation_Frame_Frac)[2]<- "p50"
  #### Fold change - Population median
  simframe_zero_conc$PredFold<- DR_up(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n) 
  ##### Quantiles
  simframe_zero_median_fold<- aggregate(PredFold~concentration, data=simframe_zero_conc, quantile, prob=0.5)
  simframe_zero_975_fold<- aggregate(PredFold~concentration, data=simframe_zero_conc, quantile, prob=0.975)
  simframe_zero_25_fold<- aggregate(PredFold~concentration, data=simframe_zero_conc, quantile, prob=0.025)
  Simulation_Frame_Zero_Fold<- simframe_zero_median_fold
  Simulation_Frame_Zero_Fold$freeconcentration <- Simulation_Frame_Zero_Fold$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_Zero_Fold$p97.5<- simframe_zero_975_fold$PredFold
  Simulation_Frame_Zero_Fold$p2.5<- simframe_zero_25_fold$PredFold
  names(Simulation_Frame_Zero_Fold)[2]<- "p50"
  #### Frac change - Population median
  simframe_zero_conc$PredFrac<- DR_up_Frac(x=simframe_zero_conc$concentration,  x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n) 
  ##### Quantiles
  simframe_zero_median_frac<- aggregate(PredFrac~concentration, data=simframe_zero_conc, quantile, prob=0.5)
  simframe_zero_975_frac<- aggregate(PredFrac~concentration, data=simframe_zero_conc, quantile, prob=0.975)
  simframe_zero_25_frac<- aggregate(PredFrac~concentration, data=simframe_zero_conc, quantile, prob=0.025)
  Simulation_Frame_Zero_Frac<- simframe_zero_median_frac
  Simulation_Frame_Zero_Frac$freeconcentration <- Simulation_Frame_Zero_Frac$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_Zero_Frac$p97.5<- simframe_zero_975_frac$PredFrac
  Simulation_Frame_Zero_Frac$p2.5<- simframe_zero_25_frac$PredFrac
  names(Simulation_Frame_Zero_Frac)[2]<- "p50"
  
  #### Fold change - Standard donor (27)
  simframe_27_conc$PredFold<- DR_up(x=simframe_27_conc$concentration, y0=1, x0=simframe_27_conc$x0, Emax=simframe_27_conc$Emax, n=simframe_27_conc$n) 
  ##### Quantiles
  simframe_27_median_fold<- aggregate(PredFold~concentration, data=simframe_27_conc, quantile, prob=0.5)
  simframe_27_975_fold<- aggregate(PredFold~concentration, data=simframe_27_conc, quantile, prob=0.975)
  simframe_27_25_fold<- aggregate(PredFold~concentration, data=simframe_27_conc, quantile, prob=0.025)
  Simulation_Frame_27_Fold<- simframe_27_median_fold
  Simulation_Frame_27_Fold$freeconcentration <- Simulation_Frame_27_Fold$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_27_Fold$p97.5<- simframe_27_975_fold$PredFold
  Simulation_Frame_27_Fold$p2.5<- simframe_27_25_fold$PredFold
  names(Simulation_Frame_27_Fold)[2]<- "p50"
  #### Frac change - Standard donor
  simframe_27_conc$PredFrac<- DR_up_Frac(x=simframe_27_conc$concentration,  x0=simframe_27_conc$x0, Emax=simframe_27_conc$Emax, n=simframe_27_conc$n) 
  ##### Quantiles
  simframe_27_median_frac<- aggregate(PredFrac~concentration, data=simframe_27_conc, quantile, prob=0.5)
  simframe_27_975_frac<- aggregate(PredFrac~concentration, data=simframe_27_conc, quantile, prob=0.975)
  simframe_27_25_frac<- aggregate(PredFrac~concentration, data=simframe_27_conc, quantile, prob=0.025)
  Simulation_Frame_27_Frac<- simframe_27_median_frac
  Simulation_Frame_27_Frac$freeconcentration <- Simulation_Frame_27_Frac$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_27_Frac$p97.5<- simframe_27_975_frac$PredFrac
  Simulation_Frame_27_Frac$p2.5<- simframe_27_25_frac$PredFrac
  names(Simulation_Frame_27_Frac)[2]<- "p50"
  
  #### Fold change - Donor 1392
  simframe_22_conc$PredFold<- DR_up(x=simframe_22_conc$concentration, y0=1, x0=simframe_22_conc$x0, Emax=simframe_22_conc$Emax, n=simframe_22_conc$n) 
  ##### Quantiles
  simframe_22_median_fold<- aggregate(PredFold~concentration, data=simframe_22_conc, quantile, prob=0.5)
  simframe_22_975_fold<- aggregate(PredFold~concentration, data=simframe_22_conc, quantile, prob=0.975)
  simframe_22_25_fold<- aggregate(PredFold~concentration, data=simframe_22_conc, quantile, prob=0.025)
  Simulation_Frame_22_Fold<- simframe_22_median_fold
  Simulation_Frame_22_Fold$freeconcentration <- Simulation_Frame_22_Fold$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_22_Fold$p97.5<- simframe_22_975_fold$PredFold
  Simulation_Frame_22_Fold$p2.5<- simframe_22_25_fold$PredFold
  names(Simulation_Frame_22_Fold)[2]<- "p50"
  #### Frac change - Donor 1392
  simframe_22_conc$PredFrac<- DR_up_Frac(x=simframe_22_conc$concentration,  x0=simframe_22_conc$x0, Emax=simframe_22_conc$Emax, n=simframe_22_conc$n) 
  ##### Quantiles
  simframe_22_median_frac<- aggregate(PredFrac~concentration, data=simframe_22_conc, quantile, prob=0.5)
  simframe_22_975_frac<- aggregate(PredFrac~concentration, data=simframe_22_conc, quantile, prob=0.975)
  simframe_22_25_frac<- aggregate(PredFrac~concentration, data=simframe_22_conc, quantile, prob=0.025)
  Simulation_Frame_22_Frac<- simframe_22_median_frac
  Simulation_Frame_22_Frac$freeconcentration <- Simulation_Frame_22_Frac$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_22_Frac$p97.5<- simframe_22_975_frac$PredFrac
  Simulation_Frame_22_Frac$p2.5<- simframe_22_25_frac$PredFrac
  names(Simulation_Frame_22_Frac)[2]<- "p50"
  
  #### Fold change - Donor 1258
  simframe_11_conc$PredFold<- DR_up(x=simframe_11_conc$concentration, y0=1, x0=simframe_11_conc$x0, Emax=simframe_11_conc$Emax, n=simframe_11_conc$n) 
  ##### Quantiles
  simframe_11_median_fold<- aggregate(PredFold~concentration, data=simframe_11_conc, quantile, prob=0.5)
  simframe_11_975_fold<- aggregate(PredFold~concentration, data=simframe_11_conc, quantile, prob=0.975)
  simframe_11_25_fold<- aggregate(PredFold~concentration, data=simframe_11_conc, quantile, prob=0.025)
  Simulation_Frame_11_Fold<- simframe_11_median_fold
  Simulation_Frame_11_Fold$freeconcentration <- Simulation_Frame_11_Fold$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_11_Fold$p97.5<- simframe_11_975_fold$PredFold
  Simulation_Frame_11_Fold$p2.5<- simframe_11_25_fold$PredFold
  names(Simulation_Frame_11_Fold)[2]<- "p50"
  #### Frac change - Donor 1258
  simframe_11_conc$PredFrac<- DR_up_Frac(x=simframe_11_conc$concentration,  x0=simframe_11_conc$x0, Emax=simframe_11_conc$Emax, n=simframe_11_conc$n) 
  ##### Quantiles
  simframe_11_median_frac<- aggregate(PredFrac~concentration, data=simframe_11_conc, quantile, prob=0.5)
  simframe_11_975_frac<- aggregate(PredFrac~concentration, data=simframe_11_conc, quantile, prob=0.975)
  simframe_11_25_frac<- aggregate(PredFrac~concentration, data=simframe_11_conc, quantile, prob=0.025)
  Simulation_Frame_11_Frac<- simframe_11_median_frac
  Simulation_Frame_11_Frac$freeconcentration <- Simulation_Frame_11_Frac$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_11_Frac$p97.5<- simframe_11_975_frac$PredFrac
  Simulation_Frame_11_Frac$p2.5<- simframe_11_25_frac$PredFrac
  names(Simulation_Frame_11_Frac)[2]<- "p50"
  
  #### Fold change - Donor 1083
  simframe_6_conc$PredFold<- DR_up(x=simframe_6_conc$concentration, y0=1, x0=simframe_6_conc$x0, Emax=simframe_6_conc$Emax, n=simframe_6_conc$n) 
  ##### Quantiles
  simframe_6_median_fold<- aggregate(PredFold~concentration, data=simframe_6_conc, quantile, prob=0.5)
  simframe_6_975_fold<- aggregate(PredFold~concentration, data=simframe_6_conc, quantile, prob=0.975)
  simframe_6_25_fold<- aggregate(PredFold~concentration, data=simframe_6_conc, quantile, prob=0.025)
  Simulation_Frame_6_Fold<- simframe_6_median_fold
  Simulation_Frame_6_Fold$freeconcentration <- Simulation_Frame_6_Fold$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_6_Fold$p97.5<- simframe_6_975_fold$PredFold
  Simulation_Frame_6_Fold$p2.5<- simframe_6_25_fold$PredFold
  names(Simulation_Frame_6_Fold)[2]<- "p50"
  #### Frac change - Donor 1083
  simframe_6_conc$PredFrac<- DR_up_Frac(x=simframe_6_conc$concentration,  x0=simframe_6_conc$x0, Emax=simframe_6_conc$Emax, n=simframe_6_conc$n) 
  ##### Quantiles
  simframe_6_median_frac<- aggregate(PredFrac~concentration, data=simframe_6_conc, quantile, prob=0.5)
  simframe_6_975_frac<- aggregate(PredFrac~concentration, data=simframe_6_conc, quantile, prob=0.975)
  simframe_6_25_frac<- aggregate(PredFrac~concentration, data=simframe_6_conc, quantile, prob=0.025)
  Simulation_Frame_6_Frac<- simframe_6_median_frac
  Simulation_Frame_6_Frac$freeconcentration <- Simulation_Frame_6_Frac$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_6_Frac$p97.5<- simframe_6_975_frac$PredFrac
  Simulation_Frame_6_Frac$p2.5<- simframe_6_25_frac$PredFrac
  names(Simulation_Frame_6_Frac)[2]<- "p50"
  
  #### Fold change - Donor 1070
  simframe_3_conc$PredFold<- DR_up(x=simframe_3_conc$concentration, y0=1, x0=simframe_3_conc$x0, Emax=simframe_3_conc$Emax, n=simframe_3_conc$n) 
  ##### Quantiles
  simframe_3_median_fold<- aggregate(PredFold~concentration, data=simframe_3_conc, quantile, prob=0.5)
  simframe_3_975_fold<- aggregate(PredFold~concentration, data=simframe_3_conc, quantile, prob=0.975)
  simframe_3_25_fold<- aggregate(PredFold~concentration, data=simframe_3_conc, quantile, prob=0.025)
  Simulation_Frame_3_Fold<- simframe_3_median_fold
  Simulation_Frame_3_Fold$freeconcentration <- Simulation_Frame_3_Fold$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_3_Fold$p97.5<- simframe_3_975_fold$PredFold
  Simulation_Frame_3_Fold$p2.5<- simframe_3_25_fold$PredFold
  names(Simulation_Frame_3_Fold)[2]<- "p50"
  #### Frac change - Donor 1070
  simframe_3_conc$PredFrac<- DR_up_Frac(x=simframe_3_conc$concentration,  x0=simframe_3_conc$x0, Emax=simframe_3_conc$Emax, n=simframe_3_conc$n) 
  ##### Quantiles
  simframe_3_median_frac<- aggregate(PredFrac~concentration, data=simframe_3_conc, quantile, prob=0.5)
  simframe_3_975_frac<- aggregate(PredFrac~concentration, data=simframe_3_conc, quantile, prob=0.975)
  simframe_3_25_frac<- aggregate(PredFrac~concentration, data=simframe_3_conc, quantile, prob=0.025)
  Simulation_Frame_3_Frac<- simframe_3_median_frac
  Simulation_Frame_3_Frac$freeconcentration <- Simulation_Frame_3_Frac$concentration*InVivoFree$FracFreeMedia
  Simulation_Frame_3_Frac$p97.5<- simframe_3_975_frac$PredFrac
  Simulation_Frame_3_Frac$p2.5<- simframe_3_25_frac$PredFrac
  names(Simulation_Frame_3_Frac)[2]<- "p50"
  
  
  tmplist$Simulation_Frame_Fold <- Simulation_Frame_Fold
  tmplist$Simulation_Frame_Frac <- Simulation_Frame_Frac
  tmplist$Simulation_Frame_Zero_Fold <- Simulation_Frame_Zero_Fold
  tmplist$Simulation_Frame_Zero_Frac <- Simulation_Frame_Zero_Frac
  tmplist$Simulation_Frame_27_Fold <- Simulation_Frame_27_Fold
  tmplist$Simulation_Frame_27_Frac <- Simulation_Frame_27_Frac
  
  tmplist$Simulation_Frame_22_Fold <- Simulation_Frame_22_Fold
  tmplist$Simulation_Frame_22_Frac <- Simulation_Frame_22_Frac
  
  tmplist$Simulation_Frame_11_Fold <- Simulation_Frame_11_Fold
  tmplist$Simulation_Frame_11_Frac <- Simulation_Frame_11_Frac
  
  tmplist$Simulation_Frame_6_Fold <- Simulation_Frame_6_Fold
  tmplist$Simulation_Frame_6_Frac <- Simulation_Frame_6_Frac
  
  tmplist$Simulation_Frame_3_Fold <- Simulation_Frame_3_Fold
  tmplist$Simulation_Frame_3_Frac <- Simulation_Frame_3_Frac
  
  data.list[[j]]<-tmplist
  lapply(data.list[[1]],head)
  rm(list = files)
}

save(data.list,chemmap,file="InVivoInVitrodat-ExperimentalFreeFrac_43.Rdata")
beep(sound=3)