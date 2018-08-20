DR_up<-function(x=1, y0=1, x0=1, Emax=1, n=1){
  
  y0*(1+((x/x0)^n)/(1+((x/x0)^n)*(1/Emax))) 
  
}

DR_up_Frac<- function(x=1, y0=1, x0=1, Emax=1, n=1){
  (y0*(1+((x/x0)^n)/(1+((x/x0)^n)*(1/Emax))))-1 
}

DR_dn<-function(x=1, y0=1, x0=1, Emax=1, n=1){
  
  y0*(1-((x/x0)^n)/(1+((x/x0)^n)/(Emax))) 
  
}

DR_dn_Frac<- function(x=1, y0=1, x0=1, Emax=1, n=1){
  (y0*(1-((x/x0)^n)/(1+((x/x0)^n)/(Emax))))-1 
}

###Master Pos function

DR_master_up<- function(chemical.name){

  setwd(paste("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\", chemical.name, "\\decay_rise_up", sep=""))
  
  ### load in CSV files, R dat file,  set up concentration frame
  
  files = list.files(pattern="*samples_*")
  
  for (l in 1:length(files)) assign(files[l], read.csv(files[l]))
  
  file.sources=list.files(pattern="*dat.R")
  
  sapply(file.sources,source,.GlobalEnv)
  
  concframe<-data.frame(x=10^(seq(-3,2,0.05)))
  
  ### extract median parameter data from csv files, create new frame with calculated parameter data, do the same when z scores and SD =0 
  
  m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))
  
  m_y0$z_score_sim<- rnorm(8000,0,1)
  
  simframe<- as.data.frame(exp(m_y0$m_y0+(m_y0$sd_y0*m_y0$z_score_sim)))
  
  rownames(simframe)<- c(1:8000)
  
  colnames(simframe)<- "sim_y0"
  
  m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))
  
  m_x0$z_score_sim<- rnorm(8000,0,1)
  
  simframe$sim_x0<-exp(m_x0$m_x0+(m_x0$sd_x0*m_x0$z_score_sim))
  
  m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))
  
  m_Emax$z_score_sim<- rnorm(8000,0,1)
  
  simframe$sim_Emax<-exp(m_Emax$m_Emax+(m_Emax$sd_Emax*m_Emax$z_score_sim))
  
  m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))
  
  m_n$z_score_sim<- rnorm(8000,0,1)
  
  simframe$sim_n<-exp(m_n$m_n+(m_n$sd_n*m_n$z_score_sim))
  
  simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), Emax=(exp(m_Emax$m_Emax)), n=(exp(m_n$m_n)))
  
  ### Load in appropriate control CSV and Parameter CSV
  
  InVivoPos<- read.csv("InVivoPred.pos.v2.csv", header=TRUE)
  
  InVivoParm<- read.csv("QTposnegparameters.v2.csv", header=TRUE)
  
  InVivo<-subset(InVivoPos, Chemical.name== chemical.name )
  
  InVivoFree<- subset(InVivoParm, Chemical== chemical.name )
  
  ### Set up working frame with concentrations, transform concentrations into freeconc
  
  simframeconc<- data.frame()
  
  simframe_zero_conc<- data.frame()
  
  for (i in 1:nrow(concframe)){
    
    simframe$concentration<-concframe$x[i]
    
    simframeconc<- rbind(simframeconc, simframe)
  }
  
  simframeconc$freeconcentration<- (simframeconc$concentration)*InVivoFree$FracFreeMedia
  
  for (i in 1:nrow(concframe)){
    
    simframe_zero$concentration<-concframe$x[i]
    
    simframe_zero_conc<- rbind(simframe_zero_conc, simframe_zero)
  }
  
  simframe_zero_conc$freeconcentration<- (simframe_zero_conc$concentration)*InVivoFree$FracFreeMedia
  
  
  ##create final frames for plotting w/ quantile values
  
  simframeconc$PredFold<- DR_up(x=simframeconc$concentration, y0=1, x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n)
  
  simframe_median_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframeconc, quantile, prob=0.5)
  
  simframe_975_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframeconc, quantile, prob=0.975)
  
  simframe_25_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframeconc, quantile, prob=0.025)
  
  Simulation_Frame_Fold<- simframe_median_fold
  
  Simulation_Frame_Fold$p97.5<- simframe_975_fold$PredFold
  
  Simulation_Frame_Fold$p2.5<- simframe_25_fold$PredFold
  
  names(Simulation_Frame_Fold)[3]<- "p50"
  
  ##########
  
  simframeconc$PredFrac<- DR_up_Frac(x=simframeconc$concentration, y0=1, x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n) 
  
  simframe_median_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframeconc, quantile, prob=0.5)
  
  simframe_975_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframeconc, quantile, prob=0.975)
  
  simframe_25_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframeconc, quantile, prob=0.025)
  
  Simulation_Frame_Frac<- simframe_median_frac
  
  Simulation_Frame_Frac$p97.5<- simframe_975_frac$PredFrac
  
  Simulation_Frame_Frac$p2.5<- simframe_25_frac$PredFrac
  
  names(Simulation_Frame_Frac)[3]<- "p50"
  
  ##########
  
  simframe_zero_conc$PredFold<- DR_up(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n)
  
  simframe_zero_median_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.5)
  
  simframe_zero_975_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.975)
  
  simframe_zero_25_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.025)
  
  Simulation_Frame_Zero_Fold<- simframe_zero_median_fold
  
  Simulation_Frame_Zero_Fold$p97.5<- simframe_zero_975_fold$PredFold
  
  Simulation_Frame_Zero_Fold$p2.5<- simframe_zero_25_fold$PredFold
  
  names(Simulation_Frame_Zero_Fold)[3]<- "p50"
  
  ######
  
  simframe_zero_conc$PredFrac<- DR_up_Frac(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n)
  
  simframe_zero_median_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.5)
  
  simframe_zero_975_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.975)
  
  simframe_zero_25_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.025)
  
  Simulation_Frame_Zero_Frac<- simframe_zero_median_frac
  
  Simulation_Frame_Zero_Frac$p97.5<- simframe_zero_975_frac$PredFrac
  
  Simulation_Frame_Zero_Frac$p2.5<- simframe_zero_25_frac$PredFrac
  
  names(Simulation_Frame_Zero_Frac)[3]<- "p50"
  
  
  ###output workind dataframes as CSV files
  
  #set wd for outputs
  
  setwd(paste("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\", chemical.name, "\\decay_rise_up", "\\invivo_invitro", sep=""))
  
  write.csv(Simulation_Frame_Fold, "InVivo_Frame_Fold_rescaled.csv")
  
  write.csv(Simulation_Frame_Zero_Fold, "Invivo_Frame_Zero_Fold_rescaled.csv")
  
  write.csv(Simulation_Frame_Frac, "Invivo_Frame_Frac_rescaled.csv")
  
  write.csv(Simulation_Frame_Zero_Frac, "Invivo_Frame_Zero_Frac_rescaled.csv")
  
  ###plot
  
  plotname_fold_log<-(paste(chemical.name, "In Vitro vs. In Vivo Fold Change Peak Decay Rise Ratio Up (Log Scale)", sep=" "))
  
  plotname_fold<-(paste(chemical.name, "In Vitro vs. In Vivo Fold Change Peak Decay Rise Ratio Up", sep=" "))
   
  plotname_frac_log<-(paste(chemical.name, "In Vitro vs. In Vivo Frac Change Peak Decay Rise Ratio Up (Log Scale)", sep=" "))
  
  plotname_frac<-(paste(chemical.name, "In Vitro vs. In Vivo Frac Change Peak Decay Rise Ratio Up (Log Scale)", sep=" "))
  
  #Fold Change plots
  
  ggplot(Simulation_Frame_Fold,aes(x=concentration, y=p50))+
    
    geom_line(data=InVivo, aes(x=xfree, y=PredFoldChange), color="blue", size=2)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p50), color="red", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line()+
    
    geom_line(aes(x=freeconcentration, y= p97.5), linetype="dotted")+
    
    geom_line(aes(x=freeconcentration, y= p2.5), linetype="dotted")+
    
    geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
    
    xlab("Free Concentration (uM)")+
    
    ylab("Decay Rise Ratio (Fold Change)")+
    
    scale_x_log10()+
    
    scale_y_log10()+
    
    ggtitle("Sotalol In Vitro vs. In Vivo Fold Change Peak Decay Rise Ratio Up (Log Scale)")+
    
    theme_minimal()
  
  ggsave("Invivo_Dose_Response_Log_Fold_rescaled.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)
  
  ggsave("InVivo_Dose_Response_Log_Fold_rescaled.png", plot=last_plot(), device="png", height= 4.25, width=4.25)
  
  ggplot(Simulation_Frame_Fold, aes(x=concentration, y=p50))+
    
    geom_line(data=InVivo, aes(x=xfree, y=PredFoldChange), color="blue", size=2)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p50), color="red", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line()+
    
    geom_line(aes(x=freeconcentration, y= p97.5), linetype="dotted")+
    
    geom_line(aes(x=freeconcentration, y= p2.5), linetype="dotted")+
    
    geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
    
    xlab("Free Concentration (uM)")+
    
    ylab("Decay Rise Ratio (Fold Change)")+
    
    scale_x_log10()+
    
    ggtitle("Sotalol In Vivo vs. In Vitro Fold Change Peak Decay Rise Ratio Up")+
    
    theme_minimal()
  
  ggsave("Invivo_Dose_Response_Fold_rescaled.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)
  
  ggsave("InVivo_Dose_Response_Fold_rescaled.png", plot=last_plot(), device="png", height= 4.25, width=4.25)
  
  # Frac Change Plots
  
  ggplot(Simulation_Frame_Frac, aes(x=concentration, y=p50))+
    
    geom_line(data=InVivo, aes(x=xfree, y=PredFracChange), color="blue", size=2)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p50), color="red", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line()+
    
    geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
    
    geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
    
    geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
    
    xlab("Log Free Concentration (uM)")+
    
    ylab("Decay Rise Ratio (Log Frac Change)")+
    
    scale_x_log10()+
    
    scale_y_log10(breaks=c(1,1e-02,1e-04,1e-06,1e-08))+
    
    ggtitle("Sotalol  In Vivo vs. In Vitro Fractional Change Peak Decay Rise Ratio Up (Log Scale)")+
    
    theme_minimal()
  
  ggsave("Invivo_Dose_Response_Log_Frac_rescaled.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)
  
  ggsave("Invivo_Dose_Response_Log_Frac_rescaled.png", plot=last_plot(), device="png", height= 4.25, width= 4.25)
  
  ggplot(Simulation_Frame_Frac, aes(x=concentration, y=p50))+
    
    geom_line(data=InVivo, aes(x=xfree, y=PredFracChange), color="blue", size=2)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p50), color="red", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line()+
    
    geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
    
    geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
    
    geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
    
    xlab("Free Concentration (uM)")+
    
    ylab("Decay Rise Ratio (Frac Change)")+
    
    scale_x_log10()+
    
    ggtitle("Sotalol  In Vivo vs. In Vitro Fractional Change Peak Decay Rise Ratio Up")+
    
    theme_minimal()
  
  ggsave("Invivo_Dose_Response_Frac_rescaled.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)
  
  ggsave("Invivo_Dose_Response_Frac_rescaled.png", plot=last_plot(), device="png", height= 4.25, width= 4.25)
  
  
  rm(list=files)
}

DR_master_dn<- function(chemical.name){

  setwd(paste("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\", chemical.name, "\\decay_rise_up", sep=""))
  
  ### load in CSV files, R dat file,  set up concentration frame
  
  files = list.files(pattern="*samples_*")
  
  for (l in 1:length(files)) assign(files[l], read.csv(files[l]))
  
  file.sources=list.files(pattern="*dat.R")
  
  sapply(file.sources,source,.GlobalEnv)
  
  concframe<-data.frame(x=10^(seq(-3,2,0.05)))
  
  ### extract median parameter data from csv files, create new frame with calculated parameter data, do the same when z scores and SD =0 
  
  m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))
  
  m_y0$z_score_sim<- rnorm(8000,0,1)
  
  simframe<- as.data.frame(exp(m_y0$m_y0+(m_y0$sd_y0*m_y0$z_score_sim)))
  
  rownames(simframe)<- c(1:8000)
  
  colnames(simframe)<- "sim_y0"
  
  m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))
  
  m_x0$z_score_sim<- rnorm(8000,0,1)
  
  simframe$sim_x0<-exp(m_x0$m_x0+(m_x0$sd_x0*m_x0$z_score_sim))
  
  m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))
  
  m_Emax$z_score_sim<- rnorm(8000,0,1)
  
  simframe$sim_Emax<-exp(m_Emax$m_Emax+(m_Emax$sd_Emax*m_Emax$z_score_sim))
  
  m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))
  
  m_n$z_score_sim<- rnorm(8000,0,1)
  
  simframe$sim_n<-exp(m_n$m_n+(m_n$sd_n*m_n$z_score_sim))
  
  simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), Emax=(exp(m_Emax$m_Emax)), n=(exp(m_n$m_n)))
  
  ### Load in appropriate control CSV and Parameter CSV
  
  InVivoPos<- read.csv("InVivoPred.pos.v2.csv", header=TRUE)
  
  InVivoParm<- read.csv("QTposnegparameters.v2.csv", header=TRUE)
  
  InVivo<-subset(InVivoPos, Chemical.name== chemical.name )
  
  InVivoFree<- subset(InVivoParm, Chemical== chemical.name )
  
  ### Set up working frame with concentrations, transform concentrations into freeconc
  
  simframeconc<- data.frame()
  
  simframe_zero_conc<- data.frame()
  
  for (i in 1:nrow(concframe)){
    
    simframe$concentration<-concframe$x[i]
    
    simframeconc<- rbind(simframeconc, simframe)
  }
  
  simframeconc$freeconcentration<- (simframeconc$concentration)*InVivoFree$FracFreeMedia
  
  for (i in 1:nrow(concframe)){
    
    simframe_zero$concentration<-concframe$x[i]
    
    simframe_zero_conc<- rbind(simframe_zero_conc, simframe_zero)
  }
  
  simframe_zero_conc$freeconcentration<- (simframe_zero_conc$concentration)*InVivoFree$FracFreeMedia
  
  
  ##create final frames for plotting w/ quantile values
  
  simframeconc$PredFold<- DR_dn(x=simframeconc$concentration, y0=1, x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n)
  
  simframe_median_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframeconc, quantile, prob=0.5)
  
  simframe_975_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframeconc, quantile, prob=0.975)
  
  simframe_25_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframeconc, quantile, prob=0.025)
  
  Simulation_Frame_Fold<- simframe_median_fold
  
  Simulation_Frame_Fold$p97.5<- simframe_975_fold$PredFold
  
  Simulation_Frame_Fold$p2.5<- simframe_25_fold$PredFold
  
  names(Simulation_Frame_Fold)[3]<- "p50"
  
  ##########
  
  simframeconc$PredFrac<- DR_dn_Frac(x=simframeconc$concentration, y0=1, x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n)
  
  simframe_median_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframeconc, quantile, prob=0.5)
  
  simframe_975_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframeconc, quantile, prob=0.975)
  
  simframe_25_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframeconc, quantile, prob=0.025)
  
  Simulation_Frame_Frac<- simframe_median_frac
  
  Simulation_Frame_Frac$p97.5<- simframe_975_frac$PredFrac
  
  Simulation_Frame_Frac$p2.5<- simframe_25_frac$PredFrac
  
  names(Simulation_Frame_Frac)[3]<- "p50"
  
  ##########
  
  simframe_zero_conc$PredFold<- DR_dn(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n)
  
  simframe_zero_median_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.5)
  
  simframe_zero_975_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.975)
  
  simframe_zero_25_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.025)
  
  Simulation_Frame_Zero_Fold<- simframe_zero_median_fold
  
  Simulation_Frame_Zero_Fold$p97.5<- simframe_zero_975_fold$PredFold
  
  Simulation_Frame_Zero_Fold$p2.5<- simframe_zero_25_fold$PredFold
  
  names(Simulation_Frame_Zero_Fold)[3]<- "p50"
  
  ######
  
  simframe_zero_conc$PredFrac<- DR_dn_Frac(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n)
  
  simframe_zero_median_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.5)
  
  simframe_zero_975_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.975)
  
  simframe_zero_25_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.025)
  
  Simulation_Frame_Zero_Frac<- simframe_zero_median_frac
  
  Simulation_Frame_Zero_Frac$p97.5<- simframe_zero_975_frac$PredFrac
  
  Simulation_Frame_Zero_Frac$p2.5<- simframe_zero_25_frac$PredFrac
  
  names(Simulation_Frame_Zero_Frac)[3]<- "p50"
  
  
  ###output workind dataframes as CSV files
  
  #set wd for outputs
  
  setwd(paste("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\", chemical.name, "\\decay_rise_up", "\\invivo_invitro", sep=""))
  
  write.csv(Simulation_Frame_Fold, "InVivo_Frame_Fold_rescaled.csv")
  
  write.csv(Simulation_Frame_Zero_Fold, "Invivo_Frame_Zero_Fold_rescaled.csv")
  
  write.csv(Simulation_Frame_Frac, "Invivo_Frame_Frac_rescaled.csv")
  
  write.csv(Simulation_Frame_Zero_Frac, "Invivo_Frame_Zero_Frac_rescaled.csv")
  
  ###plot
  
  plotname_fold_log<-(paste(chemical.name, "In Vitro vs. In Vivo Fold Change Peak Decay Rise Ratio Up (Log Scale)", sep=" "))
  
  plotname_fold<-(paste(chemical.name, "In Vitro vs. In Vivo Fold Change Peak Decay Rise Ratio Up", sep=" "))
   
  plotname_frac_log<-(paste(chemical.name, "In Vitro vs. In Vivo Frac Change Peak Decay Rise Ratio Up (Log Scale)", sep=" "))
  
  plotname_frac<-(paste(chemical.name, "In Vitro vs. In Vivo Frac Change Peak Decay Rise Ratio Up (Log Scale)", sep=" "))
  
  #Fold Change plots
  
  ggplot(Simulation_Frame_Fold,aes(x=concentration, y=p50))+
    
    geom_line(data=InVivo, aes(x=xfree, y=PredFoldChange), color="blue", size=2)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p50), color="red", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line()+
    
    geom_line(aes(x=freeconcentration, y= p97.5), linetype="dotted")+
    
    geom_line(aes(x=freeconcentration, y= p2.5), linetype="dotted")+
    
    geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
    
    xlab("Free Concentration (uM)")+
    
    ylab("Decay Rise Ratio (Fold Change)")+
    
    scale_x_log10()+
    
    scale_y_log10()+
    
    ggtitle("Sotalol In Vitro vs. In Vivo Fold Change Peak Decay Rise Ratio Up (Log Scale)")+
    
    theme_minimal()
  
  ggsave("Invivo_Dose_Response_Log_Fold_rescaled.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)
  
  ggsave("InVivo_Dose_Response_Log_Fold_rescaled.png", plot=last_plot(), device="png", height= 4.25, width=4.25)
  
  ggplot(Simulation_Frame_Fold, aes(x=concentration, y=p50))+
    
    geom_line(data=InVivo, aes(x=xfree, y=PredFoldChange), color="blue", size=2)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p50), color="red", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line()+
    
    geom_line(aes(x=freeconcentration, y= p97.5), linetype="dotted")+
    
    geom_line(aes(x=freeconcentration, y= p2.5), linetype="dotted")+
    
    geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
    
    xlab("Free Concentration (uM)")+
    
    ylab("Decay Rise Ratio (Fold Change)")+
    
    scale_x_log10()+
    
    ggtitle("Sotalol In Vivo vs. In Vitro Fold Change Peak Decay Rise Ratio Up")+
    
    theme_minimal()
  
  ggsave("Invivo_Dose_Response_Fold_rescaled.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)
  
  ggsave("InVivo_Dose_Response_Fold_rescaled.png", plot=last_plot(), device="png", height= 4.25, width=4.25)
  
  # Frac Change Plots
  
  ggplot(Simulation_Frame_Frac, aes(x=concentration, y=p50))+
    
    geom_line(data=InVivo, aes(x=xfree, y=PredFracChange), color="blue", size=2)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p50), color="red", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line()+
    
    geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
    
    geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
    
    geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
    
    xlab("Log Free Concentration (uM)")+
    
    ylab("Decay Rise Ratio (Log Frac Change)")+
    
    scale_x_log10()+
    
    scale_y_log10(breaks=c(1,1e-02,1e-04,1e-06,1e-08))+
    
    ggtitle("Sotalol  In Vivo vs. In Vitro Fractional Change Peak Decay Rise Ratio Up (Log Scale)")+
    
    theme_minimal()
  
  ggsave("Invivo_Dose_Response_Log_Frac_rescaled.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)
  
  ggsave("Invivo_Dose_Response_Log_Frac_rescaled.png", plot=last_plot(), device="png", height= 4.25, width= 4.25)
  
  ggplot(Simulation_Frame_Frac, aes(x=concentration, y=p50))+
    
    geom_line(data=InVivo, aes(x=xfree, y=PredFracChange), color="blue", size=2)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p50), color="red", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line()+
    
    geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
    
    geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
    
    geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
    
    xlab("Free Concentration (uM)")+
    
    ylab("Decay Rise Ratio (Frac Change)")+
    
    scale_x_log10()+
    
    ggtitle("Sotalol  In Vivo vs. In Vitro Fractional Change Peak Decay Rise Ratio Up")+
    
    theme_minimal()
  
  ggsave("Invivo_Dose_Response_Frac_rescaled.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)
  
  ggsave("Invivo_Dose_Response_Frac_rescaled.png", plot=last_plot(), device="png", height= 4.25, width= 4.25)
  
  
  rm(list=files)
}

DR_master_two_model_up<- function(chemical.name){

  setwd(paste("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\", chemical.name, "\\decay_rise_up", sep=""))
  
  ### load in CSV files, R dat file,  set up concentration frame
  
  files = list.files(pattern="*samples_*")
  
  for (l in 1:length(files)) assign(files[l], read.csv(files[l]))
  
  file.sources=list.files(pattern="*dat.R")
  
  sapply(file.sources,source,.GlobalEnv)
  
  concframe<-data.frame(x=10^(seq(-3,2,0.05)))
  
  ### extract median parameter data from csv files, create new frame with calculated parameter data, do the same when z scores and SD =0 
  
  m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))
  
  m_y0$z_score_sim<- rnorm(8000,0,1)
  
  simframe<- as.data.frame(exp(m_y0$m_y0+(m_y0$sd_y0*m_y0$z_score_sim)))
  
  rownames(simframe)<- c(1:8000)
  
  colnames(simframe)<- "sim_y0"
  
  m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))
  
  m_x0$z_score_sim<- rnorm(8000,0,1)
  
  simframe$sim_x0<-exp(m_x0$m_x0+(m_x0$sd_x0*m_x0$z_score_sim))
  
  m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))
  
  m_Emax$z_score_sim<- rnorm(8000,0,1)
  
  simframe$sim_Emax<-exp(m_Emax$m_Emax+(m_Emax$sd_Emax*m_Emax$z_score_sim))
  
  m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))
  
  m_n$z_score_sim<- rnorm(8000,0,1)
  
  simframe$sim_n<-exp(m_n$m_n+(m_n$sd_n*m_n$z_score_sim))
  
  simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), Emax=(exp(m_Emax$m_Emax)), n=(exp(m_n$m_n)))
  
  ### Load in appropriate control CSV and Parameter CSV
  
  InVivoPos<- read.csv("InVivoPred.pos.v2.csv", header=TRUE)
  
  InVivoParm<- read.csv("QTposnegparameters.v2.csv", header=TRUE)
  
  InVivo<-subset(InVivoPos, Chemical.name== chemical.name )
  
  InVivoFree<- subset(InVivoParm, Chemical== chemical.name )
  
  ### Set up working frame with concentrations, transform concentrations into freeconc
  
  simframeconc<- data.frame()
  
  simframe_zero_conc<- data.frame()
  
  for (i in 1:nrow(concframe)){
    
    simframe$concentration<-concframe$x[i]
    
    simframeconc<- rbind(simframeconc, simframe)
  }
  
  simframeconc$freeconcentration<- (simframeconc$concentration)*InVivoFree$FracFreeMedia
  
  for (i in 1:nrow(concframe)){
    
    simframe_zero$concentration<-concframe$x[i]
    
    simframe_zero_conc<- rbind(simframe_zero_conc, simframe_zero)
  }
  
  simframe_zero_conc$freeconcentration<- (simframe_zero_conc$concentration)*InVivoFree$FracFreeMedia
  
  
  ##create final frames for plotting w/ quantile values
  
  simframeconc$PredFold<- DR_up(x=simframeconc$concentration, y0=1, x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n)
  
  simframe_median_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframeconc, quantile, prob=0.5)
  
  simframe_975_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframeconc, quantile, prob=0.975)
  
  simframe_25_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframeconc, quantile, prob=0.025)
  
  Simulation_Frame_Fold<- simframe_median_fold
  
  Simulation_Frame_Fold$p97.5<- simframe_975_fold$PredFold
  
  Simulation_Frame_Fold$p2.5<- simframe_25_fold$PredFold
  
  names(Simulation_Frame_Fold)[3]<- "p50"
  
  ##########
  
  simframeconc$PredFrac<- DR_up_Frac(x=simframeconc$concentration, y0=1, x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n)
  
  simframe_median_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframeconc, quantile, prob=0.5)
  
  simframe_975_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframeconc, quantile, prob=0.975)
  
  simframe_25_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframeconc, quantile, prob=0.025)
  
  Simulation_Frame_Frac<- simframe_median_frac
  
  Simulation_Frame_Frac$p97.5<- simframe_975_frac$PredFrac
  
  Simulation_Frame_Frac$p2.5<- simframe_25_frac$PredFrac
  
  names(Simulation_Frame_Frac)[3]<- "p50"
  
  ##########
  
  simframe_zero_conc$PredFold<- DR_up(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n)
  
  simframe_zero_median_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.5)
  
  simframe_zero_975_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.975)
  
  simframe_zero_25_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.025)
  
  Simulation_Frame_Zero_Fold<- simframe_zero_median_fold
  
  Simulation_Frame_Zero_Fold$p97.5<- simframe_zero_975_fold$PredFold
  
  Simulation_Frame_Zero_Fold$p2.5<- simframe_zero_25_fold$PredFold
  
  names(Simulation_Frame_Zero_Fold)[3]<- "p50"
  
  ######
  
  simframe_zero_conc$PredFrac<- DR_up_Frac(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n)
  
  simframe_zero_median_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.5)
  
  simframe_zero_975_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.975)
  
  simframe_zero_25_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.025)
  
  Simulation_Frame_Zero_Frac<- simframe_zero_median_frac
  
  Simulation_Frame_Zero_Frac$p97.5<- simframe_zero_975_frac$PredFrac
  
  Simulation_Frame_Zero_Frac$p2.5<- simframe_zero_25_frac$PredFrac
  
  names(Simulation_Frame_Zero_Frac)[3]<- "p50"
  
  
  ###output workind dataframes as CSV files
  
  #set wd for outputs
  
  setwd(paste("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\", chemical.name, "\\decay_rise_up", "\\invivo_invitro", sep=""))
  
  write.csv(Simulation_Frame_Fold, "InVivo_Frame_Fold_rescaled.csv")
  
  write.csv(Simulation_Frame_Zero_Fold, "Invivo_Frame_Zero_Fold_rescaled.csv")
  
  write.csv(Simulation_Frame_Frac, "Invivo_Frame_Frac_rescaled.csv")
  
  write.csv(Simulation_Frame_Zero_Frac, "Invivo_Frame_Zero_Frac_rescaled.csv")
  
  ###plot
  
  plotname_fold_log<-(paste(chemical.name, "In Vitro vs. In Vivo Fold Change Peak Decay Rise Ratio Up (Log Scale)", sep=" "))
  
  plotname_fold<-(paste(chemical.name, "In Vitro vs. In Vivo Fold Change Peak Decay Rise Ratio Up", sep=" "))
   
  plotname_frac_log<-(paste(chemical.name, "In Vitro vs. In Vivo Frac Change Peak Decay Rise Ratio Up (Log Scale)", sep=" "))
  
  plotname_frac<-(paste(chemical.name, "In Vitro vs. In Vivo Frac Change Peak Decay Rise Ratio Up (Log Scale)", sep=" "))

	#Fold Change plots

  ggplot(Simulation_Frame_Fold,aes(x=concentration, y=p50))+
  
  geom_line(data=InVivo, aes(x=xfree,y=PredFracChange, color=Model), size=2, show.legend = FALSE)+
  
  scale_color_manual(values=c("dodgerblue4", "dodgerblue"))+
  
  geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p50), color="red", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line()+
  
  geom_line(aes(x=freeconcentration, y= p97.5), linetype="dotted")+
  
  geom_line(aes(x=freeconcentration, y= p2.5), linetype="dotted")+
  
  geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
  
  xlab("Free Concentration (uM)")+
  
  ylab("Decay Rise Ratio (Fold Change)")+
  
  scale_x_log10()+
  
  scale_y_log10()+
  
  theme(legend.position = "none")+
  
  ggtitle("N-acetylprocainamide   In Vitro vs. In Vivo Fold Change Peak Decay Rise Ratio Up (Log Scale)")+
  
  theme_minimal()

ggsave("Invivo_Dose_Response_Log_Fold.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)

ggsave("InVivo_Dose_Response_Log_Fold.png", plot=last_plot(), device="png", height= 4.25, width=4.25)

ggplot(Simulation_Frame_Fold, aes(x=concentration, y=p50))+
  
  geom_line(data=InVivo, aes(x=xfree,y=PredFracChange, color=Model), size=2, show.legend = FALSE)+
  
  scale_color_manual(values=c("dodgerblue4", "dodgerblue"))+
  
  geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p50), color="red", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line()+
  
  geom_line(aes(x=freeconcentration, y= p97.5), linetype="dotted")+
  
  geom_line(aes(x=freeconcentration, y= p2.5), linetype="dotted")+
  
  geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
  
  xlab("Free Concentration (uM)")+
  
  ylab("Decay Rise Ratio (Fold Change)")+
  
  scale_x_log10()+
  
  theme(legend.position = "none")+
  
  ggtitle("N-acetylprocainamide   In Vivo vs. In Vitro Fold Change Peak Decay Rise Ratio Up")+
  
  theme(legend.position = "none")+
  
  theme_minimal()

ggsave("Invivo_Dose_Response_Fold.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)

ggsave("InVivo_Dose_Response_Fold.png", plot=last_plot(), device="png", height= 4.25, width=4.25)

# Frac Change Plots

ggplot(Simulation_Frame_Frac, aes(x=concentration, y=p50))+
  
  geom_line(data=InVivo, aes(x=xfree,y=PredFracChange, color=Model), size=2, show.legend = FALSE)+
  
  scale_color_manual(values=c("dodgerblue4", "dodgerblue"))+
  
  geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p50), color="red", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line()+
  
  geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
  
  geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
  
  geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
  
  xlab(" Log Free Concentration (uM)")+
  
  ylab("Decay Rise Ratio ( Log Frac Change)")+
  
  scale_x_log10()+
  
  scale_y_log10(breaks=c(1, 1e-02, 1e-04, 1e-06,1e-08))+
  
  ggtitle("N-acetylprocainamide  In Vivo vs. In Vitro Fractional Change Peak Decay Rise Ratio Up (Log Scale)")+
  
  theme_minimal()

ggsave("Invivo_Dose_Response_Log_Frac.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)

ggsave("Invivo_Dose_Response_Log_Frac.png", plot=last_plot(), device="png", height= 4.25, width= 4.25)

ggplot(Simulation_Frame_Frac, aes(x=concentration, y=p50))+
  
  geom_line(data=InVivo, aes(x=xfree,y=PredFracChange, color=Model), size=2, show.legend = FALSE)+
  
  scale_color_manual(values=c("dodgerblue4", "dodgerblue"))+
  
  geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p50), color="red", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line()+
  
  geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
  
  geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
  
  geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
  
  xlab("Free Concentration (uM)")+
  
  ylab("Decay Rise Ratio (Frac Change)")+
  
  scale_x_log10()+
  
  theme(legend.position = "none")+
  
  ggtitle("N-acetylprocainamide  In Vivo vs. In Vitro Fractional Change Peak Decay Rise Ratio Up")+
  
  theme_minimal()

ggsave("Invivo_Dose_Response_Frac.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)

ggsave("Invivo_Dose_Response_Frac.png", plot=last_plot(), device="png", height= 4.25, width= 4.25)
  
  
  rm(list=files)
}

###Master Pos-neg cont function

DR_master_up_neg_cont<- function(chemical.name){

  setwd(paste("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\", chemical.name, "\\decay_rise_up", sep=""))
  
  ### load in CSV files, R dat file,  set up concentration frame
  
  files = list.files(pattern="*samples_*")
  
  for (l in 1:length(files)) assign(files[l], read.csv(files[l]))
  
  file.sources=list.files(pattern="*dat.R")
  
  sapply(file.sources,source,.GlobalEnv)
  
  concframe<-data.frame(x=10^(seq(-3,2,0.05)))
  
  ### extract median parameter data from csv files, create new frame with calculated parameter data, do the same when z scores and SD =0 
  
  m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))
  
  m_y0$z_score_sim<- rnorm(8000,0,1)
  
  simframe<- as.data.frame(exp(m_y0$m_y0+(m_y0$sd_y0*m_y0$z_score_sim)))
  
  rownames(simframe)<- c(1:8000)
  
  colnames(simframe)<- "sim_y0"
  
  m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))
  
  m_x0$z_score_sim<- rnorm(8000,0,1)
  
  simframe$sim_x0<-exp(m_x0$m_x0+(m_x0$sd_x0*m_x0$z_score_sim))
  
  m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))
  
  m_Emax$z_score_sim<- rnorm(8000,0,1)
  
  simframe$sim_Emax<-exp(m_Emax$m_Emax+(m_Emax$sd_Emax*m_Emax$z_score_sim))
  
  m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))
  
  m_n$z_score_sim<- rnorm(8000,0,1)
  
  simframe$sim_n<-exp(m_n$m_n+(m_n$sd_n*m_n$z_score_sim))
  
  simframe_zero<-data.frame(y0=(exp(m_y0$m_y0)), x0=(exp(m_x0$m_x0)), Emax=(exp(m_Emax$m_Emax)), n=(exp(m_n$m_n)))
  
  ### Load in appropriate control CSV and Parameter CSV
  
  InVivoNeg<- read.csv("InVivoPred.neg.v2.csv", header=TRUE)
  
  InVivoParm<- read.csv("QTposnegparameters.v2.csv", header=TRUE)
  
  InVivo<-subset(InVivoNeg, Chemical.name== chemical.name )
  
  InVivoFree<- subset(InVivoParm, Chemical== chemical.name )
  
  ### Set up working frame with concentrations, transform concentrations into freeconc
  
  simframeconc<- data.frame()
  
  simframe_zero_conc<- data.frame()
  
  for (i in 1:nrow(concframe)){
    
    simframe$concentration<-concframe$x[i]
    
    simframeconc<- rbind(simframeconc, simframe)
  }
  
  simframeconc$freeconcentration<- (simframeconc$concentration)*InVivoFree$FracFreeMedia
  
  for (i in 1:nrow(concframe)){
    
    simframe_zero$concentration<-concframe$x[i]
    
    simframe_zero_conc<- rbind(simframe_zero_conc, simframe_zero)
  }
  
  simframe_zero_conc$freeconcentration<- (simframe_zero_conc$concentration)*InVivoFree$FracFreeMedia
  
  
  ##create final frames for plotting w/ quantile values
  
  simframeconc$PredFold<- DR_up(x=simframeconc$concentration, y0=1, x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n)
  
  simframe_median_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframeconc, quantile, prob=0.5)
  
  simframe_975_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframeconc, quantile, prob=0.975)
  
  simframe_25_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframeconc, quantile, prob=0.025)
  
  Simulation_Frame_Fold<- simframe_median_fold
  
  Simulation_Frame_Fold$p97.5<- simframe_975_fold$PredFold
  
  Simulation_Frame_Fold$p2.5<- simframe_25_fold$PredFold
  
  names(Simulation_Frame_Fold)[3]<- "p50"
  
  ##########
  
  simframeconc$PredFrac<- DR_up_Frac(x=simframeconc$concentration, y0=1, x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n) 
  
  simframe_median_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframeconc, quantile, prob=0.5)
  
  simframe_975_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframeconc, quantile, prob=0.975)
  
  simframe_25_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframeconc, quantile, prob=0.025)
  
  Simulation_Frame_Frac<- simframe_median_frac
  
  Simulation_Frame_Frac$p97.5<- simframe_975_frac$PredFrac
  
  Simulation_Frame_Frac$p2.5<- simframe_25_frac$PredFrac
  
  names(Simulation_Frame_Frac)[3]<- "p50"
  
  ##########
  
  simframe_zero_conc$PredFold<- DR_up(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n)
  
  simframe_zero_median_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.5)
  
  simframe_zero_975_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.975)
  
  simframe_zero_25_fold<- aggregate(PredFold~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.025)
  
  Simulation_Frame_Zero_Fold<- simframe_zero_median_fold
  
  Simulation_Frame_Zero_Fold$p97.5<- simframe_zero_975_fold$PredFold
  
  Simulation_Frame_Zero_Fold$p2.5<- simframe_zero_25_fold$PredFold
  
  names(Simulation_Frame_Zero_Fold)[3]<- "p50"
  
  ######
  
  simframe_zero_conc$PredFrac<- DR_up_Frac(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n) 
  
  simframe_zero_median_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.5)
  
  simframe_zero_975_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.975)
  
  simframe_zero_25_frac<- aggregate(PredFrac~freeconcentration+concentration, data=simframe_zero_conc, quantile, prob=0.025)
  
  Simulation_Frame_Zero_Frac<- simframe_zero_median_frac
  
  Simulation_Frame_Zero_Frac$p97.5<- simframe_zero_975_frac$PredFrac
  
  Simulation_Frame_Zero_Frac$p2.5<- simframe_zero_25_frac$PredFrac
  
  names(Simulation_Frame_Zero_Frac)[3]<- "p50"
  
  
  ###output workind dataframes as CSV files
  
  #set wd for outputs
  
  setwd(paste("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\", chemical.name, "\\decay_rise_up", "\\invivo_invitro", sep=""))
  
  write.csv(Simulation_Frame_Fold, "InVivo_Frame_Fold_rescaled.csv")
  
  write.csv(Simulation_Frame_Zero_Fold, "Invivo_Frame_Zero_Fold_rescaled.csv")
  
  write.csv(Simulation_Frame_Frac, "Invivo_Frame_Frac_rescaled.csv")
  
  write.csv(Simulation_Frame_Zero_Frac, "Invivo_Frame_Zero_Frac_rescaled.csv")
  
  ###plot
  
  plotname_fold_log<-(paste(chemical.name, "In Vitro vs. In Vivo Fold Change Peak Decay Rise Ratio Up (Log Scale)", sep=" "))
  
  plotname_fold<-(paste(chemical.name, "In Vitro vs. In Vivo Fold Change Peak Decay Rise Ratio Up", sep=" "))
   
  plotname_frac_log<-(paste(chemical.name, "In Vitro vs. In Vivo Frac Change Peak Decay Rise Ratio Up (Log Scale)", sep=" "))
  
  plotname_frac<-(paste(chemical.name, "In Vitro vs. In Vivo Frac Change Peak Decay Rise Ratio Up (Log Scale)", sep=" "))
  
  #Fold Change plots
  
  ggplot(Simulation_Frame_Fold,aes(x=concentration, y=p50))+
    
    geom_line(data=InVivo, aes(x=xfree, y=PredFoldChange), color="blue", size=2)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p50), color="red", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line()+
    
    geom_line(aes(x=freeconcentration, y= p97.5), linetype="dotted")+
    
    geom_line(aes(x=freeconcentration, y= p2.5), linetype="dotted")+
    
    geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
    
    xlab("Free Concentration (uM)")+
    
    ylab("Decay Rise Ratio (Fold Change)")+
    
    scale_x_log10()+
    
    scale_y_log10()+
    
    ggtitle("Sotalol In Vitro vs. In Vivo Fold Change Peak Decay Rise Ratio Up (Log Scale)")+
    
    theme_minimal()
  
  ggsave("Invivo_Dose_Response_Log_Fold_rescaled.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)
  
  ggsave("InVivo_Dose_Response_Log_Fold_rescaled.png", plot=last_plot(), device="png", height= 4.25, width=4.25)
  
  ggplot(Simulation_Frame_Fold, aes(x=concentration, y=p50))+
    
    geom_line(data=InVivo, aes(x=xfree, y=PredFoldChange), color="blue", size=2)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p50), color="red", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Fold, aes(x=freeconcentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line()+
    
    geom_line(aes(x=freeconcentration, y= p97.5), linetype="dotted")+
    
    geom_line(aes(x=freeconcentration, y= p2.5), linetype="dotted")+
    
    geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
    
    xlab("Free Concentration (uM)")+
    
    ylab("Decay Rise Ratio (Fold Change)")+
    
    scale_x_log10()+
    
    ggtitle("Sotalol In Vivo vs. In Vitro Fold Change Peak Decay Rise Ratio Up")+
    
    theme_minimal()
  
  ggsave("Invivo_Dose_Response_Fold_rescaled.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)
  
  ggsave("InVivo_Dose_Response_Fold_rescaled.png", plot=last_plot(), device="png", height= 4.25, width=4.25)
  
  # Frac Change Plots
  
  ggplot(Simulation_Frame_Frac, aes(x=concentration, y=p50))+
    
    geom_line(data=InVivo, aes(x=xfree, y=PredFracChange), color="blue", size=2)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p50), color="red", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line()+
    
    geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
    
    geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
    
    geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
    
    xlab("Log Free Concentration (uM)")+
    
    ylab("Decay Rise Ratio (Log Frac Change)")+
    
    scale_x_log10()+
    
    scale_y_log10(breaks=c(1,1e-02,1e-04,1e-06,1e-08))+
    
    ggtitle("Sotalol  In Vivo vs. In Vitro Fractional Change Peak Decay Rise Ratio Up (Log Scale)")+
    
    theme_minimal()
  
  ggsave("Invivo_Dose_Response_Log_Frac_rescaled.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)
  
  ggsave("Invivo_Dose_Response_Log_Frac_rescaled.png", plot=last_plot(), device="png", height= 4.25, width= 4.25)
  
  ggplot(Simulation_Frame_Frac, aes(x=concentration, y=p50))+
    
    geom_line(data=InVivo, aes(x=xfree, y=PredFracChange), color="blue", size=2)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p50), color="red", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
    
    geom_line()+
    
    geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
    
    geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
    
    geom_vline(xintercept=(0.1*InVivoFree$FracFreeMedia), color="grey", linetype="dotted", size=2)+
    
    xlab("Free Concentration (uM)")+
    
    ylab("Decay Rise Ratio (Frac Change)")+
    
    scale_x_log10()+
    
    ggtitle("Sotalol  In Vivo vs. In Vitro Fractional Change Peak Decay Rise Ratio Up")+
    
    theme_minimal()
  
  ggsave("Invivo_Dose_Response_Frac_rescaled.pdf", plot=last_plot(), device="pdf", height= 7.5, width= 10)
  
  ggsave("Invivo_Dose_Response_Frac_rescaled.png", plot=last_plot(), device="png", height= 4.25, width= 4.25)
  
  
  rm(list=files)
}

