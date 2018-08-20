###set WD and load packages

setwd("C:\\Users\\ablanchette\\Desktop\\Test_Env")

library(ggplot2)

### load in CSV files, set up concentration frame

files = list.files(pattern="*samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

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

### set up working frames with concentrations specified in concframe

simframeconc<- data.frame()

simframe_zero_conc<- data.frame()

for (i in 1:nrow(concframe)){
  
  simframe$concentration<-concframe$x[i]
  
  simframeconc<- rbind(simframeconc, simframe)
}

for (i in 1:nrow(concframe)){
  
  simframe_zero$concentration<-concframe$x[i]
  
  simframe_zero_conc<- rbind(simframe_zero_conc, simframe_zero)
}

##set up DR Functions and create final frames for plotting w/ quantile values

DR_up<-function(x=1, y0=1, x0=1, Emax=1, n=1){
  
  y0*(1+((x/x0)^n)/(1+((x/x0)^n)*(1/Emax))) 
  
}

DR_up_Frac<- function(x=1, y0=1, x0=1, Emax=1, n=1){
  (y0*(1+((x/x0)^n)/(1+((x/x0)^n)*(1/Emax))))-1 
}

simframeconc$PredFold<- DR_up(x=simframeconc$concentration, y0=1, x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n) 

simframe_median_fold<- aggregate(PredFold~concentration, data=simframeconc, quantile, prob=0.5)

simframe_975_fold<- aggregate(PredFold~concentration, data=simframeconc, quantile, prob=0.975)

simframe_25_fold<- aggregate(PredFold~concentration, data=simframeconc, quantile, prob=0.025)

Simulation_Frame_Fold<- simframe_median_fold

Simulation_Frame_Fold$p97.5<- simframe_975_fold$PredFold

Simulation_Frame_Fold$p2.5<- simframe_25_fold$PredFold

names(Simulation_Frame_Fold)[2]<- "p50"

##########

simframeconc$PredFrac<- DR_up_Frac(x=simframeconc$concentration, y0=1, x0=simframeconc$sim_x0, Emax=simframeconc$sim_Emax, n=simframeconc$sim_n) 

simframe_median_frac<- aggregate(PredFrac~concentration, data=simframeconc, quantile, prob=0.5)

simframe_975_frac<- aggregate(PredFrac~concentration, data=simframeconc, quantile, prob=0.975)

simframe_25_frac<- aggregate(PredFrac~concentration, data=simframeconc, quantile, prob=0.025)

Simulation_Frame_Frac<- simframe_median_frac

Simulation_Frame_Frac$p97.5<- simframe_975_frac$PredFrac

Simulation_Frame_Frac$p2.5<- simframe_25_frac$PredFrac

names(Simulation_Frame_Frac)[2]<- "p50"

##########

simframe_zero_conc$PredFold<- DR_up(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n) 

simframe_zero_median_fold<- aggregate(PredFold~concentration, data=simframe_zero_conc, quantile, prob=0.5)

simframe_zero_975_fold<- aggregate(PredFold~concentration, data=simframe_zero_conc, quantile, prob=0.975)

simframe_zero_25_fold<- aggregate(PredFold~concentration, data=simframe_zero_conc, quantile, prob=0.025)

Simulation_Frame_Zero_Fold<- simframe_zero_median_fold

Simulation_Frame_Zero_Fold$p97.5<- simframe_zero_975_fold$PredFold

Simulation_Frame_Zero_Fold$p2.5<- simframe_zero_25_fold$PredFold

names(Simulation_Frame_Zero_Fold)[2]<- "p50"

######

simframe_zero_conc$PredFrac<- DR_up_Frac(x=simframe_zero_conc$concentration, y0=1, x0=simframe_zero_conc$x0, Emax=simframe_zero_conc$Emax, n=simframe_zero_conc$n) 

simframe_zero_median_frac<- aggregate(PredFrac~concentration, data=simframe_zero_conc, quantile, prob=0.5)

simframe_zero_975_frac<- aggregate(PredFrac~concentration, data=simframe_zero_conc, quantile, prob=0.975)

simframe_zero_25_frac<- aggregate(PredFrac~concentration, data=simframe_zero_conc, quantile, prob=0.025)

Simulation_Frame_Zero_Frac<- simframe_zero_median_frac

Simulation_Frame_Zero_Frac$p97.5<- simframe_zero_975_frac$PredFrac

Simulation_Frame_Zero_Frac$p2.5<- simframe_zero_25_frac$PredFrac

names(Simulation_Frame_Zero_Frac)[2]<- "p50"


###output workind dataframes as CSV files

write.csv(Simulation_Frame_Fold, "Simulation_Frame_Fold.csv")

write.csv(Simulation_Frame_Zero_Fold, "Simulation_Frame_Zero_Fold.csv")

write.csv(Simulation_Frame_Frac, "Simulation_Frame_Frac.csv")

write.csv(Simulation_Frame_Zero_Frac, "Simulation_Frame_Zero_Frac.csv")

write.csv(simframeconc, "simframeconc.csv")

write.csv(simframe_zero_conc, "simframe_zero_conc.csv")

###plot

#Fold Change plots

ggplot(Simulation_Frame_Fold, aes(x=concentration, y=p50))+
  
  geom_line(data=Simulation_Frame_Zero_Fold, aes(x=concentration, y=p50), color="red", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Fold, aes(x=concentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Fold, aes(x=concentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line()+
  
  geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
  
  geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
  
  xlab("Concentration (uM)")+
  
  ylab("Decay Rise Ratio (Fold Change)")+
  
  scale_x_log10()+
  
  scale_y_log10()+
  
  ggtitle("Sotalol Fold Change Peak Decay Rise Ratio Up (Log Scale)")+
  
  theme_minimal()

ggsave("Simulated_Dose_Response_Log_Fold.pdf", plot=last_plot(), device="pdf", height= 10, width= 7.5)

ggplot(Simulation_Frame_Fold, aes(x=concentration, y=p50))+
  
  geom_line(data=Simulation_Frame_Zero_Fold, aes(x=concentration, y=p50), color="red", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Fold, aes(x=concentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Fold, aes(x=concentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line()+
  
  geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
  
  geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
  
  xlab("Concentration (uM)")+
  
  ylab("Decay Rise Ratio (Fold Change)")+
  
  scale_x_log10()+
  
  ggtitle("Sotalol Fold Change Peak Decay Rise Ratio Up")+
  
  theme_minimal()

ggsave("Simulated_Dose_Response_Fold.pdf", plot=last_plot(), device="pdf", height= 10, width= 7.5)

# Frac Change Plots

ggplot(Simulation_Frame_Frac, aes(x=concentration, y=p50))+
  
  geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p50), color="red", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line()+
  
  geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
  
  geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
  
  xlab("Concentration (uM)")+
  
  ylab("Decay Rise Ratio (Frac Change)")+
  
  scale_x_log10()+
  
  scale_y_log10()+
  
  ggtitle("Sotalol Fractional Change Peak Decay Rise Ratio Up (Log Scale)")+
  
  theme_minimal()

ggsave("Simulated_Dose_Response_Log_Frac.pdf", plot=last_plot(), device="pdf", height= 10, width= 7.5)

ggplot(Simulation_Frame_Frac, aes(x=concentration, y=p50))+
  
  geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p50), color="red", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p97.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line(data=Simulation_Frame_Zero_Frac, aes(x=concentration, y=p2.5), color="red", linetype="dashed", size=1.5)+
  
  geom_line()+
  
  geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
  
  geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
  
  xlab("Concentration (uM)")+
  
  ylab("Decay Rise Ratio (Frac Change)")+
  
  scale_x_log10()+
  
  ggtitle("Sotalol Fractional Change Peak Decay Rise Ratio Up")+
  
  theme_minimal()

ggsave("Simulated_Dose_Response_Frac.pdf", plot=last_plot(), device="pdf", height= 10, width= 7.5)
rm(list=ls())

