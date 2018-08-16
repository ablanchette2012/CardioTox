##load packages

library(reshape2)

library(ggplot2)

DR_up<-function(x=1, y0=1, x0=1, Emax=1, n=1){
  
   y0*(1+((x/x0)^n)/(1+((x/x0)^n)*(1/Emax))) 
  
}

obs_vs_pred_func<-function(chemical.name) {
  
  setwd(paste("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_CO\\", chemical.name, "\\decay_rise_up", sep=""))
  
  file.sources=list.files(pattern="*dat.R")
  
  sapply(file.sources,source,.GlobalEnv)
  
  y_obs<-data.frame(cell,x,ys)
  
  colnames(y_obs)<- c( "Individual","Dose", "Response")
  
  ##call in and process output .csv files into seperate frames
  
  files = list.files(pattern="*_samples_*")
  
  for (l in 1:length(files)) assign(files[l], read.csv(files[l]))
  
  y0<-as.data.frame(rbind((eval(parse(text=files[1])))[4004, 125:151], (eval(parse(text=files[2])))[4004,125:151], (eval(parse(text=files[3])))[4004, 125:151], (eval(parse(text=files[4])))[4004, 125:151]))
  
  colnames(y0)<- c(1:27)
  
  x0<-as.data.frame(rbind((eval(parse(text=files[1])))[4004, 152:178], (eval(parse(text=files[2])))[4004, 152:178], (eval(parse(text=files[3])))[4004, 152:178], (eval(parse(text=files[4])))[4004, 152:178]))
  
  colnames(x0)<- c(1:27)
  
  Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[4004, 179:205], (eval(parse(text=files[2])))[4004, 179:205], (eval(parse(text=files[3])))[4004, 179:205], (eval(parse(text=files[4])))[4004, 179:205]))
  
  colnames(Emax)<- c(1:27)
  
  n<-as.data.frame(rbind((eval(parse(text=files[1])))[4004, 206:232], (eval(parse(text=files[2])))[4004, 206:232], (eval(parse(text=files[3])))[4004, 206:232], (eval(parse(text=files[4])))[4004, 206:232]))
  
  colnames(n)<- c(1:27)
  
  ##set up DR function and create Predicted column in y_obs frame
  
  y_obs$Predicted<- t(DR_up(x=y_obs$Dose, x0=x0[1,y_obs$Individual], y0=y0[1,y_obs$Individual], n=n[1,y_obs$Individual], Emax=Emax[1,y_obs$Individual]))
  
  y_obs$Predicted<-(y_obs$Predicted)*scale_factor
  
  y_obs$Response<-(y_obs$Response*scale_factor)
  
  ##Plot
  
  ggplot()+
    
    geom_point(data=y_obs, aes(x=Response, y=Predicted, color=Individual))+
    
    scale_x_log10(limits=c(1,100))+
    
    scale_y_log10(limits=c(1,100)) +
    
    annotation_logticks() +
    
    theme(panel.grid.minor = element_blank())+
    
    geom_abline(intercept=0, slope=1, color="grey60")+
    
    geom_abline(intercept=0.5, slope=1, color="grey60", linetype="dashed")+
    
    geom_abline(intercept=-0.5, slope=1, color="grey60", linetype="dashed")
  
  y_obs$chemical.name<-chemical.name
  
  return(y_obs)
  
}

chem.list<-(c("Cisapride", "Citalopram", "Disopyramide", 'Dofetilide', "Moxifloxacin", "N-acetylprocainamide", "Quinidine", "Sematilide", "Sotalol", "Vernacalant", "Cabazitaxel", "Lamotrigine", "Mifepristone"))

Obs_vs_Pred_All<-data.frame()

for (chem in chem.list){

Obs_vs_Pred_All<-rbind(Obs_vs_Pred_All, obs_vs_pred_func(chem))
}

Obs_vs_Pred_All_Cont<-subset(Obs_vs_Pred_All, Dose==0 & chemical.name=="Cisapride")

Obs_vs_Pred_All_Treatment<-subset(Obs_vs_Pred_All, Dose>0)

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_CO\\All")

write.csv(Obs_vs_Pred_All, "Observed_Vs_Predicted_All_Chems.csv")

### Plot

##Plot All

ggplot(Obs_vs_Pred_All_Treatment, aes(x= Response, y=Predicted, color=Individual, size=chemical.name))+
  
  geom_point(data=Obs_vs_Pred_All_Treatment, aes(x=Response, y=Predicted, color=Individual))+
  
  scale_x_log10(limits=c(1,100))+
  
  scale_y_log10(limits=c(1,100)) +
  
  annotation_logticks() +
  
  theme(panel.grid.minor = element_blank())+
  
  geom_abline(intercept=0, slope=1, color="grey60")+
  
  geom_abline(intercept=0.5, slope=1, color="grey60", linetype="dashed")+
  
  geom_abline(intercept=-0.5, slope=1, color="grey60", linetype="dashed")


  #ggtitle("Observed vs. Predicted Response for All Chemicals")

ggsave("Observed_vs_Predicted_All_Chems_CO_Treatment.pdf", plot=last_plot(), device="pdf", height= 5.5, width=7.5)

ggsave("Observed_vs_Predicted_All_Chems_CO.png", plot=last_plot(), device="png", height= 4.25, width=4.25)

##Plot Controls

ggplot(Obs_vs_Pred_All_Cont, aes(x= Response, y=Predicted, color=Individual, size=Chemical.name))+
  
  geom_point()+
  
  scale_x_log10(limits=c(-10,50))+
  
  scale_y_log10(limits=c(-10,50))+
  
  labs(x="Observed Response", y="Predicted Response")+
  
  geom_abline(intercept=0, slope=1, color="grey60")+
  
  geom_abline(intercept=0.5, slope=1, color="grey60", linetype="dashed")+
  
  geom_abline(intercept=-0.5, slope=1, color="grey60", linetype="dashed")

#ggtitle("Observed vs. Predicted Response for All Chemicals")

ggsave("Observed_vs_Predicted_All_Chems_CO_Cont.pdf", plot=last_plot(), device="pdf", height= 5.5, width=7.5)

ggsave("Observed_vs_Predicted_All_Chems_CO_Cont.png", plot=last_plot(), device="png", height= 4.25, width=4.25)


