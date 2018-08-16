setwd("C:\\Users\\ablanchette\\Desktop\\Test_Env\\decay_rise_up")

library(ggplot2)

library(reshape2)

file.sources=list.files(pattern="*dat.R")

sapply(file.sources,source,.GlobalEnv)

DR_Frame<- data.frame(cell, x, ys)

colnames(DR_Frame)<- c("Individual", "Dose", "Response")

DR_Frame$Dose[DR_Frame$Dose==0]<- 0.001

files = list.files(pattern="*samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

concframe<-data.frame(x=10^(seq(-3,2,0.05)))

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

parmframeconc<- data.frame()

for (i in 1:nrow(concframe)){
  
  parmframe$concentration<-concframe$x[i]
  
  parmframeconc<- rbind(parmframeconc, parmframe)
  }

##set up DR function and create Predicted column in y_obs frame

DR_dn<-function(x=1, y0=1, x0=1, Emax=1, n=1){
  
  
  y0*(1-((x/x0)^n)/(1+((x/x0)^n)/(Emax))) 
  
}

parmframeconc$Pred<- (DR_dn(x=parmframeconc$concentration, y0=parmframeconc$y0, x0=parmframeconc$x0, Emax=parmframeconc$Emax, n=parmframeconc$n))*(scale_factor)

predframe_median<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.5)

predframe_975<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.975)

predframe_25<- aggregate(Pred~concentration+ID, data=parmframeconc, quantile, prob=0.025)

Prediction_Frame<- predframe_median

Prediction_Frame$p97.5<- predframe_975$Pred

Prediction_Frame$p2.5<- predframe_25$Pred

names(Prediction_Frame)[3]<- "p50"

names(Prediction_Frame)[2]<- "Individual"

write.csv(DR_Frame, "DR_Frame.csv")

write.csv(Prediction_Frame, "Prediction_Frame.csv")

#write.csv(parmframeconc, "parmframeconc.csv")

ggplot(Prediction_Frame, aes(x=concentration, y=p50))+
  
  geom_line(color="red")+
  
  geom_line(aes(x=concentration, y= p97.5), linetype="dotted")+
  
  geom_line(aes(x=concentration, y= p2.5), linetype="dotted")+
  
  geom_point(data=DR_Frame, aes(x=Dose, y=Response))+
  
  xlab("Concentration (uM)")+
  
  ylab("Decay Rise Ratio")+
  
  scale_x_log10()+
  
  ggtitle("Test Peak Decay Rise Ratio Up")+
  
  theme_minimal()+
  
  facet_wrap(~Individual, nrow=9, ncol=3, scales = "free_y")
  
ggsave("Individual_Dose_Response_rescaled.pdf", plot=last_plot(), device="pdf", height= 10, width= 7.5)

rm(list=ls())
