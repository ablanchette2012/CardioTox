## set wd and load in packages

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\All")

library(ggplot2)

##read in and format in vivo data

Invivo<- read.csv("InvivoPred.all.csv")

QT_Parameters<- read.csv("QTposnegparameters.csv")

##set up fold and frac equations

DR_up<-function(x=1, y0=1, x0=1, Emax=1, n=1){
  
  y0*(1+((x/x0)^n)/(1+((x/x0)^n)*(1/Emax))) 
  
}

DR_up_Frac<- function(x=1, y0=1, x0=1, Emax=1, n=1){
  (y0*(1+((x/x0)^n)/(1+((x/x0)^n)*(1/Emax))))-1 
}

################### Begin Chems #############################

### Cisapride ###

##set up frame with parameters

file.sources=list.files(pattern="*_17_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_17_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))


Cisapride<- as.data.frame(exp(m_y0$m_y0))

rownames(Cisapride)<- c(1:8000)

colnames(Cisapride)<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))


Cisapride$x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))


Cisapride$Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))


Cisapride$n<-exp(m_n$m_n)

## Calc out Cisapride response at xfree and create Quantile frames

Cisapride$Pred_Fold<- (DR_up(y0=1, x0=Cisapride$x0, Emax=Cisapride$Emax, n=Cisapride$n, x=((Invivo$xfree[1])/QT_Parameters$FracFreeMedia[1])))*(scale_factor)

Cisapride$Pred_Frac<- (DR_up_Frac(y0=1, x0=Cisapride$x0, Emax=Cisapride$Emax, n=Cisapride$n, x=((Invivo$xfree[1])/QT_Parameters$FracFreeMedia[1])))*(scale_factor)

Quantiles_Fold<- data.frame(p50= median(Cisapride$Pred_Fold))

Quantiles_Fold$p975<- quantile(Cisapride$Pred_Fold, prob=0.975)

Quantiles_Fold$p25<- quantile(Cisapride$Pred_Fold, prob=0.025)

Quantiles_Frac<- data.frame(p50=median(Cisapride$Pred_Frac))

Quantiles_Frac$p975<- quantile(Cisapride$Pred_Frac, prob=0.975)

Quantiles_Frac$p25<-quantile(Cisapride$Pred_Frac, prob=0.025)

rm(list=files, Cisapride)

### Citalopram ###

##set up frame with parameters

file.sources=list.files(pattern="*_5_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_5_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Citalopram<- as.data.frame(exp(m_y0$m_y0))

rownames(Citalopram)<- c(1:8000)

colnames(Citalopram)<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Citalopram$x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Citalopram$Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Citalopram$n<-exp(m_n$m_n)

## Calc out Citalopram response at xfree and create Quantile frames

Citalopram$Pred_Fold<- (DR_up(y0=1, x0=Citalopram$x0, Emax=Citalopram$Emax, n=Citalopram$n, x=((Invivo$xfree[2])/QT_Parameters$FracFreeMedia[2])))*(scale_factor)

Citalopram$Pred_Frac<- (DR_up_Frac(y0=1, x0=Citalopram$x0, Emax=Citalopram$Emax, n=Citalopram$n, x=((Invivo$xfree[2])/QT_Parameters$FracFreeMedia[2])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Citalopram$Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Citalopram$Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Citalopram$Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Citalopram$Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Citalopram$Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Citalopram$Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Citalopram)

### Disopyramide Linear ###

##set up frame with parameters

file.sources=list.files(pattern="*_6_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_6_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Disopyramide<- as.data.frame(exp(m_y0$m_y0))

rownames(Disopyramide)<- c(1:8000)

colnames(Disopyramide)<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Disopyramide$x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Disopyramide$Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Disopyramide$n<-exp(m_n$m_n)

## Calc out Disopyramide response at xfree and create Quantile frames

Disopyramide$Pred_Fold<- (DR_up(y0=1, x0=Disopyramide$x0, Emax=Disopyramide$Emax, n=Disopyramide$n, x=((Invivo$xfree[3])/QT_Parameters$FracFreeMedia[3])))*(scale_factor)

Disopyramide$Pred_Frac<- (DR_up_Frac(y0=1, x0=Disopyramide$x0, Emax=Disopyramide$Emax, n=Disopyramide$n, x=((Invivo$xfree[3])/QT_Parameters$FracFreeMedia[3])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Disopyramide$Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Disopyramide$Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Disopyramide$Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Disopyramide$Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Disopyramide$Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Disopyramide$Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Disopyramide)

### Disopyramide Hill ###

##set up frame with parameters

file.sources=list.files(pattern="*_6_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_6_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Disopyramide<- as.data.frame(exp(m_y0$m_y0))

rownames(Disopyramide)<- c(1:8000)

colnames(Disopyramide)<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Disopyramide$x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Disopyramide$Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Disopyramide$n<-exp(m_n$m_n)

## Calc out Disopyramide response at xfree and create Quantile frames

Disopyramide$Pred_Fold<- (DR_up(y0=1, x0=Disopyramide$x0, Emax=Disopyramide$Emax, n=Disopyramide$n, x=((Invivo$xfree[4])/QT_Parameters$FracFreeMedia[3])))*(scale_factor)

Disopyramide$Pred_Frac<- (DR_up_Frac(y0=1, x0=Disopyramide$x0, Emax=Disopyramide$Emax, n=Disopyramide$n, x=((Invivo$xfree[4])/QT_Parameters$FracFreeMedia[3])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Disopyramide$Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Disopyramide$Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Disopyramide$Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Disopyramide$Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Disopyramide$Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Disopyramide$Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Disopyramide)


### Dofetilide ###

##set up frame with parameters

file.sources=list.files(pattern="*_15_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_15_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Dofetilide<- as.data.frame(exp(m_y0$m_y0))

rownames(Dofetilide)<- c(1:8000)

colnames(Dofetilide)<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Dofetilide$x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Dofetilide$Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Dofetilide$n<-exp(m_n$m_n)

## Calc out Dofetilide response at xfree and create Quantile frames

Dofetilide$Pred_Fold<- (DR_up(y0=1, x0=Dofetilide$x0, Emax=Dofetilide$Emax, n=Dofetilide$n, x=((Invivo$xfree[5])/QT_Parameters$FracFreeMedia[4])))*(scale_factor)

Dofetilide$Pred_Frac<- (DR_up_Frac(y0=1, x0=Dofetilide$x0, Emax=Dofetilide$Emax, n=Dofetilide$n, x=((Invivo$xfree[5])/QT_Parameters$FracFreeMedia[4])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Dofetilide$Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Dofetilide$Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Dofetilide$Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Dofetilide$Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Dofetilide$Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Dofetilide$Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Dofetilide)

### Moxifloxacin ###

##set up frame with parameters

file.sources=list.files(pattern="*_37_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_37_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Moxifloxacin<- as.data.frame(exp(m_y0$m_y0))

rownames(Moxifloxacin)<- c(1:8000)

colnames(Moxifloxacin)<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Moxifloxacin$x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Moxifloxacin$Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Moxifloxacin$n<-exp(m_n$m_n)

## Calc out Moxifloxacin response at xfree and create Quantile frames

Moxifloxacin$Pred_Fold<- (DR_up(y0=1, x0=Moxifloxacin$x0, Emax=Moxifloxacin$Emax, n=Moxifloxacin$n, x=((Invivo$xfree[6])/QT_Parameters$FracFreeMedia[5])))*(scale_factor)

Moxifloxacin$Pred_Frac<- (DR_up_Frac(y0=1, x0=Moxifloxacin$x0, Emax=Moxifloxacin$Emax, n=Moxifloxacin$n, x=((Invivo$xfree[6])/QT_Parameters$FracFreeMedia[5])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Moxifloxacin$Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Moxifloxacin$Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Moxifloxacin$Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Moxifloxacin$Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Moxifloxacin$Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Moxifloxacin$Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Moxifloxacin)

### N_acetylprocainamide Linear ###

##set up frame with parameters

file.sources=list.files(pattern="*_2_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_2_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

N_acetylprocainamide<- as.data.frame(exp(m_y0$m_y0))

rownames(N_acetylprocainamide)<- c(1:8000)

colnames(N_acetylprocainamide)<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

N_acetylprocainamide$x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

N_acetylprocainamide$Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

N_acetylprocainamide$n<-exp(m_n$m_n)

## Calc out N_acetylprocainamide response at xfree and create Quantile frames

N_acetylprocainamide$Pred_Fold<- (DR_up(y0=1, x0=N_acetylprocainamide$x0, Emax=N_acetylprocainamide$Emax, n=N_acetylprocainamide$n, x=((Invivo$xfree[7])/QT_Parameters$FracFreeMedia[6])))*(scale_factor)

N_acetylprocainamide$Pred_Frac<- (DR_up_Frac(y0=1, x0=N_acetylprocainamide$x0, Emax=N_acetylprocainamide$Emax, n=N_acetylprocainamide$n, x=((Invivo$xfree[7])/QT_Parameters$FracFreeMedia[6])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(N_acetylprocainamide$Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(N_acetylprocainamide$Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(N_acetylprocainamide$Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(N_acetylprocainamide$Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(N_acetylprocainamide$Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(N_acetylprocainamide$Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, N_acetylprocainamide)

### N_acetylprocainamide Hill ###

##set up frame with parameters

file.sources=list.files(pattern="*_2_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_2_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

N_acetylprocainamide<- as.data.frame(exp(m_y0$m_y0))

rownames(N_acetylprocainamide)<- c(1:8000)

colnames(N_acetylprocainamide)<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

N_acetylprocainamide$x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

N_acetylprocainamide$Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

N_acetylprocainamide$n<-exp(m_n$m_n)

## Calc out N_acetylprocainamide response at xfree and create Quantile frames

N_acetylprocainamide$Pred_Fold<- (DR_up(y0=1, x0=N_acetylprocainamide$x0, Emax=N_acetylprocainamide$Emax, n=N_acetylprocainamide$n, x=((Invivo$xfree[8])/QT_Parameters$FracFreeMedia[6])))*(scale_factor)

N_acetylprocainamide$Pred_Frac<- (DR_up_Frac(y0=1, x0=N_acetylprocainamide$x0, Emax=N_acetylprocainamide$Emax, n=N_acetylprocainamide$n, x=((Invivo$xfree[8])/QT_Parameters$FracFreeMedia[6])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(N_acetylprocainamide$Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(N_acetylprocainamide$Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(N_acetylprocainamide$Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(N_acetylprocainamide$Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(N_acetylprocainamide$Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(N_acetylprocainamide$Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, N_acetylprocainamide)

### Quinidine ###

##set up frame with parameters

file.sources=list.files(pattern="*_34_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_34_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Quinidine<- as.data.frame(exp(m_y0$m_y0))

rownames(Quinidine)<- c(1:8000)

colnames(Quinidine)<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Quinidine$x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Quinidine$Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Quinidine$n<-exp(m_n$m_n)

## Calc out Quinidine response at xfree and create Quantile frames

Quinidine$Pred_Fold<- (DR_up(y0=1, x0=Quinidine$x0, Emax=Quinidine$Emax, n=Quinidine$n, x=((Invivo$xfree[9])/QT_Parameters$FracFreeMedia[7])))*(scale_factor)

Quinidine$Pred_Frac<- (DR_up_Frac(y0=1, x0=Quinidine$x0, Emax=Quinidine$Emax, n=Quinidine$n, x=((Invivo$xfree[9])/QT_Parameters$FracFreeMedia[7])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Quinidine$Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Quinidine$Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Quinidine$Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Quinidine$Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Quinidine$Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Quinidine$Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Quinidine)

### Sematilide Hill  ###

##set up frame with parameters

file.sources=list.files(pattern="*_35_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_35_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Sematilide <- as.data.frame(exp(m_y0$m_y0))

rownames(Sematilide )<- c(1:8000)

colnames(Sematilide )<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Sematilide $x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Sematilide $Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Sematilide $n<-exp(m_n$m_n)

## Calc out Sematilide  response at xfree and create Quantile frames

Sematilide $Pred_Fold<- (DR_up(y0=1, x0=Sematilide $x0, Emax=Sematilide $Emax, n=Sematilide $n, x=((Invivo$xfree[10])/QT_Parameters$FracFreeMedia[8])))*(scale_factor)

Sematilide $Pred_Frac<- (DR_up_Frac(y0=1, x0=Sematilide $x0, Emax=Sematilide $Emax, n=Sematilide $n, x=((Invivo$xfree[10])/QT_Parameters$FracFreeMedia[8])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Sematilide $Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Sematilide $Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Sematilide $Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Sematilide $Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Sematilide $Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Sematilide $Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Sematilide )

### Sematilide Linear  ###

##set up frame with parameters

file.sources=list.files(pattern="*_35_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_35_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Sematilide <- as.data.frame(exp(m_y0$m_y0))

rownames(Sematilide )<- c(1:8000)

colnames(Sematilide )<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Sematilide $x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Sematilide $Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Sematilide $n<-exp(m_n$m_n)

## Calc out Sematilide  response at xfree and create Quantile frames

Sematilide $Pred_Fold<- (DR_up(y0=1, x0=Sematilide $x0, Emax=Sematilide $Emax, n=Sematilide $n, x=((Invivo$xfree[11])/QT_Parameters$FracFreeMedia[8])))*(scale_factor)

Sematilide $Pred_Frac<- (DR_up_Frac(y0=1, x0=Sematilide $x0, Emax=Sematilide $Emax, n=Sematilide $n, x=((Invivo$xfree[11])/QT_Parameters$FracFreeMedia[8])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Sematilide $Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Sematilide $Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Sematilide $Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Sematilide $Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Sematilide $Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Sematilide $Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Sematilide )

### Sotalol  ###

##set up frame with parameters

file.sources=list.files(pattern="*_29_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_29_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Sotalol <- as.data.frame(exp(m_y0$m_y0))

rownames(Sotalol )<- c(1:8000)

colnames(Sotalol )<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Sotalol $x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Sotalol $Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Sotalol $n<-exp(m_n$m_n)

## Calc out Sotalol  response at xfree and create Quantile frames

Sotalol $Pred_Fold<- (DR_up(y0=1, x0=Sotalol $x0, Emax=Sotalol $Emax, n=Sotalol $n, x=((Invivo$xfree[12])/QT_Parameters$FracFreeMedia[9])))*(scale_factor)

Sotalol $Pred_Frac<- (DR_up_Frac(y0=1, x0=Sotalol $x0, Emax=Sotalol $Emax, n=Sotalol $n, x=((Invivo$xfree[12])/QT_Parameters$FracFreeMedia[9])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Sotalol $Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Sotalol $Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Sotalol $Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Sotalol $Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Sotalol $Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Sotalol $Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Sotalol )

### Vernacalant     ###

##set up frame with parameters

file.sources=list.files(pattern="*_12_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_12_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Vernacalant <- as.data.frame(exp(m_y0$m_y0))

rownames(Vernacalant)<- c(1:8000)

colnames(Vernacalant)<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Vernacalant$x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Vernacalant$Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Vernacalant$n<-exp(m_n$m_n)

## Calc out Vernacalant  response at xfree and create Quantile frames

Vernacalant $Pred_Fold<- (DR_up(y0=1, x0=Vernacalant $x0, Emax=Vernacalant $Emax, n=Vernacalant $n, x=((Invivo$xfree[13])/QT_Parameters$FracFreeMedia[10])))*(scale_factor)

Vernacalant $Pred_Frac<- (DR_up_Frac(y0=1, x0=Vernacalant $x0, Emax=Vernacalant $Emax, n=Vernacalant $n, x=((Invivo$xfree[13])/QT_Parameters$FracFreeMedia[10])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Vernacalant $Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Vernacalant $Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Vernacalant $Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Vernacalant $Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Vernacalant $Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Vernacalant $Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Vernacalant )

### Cabazitaxel    ###

##set up frame with parameters

file.sources=list.files(pattern="*_21_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_21_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Cabazitaxel<- as.data.frame(exp(m_y0$m_y0))

rownames(Cabazitaxel)<- c(1:8000)

colnames(Cabazitaxel)<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Cabazitaxel$x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Cabazitaxel$Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Cabazitaxel$n<-exp(m_n$m_n)

## Calc out Cabazitaxel   response at xfree and create Quantile frames

Cabazitaxel  $Pred_Fold<- (DR_up(y0=1, x0=Cabazitaxel  $x0, Emax=Cabazitaxel  $Emax, n=Cabazitaxel  $n, x=((Invivo$xfree[14])/QT_Parameters$FracFreeMedia[11])))*(scale_factor)

Cabazitaxel  $Pred_Frac<- (DR_up_Frac(y0=1, x0=Cabazitaxel  $x0, Emax=Cabazitaxel  $Emax, n=Cabazitaxel  $n, x=((Invivo$xfree[14])/QT_Parameters$FracFreeMedia[11])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Cabazitaxel  $Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Cabazitaxel  $Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Cabazitaxel  $Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Cabazitaxel  $Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Cabazitaxel  $Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Cabazitaxel  $Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Cabazitaxel  )

### Lamotrigine     ###

##set up frame with parameters

file.sources=list.files(pattern="*_14_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_14_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Lamotrigine   <- as.data.frame(exp(m_y0$m_y0))

rownames(Lamotrigine   )<- c(1:8000)

colnames(Lamotrigine   )<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Lamotrigine   $x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Lamotrigine   $Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Lamotrigine   $n<-exp(m_n$m_n)

## Calc out Lamotrigine    response at xfree and create Quantile frames

Lamotrigine   $Pred_Fold<- (DR_up(y0=1, x0=Lamotrigine   $x0, Emax=Lamotrigine   $Emax, n=Lamotrigine   $n, x=((Invivo$xfree[15])/QT_Parameters$FracFreeMedia[12])))*(scale_factor)

Lamotrigine   $Pred_Frac<- (DR_up_Frac(y0=1, x0=Lamotrigine   $x0, Emax=Lamotrigine   $Emax, n=Lamotrigine   $n, x=((Invivo$xfree[15])/QT_Parameters$FracFreeMedia[12])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Lamotrigine   $Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Lamotrigine   $Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Lamotrigine   $Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Lamotrigine   $Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Lamotrigine   $Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Lamotrigine   $Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Lamotrigine   )

### Mifepristone     ###

##set up frame with parameters

file.sources=list.files(pattern="*_1_dat.R")

sapply(file.sources,source,.GlobalEnv)

files = list.files(pattern="*_1_samples_*")

for (l in 1:length(files)) assign(files[l], read.csv(files[l]))

m_y0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004,c(8,12)],(eval(parse(text=files[2])))[2005:4004, c(8,12)], (eval(parse(text=files[3])))[2005:4004,c(8,12)], (eval(parse(text=files[4])))[2005:4004,c(8,12)]))

Mifepristone   <- as.data.frame(exp(m_y0$m_y0))

rownames(Mifepristone   )<- c(1:8000)

colnames(Mifepristone   )<- "y0"

m_x0<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(9,13)], (eval(parse(text=files[2])))[2005:4004, c(9,13)], (eval(parse(text=files[3])))[2005:4004, c(9,13)], (eval(parse(text=files[4])))[2005:4004, c(9,13)]))

Mifepristone   $x0<-exp(m_x0$m_x0)

m_Emax<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(10,14)], (eval(parse(text=files[2])))[2005:4004,c(10,14)], (eval(parse(text=files[3])))[2005:4004, c(10,14)], (eval(parse(text=files[4])))[2005:4004, c(10,14)]))

Mifepristone$Emax<-exp(m_Emax$m_Emax)

m_n<-as.data.frame(rbind((eval(parse(text=files[1])))[2005:4004, c(11,15)], (eval(parse(text=files[2])))[2005:4004, c(11,15)], (eval(parse(text=files[3])))[2005:4004, c(11,15)], (eval(parse(text=files[4])))[2005:4004, c(11,15)]))

Mifepristone   $n<-exp(m_n$m_n)

## Calc out Mifepristone    response at xfree and create Quantile frames

Mifepristone   $Pred_Fold<- (DR_up(y0=1, x0=Mifepristone   $x0, Emax=Mifepristone   $Emax, n=Mifepristone   $n, x=((Invivo$xfree[16])/QT_Parameters$FracFreeMedia[13])))*(scale_factor)

Mifepristone   $Pred_Frac<- (DR_up_Frac(y0=1, x0=Mifepristone   $x0, Emax=Mifepristone   $Emax, n=Mifepristone   $n, x=((Invivo$xfree[16])/QT_Parameters$FracFreeMedia[13])))*(scale_factor)

Temp_Quantiles_Fold<- data.frame(p50= median(Mifepristone   $Pred_Fold))

Temp_Quantiles_Fold$p975<- quantile(Mifepristone   $Pred_Fold, prob=0.975)

Temp_Quantiles_Fold$p25<- quantile(Mifepristone   $Pred_Fold, prob=0.025)

Temp_Quantiles_Frac<- data.frame(p50=median(Mifepristone   $Pred_Frac))

Temp_Quantiles_Frac$p975<- quantile(Mifepristone   $Pred_Frac, prob=0.975)

Temp_Quantiles_Frac$p25<-quantile(Mifepristone   $Pred_Frac, prob=0.025)

Quantiles_Fold<-rbind(Quantiles_Fold, Temp_Quantiles_Fold)

Quantiles_Frac<- rbind(Quantiles_Frac, Temp_Quantiles_Frac)

rm(list=files, Mifepristone   )

#################### End Chems #############################

## Combine frames and export as cvm

colnames(Quantiles_Fold)<- c("Fold_p50", "Fold_p97.5", "Fold_p2.5")

colnames(Quantiles_Frac)<- c("Frac_p50", "Frac_p97.5", "Frac_p2.5")

Invivo<-cbind(Invivo, Quantiles_Fold, Quantiles_Frac)

write.csv(Invivo, "In_Vivo_Frame_rescaled.csv")

##plot

ggplot(subset(Invivo, PredFracChange>0), aes(x=PredFracChange, y=Frac_p50))+
  
  geom_point(color="red", size=3)+
  
  geom_errorbar(aes(ymin= Frac_p2.5, ymax= Frac_p97.5), color="red")+
  
  scale_x_log10(limits=c(3e-4,2))+
  
  scale_y_log10(limits=c(3e-4,2))+
  
  geom_abline(intercept=0, slope=1, color="grey60")+
  
  geom_abline(intercept=0.5, slope=1, color="grey60", linetype="dashed")+
  
  geom_abline(intercept=-0.5, slope=1, color="grey60", linetype="dashed")+
  
  ggtitle("Frac Change @ Cmax")+
  
  labs(x="log observed", y="log predicted")

ggsave("Observed_vs_Predicted_Cmax_Frac_rescaled.pdf", plot=last_plot(), device="pdf", height=7.5, width = 7.5)

ggsave("Observed_vs_Predicted_Cmax_Frac_rescaled.png", plot=last_plot(), device="png", height=4.25, width = 4.25)



  