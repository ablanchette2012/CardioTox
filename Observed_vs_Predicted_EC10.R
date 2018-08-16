 
###set WD and load packages

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\All")

library(ggplot2)

InVivoPos<-read.csv("InVivoPred.pos.csv")

InVivoAll<- read.csv("InVivoPred.all.csv")

##### generate EC10, EC05, EC01 values for in vivo and in vitro data #####

##Cisapride

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\Cisapride\\decay_rise_up\\invivo_invitro")

Cisapride<-subset(InVivoPos, Chemical.name=="cisapride")

VivoEC<- data.frame(EC10=(approx(x=Cisapride$PredFracChange, y=Cisapride$xfree, xout=0.1)),EC05=(approx(x=Cisapride$PredFracChange, y=Cisapride$xfree, xout=0.05)), EC01=(approx(x=Cisapride$PredFracChange, y=Cisapride$xfree, xout=0.01)))

Cisapride_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

VitroEC<- data.frame(EC10_p50=(approx(x=Cisapride_vitro$p50, y=Cisapride_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=Cisapride_vitro$p97.5,y=Cisapride_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=Cisapride_vitro$p2.5, y=Cisapride_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=Cisapride_vitro$p50, y=Cisapride_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=Cisapride_vitro$p97.5,y=Cisapride_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=Cisapride_vitro$p2.5, y=Cisapride_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=Cisapride_vitro$p50, y=Cisapride_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=Cisapride_vitro$p97.5,y=Cisapride_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=Cisapride_vitro$p2.5, y=Cisapride_vitro$freeconcentration, xout=0.01)) )

rm(Cisapride, Cisapride_vitro)

##Citalopram

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\Citalopram\\decay_rise_up\\invivo_invitro")

Citalopram<-subset(InVivoPos, Chemical.name=="citalopram hydrobromide")

TempVivo<- data.frame(EC10=(approx(x=Citalopram$PredFracChange, y=Citalopram$xfree, xout=0.1)), EC05=(approx(x=Citalopram$PredFracChange, y=Citalopram$xfree, xout=0.05)), EC01=(approx(x=Citalopram$PredFracChange, y=Citalopram$xfree, xout=0.01)))

VivoEC<- rbind(VivoEC,TempVivo)

Citalopram_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

TempVitro<- data.frame(EC10_p50=(approx(x=Citalopram_vitro$p50, y=Citalopram_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=Citalopram_vitro$p97.5,y=Citalopram_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=Citalopram_vitro$p2.5, y=Citalopram_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=Citalopram_vitro$p50, y=Citalopram_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=Citalopram_vitro$p97.5,y=Citalopram_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=Citalopram_vitro$p2.5, y=Citalopram_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=Citalopram_vitro$p50, y=Citalopram_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=Citalopram_vitro$p97.5,y=Citalopram_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=Citalopram_vitro$p2.5, y=Citalopram_vitro$freeconcentration, xout=0.01)) )

VitroEC<- rbind(VitroEC, TempVitro)

rm(Citalopram, Citalopram_vitro)

##Disopyramide Linear

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\Disopyramide\\decay_rise_up\\invivo_invitro")

Disopyramide<-subset(InVivoPos, Chemical.name=="disopyramide phosphate" & Model == "in vivo Linear")

TempVivo<- data.frame(EC10=(approx(x=Disopyramide $PredFracChange, y=Disopyramide $xfree, xout=0.1)),EC05=(approx(x=Disopyramide$PredFracChange, y=Disopyramide$xfree, xout=0.05)), EC01=(approx(x=Disopyramide $PredFracChange, y=Disopyramide $xfree, xout=0.01)))

VivoEC<- rbind(VivoEC,TempVivo)

Disopyramide_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

TempVitro<- data.frame(EC10_p50=(approx(x=Disopyramide_vitro$p50, y=Disopyramide_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=Disopyramide_vitro$p97.5,y=Disopyramide_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=Disopyramide_vitro$p2.5, y=Disopyramide_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=Disopyramide_vitro$p50, y=Disopyramide_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=Disopyramide_vitro$p97.5,y=Disopyramide_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=Disopyramide_vitro$p2.5, y=Disopyramide_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=Disopyramide_vitro$p50, y=Disopyramide_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=Disopyramide_vitro$p97.5,y=Disopyramide_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=Disopyramide_vitro$p2.5, y=Disopyramide_vitro$freeconcentration, xout=0.01)) )

VitroEC<- rbind(VitroEC, TempVitro)

rm(Disopyramide, Disopyramide_vitro)


##Disopyramide Hill

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\Disopyramide\\decay_rise_up\\invivo_invitro")

Disopyramide<-subset(InVivoPos, Chemical.name=="disopyramide phosphate" & Model == "in vivo Hill")

TempVivo<- data.frame(EC10=(approx(x=Disopyramide $PredFracChange, y=Disopyramide $xfree, xout=0.1)),EC05=(approx(x=Disopyramide$PredFracChange, y=Disopyramide$xfree, xout=0.05)), EC01=(approx(x=Disopyramide $PredFracChange, y=Disopyramide $xfree, xout=0.01)))

VivoEC<- rbind(VivoEC,TempVivo)

Disopyramide_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

TempVitro<- data.frame(EC10_p50=(approx(x=Disopyramide_vitro$p50, y=Disopyramide_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=Disopyramide_vitro$p97.5,y=Disopyramide_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=Disopyramide_vitro$p2.5, y=Disopyramide_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=Disopyramide_vitro$p50, y=Disopyramide_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=Disopyramide_vitro$p97.5,y=Disopyramide_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=Disopyramide_vitro$p2.5, y=Disopyramide_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=Disopyramide_vitro$p50, y=Disopyramide_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=Disopyramide_vitro$p97.5,y=Disopyramide_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=Disopyramide_vitro$p2.5, y=Disopyramide_vitro$freeconcentration, xout=0.01)) )

VitroEC<- rbind(VitroEC, TempVitro)

rm(Disopyramide, Disopyramide_vitro)

##Dofetilide 

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\Dofetilide\\decay_rise_up\\invivo_invitro")

Dofetilide<-subset(InVivoPos, Chemical.name=="dofetilide")

TempVivo<- data.frame(EC10=(approx(x=Dofetilide $PredFracChange, y=Dofetilide $xfree, xout=0.1)),EC05=(approx(x=Dofetilide$PredFracChange, y=Dofetilide$xfree, xout=0.05)), EC01=(approx(x=Dofetilide $PredFracChange, y=Dofetilide $xfree, xout=0.01)))

VivoEC<- rbind(VivoEC,TempVivo)

Dofetilide_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

TempVitro<- data.frame(EC10_p50=(approx(x=Dofetilide_vitro$p50, y=Dofetilide_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=Dofetilide_vitro$p97.5,y=Dofetilide_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=Dofetilide_vitro$p2.5, y=Dofetilide_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=Dofetilide_vitro$p50, y=Dofetilide_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=Dofetilide_vitro$p97.5,y=Dofetilide_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=Dofetilide_vitro$p2.5, y=Dofetilide_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=Dofetilide_vitro$p50, y=Dofetilide_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=Dofetilide_vitro$p97.5,y=Dofetilide_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=Dofetilide_vitro$p2.5, y=Dofetilide_vitro$freeconcentration, xout=0.01)) )

VitroEC<- rbind(VitroEC, TempVitro)

rm(Dofetilide, Dofetilide_vitro)

##Moxifloxacin

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\Moxifloxacin\\decay_rise_up\\invivo_invitro")

Moxifloxacin<-subset(InVivoPos, Chemical.name=="moxifloxacin hydrochloride")

TempVivo<- data.frame(EC10=(approx(x=Moxifloxacin$PredFracChange, y=Moxifloxacin$xfree, xout=0.1)),EC05=(approx(x=Moxifloxacin$PredFracChange, y=Moxifloxacin$xfree, xout=0.05)), EC01=(approx(x=Moxifloxacin$PredFracChange, y=Moxifloxacin$xfree, xout=0.01)))

VivoEC<- rbind(VivoEC,TempVivo)

Moxifloxacin_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

TempVitro<- data.frame(EC10_p50=(approx(x=Moxifloxacin_vitro$p50, y=Moxifloxacin_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=Moxifloxacin_vitro$p97.5,y=Moxifloxacin_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=Moxifloxacin_vitro$p2.5, y=Moxifloxacin_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=Moxifloxacin_vitro$p50, y=Moxifloxacin_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=Moxifloxacin_vitro$p97.5,y=Moxifloxacin_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=Moxifloxacin_vitro$p2.5, y=Moxifloxacin_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=Moxifloxacin_vitro$p50, y=Moxifloxacin_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=Moxifloxacin_vitro$p97.5,y=Moxifloxacin_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=Moxifloxacin_vitro$p2.5, y=Moxifloxacin_vitro$freeconcentration, xout=0.01)) )

VitroEC<- rbind(VitroEC, TempVitro)

rm(Moxifloxacin, Moxifloxacin_vitro)

##N-acetylprocainamide Linear

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\N-acetylprocainamide\\decay_rise_up\\invivo_invitro")

N_acetylprocainamide<-subset(InVivoPos, Chemical.name=="N-acetylprocainamide" & Model == "in vivo Linear")

TempVivo<- data.frame(EC10=(approx(x=N_acetylprocainamide$PredFracChange, y=N_acetylprocainamide$xfree, xout=0.1)),EC05=(approx(x=N_acetylprocainamide$PredFracChange, y=N_acetylprocainamide$xfree, xout=0.05)), EC01=(approx(x=N_acetylprocainamide$PredFracChange, y=N_acetylprocainamide$xfree, xout=0.01)))

VivoEC<- rbind(VivoEC,TempVivo)

N_acetylprocainamide_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

TempVitro<- data.frame(EC10_p50=(approx(x=N_acetylprocainamide_vitro$p50, y=N_acetylprocainamide_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=N_acetylprocainamide_vitro$p97.5,y=N_acetylprocainamide_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=N_acetylprocainamide_vitro$p2.5, y=N_acetylprocainamide_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=N_acetylprocainamide_vitro$p50, y=N_acetylprocainamide_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=N_acetylprocainamide_vitro$p97.5,y=N_acetylprocainamide_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=N_acetylprocainamide_vitro$p2.5, y=N_acetylprocainamide_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=N_acetylprocainamide_vitro$p50, y=N_acetylprocainamide_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=N_acetylprocainamide_vitro$p97.5,y=N_acetylprocainamide_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=N_acetylprocainamide_vitro$p2.5, y=N_acetylprocainamide_vitro$freeconcentration, xout=0.01)) )

VitroEC<- rbind(VitroEC, TempVitro)

rm(N_acetylprocainamide, N_acetylprocainamide_vitro)

##N-acetylprocainamide Hill

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\N-acetylprocainamide\\decay_rise_up\\invivo_invitro")

N_acetylprocainamide<-subset(InVivoPos, Chemical.name=="N-acetylprocainamide" & Model == "in vivo Hill")

TempVivo<- data.frame(EC10=(approx(x=N_acetylprocainamide$PredFracChange, y=N_acetylprocainamide$xfree, xout=0.1)),EC05=(approx(x=N_acetylprocainamide$PredFracChange, y=N_acetylprocainamide$xfree, xout=0.05)), EC01=(approx(x=N_acetylprocainamide$PredFracChange, y=N_acetylprocainamide$xfree, xout=0.01)))

VivoEC<- rbind(VivoEC,TempVivo)

N_acetylprocainamide_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

TempVitro<- data.frame(EC10_p50=(approx(x=N_acetylprocainamide_vitro$p50, y=N_acetylprocainamide_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=N_acetylprocainamide_vitro$p97.5,y=N_acetylprocainamide_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=N_acetylprocainamide_vitro$p2.5, y=N_acetylprocainamide_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=N_acetylprocainamide_vitro$p50, y=N_acetylprocainamide_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=N_acetylprocainamide_vitro$p97.5,y=N_acetylprocainamide_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=N_acetylprocainamide_vitro$p2.5, y=N_acetylprocainamide_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=N_acetylprocainamide_vitro$p50, y=N_acetylprocainamide_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=N_acetylprocainamide_vitro$p97.5,y=N_acetylprocainamide_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=N_acetylprocainamide_vitro$p2.5, y=N_acetylprocainamide_vitro$freeconcentration, xout=0.01)) )

VitroEC<- rbind(VitroEC, TempVitro)

rm(N_acetylprocainamide, N_acetylprocainamide_vitro)

##Quinidine

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\Quinidine\\decay_rise_up\\invivo_invitro")

Quinidine<-subset(InVivoPos, Chemical.name=="quinidine sulfate")

TempVivo<- data.frame(EC10=(approx(x=Quinidine$PredFracChange, y=Quinidine$xfree, xout=0.1)),EC05=(approx(x=Quinidine$PredFracChange, y=Quinidine$xfree, xout=0.05)), EC01=(approx(x=Quinidine$PredFracChange, y=Quinidine$xfree, xout=0.01)))

VivoEC<- rbind(VivoEC,TempVivo)

Quinidine_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

TempVitro<- data.frame(EC10_p50=(approx(x=Quinidine_vitro$p50, y=Quinidine_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=Quinidine_vitro$p97.5,y=Quinidine_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=Quinidine_vitro$p2.5, y=Quinidine_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=Quinidine_vitro$p50, y=Quinidine_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=Quinidine_vitro$p97.5,y=Quinidine_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=Quinidine_vitro$p2.5, y=Quinidine_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=Quinidine_vitro$p50, y=Quinidine_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=Quinidine_vitro$p97.5,y=Quinidine_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=Quinidine_vitro$p2.5, y=Quinidine_vitro$freeconcentration, xout=0.01)) )

VitroEC<- rbind(VitroEC, TempVitro)

rm(Quinidine, Quinidine_vitro)

##Sematilide Linear

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\Sematilide\\decay_rise_up\\invivo_invitro")

Sematilide<-subset(InVivoPos, Chemical.name=="sematilide" & Model=="in vivo Linear")

TempVivo<- data.frame(EC10=(approx(x=Sematilide$PredFracChange, y=Sematilide$xfree, xout=0.1)),EC05=(approx(x=Sematilide$PredFracChange, y=Sematilide$xfree, xout=0.05)), EC01=(approx(x=Sematilide$PredFracChange, y=Sematilide$xfree, xout=0.01)))

VivoEC<- rbind(VivoEC,TempVivo)

Sematilide_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

TempVitro<- data.frame(EC10_p50=(approx(x=Sematilide_vitro$p50, y=Sematilide_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=Sematilide_vitro$p97.5,y=Sematilide_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=Sematilide_vitro$p2.5, y=Sematilide_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=Sematilide_vitro$p50, y=Sematilide_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=Sematilide_vitro$p97.5,y=Sematilide_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=Sematilide_vitro$p2.5, y=Sematilide_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=Sematilide_vitro$p50, y=Sematilide_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=Sematilide_vitro$p97.5,y=Sematilide_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=Sematilide_vitro$p2.5, y=Sematilide_vitro$freeconcentration, xout=0.01)) )

VitroEC<- rbind(VitroEC, TempVitro)

rm(Sematilide, Sematilide_vitro)

##Sematilide Hill

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\Sematilide\\decay_rise_up\\invivo_invitro")

Sematilide<-subset(InVivoPos, Chemical.name=="sematilide" & Model=="in vivo Hill")

TempVivo<- data.frame(EC10=(approx(x=Sematilide$PredFracChange, y=Sematilide$xfree, xout=0.1)),EC05=(approx(x=Sematilide$PredFracChange, y=Sematilide$xfree, xout=0.05)), EC01=(approx(x=Sematilide$PredFracChange, y=Sematilide$xfree, xout=0.01)))

VivoEC<- rbind(VivoEC,TempVivo)

Sematilide_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

TempVitro<- data.frame(EC10_p50=(approx(x=Sematilide_vitro$p50, y=Sematilide_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=Sematilide_vitro$p97.5,y=Sematilide_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=Sematilide_vitro$p2.5, y=Sematilide_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=Sematilide_vitro$p50, y=Sematilide_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=Sematilide_vitro$p97.5,y=Sematilide_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=Sematilide_vitro$p2.5, y=Sematilide_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=Sematilide_vitro$p50, y=Sematilide_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=Sematilide_vitro$p97.5,y=Sematilide_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=Sematilide_vitro$p2.5, y=Sematilide_vitro$freeconcentration, xout=0.01)) )

VitroEC<- rbind(VitroEC, TempVitro)

rm(Sematilide, Sematilide_vitro)

##Sotalol

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\Sotalol\\decay_rise_up\\invivo_invitro")

Sotalol<-subset(InVivoPos, Chemical.name=="sotalol")

TempVivo<- data.frame(EC10=(approx(x=Sotalol$PredFracChange, y=Sotalol$xfree, xout=0.1)),EC05=(approx(x=Sotalol$PredFracChange, y=Sotalol$xfree, xout=0.05)), EC01=(approx(x=Sotalol$PredFracChange, y=Sotalol$xfree, xout=0.01)))

VivoEC<- rbind(VivoEC,TempVivo)

Sotalol_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

TempVitro<- data.frame(EC10_p50=(approx(x=Sotalol_vitro$p50, y=Sotalol_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=Sotalol_vitro$p97.5,y=Sotalol_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=Sotalol_vitro$p2.5, y=Sotalol_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=Sotalol_vitro$p50, y=Sotalol_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=Sotalol_vitro$p97.5,y=Sotalol_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=Sotalol_vitro$p2.5, y=Sotalol_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=Sotalol_vitro$p50, y=Sotalol_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=Sotalol_vitro$p97.5,y=Sotalol_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=Sotalol_vitro$p2.5, y=Sotalol_vitro$freeconcentration, xout=0.01)) )

VitroEC<- rbind(VitroEC, TempVitro)

rm(Sotalol, Sotalol_vitro)

##Vernacalant

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\Vernacalant\\decay_rise_up\\invivo_invitro")

Vernacalant<-subset(InVivoPos, Chemical.name=="vernacalant")

TempVivo<- data.frame(EC10=(approx(x=Vernacalant$PredFracChange, y=Vernacalant$xfree, xout=0.1)),EC05=(approx(x=Vernacalant$PredFracChange, y=Vernacalant$xfree, xout=0.05)), EC01=(approx(x=Vernacalant$PredFracChange, y=Vernacalant$xfree, xout=0.01)))

VivoEC<- rbind(VivoEC,TempVivo)

Vernacalant_vitro<-read.csv("Invivo_Frame_Zero_Frac.csv")

TempVitro<- data.frame(EC10_p50=(approx(x=Vernacalant_vitro$p50, y=Vernacalant_vitro$freeconcentration, xout=0.1)), EC10_p97.5=(approx(x=Vernacalant_vitro$p97.5,y=Vernacalant_vitro$freeconcentration, xout=0.1)), EC10_p2.5=(approx(x=Vernacalant_vitro$p2.5, y=Vernacalant_vitro$freeconcentration, xout=0.1)),EC05_p50=(approx(x=Vernacalant_vitro$p50, y=Vernacalant_vitro$freeconcentration, xout=0.05)), EC05_p97.5=(approx(x=Vernacalant_vitro$p97.5,y=Vernacalant_vitro$freeconcentration, xout=0.05)), EC05_p2.5=(approx(x=Vernacalant_vitro$p2.5, y=Vernacalant_vitro$freeconcentration, xout=0.05)), EC01_p50=(approx(x=Vernacalant_vitro$p50, y=Vernacalant_vitro$freeconcentration, xout=0.01)), EC01_p97.5=(approx(x=Vernacalant_vitro$p97.5,y=Vernacalant_vitro$freeconcentration, xout=0.01)), EC01_p2.5=(approx(x=Vernacalant_vitro$p2.5, y=Vernacalant_vitro$freeconcentration, xout=0.01)) )

VitroEC<- rbind(VitroEC, TempVitro)

rm(Vernacalant, Vernacalant_vitro)

##### End Calcs ##### 

##Set up ECFrame

setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template\\All")

ECFrame<- data.frame(Chemical.name=InVivoAll$Chemical.name[1:13], Model=InVivoAll$Model[1:13], Invivo_EC10=VivoEC$EC10.y, Invivo_EC05=VivoEC$EC05.y, Invivo_EC01=VivoEC$EC01.y, Invitro_EC10_p50=VitroEC$EC10_p50.y, Invitro_EC10_p97.5=VitroEC$EC10_p97.5.y, Invitro_EC10_p2.5=VitroEC$EC10_p2.5.y, Invitro_EC05_p50=VitroEC$EC05_p50.y, Invitro_EC05_p97.5=VitroEC$EC05_p97.5.y, Invitro_EC05_p2.5=VitroEC$EC05_p2.5.y, Invitro_EC01_p50=VitroEC$EC01_p50.y, Invitro_EC01_p97.5=VitroEC$EC01_p97.5.y, Invitro_EC01_p2.5=VitroEC$EC01_p2.5.y)


### ouput ECFrame as csv and generate plots

ECFrame[is.na(ECFrame)]<- 0.0001

write.csv(ECFrame, "ECFrame.csv")

##EC10

ggplot(ECFrame, aes(x=Invivo_EC10, y=Invitro_EC10_p50))+
  
  geom_point(size=3, color="red")+
  
  geom_errorbar(aes(ymin=Invitro_EC10_p2.5, ymax=Invitro_EC10_p97.5), color="red")+
  
  scale_x_log10(limits=c(3e-04,100))+
  
  scale_y_log10(limits=c(3e-04,100))+

  geom_abline(intercept=0, slope=1, color="grey60")+
  
  geom_abline(intercept=0.5, slope=1, color="grey60", linetype="dashed")+
  
  geom_abline(intercept=-0.5, slope=1, color="grey60", linetype="dashed")+
  
  ggtitle("Observed vs. Predicted EC10")+
  
  labs(x="Observed EC10", y="Predicted EC10")

ggsave("Observed_vs_Predicted_EC10.pdf", plot = last_plot(), device="pdf", height=7.5, width=7.5)

ggsave("Observed_vs_Predicted_EC10.png", plot = last_plot(), device="png", height=4.25, width = 4.25)

##EC05

ggplot(ECFrame, aes(x=Invivo_EC05, y=Invitro_EC05_p50))+
  
  geom_point(size=3, color="red")+
  
  geom_errorbar(aes(ymin=Invitro_EC05_p2.5, ymax=Invitro_EC05_p97.5), color="red")+
  
  scale_x_log10(limits=c(3e-04,100))+
  
  scale_y_log10(limits=c(3e-04,100))+
  
  geom_abline(intercept=0, slope=1, color="grey60")+
  
  geom_abline(intercept=0.5, slope=1, color="grey60", linetype="dashed")+
  
  geom_abline(intercept=-0.5, slope=1, color="grey60", linetype="dashed")+
  
  ggtitle("Observed vs. Predicted EC05")+
  
  labs(x="Observed EC05", y="Predicted EC05")

ggsave("Observed_vs_Predicted_EC05.pdf", plot = last_plot(), device="pdf", height=7.5, width=7.5)

ggsave("Observed_vs_Predicted_EC05.png", plot = last_plot(), device="png", height=4.25, width = 4.25)

##EC01

ggplot(ECFrame, aes(x=Invivo_EC01, y=Invitro_EC01_p50))+
  
  geom_point(size=3, color="red")+
  
  geom_errorbar(aes(ymin=Invitro_EC01_p2.5, ymax=Invitro_EC01_p97.5), color="red")+
  
  scale_x_log10(limits=c(3e-04,100))+
  
  scale_y_log10(limits=c(3e-04,100))+
  
  geom_abline(intercept=0, slope=1, color="grey60")+
  
  geom_abline(intercept=0.5, slope=1, color="grey60", linetype="dashed")+
  
  geom_abline(intercept=-0.5, slope=1, color="grey60", linetype="dashed")+
  
  ggtitle("Observed vs. Predicted EC01")+
  
  labs(x="Observed EC01", y="Predicted EC01")

ggsave("Observed_vs_Predicted_EC01.pdf", plot = last_plot(), device="pdf", height=7.5, width=7.5)

ggsave("Observed_vs_Predicted_EC01.png", plot = last_plot(), device="png", height=4.25, width = 4.25)
