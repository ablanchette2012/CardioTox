setwd("C:\\Users\\ablanchette\\Documents\\Scripts\\QTIVIVEposnegdrugs")

# ProbGEQ10ms.R
library(ggplot2)
library(viridis)
# created by source("MakeInVivoInVitrodata-ExperimentalFreeFracSamples.R")
load("InVivoInVitrodat-ExperimentalFreeFracSamples.Rdata")
celllines<-as.character(read.csv("Cell.lines.csv",as.is=TRUE)[,2])
names(celllines)<-1:length(celllines)
chemmap<-read.csv("IVIVEchems.csv",as.is=TRUE)

add_qtc_func <-function(simframe_in,qtc0.m=421.5,qtc0.sd=0 
                        # For random individual, use (458.5-421.5)/qnorm(0.95) 
                        ) {
  n <- nrow(simframe_in)
  qtc0<-qtc0.m + qtc0.sd*rnorm(n)
  simframe_out<-simframe_in
  simframe_out$qtc <- simframe_out$PredFrac * qtc0
  simframe_out
}

qtc_geq_func <- function(qtc,qtcstar=10) {
  p_geq <- sum(qtc>=qtcstar)/length(qtc)
  p_geq
}

probgeq10.df <- data.frame()
for (j in 1:length(data.samples.list)) {
  # Population Median  
  simframe_zero_conc<-data.samples.list[[j]]$simframe_zero_conc
  simframe_out <- add_qtc_func(simframe_zero_conc)
  probgeq10 <- aggregate(qtc~freeconcentration, data=simframe_out, qtc_geq_func)
  probgeq10$safemin<- -0.05
  probgeq10$safemax<- 0.05
  probgeq10$safemin[probgeq10$qtc > 0.05] <- NA
  probgeq10$safemax[probgeq10$qtc > 0.05] <- NA
  probgeq10$InVitroconcentration <- probgeq10$freeconcentration / data.samples.list[[j]]$InVivoFree$FracFreeMedia
  probgeq10$InVivoconcentration <- probgeq10$freeconcentration / data.samples.list[[j]]$InVivoFree$FracFreePlasma
  probgeq10$Chemical.num <- data.samples.list[[j]]$chemical.num 
  probgeq10$Chemical.name <- data.samples.list[[j]]$chemical.name
  probgeq10$InVitroModel <- "Population median"
  probgeq10$CmaxFreeInVivo.x<-max(data.samples.list[[j]]$InVivo$xfree)
  probgeq10$CmaxInVivo.x<-max(data.samples.list[[j]]$InVivo$xfree)/data.samples.list[[j]]$InVivoFree$FracFreePlasma
  probgeq10$CmaxInVivo.y<-1
  probgeq10.df <- rbind(probgeq10.df,probgeq10)
  # Standard donor
  simframe_27_conc<-data.samples.list[[j]]$simframe_27_conc
  simframe_out <- add_qtc_func(simframe_27_conc)
  probgeq10 <- aggregate(qtc~freeconcentration, data=simframe_out, qtc_geq_func)
  probgeq10$safemin<- -0.025
  probgeq10$safemax<- 0.05
  probgeq10$safemin[probgeq10$qtc > 0.05] <- NA
  probgeq10$safemax[probgeq10$qtc > 0.05] <- NA
  probgeq10$InVitroconcentration <- probgeq10$freeconcentration / data.samples.list[[j]]$InVivoFree$FracFreeMedia
  probgeq10$InVivoconcentration <- probgeq10$freeconcentration / data.samples.list[[j]]$InVivoFree$FracFreePlasma
  probgeq10$Chemical.num <- data.samples.list[[j]]$chemical.num 
  probgeq10$Chemical.name <- data.samples.list[[j]]$chemical.name
  probgeq10$InVitroModel <- "Standard donor (1434)"
  probgeq10$CmaxFreeInVivo.x<-max(data.samples.list[[j]]$InVivo$xfree)
  probgeq10$CmaxInVivo.x<-max(data.samples.list[[j]]$InVivo$xfree)/data.samples.list[[j]]$InVivoFree$FracFreePlasma
  probgeq10$CmaxInVivo.y<-1
  probgeq10.df <- rbind(probgeq10.df,probgeq10)
  # Random individual donor
  simframeconc<-data.samples.list[[j]]$simframeconc
  simframe_out <- add_qtc_func(simframeconc,qtc0.sd=(458.5-421.5)/qnorm(0.95))
  probgeq10 <- aggregate(qtc~freeconcentration, data=simframe_out, qtc_geq_func)
  probgeq10$safemin<- 0
  probgeq10$safemax<- 0.05
  probgeq10$safemin[probgeq10$qtc > 0.05] <- NA
  probgeq10$safemax[probgeq10$qtc > 0.05] <- NA
  probgeq10$InVitroconcentration <- probgeq10$freeconcentration / data.samples.list[[j]]$InVivoFree$FracFreeMedia
  probgeq10$InVivoconcentration <- probgeq10$freeconcentration / data.samples.list[[j]]$InVivoFree$FracFreePlasma
  probgeq10$Chemical.num <- data.samples.list[[j]]$chemical.num 
  probgeq10$Chemical.name <- data.samples.list[[j]]$chemical.name
  probgeq10$InVitroModel <- "Random individual donor"
  probgeq10$CmaxFreeInVivo.x<-max(data.samples.list[[j]]$InVivo$xfree)
  probgeq10$CmaxInVivo.x<-max(data.samples.list[[j]]$InVivo$xfree)/data.samples.list[[j]]$InVivoFree$FracFreePlasma
  probgeq10$CmaxInVivo.y<-1
  probgeq10.df <- rbind(probgeq10.df,probgeq10)
}
names(probgeq10.df)[2]<-"PrQTc10"

probgeq10.df$Chemical.name <- factor(probgeq10.df$Chemical.name,
                                 levels=c(chemmap$Chemical.name))
probgeq10.df$InVitroModel <- factor(probgeq10.df$InVitroModel,
                                    levels=c("Population median",
                                             "Standard donor (1434)",
                                             "Random individual donor"))

cols=c("Standard donor (1434)"=magma(5)[2],
       "Population median"=magma(5)[4],
       "Random individual donor"=magma(5)[5])

ltys<-c("Population median"=1,
        "Standard donor (1434)"=2,
        "Random individual donor"=3)

p<-ggplot(probgeq10.df)+
  scale_x_log10()+
  scale_fill_manual(name="Below Threshold",values=cols)+
  scale_linetype_manual(name="In Vitro Prediction",values=ltys)+
  scale_shape_manual(name="",values=c("In Vivo Cmax"=25))+
  annotation_logticks(color="grey50",side="t")+
  geom_ribbon(aes(x=InVivoconcentration,ymin=safemin,ymax=safemax,fill=InVitroModel))+
  geom_line(aes(x=InVivoconcentration,y=PrQTc10,linetype=InVitroModel))+
  geom_point(aes(x=CmaxInVivo.x,y=CmaxInVivo.y,shape="In Vivo Cmax"),size=6,fill="black")+
  xlab("Plasma Concentration (uM)")+
  ylab("Probability that Change QTc >= 10 ms")+
  coord_cartesian(expand=FALSE)+
  theme_bw()+
  facet_wrap(~Chemical.name)
print(p)
ggsave("ProbDeltaQTc10.pdf",p,height=7,width=9,dpi=600)
write.csv(probgeq10.df,file="ProbDeltaQTc10.csv",row.names = FALSE)
