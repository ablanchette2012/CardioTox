Cmax.43<- read.csv("Cmax.43.csv", header = TRUE, as.is = TRUE)

Cmax.27<- read.csv("Cmax.27.csv", header = T, as.is = T)

cols<- c("Vivo", "Vitro")

compframe.43<- data.frame(Cmax.43$InVivo.FracChangeCmaxFree*100, Cmax.43$InVitro.FracChangeCmaxFree_p50.y*100)

colnames(compframe.43)<- cols

compframe.27<- data.frame(Cmax.27$InVivo.FracChangeCmaxFree*100, Cmax.27$InVitro.FracChangeCmaxFree_p50.y*100)

colnames(compframe.27)<- cols

Corr.43<- cor.test(compframe.43$Vivo, compframe.43$Vitro, method = "pearson", use = "complete.obs")

r2.43<- cor(compframe.43$Vivo, compframe.43$Vitro, use = "complete.obs")^2

Corr.27<- cor.test(compframe.27$Vivo, compframe.27$Vitro, method = "pearson", use = "complete.obs")
