library(ggplot2)

invivodf.pos <- read.csv("InVivoPred-ExperimentalFreeFrac.pos.csv")
invivodf.neg <- read.csv("InVivoPred-ExperimentalFreeFrac.neg.csv")

label_log10 <- function(x) parse(text = paste0('10^', log(x, 10)))

pdf("FoldChange.invivo-ExperimentalFreeFrac.pos.pdf",width=6,height=8)
print(ggplot(data=invivodf.pos,aes(x=xfree,y=PredFoldChange,linetype=Model)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA,size=1.5))+
#  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_line(size=1) + facet_wrap(~Chemical.name,scales="free_y",nrow=5) +
  scale_x_continuous("Free Concentration (uM)",
                     trans = 'log10',breaks=10^(-3:2),
                     labels = label_log10) +
  ylab("Fold-change in QTc"))
dev.off()

pdf("FoldChange.invivo-ExperimentalFreeFrac.neg.pdf",width=6,height=8*2/5)
print(ggplot(data=invivodf.neg,aes(x=xfree,y=PredFoldChange,linetype=Model)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA,size=1.5))+
#  theme(panel.background = element_rect(fill = "white", colour = "black")) +
  geom_line(size=1.5) + facet_wrap(~Chemical.name,nrow=2) +
  scale_x_continuous("Free Concentration (uM)",
                     trans = 'log10',breaks=10^(-3:2),
                     labels = label_log10) +
  ylab("Fold-change in QTc"))
dev.off()


pdf("FracChange.invivo-ExperimentalFreeFrac.pos.log.pdf",width=6,height=8)
print(ggplot(data=invivodf.pos,aes(x=xfree,y=PredFracChange,linetype=Model)) +
  theme_classic()+
  theme(panel.border = element_rect(color = "black", fill = NA,size=1.5))+
  geom_line(size=1) + facet_wrap(~Chemical.name,nrow=5) +
  scale_x_continuous("Free Concentration (uM)",
                     trans = 'log10',breaks=10^(-3:2),
                     labels = label_log10) +
  scale_y_continuous("Fractional change in QTc",
                     trans = 'log10',breaks=10^(-5:0),
                     labels = label_log10) )
dev.off()
