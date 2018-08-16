## Read in parms 

Parms<- read.csv("Parms.csv")

Parms.new<-subset(Parms, Chemical.Concentration==0)

ChemNumList<-Parms$Chemical.Number[!duplicated(Parms$Chemical.Number)]

indivlist<-Parms$Cell.line[!duplicated(Parms$Cell.line)]

for (Chemnum in ChemNumList){
  
  print(Chemnum)
  
  for (indiv in indivlist){
    
    print (indiv)
    
    parms.tmp<-subset(Parms, Chemical.Number==Chemnum & Cell.line==indiv)
    
    if (sum(parms.tmp$notch.frac, na.rm =TRUE)==0){
      
      Parms.new<-rbind(Parms.new, parms.tmp)
      
      }else{
      
          tmp.sub<-subset(parms.tmp, notch.frac>0)
      
          minconc<-min(tmp.sub$Chemical.Concentration)
      
          parms.tmp<-subset(parms.tmp, Chemical.Concentration<=minconc)
      
          Parms.new<-rbind(Parms.new, parms.tmp)
      }
  }
}
write.csv(Parms.new, "Nonblinded.peak.parms.sum.CO.csv")