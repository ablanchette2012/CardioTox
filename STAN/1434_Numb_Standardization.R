library(tidyr)

#read in peak parms and chemical.dat files
alldat<-read.table("Nonblinded.peak.parms.sum.dat",header=TRUE,as.is=TRUE,sep="\t")

chemicals<-read.table("chemicals.dat", header=FALSE, as.is=TRUE, sep = "\t")

colnames(chemicals)<- c("Chemical.Number", "Chemical.Name")

chemicals_1434<- read.table("chemicals_1434.dat", header=FALSE, as.is= TRUE, sep = "\t")

colnames(chemicals_1434)<- c("Chemical.Number", "Chemical.Name")

#subset alldat to isolate 1434 rows

chems1434<-subset(alldat, Cell.line==1434)

chems1434<- drop_na(chems1434)

chems1434_num<- (chems1434$Chemical.Number)

#match chemical numbers with 1434 names

chemxnames_1434<-data.frame()

for (chemnum in chems1434_num){
  
  chem_name<- data.frame(as.character(chemicals_1434[match(chemnum, chemicals_1434$Chemical.Number),2]))
  
  chemxnames_1434<-rbind(chemxnames_1434, chem_name)
}

chemxnames_1434<-data.frame(chems1434_num, chemxnames_1434)

colnames(chemxnames_1434)<-c("1434_Num", "1434_Name")

#match 1434 names with regular key numbers

chemxnames<- data.frame()

for (chemname in chemxnames_1434$`1434_Name`) {
  
  chem_num<-data.frame(chemicals[match(chemname, chemicals$Chemical.Name),1])

  chemxnames<- rbind(chemxnames, chem_num)
  
}

colnames(chemxnames)<- "Std_Num"

tempframe<-data.frame()

for (num in chemxnames$Std_Num) {
  
  chem_name<- data.frame(as.character(chemicals[match(num, chemicals$Chemical.Number),2]))
  
  tempframe<- rbind(tempframe, chem_name)
}

chemxnames<-cbind(chemxnames, tempframe)

colnames(chemxnames)<-c("Std_Num", "Std_Name")

chemxnames<-cbind(chemxnames_1434,chemxnames)

write.csv(chemxnames, "Standard_and_1434_Names_Numbs.csv")

chems1434_new<-chems1434

chems1434_new$Chemical.Number<-chemxnames$Std_Num

for (i in 1:nrow(chems1434_new)) {
  
  rownum<- match(chems1434_new[i,8], alldat$baseline)
  
  print(rownum)
  
  alldat[rownum, 5]<- chems1434_new[i, 5]
  
}

write.table(alldat, "Nonblinded.peak.parms.stdnum.sum.dat")


