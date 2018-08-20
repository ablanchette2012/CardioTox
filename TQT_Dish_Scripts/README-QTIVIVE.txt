Each of the "*" files should be in this archive.  The "o" files are ones that are created.

Place all the files in this archive in its own directory, and then copy/move the "All" folder as a subfolder in that directory.  The only possible addition is that you might need to copy IVIVEchems.csv to the "All" subfolder.

* MakeInVitrodataprediction.R - takes modeling results from "All" subdirectory and creates data frames for observed and predicted that are saved in the following Rdata file
  o InVitrodatandprediction.Rdata
* MakeInVitrodatapredictionSmall.R - same as above, but smaller Rdata file
  o InVitrodatandpredictionSmall.Rdata
* Cell.lines.csv - mapping between cell line # and donor ID
* In Vitro Obs Pred_V2.Rmd - in vitro only observed vs. modeled - uses Rdata files above, as well as Cell.lines.csv
* QTposnegparameters-Experimental.csv - parameters for In Vivo models using experimental Free Fraction
* QTposnegparameters+Table1-Experimental.xlsx - corresponding Excel file
* MakeInVivoPred-ExperimentalFreeFrac.R - output is in vivo conc-resp models, converted to free concentration using Experimental values.  Creates:
  o InVivoPred-ExperimentalFreeFrac.pos.csv
  o InVivoPred-ExperimentalFreeFrac.neg.csv
* PlotInVivo-ExperimentalFreeFrac.R - plots of in vivo conc-resp using above CSV files.  Creates 3 PDF files
  o FoldChange.invivo-ExperimentalFreeFrac.pos.pdf
  o FoldChange.invivo-ExperimentalFreeFrac.neg.pdf
  o FracChange.invivo-ExperimentalFreeFrac.pos.log.pdf
* QTposnegparameters-Literature.csv - parameters for In Vivo models using experimental Free Fraction
* QTposnegparameters+Table1-LiteratureFreeFrac.xlsx - corresponding Excel file
* MakeInVivoPred-LiteratureFreeFrac.R - output is in vivo conc-resp models, converted to free concentration using Literature values (including Armitage model for in vitro).  Creates:
  o InVivoPred-LiteratureFreeFrac.pos.csv
  o InVivoPred-LiteratureFreeFrac.neg.csv
* PlotInVivo-LiteratureFreeFrac.R - plots of in vivo conc-resp using above CSV files.   Creates 3 PDF files:
  o FoldChange.invivo-LiteratureFreeFrac.pos.pdf
  o FoldChange.invivo-LiteratureFreeFrac.neg.pdf
  o FracChange.invivo-LiteratureFreeFrac.pos.log.pdf
* IVIVEchems.csv - mapping of chemical number to names
* MakeInVivoInVitrodata-ExperimentalFreeFrac.R - makes a data frame with both in vivo and in vitro model predictions - uses previously created CSV files, as well as files from the "All" subdirectory.  Saves data frames in 
  o InVivoInVitrodat-ExperimentalFreeFrac.Rdata 
* MakeInVivoInVitrodata-LiteratureFreeFrac.R - makes a data frame with both in vivo and in vitro model predictions (uses output from ...Experimental... to save time) Saves data frames in 
  o InVivoInVitrodat-LiteratureFreeFrac.Rdata
* In vitro Dose-Response Results_v4_ExperimentalFreeFrac.Rmd - plots comparing in vivo and in vitro.  Uses InVivoInVitrodat-ExperimentalFreeFrac.Rdata file and outputs several files:
  o IVIVE-ExperimentalFreeFrac.csv - table comparing in vivo and in vitro
  o IVIVE-ExperimentalFreeFracTable.csv - more formatted table
  o _ChemicalName_.InVitroInVivoConcResp-ExperimentalFreeFrac.pdf
  o InVitroInVivo-ExperimentalFreeFrac.PopMed.ECmax.pdf - Pop median results
  o InVitroInVivo-ExperimentalFreeFrac.StdDonor.ECmax.pdf - Std donor results
  o InVitroInVivo-ExperimentalFreeFrac.ECmax.pdf - both results side by side
  o InVitroInVivo-ExperimentalFreeFrac.PopMed.EC10.pdf - Pop median results
  o InVitroInVivo-ExperimentalFreeFrac.StdDonor.EC10.pdf - Std donor results
  o InVitroInVivo-ExperimentalFreeFrac.EC10.pdf - both results side by side
  o InVitroInVivo-ExperimentalFreeFrac.PopMed.EC05.pdf - Pop median results
  o InVitroInVivo-ExperimentalFreeFrac.StdDonor.EC05.pdf - Std donor results
  o InVitroInVivo-ExperimentalFreeFrac.EC05.pdf - both results side by side
  o InVitroInVivo-ExperimentalFreeFrac.PopMed.EC01.pdf - Pop median results
  o InVitroInVivo-ExperimentalFreeFrac.StdDonor.EC01.pdf - Std donor results
  o InVitroInVivo-ExperimentalFreeFrac.EC01.pdf - both results side by side
* In vitro Dose-Response Results_v4_LiteratureFreeFrac.Rmd - plots comparing in vivo and in vitro - same as previous but using literature values for free fraction
* PlotTraces.Rmd - plot Ca flux traces - this you will have to modify the directories to get to the right files
* MakeInVivoInVitrodata-ExperimentalFreeFrac.Samples.R - makes a data frame with both in vivo and in vitro model predictions but with all in vitro samples.  Saves data frames in 
  o InVivoInVitrodat-ExperimentalFreeFracSamples.Rdata
* ProbGEQ10ms.R - calculates the probability that the regulatory threshold of 10 msec is breached.  Uses InVivoInVitrodat-ExperimentalFreeFracSamples.Rdata and several previously CSV files.  Creates:
  o ProbDeltaQTc10.csv - data frame with probability as a function of concentration
  o ProbDeltaQTc10.pdf - PDF plot
