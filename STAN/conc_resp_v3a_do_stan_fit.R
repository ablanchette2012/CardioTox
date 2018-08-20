# Redo scripts for this conc_response
# nameprefix <- "dose_response_test"
# modelname<-"conc_resp_v3a"
setwd("C:\\Users\\ablanchette\\Documents\\dose_response_template\\dose_response_template")
library(parallel)
library(rstan)
library(MASS)
library(Hmisc) ##getting error-->no package named Hmisc
library(reshape2)
library(lattice)
######
##retrieve chemical list, but why retreieve both chemicals.dat and chemicals_1434? 
get.chem.map <- function(filepath=".") {
  chem.map<-read.table(paste(filepath,"/chemicals.dat",sep=""),sep="\t",header=TRUE,row.names=1,comment.char="",as.is=TRUE,encoding="UTF-8")
  chem.map<-cbind(chem.map,read.table(paste(filepath,"/chemicals_1434.dat",sep=""),sep="\t",header=TRUE,row.names=1,comment.char="",as.is=TRUE,encoding="UTF-8"))
  names(chem.map)[2]<-paste(names(chem.map)[2],"1434",sep=".")
  chem.map
}
##create/organize data for STAN to utilize in modeling
make_stan_dat <- function(dat,chemnum, parmcol, parmname="", parmnamenorm="",
                          quants = c(0.01,0.025,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.975,0.99)) {
  dat$cell.f<-factor(dat$Cell.line)
  chemmap<-get.chem.map()
  dat$y <- dat[,parmcol];
  chemname<-chemmap[as.character(chemnum),1]
  chemnum.1434<-as.numeric(rownames(chemmap)[match(chemname,chemmap[,2])])
  dat.tmp<-rbind(subset(dat,Cell.line != 1434 & (Chemical.Number == chemnum | is.na(Chemical.Number))),
                 subset(dat,Cell.line == 1434 & (Chemical.Number == chemnum.1434 | is.na(Chemical.Number))))
  dat.tmp<-subset(dat.tmp,!is.na(y))
  ###
  ### For all other chemicals, need to do additional subset for same plate-controls
  ###
  fileprefix<-paste("standat_v3a_do_stan_fit_",parmnamenorm,"_",chemnum,sep="")
  cat("##############",fileprefix,"\n")
  stan_dat <- list(
    scale_factor = median(dat$y,na.rm=TRUE), # dat not dat.tmp
    Ni = length(unique(dat$cell.f)), # dat not dat.tmp
    Nj = dim(dat.tmp)[1],
    x = dat.tmp$Chemical.Concentration,
    ys = dat.tmp$y/median(dat$y,na.rm=TRUE),
    cell = as.numeric(dat.tmp$cell.f),
    quants = quants,
    Nquants = length(quants)
  )
  with(stan_dat, {
    stan_rdump(names(stan_dat),file=paste(fileprefix,"_dat.R",sep=""))
  })
  list(stan_dat=stan_dat,fileprefix=fileprefix,chemname=chemname,parmcol=parmcol,
       parmname=parmname,parmnamenorm=parmnamenorm)
}
#specify the model to use based on direction of data
do_stan_fit <- function(stan_dat, fileprefix="",direction = 1,
                        verbose=FALSE,seed=271828,iter=4000,
                        ...) {
  if (direction == 1) {
    stanmodel <- 'conc_resp_v3a_up.stan' 
  } else {
    stanmodel <- 'conc_resp_v3a_dn.stan' 
  }
  time.start<-proc.time()
  stan_fit <- stan(file=stanmodel,data=stan_dat,verbose=verbose,seed=seed,iter=iter,
                   sample_file=paste(fileprefix,"_samples",sep=""),...);
  time.end<-proc.time()
  print(time.end-time.start)
  stan_samples<-extract(stan_fit)
  list(stan_dat=stan_dat,stan_fit=stan_fit,stan_samples=stan_samples,direction=direction)
}

get_stan_fit_csv <- function(stan_dat,fileprefix,direction) {
  stan_fit <- read_stan_csv(dir(pattern=paste(fileprefix,"_samples_","[[:digit:]]",sep=""), ##should be creating files in working directory?
                                full.names = TRUE))
#  stan_fit[[1]] <- NULL # get rid of "energy"
  stan_samples<-extract(stan_fit) ##doesn't work if stan_fit doesn't work, don't understand the extraxt command
#  stan_samples[[1]]<-NULL # get rid of "energy"
  list(stan_dat=stan_dat,stan_fit=stan_fit,stan_samples=stan_samples,direction=direction)
}

plot_stan_fit_rhat <- function(stan_fit,stan_samples,fileprefix,dopdf=FALSE) {    ##not plotting, outputting pdf file-->is it looking for outputted csv files from previous blocks?
  if (dopdf) pdf(paste(fileprefix,"rhat","pdf",sep="."),height=10.5,width=8)
  print(plot(stan_fit,plotfun="stan_rhat",par=names(stan_samples)))
  if (dopdf) dev.off()
}

plot_stan_fit_parms <- function(stan_fit,fileprefix, dopdf=FALSE) { #not working because not outputting stanfit csv?
  if (dopdf) pdf(paste(fileprefix,"parms","pdf",sep="."),height=10.5,width=8)
  print(stan_plot(stan_fit,
                  par=c("m_y0","m_x0","m_Emax","m_n","sd_y0","sd_x0","sd_Emax","sd_n","sigma_y"))+ 
          ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="z_y0")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="z_x0")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="z_Emax")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="z_n")+ ggtitle(fileprefix))
  if (dopdf) dev.off()
}

plot_stan_fit_pop <- function(stan_fit,fileprefix, dopdf=FALSE) { ##same
  if (dopdf) pdf(paste(fileprefix,"pop","pdf",sep="."),height=10.5,width=8)
  print(stan_plot(stan_fit,par="y0")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="x0")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="Emax")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="n")+ ggtitle(fileprefix))
  print(stan_plot(stan_fit,par="ec10")+ ggtitle(fileprefix))
  if (dopdf) dev.off()
}

plot_stan_fit_iter <- function(stan_samples,stan_dat,direction,chemname,parmnamenorm, ##not working if stan_samples isn't
                               fileprefix="",dopdf=FALSE,niter=50) {
  if (dopdf) pdf(paste(fileprefix,"iter","pdf",sep="."),height=10.5,width=8)
  par(mfrow=c(9,3),mar=c(0,2,0.5,1),oma=c(4,3,3,1))
  x50 <- stan_samples$x0^stan_samples$n*stan_samples$Emax
  xmin=max(1e-6,min(0.001,10^floor(log10(quantile(x50,0.1,na.rm=TRUE))),na.rm=TRUE))
  #yrange=range(stan_dat$scale_factor*stan_dat$ys)
  for (i in 1:stan_dat$Ni) {
    xtmp <- stan_dat$x[stan_dat$cell == i]+xmin
    ytmp <- stan_dat$scale_factor*stan_dat$ys[stan_dat$cell == i]
    yrange=range(ytmp,na.rm=TRUE)
    xx<- 10^(((10*(ceiling(log10(xmin)))):30)/10)
    lastiter<-dim(stan_samples$y0)[1]
    for (iter in 0:(niter-1)) {
      y0<-stan_samples$y0[lastiter-iter,i]
      x0<-stan_samples$x0[lastiter-iter,i]
      Emax<-stan_samples$Emax[lastiter-iter,i]
      n<-stan_samples$n[lastiter-iter,i]
      if (iter == 0) {
        plot(xx,y0*stan_dat$scale_factor*(1 + direction * (xx / x0)^n / (1 + (xx / x0)^n/Emax)),
             ylim=yrange, #c(yrange[1]-0.25*(yrange[2]-yrange[1]),yrange[2]),
             type="l",log="x",col=i,xlab="",ylab="",axes=FALSE,lwd=0.5);
        axis(2,line=0)
        axis(1,at=10^(log10(xmin):2),lab=c("C",log10(xmin*10):2),line=0,outer=TRUE)
        box()
      } else {
        lines(xx,y0*stan_dat$scale_factor*(1 + direction * (xx / x0)^n / (1 + (xx / x0)^n/Emax)),col=i,lwd=0.5)
      }
    }
    xec10<-10^seq(max(log10(xmin),quantile(log10(stan_samples$ec10[,i]),prob=c(0.05))),
                  quantile(log10(stan_samples$ec10[,i]),prob=c(0.95)),
                  length.out=100)
    yec10<-stan_dat$scale_factor*rep(median(stan_samples$y0[,i])*(1+0.1*direction),100)
    points(xec10,yec10,pch=15)
    points(xec10,yec10,pch=15,col="grey",cex=0.5)
    points(xtmp,ytmp,pch=21,bg="white")
    text(xmin,yrange[2],lab=i,adj=c(0,1))
  }
  mtext(paste(chemname,parmnamenorm),line=1,outer=TRUE)
  mtext("Log10 Concentration",side=1,line=2.2,outer=TRUE)
  mtext(parmnamenorm,side=2,line=1,outer=TRUE)
  if (dopdf) dev.off()
}


###### Set up Rstan
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

###### Clean up data
alldat<-read.table("Nonblinded.peak.parms.sum.dat",header=TRUE,as.is=TRUE,sep="\t")
# Address quiescent cells - set most parameters to NA
quiescent <- alldat$amplitude == 0 | alldat$n.peaks == 0
alldat[quiescent,c(11:15,17:18,20:(dim(alldat)[2]))] <- NA

# excluding 3 cell lines
dat.clean<-subset(alldat,(Cell.line != 1384 & Cell.line != 1236 & Cell.line !=1066))

write.table(dat.clean,"Nonblinded.peak.parms.clean.dat",sep="\t",row.names=FALSE)

###### Read in clean data
dat<-read.table("Nonblinded.peak.parms.clean.dat",header=TRUE,as.is=TRUE,sep="\t")

parmcolvec<-c("peak.amp.avg","BPM","BPM","peak.amp.cv","peak.spacing.cv","decay.rise.ratio.avg")
parmnamevec<-c("Peak Amplitude (RFU)","BPM","BPM","Peak Amplitude CV","Peak Spacing CV","Peak Decay-Rise Ratio")
parmnamenormvec<-c("Peak_Amplitude_dn","BPM_up","BPM_dn",
                  "Peak_Amplitude_CV_up","Peak_Spacing_CV_up",
                  "Peak_Decay_Rise_Ratio_up")
# set direction (pos or negative) for each parameter
dirvec<-c(-1,1,-1,1,1,1)

# Iterate through each of the 6 peak parameters
for (k in 1:6) {
  parmcol<-parmcolvec[k]
  parmname<-parmnamevec[k]
  parmnamenorm<-parmnamenormvec[k]
  direction<-dirvec[k]
  # Iterate through 3 chemicals
  for (chemnum in c(5,6,15,37,2,34,35,28,12)) {
    stan_dat_results<-make_stan_dat(dat=dat,chemnum=chemnum,
                                    parmcol=parmcol,parmname=parmname,parmnamenorm=parmnamenorm)
    stan_fit_results<-do_stan_fit(stan_dat_results$stan_dat,
                                  stan_dat_results$fileprefix,direction=direction)
    # If Rhat is > 1.2 for some parameters, may need to increase number of iterations in "do_stan_fit"
    plot_stan_fit_rhat(stan_fit_results$stan_fit,stan_fit_results$stan_samples,
                       stan_dat_results$fileprefix,dopdf=TRUE)
    stan_fit_results<-get_stan_fit_csv(stan_dat_results$stan_dat,  ##should get_stan_fit_csv be a .csv? Doesn't seem to fix the error
                                       stan_dat_results$fileprefix,direction=direction)
    plot_stan_fit_parms(stan_fit_results$stan_fit,
                        stan_dat_results$fileprefix,dopdf=TRUE)
    plot_stan_fit_pop(stan_fit_results$stan_fit,
                      stan_dat_results$fileprefix,dopdf=TRUE)
    plot_stan_fit_iter(stan_fit_results$stan_samples,stan_dat_results$stan_dat,
                       stan_fit_results$direction,stan_dat_results$chemname,
                       stan_dat_results$parmnamenorm,stan_dat_results$fileprefix,dopdf=TRUE)
  }
}

















