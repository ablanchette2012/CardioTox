get.conc.resp<-function(dat,yname,ydirection,controlfactor=400) {
  dat$y<-dat[,yname]
  controldat<-subset(dat,is.na(Chemical.Number))
  bval<-median(controldat$y)
  if (ydirection==1) { # positive direction - ln(fold-change)
    dat$resp <- log(dat$y/bval)
  } else { # negative direction - percent change
    pval<-0 # maximal effect is at y=0
    dat$resp <- 100*(dat$y-bval)/(pval-bval)
  }
  dat$logc <- log10(dat$Chemical.Concentration)
  minconc<-min(dat$logc[dat$Chemical.Concentration!=0])
  dat$logc[dat$Chemical.Concentration==0] <- minconc-log10(controlfactor)
  dat
}

dotcplfit<-function(dat.y,ynamenow,ydirectionnow) {
  if (ydirectionnow==1) {
    ydirtext<-"up"
  } else {
    ydirtext<-"down"
  }
  dat.controls<-subset(dat.y,Chemical.Concentration==0)
  bmad<-median(abs(dat.controls$resp-median(dat.controls$resp)))
  dat.treated<-subset(dat.y,Chemical.Number==n)
  logc<-c(dat.controls$logc,dat.treated$logc)
  xlim<-range(c(logc,2))
  resp<-c(dat.controls$resp,dat.treated$resp)
  if (ydirectionnow==1) {
    ylab<-"ln(fold-increase))"
    logc<-logc[!is.infinite(resp) & !is.na(resp)] # remove infinite and NA log-fold-changes
    resp<-resp[!is.infinite(resp) & !is.na(resp)]
    ylim<-range(c(0,3*bmad,-3*bmad,resp),na.rm=TRUE)
    bmr10<-log(1.1)
  } else {
    ylab<-"% Decrease"
    ylim<-range(c(resp,3*bmad,-3*bmad,0,100),na.rm=TRUE)
    bmr10<-10
  }
  plot(resp~logc,main="",xlab="log10(Concentration)",ylab=ylab,xlim=xlim,ylim=ylim)
  mtext(paste("Chemical",n,ynamenow,ydirtext),line=2.,cex=0.7)
  parms<-tcplFit(logc=logc, resp=resp, bmad=bmad)
  if (!is.na(parms$gnls_aic)) tcplAddModel(pars = parms, modl = "gnls",lwd=3,lty=2,col="maroon")
  if (!is.na(parms$hill_aic)) tcplAddModel(pars = parms, modl = "hill",lwd=3,lty=2,col="purple")
  if (!is.na(parms$cnst_aic)) tcplAddModel(pars = parms, modl = "cnst",lwd=3,lty=2,col="orange")
  rect(-10,-3*bmad,10,3*bmad,density=10,col="grey")
  bestmodel<-c("cnst","hill","gnls")[which.min(c(parms$cnst_aic,parms$hill_aic,parms$gnls_aic))]
  if (length(bestmodel)==0) bestmodel<-"cnst"
  if (bestmodel=="cnst") {
    ac50<-NA; ec10<-NA; ecbmad<-NA; tp<-0; gw<-0
  } else if (bestmodel=="hill") {
    ac50<-10^parms$hill_ga;
    ec10<-ac50*(parms$hill_tp/bmr10-1)^(-1/parms$hill_gw)
    ecbmad<-ac50*(parms$hill_tp/bmad-1)^(-1/parms$hill_gw)
    tp<-parms$hill_tp
    gw<-parms$hill_gw
  } else {
    ac50<-10^parms$gnls_ga;
    ec10<-ac50*(parms$gnls_tp/bmr10-1)^(-1/parms$gnls_gw)
    ecbmad<-ac50*(parms$gnls_tp/bmad-1)^(-1/parms$gnls_gw)
    tp<-parms$gnls_tp
    gw<-parms$gnls_gw
  }
  abline(v=log10(ac50),col="red")
  abline(h=tp/2,col="red",lty=2)
  abline(v=log10(ec10),col="blue")
  abline(h=bmr10,col="blue",lty=2)
  abline(v=log10(ecbmad),col="grey")
  abline(h=bmad,col="grey",lty=2)
  points(dat.treated$resp~dat.treated$logc,pch=16)
  results.tmp<-data.frame(Chemical.Number=n,yname=ynamenow,direction=ydirtext,
                          bmad,bestmodel,hillcoeff=gw,ac50,ec10,ecbmad)
  results.label<-paste("BMAD:",signif(bmad,3)," Best Model:",bestmodel,
                       " Hill coeff:",signif(gw,3),"\n",
                       " AC50:",signif(ac50,3),
                       " EC10:",signif(ec10,3),
                       " ECBMAD:",signif(ecbmad,3))
  mtext(results.label,line=0.05,cex=0.5)
  results.tmp
}
