# CalciumFlux_peak_processing_functions 
# 15-Nov-2016
library(stringr)
## EXAMPLE:
##> seq1.file.list<- get.seq1.file.list()
##> n.files<-length(seq1.file.list)
##> for (n in 1:n.files) {
##> 	fname<-seq1.file.list[n]
##> 	dat<-get.oneplate.dat(fname)
##> 	peak.parms.vec<-get.onechem.peak.parms(dat,chemical.num=NA)
##> 	peak.parms.sum<-summarize.peak.parms(peak.parms.vec)
##> 	multiplot.check.peak.parms(dat,chemical.num=NA,lab="")
##> }

### GET FILE NAMES
get.seq1.file.list <- function(filepath="seq1") {
	seq1.file.list<-paste(filepath,"/",list.files(filepath,pattern="seq1"),sep="")
}

## GETTING PLATE DESIGN 
get.plate.design<- function(filepath=".",fnameprefix="PlateDesign",fext="dat",nplates=1,removeempty=TRUE) {
	plate.map<-list()
	plate.map[[1]]<-numeric(0)
	plate.map[[2]]<-character(0)
	plate.map[[3]]<-numeric(0)
	plate.map[[4]]<-numeric(0)
	plate.map[[5]]<-character(0)
	names(plate.map)[1]<-"Plate"
	names(plate.map)[2]<-"rowcol"
	for (j in 1:nplates) { # loop through plates 1 to nplates
		fname<-paste(filepath,"/",fnameprefix,j,".",fext,sep="")
		for (i in 1:3) { # There are 3 parameters in the plate design
			nskip <- (i-1)*18
			# name of parameter
			data.parm.tmp<-str_trim(read.delim(fname,skip=nskip,nrows=1,header=FALSE,sep="=",as.is=TRUE,encoding="UTF-8")[[1]])
			# data on that parameter
			dat.tmp<-read.table(fname,skip=nskip+1,header=TRUE,row.names=1,sep="\t",as.is=TRUE,encoding="UTF-8",nrows=16)
			if (i == 1) { # row and column in this plate
				rowcol.tmp<-c( outer( rownames(dat.tmp), 1:24,FUN=paste ,sep=""))
				# number of data points on plate
				n.tmp<-length(rowcol.tmp)
				plate.map[["Plate"]]<-c(plate.map[["Plate"]],rep(j,n.tmp))
				plate.map[["rowcol"]]<-c(plate.map[["rowcol"]],rowcol.tmp)
			}
			if (j == 1) names(plate.map)[i+2]<-data.parm.tmp
			plate.map[[data.parm.tmp]]<-c(plate.map[[data.parm.tmp]],as.vector(as.matrix(dat.tmp)))
		}
	}
	plate.map.frame<-as.data.frame(plate.map,stringsAsFactors=FALSE)
	# Remove empty wells
	if (removeempty) {
	  keep.indx<- plate.map.frame[,5]!=""
	  keep.indx[is.na(keep.indx)]<-TRUE
	  plate.map.frame<-plate.map.frame[keep.indx,]
	}
	plate.map.frame
}

### GETTING CHEMICALS
get.chem.map <- function(filepath=".",fname="chemicals.dat") {
	chem.map<-read.table(paste(filepath,"/",fname,sep=""),sep="\t",header=TRUE,row.names=1,comment.char="",as.is=TRUE,encoding="UTF-8")
	chem.map
}

peak.parms.diagnostic.plots <- function(well,peak.parms.vec,peak.parms.sum,
                                       lower=0.1,upper=0.8,lab="") {
  j <- match(well,peak.parms.sum$Well)
  yrange<-numeric(0)
  for (n in 1:dim(peak.parms.sum)[1]) {
    yrange <- range(yrange,peak.parms.vec[[j]]$y)
  }
  breaks <- (floor(yrange[1]/10):ceiling(yrange[2]/10))*10
  with(peak.parms.vec[[j]], {
    with(peak.parms.sum[j,], {
      if (amp != 0) lower.tmp <- max(lower,(baseline-min(y))/amp) else lower.tmp <- lower # defining beginning of peak
      ### Find baseline and amplitude - use density (first and last mode)
      ### Plots
      # whole thing
      plot(t,y,axes=FALSE,xlab="",ylab="RFU",ylim=c(0,yrange[2]),cex=0.5)
      lines(t,y)
      points(peak.times,peak.max,pch=15,col="red",cex=0.75)
      lines(peak.times,peak.max,col="red")
      box()
      axis(2)
      abline(h=baseline,lty=2,col="red")
      abline(h=height,lty=2,col="red")
      abline(h=baseline+lower.tmp*amp,col="red",lty=3)
      abline(h=baseline+upper*amp,col="red",lty=3)
      mtext(paste(" PkA =",signif(peak.amp.avg,4),
                  #                      " PkFrq =",signif(peak.freq.avg,3),
                  " BPM =",signif(BPM,3),
                  " PkSpc =",signif(peak.spacing.avg,3),sep=""),
            side=1,line=-1.9,cex=0.7)
      mtext(paste(" PkACV =",signif(peak.amp.cv,4),
                  #                      " PkFrqCV =",signif(peak.freq.cv,3),
                  " PkSpcCV =",signif(peak.spacing.cv,3),sep=""),
            side=1,line=-1.05,cex=0.7)
      axis(1); mtext("t (sec)",side=1,line=2,cex=0.7);
      # 3 peaks
      xrange <- min(t)+c(0,min(100,max(3*peak.spacing.avg,5,na.rm=TRUE)))
      plot(t,y,axes=FALSE,xlab="",ylab="RFU",ylim=c(0,yrange[2]),cex=0.5,xlim=xrange)
      lines(t,y)
      points(peak.times,peak.max,pch=15,col="red",cex=0.75)
      lines(peak.times,peak.max,col="red")
      box()
      axis(2)
      axis(1,lab=FALSE)
      axis(1); mtext("t (sec)",side=1,line=2,cex=0.7);
      abline(h=baseline,lty=2,col="red")
      abline(h=height,lty=2,col="red")
      abline(h=baseline+lower.tmp*amp,col="red",lty=3)
      if (notch.frac>0) abline(h=baseline+notch.frac*amp,col="blue",lty=1)
      abline(h=baseline+upper*amp,col="red",lty=3)
      points(c(peak.times[1]-rise.times[1],peak.times[1]+decay.times[1]),
             c(baseline,baseline),col="red",pch=c(2,6))
      points(peak.times[1],baseline,pch=3,col="red")
      lines(c(peak.times[1]-rise.times[1],peak.times[1]+decay.times[1]),
            c(baseline,baseline),col="red",lwd=2)
      mtext(paste("PkRsT =",signif(rise.times.avg,3),
                  "PkDecT =",signif(decay.times.avg,3)),
            side=1,line=-1.9,cex=0.7)
      mtext(paste("PkWd =",signif(peak.width.avg,3),
                  "PkDecRsRatio =",signif(decay.rise.ratio.avg,3)),
            side=1,line=-1.05,cex=0.7)
      title(lab,cex.main=1)
      # finding baseline and amplitude
      hist(y,breaks=breaks,main="",axes=FALSE,freq=FALSE,xlab="",ylab="")
      lines(densy,col="red")
      points(c(densy$x[densy$x==baseline],densy$x[densy$x==height]),
             c(densy$y[densy$x==baseline],densy$y[densy$x==height]),pch=15,col="red")
      box()
      axis(1); 
      mtext("RFU",side=1,line=2,cex=0.7);
      mtext("Density",side=2,line=1,cex=0.7)
      abline(v=baseline,col="red",lty=2)
      abline(v=baseline,col="red",lty=2)
      abline(v=baseline+lower.tmp*amp,col="red",lty=3)
      if (notch.frac>0) abline(v=baseline+notch.frac*amp,col="blue",lty=1)
      abline(v=baseline+upper*amp,col="red",lty=3)
      if (notchpeak > 0) {
        points(notchpeak,densy$y[densy$x==notchpeak][1],pch=17,col="red",cex=2)
      }
      mtext(paste("sd:",signif(sd(y,na.rm=TRUE),3)," baseline:",signif(baseline,3),
                  " height:",signif(baseline+amp,3),"\nnotch.frac:",signif(notch.frac,3)),
                  side=3,line=-2.2,cex=0.7)
    }
    ) # with
  }
  ) # with
}

### FUNCTION TO DO DIAGNOSTIC PLOTS OF PEAK PROCESSING ALGORITHMS
multiplot.check.peak.parms <- function(dat,chemical.num=NA,lab="",lower=0.1,upper=0.8,
                                       sdlim=1,getmedia=FALSE,getsolvent=TRUE,nheader=8,savedat=FALSE) {
  nmax<-12;
  chem.map<-get.chem.map()
  peak.parms.vec<-get.onechem.peak.parms(dat,chemical.num=chemical.num,lower=lower,upper=upper,
                                         sdlim=sdlim,getmedia=getmedia,getsolvent=getsolvent)
  peak.parms.sum<-summarize.peak.parms(peak.parms.vec)
  n <- length(peak.parms.sum[[1]])
  yrange<-numeric(0)
  for (j in 1:n) {
    yrange <- range(yrange,peak.parms.vec[[j]]$y)
  }
  breaks <- (floor(yrange[1]/10):ceiling(yrange[2]/10))*10
  if (n <= nmax) {
    npages<-1
    nperpage<-n
  } else {
    npages<-ceiling(n/nmax)
    nperpage<-nmax
  }
  jstart<-1
  par(mfrow=c(nperpage,3),mar=c(0,4,0,1),oma=c(6,1,6,1));
  for (pagenow in 1:npages) {
    jstart <- 1+(pagenow-1)*nmax
    jend <- min(n,jstart+nmax-1)
    for (j in jstart:jend) {
      with(peak.parms.vec[[j]], {
        with(peak.parms.sum[j,], {
          if (amp != 0) lower.tmp <- max(lower,(baseline-min(y))/amp) else lower.tmp <- lower # defining beginning of peak
          ### Find baseline and amplitude - use density (first and last mode)
          ### Plots
          # whole thing
          plot(t,y,axes=FALSE,ylab="RFU",ylim=c(0,yrange[2]),cex=0.5)
          lines(t,y)
          points(peak.times,peak.max,pch=15,col="red",cex=0.75)
          lines(peak.times,peak.max,col="red")
          box()
          axis(2)
          abline(h=baseline,lty=2,col="red")
          abline(h=height,lty=2,col="red")
          abline(h=baseline+lower.tmp*amp,col="red",lty=3)
          abline(h=baseline+upper*amp,col="red",lty=3)
          mtext(paste("Well ",Well,
                      " PkA =",signif(peak.amp.avg,4),
#                      " PkFrq =",signif(peak.freq.avg,3),
                        " BPM =",signif(BPM,3),
                      " PkSpc =",signif(peak.spacing.avg,3),sep=""),
                side=1,line=-1.9,cex=0.7)
          mtext(paste("Conc =",Chemical.Concentration,
                      " PkACV =",signif(peak.amp.cv,4),
#                      " PkFrqCV =",signif(peak.freq.cv,3),
                      " PkSpcCV =",signif(peak.spacing.cv,3),sep=""),
                side=1,line=-1.05,cex=0.7)
          if (j == jend) { axis(1); mtext("t (sec)",side=1,line=2);}
          # 3 peaks
          xrange <- c(0,min(100,max(3*peak.spacing.avg,5,na.rm=TRUE)))
          plot(t,y,axes=FALSE,ylab="RFU",ylim=c(0,yrange[2]),cex=0.5,xlim=xrange)
          lines(t,y)
          points(peak.times,peak.max,pch=15,col="red",cex=0.75)
          lines(peak.times,peak.max,col="red")
          box()
          axis(2)
          axis(1,lab=FALSE)
          if (j == jend) { axis(1); mtext("t (sec)",side=1,line=2);}
          abline(h=baseline,lty=2,col="red")
          abline(h=height,lty=2,col="red")
          abline(h=baseline+lower.tmp*amp,col="red",lty=3)
          if (notch.frac>0) abline(h=baseline+notch.frac*amp,col="blue",lty=1)
          abline(h=baseline+upper*amp,col="red",lty=3)
          points(c(peak.times[1]-rise.times[1],peak.times[1]+decay.times[1]),
                 c(baseline,baseline),col="red",pch=c(2,6))
          points(peak.times[1],baseline,pch=3,col="red")
          lines(c(peak.times[1]-rise.times[1],peak.times[1]+decay.times[1]),
                c(baseline,baseline),col="red",lwd=2)
          mtext(paste("PkRsT =",signif(rise.times.avg,3),
                 "PkDecT =",signif(decay.times.avg,3)),
                side=1,line=-1.9,cex=0.7)
          mtext(paste("PkWd =",signif(peak.width.avg,3),
                      "PkDecRsRatio =",signif(decay.rise.ratio.avg,3)),
                side=1,line=-1.05,cex=0.7)
          # finding baseline and amplitude
          hist(y,breaks=breaks,main="",axes=FALSE,freq=FALSE)
          lines(densy,col="red")
          points(c(densy$x[densy$x==baseline],densy$x[densy$x==height]),
                 c(densy$y[densy$x==baseline],densy$y[densy$x==height]),pch=15,col="red")
          box()
          if (j == jend) { axis(1); mtext("RFU",side=1,line=2);}
          abline(v=baseline,col="red",lty=2)
          abline(v=baseline,col="red",lty=2)
          abline(v=baseline+lower.tmp*amp,col="red",lty=3)
          if (notch.frac>0) abline(v=baseline+notch.frac*amp,col="blue",lty=1)
          abline(v=baseline+upper*amp,col="red",lty=3)
          if (notchpeak > 0) {
            points(notchpeak,densy$y[densy$x==notchpeak][1],pch=17,col="red",cex=2)
          }
          mtext(paste("sd:",signif(sd(y,na.rm=TRUE),3)," baseline:",signif(baseline,3),
                      " height:",signif(baseline+amp,3)," notch.frac:",signif(notch.frac,3)),side=3,line=-1.2,cex=0.7)
          if (j == jstart) {
            if (lab == "") {
              lab.tmp <- paste("Cell line",Cell.line,
                               #"Time point",
                               timepoint,
                               "/ Plate",Plate)
              if (!is.na(chemical.num)) {
                chem.name<-chem.map[paste(chemical.num),1] 
                lab.tmp <- paste(lab.tmp,"/",chem.name)
              } else {
                lab.tmp <- paste(lab.tmp,"/",Solvent)
              }
            }
            title(lab.tmp,outer=TRUE)
          }
        }
        ) # with
      }
      ) # with
    }
  }
  if (savedat) {
    outfilename<-paste("peak.parms.sum",peak.parms.sum$Cell.line[1],peak.parms.sum$Plate[1],
                       chemical.num,"L",lower,"U",upper,"dat",sep=".")
    write.table(peak.parms.sum,file=outfilename,row.names=FALSE,sep="\t")
  }
  peak.parms.sum
}

find.one.plate<-function(cellline,plate=1,t="90min") {
	seq1.file.list<- get.seq1.file.list()
	file.array<-t(array(unlist(strsplit(seq1.file.list,"_")),dim=c(5,length(seq1.file.list))))
	fileidstring<-paste(cellline,"_P",plate,"_",t,sep="")
	print(fileidstring)
	file.num<-grep(fileidstring,seq1.file.list)
	if (length(file.num) > 1) {
		print("ambiguous files")
		return("ambiguous files")
	} else if (length(file.num) == 0) {
		print("file not found")
		return("file not found")
	} else {
		fname <- seq1.file.list[file.num]
		print(fname)
		get.oneplate.dat(fname)
	}
}

plot.one.well<-function(cellline,plate=1,well="B2",t="90min",lab="",lower=0.1,upper=0.8,sdlim=1,
                        samplingrate=8,nheader=8,nsec=10,...) {
  par(mfrow=c(2,1),mar=c(3,4,.1,.1),oma=c(0,0,2,0));
  dat<-find.one.plate(cellline,plate=plate,t=t)
  dat<-trim.oneplate.dat(dat,samplingrate=samplingrate,nheader=nheader,nsec=nsec)
  well.indx<-grep(well,dat$Well)
  well.parms.vec<-get.peak.parms(dat[well.indx,],lower=lower,upper=upper,sdlim=sdlim,...)
  well.parms.sum<-summarize.peak.parms(well.parms.vec)
  with(well.parms.vec[[1]], {
    with(well.parms.sum[1,], {
      if (amp != 0) lower.tmp <- max(lower,(baseline-min(y))/amp) else lower.tmp <- lower # defining beginning of peak
      plot(t,y,axes=FALSE,ylab="RFU",xlab="",ylim=c(0,max(y)),cex=0.5)
      lines(t,y)
      points(peak.times,peak.max,pch=15,col="red",cex=0.75)
      lines(peak.times,peak.max,col="red")
      box()
      axis(2)
      abline(h=baseline,lty=2,col="red")
      abline(h=height,lty=2,col="red")
      abline(h=baseline+lower.tmp*amp,col="red",lty=3)
      abline(h=baseline+upper*amp,col="red",lty=3)
      if (notch.frac >0) abline(h=baseline+notch.frac*amp,col="blue")
      mtext(paste("PkAmp =",signif(amp,4),
                  "; pk.amp.m=",signif(peak.amp.avg,4),
                  "; pk.amp.sd=",signif(peak.amp.sd,4),
                  "; pk.amp.cv=",signif(peak.amp.cv,3),
                  sep=""),
            side=1,line=-2.85)#,cex=0.7)
      mtext(paste("BPM =",signif(BPM,3),"; pk.obs=",n.peaks,"; pk.exp=", round(n.peaks.expected,1),
                  "; pk.exp/pk.obs=",signif(n.peaks.expected/n.peaks,3),
                  sep=""),
            side=1,line=-1.9)#,cex=0.7)
      mtext(paste("PkSp.m =",signif(peak.spacing.avg,3),
                  "; PkSp.sd =",signif(peak.spacing.sd,3),
                  "; PkSp.cv =",signif(peak.spacing.cv,3),
                  sep=""),
            side=1,line=-1.05)#,cex=0.7)
      axis(1); 
      mtext("t (sec)",side=1,line=2)
      xrange <- c(0,min(100,max(3*peak.spacing.avg,5,na.rm=TRUE)))
      plot(t,y,axes=FALSE,ylab="RFU",xlab="",ylim=c(0,max(y)),cex=0.5,xlim=xrange)
      lines(t,y)
      points(peak.times,peak.max,pch=15,col="red",cex=0.75)
      lines(peak.times,peak.max,col="red")
      box()
      axis(2)
      axis(1,lab=FALSE)
      axis(1); mtext("t (sec)",side=1,line=2);
      abline(h=baseline,lty=2,col="red")
      abline(h=height,lty=2,col="red")
      abline(h=baseline+lower.tmp*amp,col="red",lty=3)
      abline(h=baseline+upper*amp,col="red",lty=3)
      if (notch.frac >0) abline(h=baseline+notch.frac*amp,col="blue")
      points(c(peak.times[1]-rise.times[1],peak.times[1]+decay.times[1]),
             c(baseline,baseline),col="red",pch=c(2,6))
      points(peak.times[1],baseline,pch=3,col="red")
      lines(c(peak.times[1]-rise.times[1],peak.times[1]+decay.times[1]),
            c(baseline,baseline),col="red",lwd=2)
      mtext(paste("PkSpc.m =",signif(peak.spacing.avg,3),
                  "; PkSpc.sd =",signif(peak.spacing.sd,3),
                  "; PkSpc.cv =",signif(peak.spacing.cv,3),
                  sep=""),
            side=1,line=-2.85)#,cex=0.7)
      mtext(paste("PkRsT.m =",signif(rise.times.avg,3),
                  "; PkDecT.m =",signif(decay.times.avg,3),
                  "; PkWd.m =",signif(peak.width.avg,3)
      ),
      side=1,line=-1.9)#,cex=0.7)
      mtext(paste("notch.frac:",signif(notch.frac,3),
                  "; (PkDecT/PkRsT).m =",signif(decay.rise.ratio.avg,3),
                  "; (PkDecT-PkRsT).m =",signif(decay.times.avg-rise.times.avg,3)
      ),
      side=1,line=-1.05)#,cex=0.7)
      if (lab == "") {
        chem.map<-get.chem.map()
        lab.tmp <- paste("Cell line",cellline,
                         #"Time point",
                         timepoint,
                         "/ Plate",plate,"/ Well ",Well, "/ Conc =",Chemical.Concentration)
        if (!is.na(Chemical.Number)) {
          chem.name<-chem.map[paste(Chemical.Number),1] 
          lab.tmp <- paste(lab.tmp,"/",chem.name)
        } else {
          lab.tmp <- paste(lab.tmp,"/",Solvent)
        }
      }
      title(lab.tmp,outer=TRUE,cex.main=1)
    })
  })
}

plot.one.well.subplot<-function(cellline,plate=1,well="B2",t="90min",xrange.sub=NA,yrange.sub=NA,
                                lab="",lower=0.1,upper=0.8,sdlim=1,
                                samplingrate=8,nheader=8,nsec=10,...) {
  dat<-find.one.plate(cellline,plate=plate,t=t)
  dat<-trim.oneplate.dat(dat,samplingrate=samplingrate,nheader=nheader,nsec=nsec)
  well.indx<-grep(well,dat$Well)
  well.parms.vec<-get.peak.parms(dat[well.indx,],lower=lower,upper=upper,sdlim=sdlim,...)
  well.parms.sum<-summarize.peak.parms(well.parms.vec)
  with(well.parms.vec[[1]], {
    with(well.parms.sum[1,], {
      if (is.na(xrange.sub[1])) xrange.sub <- c(0,min(100,max(1.5*peak.spacing.avg,5,na.rm=TRUE)))
      if (is.na(yrange.sub[1])) yrange.sub <- range(y[t<max(xrange.sub)])
      #      if (amp != 0) lower.tmp <- max(lower,(baseline-min(y))/amp) else lower.tmp <- lower # defining beginning of peak
      plot(t,y,axes=FALSE,ylab="RFU",xlab="",ylim=c(-490,yrange.sub[2]),
           xaxs="i",yaxs="i",col="red",lwd=2,type="l")
      #      lines(t,y)
      #      points(peak.times,peak.max,pch=15,col="red",cex=0.75)
      #      lines(peak.times,peak.max,col="red")
      box()
      axis(2)
      #      abline(h=baseline,lty=2,col="red")
      #      abline(h=height,lty=2,col="red")
      #      abline(h=baseline+lower.tmp*amp,col="red",lty=3)
      #      abline(h=baseline+upper*amp,col="red",lty=3)
      #      if (notch.frac >0) abline(h=baseline+notch.frac*amp,col="blue")
      #      mtext(paste("PkAmp =",signif(amp,4),
      #                  "; pk.amp.m=",signif(peak.amp.avg,4),
      #                  "; pk.amp.sd=",signif(peak.amp.sd,4),
      #                  "; pk.amp.cv=",signif(peak.amp.cv,3),
      #                  sep=""),
      #            side=1,line=-2.85)#,cex=0.7)
      #      mtext(paste("pk.obs=",n.peaks,"; pk.exp=", round(n.peaks.expected,1),
      #                  "; pk.exp/pk.obs=",signif(n.peaks.expected/n.peaks,3),
      #                  sep=""),
      #            side=1,line=-1.9)#,cex=0.7)
      #      mtext(paste("PkFrq.m =",signif(peak.freq.avg,3),
      #                  "; PkFrq.sd =",signif(peak.freq.sd,3),
      #                  "; PkFrq.cv =",signif(peak.freq.cv,3),
      #                  sep=""),
      #            side=1,line=-1.05)#,cex=0.7)
      axis(1); 
      mtext("t (sec)",side=1,line=2)
      subplot(plot(t,y,ylab="",xlab="",ylim=yrange.sub,xlim=xrange.sub,lwd=3,col="red",type="l")
              ,5,-250,size=c(1.5,1),vadj=0,hadj=0
              ,pars=list(cex=0.8,tcl=-0.2,lab=c(5,2,2),mgp=c(0,.2,0)))
      #      lines(t,y)
      #      points(peak.times,peak.max,pch=15,col="red",cex=0.75)
      #      lines(peak.times,peak.max,col="red")
      box()
      axis(2,lab=FALSE)
      axis(1,lab=FALSE)
      axis(1); mtext("t (sec)",side=1,line=2);
      #      abline(h=baseline,lty=2,col="red")
      #      abline(h=height,lty=2,col="red")
      #      abline(h=baseline+lower.tmp*amp,col="red",lty=3)
      #      abline(h=baseline+upper*amp,col="red",lty=3)
      #      if (notch.frac >0) abline(h=baseline+notch.frac*amp,col="blue")
      #      points(c(peak.times[1]-rise.times[1],peak.times[1]+decay.times[1]),
      #             c(baseline,baseline),col="red",pch=c(2,6))
      #      points(peak.times[1],baseline,pch=3,col="red")
      #      lines(c(peak.times[1]-rise.times[1],peak.times[1]+decay.times[1]),
      #            c(baseline,baseline),col="red",lwd=2)
      #      mtext(paste("PkSpc.m =",signif(peak.spacing.avg,3),
      #                  "; PkSpc.sd =",signif(peak.spacing.sd,3),
      #                  "; PkSpc.cv =",signif(peak.spacing.cv,3),
      #                  sep=""),
      #            side=1,line=-2.85)#,cex=0.7)
      #      mtext(paste("PkRsT.m =",signif(rise.times.avg,3),
      #                  "; PkDecT.m =",signif(decay.times.avg,3),
      #                  "; PkWd.m =",signif(peak.width.avg,3)
      #      ),
      #      side=1,line=-1.9)#,cex=0.7)
      #      mtext(paste("notch.frac:",signif(notch.frac,3),
      #                  "; (PkDecT/PkRsT).m =",signif(decay.rise.ratio.avg,3),
      #                  "; (PkDecT-PkRsT).m =",signif(decay.times.avg-rise.times.avg,3)
      #      ),
      #      side=1,line=-1.05)#,cex=0.7)
      if (lab == "") {
        chem.map<-get.chem.map()
        if (!is.na(Chemical.Number)) {
          chem.name<-chem.map[paste(Chemical.Number),1]
          lab.tmp <- paste(chem.name,"Conc =",Chemical.Concentration)
        } else {
          lab.tmp <- paste(Solvent)
        }
      }
      title(lab.tmp)
    })
  })
}

plot.one.well.subplot.compare<-function(peak.parms.sum.tmp,printparms=FALSE,
                                        xrange.sub=NA,yrange.sub=NA,
                                        lab="",lower=0.1,upper=0.8,sdlim=1,
                                        samplingrate=8,nheader=8,nsec=10,...) {
  cellline<-peak.parms.sum.tmp$Cell.line
  plate<-peak.parms.sum.tmp$Plate
  well<-peak.parms.sum.tmp$Well
  t<-peak.parms.sum.tmp$Timepoint
  
  dat<-find.one.plate(cellline,plate=plate,t=t)
  dat<-trim.oneplate.dat(dat,samplingrate=samplingrate,nheader=nheader,nsec=nsec)
  well.indx<-grep(well,dat$Well)
  well.parms.vec<-get.peak.parms(dat[well.indx,],lower=lower,upper=upper,sdlim=sdlim,...)
  well.parms.sum<-summarize.peak.parms(well.parms.vec)
  with(well.parms.vec[[1]], {
    with(well.parms.sum[1,], {
      if (is.na(xrange.sub[1])) xrange.sub <- c(max(0,peak.times[1]-1),min(100,peak.times[1]+max(1.5*peak.spacing.avg,3.5,na.rm=TRUE)))
      if (is.na(yrange.sub[1])) yrange.sub <- range(y)#[t<max(xrange.sub)])
      #      if (amp != 0) lower.tmp <- max(lower,(baseline-min(y))/amp) else lower.tmp <- lower # defining beginning of peak
      plot(t,y,axes=FALSE,ylab="RFU",xlab="",ylim=c(-490,yrange.sub[2]*1.05),
           xaxs="i",yaxs="i",col="red",lwd=2,type="l")
      #      lines(t,y)
      #      points(peak.times,peak.max,pch=15,col="red",cex=0.75)
      #      lines(peak.times,peak.max,col="red")
      box()
      axis(2)
      abline(h=baseline,lwd=2)
      abline(h=(baseline+peak.amp.avg),lty=2,lwd=2)
      abline(h=peak.parms.sum.tmp$AvrgPkB,col="grey",lwd=2)
      abline(h=(peak.parms.sum.tmp$AvrgPkB+peak.parms.sum.tmp$AvrgPkA),lwd=2,col="grey",lty=2)
      if (printparms) {
        text(c(55,70),550,c("MolDev","In-House"),adj=c(0,0))
        text(c(35,55,70),350,c("Amplitude",signif(peak.parms.sum.tmp$AvrgPkA,3),signif(peak.amp.avg,3)),adj=c(0,0))
        lines(c(55,60),c(325,325),lwd=2,col="grey",lty=2)
        lines(c(70,75),c(325,325),lwd=2,lty=2)
        text(c(35,55,70),150,c("Baseline",signif(peak.parms.sum.tmp$AvrgPkB,3),signif(baseline,3)),adj=c(0,0))
        lines(c(55,60),c(125,125),lwd=2,col="grey")
        lines(c(70,75),c(125,125),lwd=2)
        text(c(35,55,70),-50,c("BPM",signif(peak.parms.sum.tmp$PkFrBPM,3),signif(BPM,3)),adj=c(0,0))
        text(c(35,55,70),-250,c("Width",signif(peak.parms.sum.tmp$AvPW10A,3),signif(peak.width.avg,3)),adj=c(0,0))
      }
      #      abline(h=baseline+lower.tmp*amp,col="red",lty=3)
      #      abline(h=baseline+upper*amp,col="red",lty=3)
      #      if (notch.frac >0) abline(h=baseline+notch.frac*amp,col="blue")
      #      mtext(paste("PkAmp =",signif(amp,4),
      #                  "; pk.amp.m=",signif(peak.amp.avg,4),
      #                  "; pk.amp.sd=",signif(peak.amp.sd,4),
      #                  "; pk.amp.cv=",signif(peak.amp.cv,3),
      #                  sep=""),
      #            side=1,line=-2.85)#,cex=0.7)
      #      mtext(paste("pk.obs=",n.peaks,"; pk.exp=", round(n.peaks.expected,1),
      #                  "; pk.exp/pk.obs=",signif(n.peaks.expected/n.peaks,3),
      #                  sep=""),
      #            side=1,line=-1.9)#,cex=0.7)
      #      mtext(paste("PkFrq.m =",signif(peak.freq.avg,3),
      #                  "; PkFrq.sd =",signif(peak.freq.sd,3),
      #                  "; PkFrq.cv =",signif(peak.freq.cv,3),
      #                  sep=""),
      #            side=1,line=-1.05)#,cex=0.7)
      axis(1); 
      mtext("t (sec)",side=1,line=2)
      subplotfunc <- function (t,y,xrange.sub,yrange.sub,baseline,peak.amp.avg,peak.times,rise.times,decay.times,
                               peak.parms.sum.tmp) {
        plot(t,y,ylab="",xlab="",ylim=yrange.sub,xlim=xrange.sub,lwd=3,col="red",type="l")
        abline(h=baseline,lwd=2)
        abline(h=(baseline+peak.amp.avg),lty=2,lwd=2)
        abline(h=peak.parms.sum.tmp$AvrgPkB,col="grey",lwd=2)
        abline(h=(peak.parms.sum.tmp$AvrgPkB+peak.parms.sum.tmp$AvrgPkA),lwd=2,col="grey",lty=2)
        for (j in 1:length(peak.times)) {
          points(c(peak.times[j]-rise.times[j],peak.times[j]+decay.times[j]),
                 c(baseline,baseline),pch=c(2,6),lwd=2)
          points(peak.times[j],baseline+peak.amp.avg,pch=15,lwd=2)
          lines(c(peak.times[j]-rise.times[j],peak.times[j]+decay.times[j]),
                c(baseline,baseline),lwd=2)
        }
      }
      subplot(subplotfunc(t,y,xrange.sub,yrange.sub,baseline,peak.amp.avg,peak.times,rise.times,decay.times,
                          peak.parms.sum.tmp)
              ,5,-250,size=c(1.5,1),vadj=0,hadj=0
              ,pars=list(cex=0.8,tcl=-0.2,lab=c(5,2,2),mgp=c(0,.2,0)))
      #      lines(t,y)
      #      points(peak.times,peak.max,pch=15,col="red",cex=0.75)
      #      lines(peak.times,peak.max,col="red")
      #      box()
      #      axis(2,lab=FALSE)
      #      axis(1,lab=FALSE)
      #      axis(1); mtext("t (sec)",side=1,line=2);
      #      if (notch.frac >0) abline(h=baseline+notch.frac*amp,col="blue")
      #      mtext(paste("PkSpc.m =",signif(peak.spacing.avg,3),
      #                  "; PkSpc.sd =",signif(peak.spacing.sd,3),
      #                  "; PkSpc.cv =",signif(peak.spacing.cv,3),
      #                  sep=""),
      #            side=1,line=-2.85)#,cex=0.7)
      #      mtext(paste("PkRsT.m =",signif(rise.times.avg,3),
      #                  "; PkDecT.m =",signif(decay.times.avg,3),
      #                  "; PkWd.m =",signif(peak.width.avg,3)
      #      ),
      #      side=1,line=-1.9)#,cex=0.7)
      #      mtext(paste("notch.frac:",signif(notch.frac,3),
      #                  "; (PkDecT/PkRsT).m =",signif(decay.rise.ratio.avg,3),
      #                  "; (PkDecT-PkRsT).m =",signif(decay.times.avg-rise.times.avg,3)
      #      ),
      #      side=1,line=-1.05)#,cex=0.7)
      if (lab == "") {
        chem.map<-get.chem.map()
        if (!is.na(Chemical.Number)) {
          chem.name<-chem.map[paste(Chemical.Number),1]
          lab.tmp <- paste(chem.name,"Conc =",Chemical.Concentration)
        } else {
          lab.tmp <- paste(Solvent)
        }
      }
      title(lab.tmp)
    })
  })
}

## GETTING PLATE DATA 
get.plate.dat<- function(fname="plate1.dat",nparm=1,removeempty=TRUE) {
  plate.dat<-list()
  plate.dat[[1]]<-character(0)
  names(plate.dat)[1]<-"rowcol"
  for (i in 1:nparm) { # There are nparm parameters in the plate data
    nskip <- (i-1)*19
    # name of parameter
    data.parm.tmp<-str_trim(read.delim(fname,skip=(nskip+1),nrows=1,header=FALSE,sep="\t",as.is=TRUE,encoding="UTF-8")[2])
    # data on that parameter
    dat.tmp<-read.table(fname,skip=nskip+2,header=TRUE,row.names=1,sep="\t",as.is=TRUE,encoding="UTF-8",nrows=16)
    if (i == 1) { # row and column in this plate
      rowcol.tmp<-c( outer( rownames(dat.tmp), 1:24,FUN=paste ,sep=""))
      # number of data points on plate
      n.tmp<-length(rowcol.tmp)
      plate.dat[["rowcol"]]<-c(plate.dat[["rowcol"]],rowcol.tmp)
    }
    plate.dat[[data.parm.tmp]]<-c(plate.dat[[data.parm.tmp]],as.vector(as.matrix(dat.tmp)))
  }
  plate.dat.frame<-as.data.frame(plate.dat,stringsAsFactors=FALSE)
  # Remove empty wells -- all NA
  if (removeempty) {
    keep.indx<- apply(!is.na(as.matrix(as.data.frame(plate.dat,stringsAsFactors=FALSE)[,-1])),1,sum) != 0
    plate.dat.frame<-plate.dat.frame[keep.indx,]
  }
  plate.dat.frame
}