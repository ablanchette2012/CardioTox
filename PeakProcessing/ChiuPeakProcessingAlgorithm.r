# ChiuPeakProcessingAlgorithm 
library(stringr)
### GETTING DATA FRAME ONE PLATE
get.seq1.dat <- function(fname,printname=FALSE) {
  if (printname) print(fname)
  dat.tmp<-read.table(fname,as.is=TRUE,header=TRUE,sep="\t")
  dat.tmp<-dat.tmp[,-dim(dat.tmp)[2]] # remove last column (blank)
  dat.tmp<-dat.tmp[,-(1:4)] # remove first four columns
  dat.tmp
}

### FUNCTION TO FIND PEAK AMPLITUDE (INCLUDING IF THERE IS A "NOTCH")
## Peak baseline is defined as the first "mode" in the density distribution of y
## Peak height is defined as the last "mode" in the density distribution of y
## Peak amplitude is peak height - peak baseline
## The "notch" is defined as a "mode" that 
##   (a) occurs between "lower" and "upper" AND
##   (b) has a density at the mode that is greater than "notch.min" times the density of the peak "height"
## If the "baseline" is noisy, "lower" might be made higher than the input parameter 
##   (default = 0.1), and defined as twice the distance from the min(y) to the first mode
find.peak.amplitude <- function(y,bw=20, lower = 0.1, upper=0.8, notch.min = 1) {
  densy <- density(y,bw=bw)
  peakbaseline.indx<-match(-1,sign(diff(densy$y)))  # first time density distribution changes sign
  peakbaseline<-densy$x[peakbaseline.indx]
  peakheight.revindx<-match(-1,sign(diff(rev(densy$y)))) # last time density distribution changes sign
  peakheight<-rev(densy$x)[peakheight.revindx]
  peakamp<-peakheight-peakbaseline # height - baseline
  if (peakamp != 0) lower <- max(lower,(peakbaseline-min(y))/peakamp) # redefine defining beginning of peak
  if (peakheight <= peakbaseline) { # no peaks
    notchpeak <- 0 
  } else { # find if there is a "notch"
    # keep only y values between "lower" and "upper"
    tmp.indx<-densy$x > (peakbaseline+peakamp*lower) & densy$x < (peakbaseline+peakamp*upper);
    if (sum(tmp.indx)<3) { # too few samples between baseline and peak
      notchpeak <- 0
    } else {
      densy.y.tmp<-densy$y[tmp.indx]
      # remove declining density at beginning and end (these are parts of the baseline and peak)
      while (length(densy.y.tmp) >= 3 & densy.y.tmp[2] < densy.y.tmp[1]) densy.y.tmp <- densy.y.tmp[-1]
      while (length(densy.y.tmp) >= 3 & 
             densy.y.tmp[length(densy.y.tmp)-1] < densy.y.tmp[length(densy.y.tmp)]) densy.y.tmp <- densy.y.tmp[-length(densy.y.tmp)]
      if (length(densy.y.tmp) < 3) { # too few samples between baseline and peak
        notchpeak <- 0
      } else { # check the maximum density remaining (mode)
        if (max(densy.y.tmp) < notch.min*max(densy$y[densy$x > (peakbaseline+peakamp*upper)])) { # notch is too small
          notchpeak <- 0
        } else { # notch location (only records the first max) - point is whether it exists or not
          notchpeakdens <- max(densy.y.tmp) 
          notchpeak <- densy$x[densy$y == notchpeakdens][1] # value of y where notch is
        }
      }
    }
  }
  if (peakamp > 0 & notchpeak > 0) {
    notchpeakfrac <- (notchpeak - peakbaseline)/ peakamp # value of y as fraction of amplitude
  } else {
    notchpeakfrac <- 0; 
  }
  list(densy = densy,
       baseline = peakbaseline, 
       height=peakheight,
       amp=peakamp,
       notchpeak=notchpeak,
       notchpeakfrac=notchpeakfrac)
}

### FUNCTION TO FIND OTHER PEAK STATISTICS
##  A peak "rise" is defined as the data rising from below "lower" to above "upper"
##  A peak "decay" is defined as the data falling from above "upper" to below "lower"
##  "lower" might be made higher than the input parameter (default 0.1) if
##  the "baseline" is noisy
##  Assume there are no peaks if the sd(y) is less than "sdlim" or if peakamp = 0
##  Identifies for each peak:
##    peak maximum (must be above upper)
##    peak times (times where max occurs)
##    rise time (time from below "lower" to "max"
##    decay time (time from "max" to below "lower"
##    peak spaces (time between "max" values)
##    peak frequency in Hz (1/peak spaces)
##    number of total peak 
find.peak.stats <- function(t,y,peak.amp,lower=0.1,upper=0.8,sdlim=1) {
  peakbaseline<-peak.amp$baseline
  peakamp<-peak.amp$amp
  peakheight<-peak.amp$height
  densy<-peak.amp$densy
  if (peakamp != 0) lower <- max(lower,(peakbaseline-min(y))/peakamp) # defining beginning of peak
  if (peakamp != 0 & sd(y) > sdlim & (lower < upper)) {
    peaks.baseline.indx<-!(y >= peakbaseline+peakamp*lower & y <=peakbaseline+peakamp*upper)
    # remove values between lower and upper
    y.tmp<-y[peaks.baseline.indx]
    t.tmp<-t[peaks.baseline.indx]
    n.tmp<-length(y.tmp)
    # time between points
    dy.peaks.baseline <- diff(y.tmp)
    dt.peaks.baseline <- diff(t.tmp)
    # indices of all the rises and decays (travels further than upper-lower)
    dy.rise.indx <- (dy.peaks.baseline >= (upper-lower)*peakamp) & 
      (y.tmp[-1] > (peakbaseline+peakamp*upper)) &
      (y.tmp[-(length(y.tmp))]) < (peakbaseline+peakamp*lower)
    dy.decay.indx <- (dy.peaks.baseline <= -(upper-lower)*peakamp) &
      (y.tmp[-(length(y.tmp))] > (peakbaseline+peakamp*upper)) &
      (y.tmp[-1]) < (peakbaseline+peakamp*lower)
    # indices of all the "starts" and "ends" of peaks
    peak.starts.indx <- (1:n.tmp)[c(dy.rise.indx,FALSE)]
    peak.ends.indx <- (1:n.tmp)[c(FALSE,dy.decay.indx)]
    #		print(paste("raw",length(peak.starts.indx),length(peak.ends.indx)))
    #	  print(peak.starts.indx)
    #		print(peak.ends.indx)
    # make sure have a "whole peak" (start followed by end)
    if (length(peak.starts.indx) == 0 | length(peak.ends.indx) ==0) {
      n.peaks <- 0;
      peak.widths <- 0;
      peak.max <- peakbaseline;
      peak.times <- 0;
      rise.times <- 0;
      decay.times <- 0;
      peak.spaces <- 0;
      peak.freqs <- 0;
    } else {
      if (t.tmp[peak.ends.indx[1]] <= t.tmp[peak.starts.indx[1]]) peak.ends.indx <- peak.ends.indx[-1];
      if (max(t.tmp[peak.starts.indx]) >= max(t.tmp[peak.ends.indx])) peak.starts.indx <- peak.starts.indx[-length(peak.starts.indx)];
      # peak statistics
      n.peaks <- length(peak.starts.indx)
      peak.widths <- numeric(n.peaks)
      peak.max <- numeric(n.peaks)
      peak.max.indx <- numeric(n.peaks)
      peak.times <- numeric(n.peaks)
      rise.times <- numeric(n.peaks)
      decay.times <- numeric(n.peaks)
      for (k in 1:n.peaks) { # for each peak, get width, value at max, time of max, rise time, decay time
        peak.widths[k] <- t.tmp[peak.ends.indx[k]] - t.tmp[peak.starts.indx[k]]
        peak.max[k] <- max(y.tmp[peak.starts.indx[k]:peak.ends.indx[k]])
        peak.max.indx[k] <- ((peak.starts.indx[k]:peak.ends.indx[k])[y.tmp[peak.starts.indx[k]:peak.ends.indx[k]]==peak.max[k]])[1]
        peak.times[k] <- t.tmp[peak.max.indx[k]]
        rise.times[k] <- t.tmp[peak.max.indx[k]]-t.tmp[peak.starts.indx[k]]
        decay.times[k] <- t.tmp[peak.ends.indx[k]]-t.tmp[peak.max.indx[k]]
      }
      peak.spaces <- diff(peak.times)
      peak.freqs <- 1/peak.spaces
    }
  } else {
    n.peaks <- 0;
    peak.widths <- 0;
    peak.max <- peakbaseline;
    peak.times <- 0;
    rise.times <- 0;
    decay.times <- 0;
    peak.spaces <- 0;
    peak.freqs <- 0;
  }
  list(n.peaks = n.peaks,
       peak.widths = peak.widths,
       peak.max = peak.max,
       peak.times = peak.times,
       rise.times = rise.times,
       decay.times = decay.times,
       peak.spaces = peak.spaces,
       peak.freqs = peak.freqs
  )
}

get.peak.parms <- function(dat.in,bw=20,lower=0.1,upper=0.8,sdlim=1,nheader=1,
                           trimsec=0) {
  n <- dim(dat.in)[1]
  t<-as.numeric(sub("X","",colnames(dat.in)[-(1:nheader)]))
  t.keep<-t>=trimsec
  t<-t[t.keep]
  dat<-dat.in[,c(rep(TRUE,nheader),t.keep)]
  peak.parms.vec <- vector(mode="list",length=n)
  for (j in 1:n) {
    y <-as.numeric(dat[j,-(1:nheader)])
    ### Find baseline and amplitude - use density (first and last mode)
    peak.amp <- find.peak.amplitude(y,bw=bw,lower=lower,upper=upper)
    peak.stats <- find.peak.stats(t,y,peak.amp,lower=lower,upper=upper,sdlim=sdlim)
    peak.amp$n.peaks.expected <- max(t)*mean(peak.stats$peak.freqs)
    peak.amp$t<-t
    peak.amp$y<-y
    peak.parms.vec[[j]] <- c(dat[j,1:nheader],peak.amp,peak.stats)
    names(peak.parms.vec[[j]])[1]<-"Well"
  }
  peak.parms.vec
}

summarize.peak.parms <- function(peak.parms.vec) {
  peak.parms.vec.summary.names <- c(
    "Well",
    "baseline",
    "height",
    "amplitude",
    "notch.amp",
    "notch.frac",
    "n.peaks",
    "n.peaks.expected",
    "BPM",
    "peak.max.avg",
    "peak.max.sd",
    "peak.max.cv",
    "peak.amp.avg",
    "peak.amp.sd",
    "peak.amp.cv",
    "peak.spacing.avg",
    "peak.spacing.sd",
    "peak.spacing.cv",
    "peak.freq.avg",
    "peak.freq.sd",
    "peak.freq.cv",
    "peak.width.avg",
    "peak.width.sd",
    "peak.width.cv",
    "rise.times.avg",
    "rise.times.sd",
    "rise.times.cv",
    "decay.times.avg",
    "decay.times.sd",
    "decay.times.cv",
    "decay.rise.ratio.avg",
    "decay.rise.ratio.sd",
    "decay.rise.ratio.cv"
  )
  peak.parms.vec.summary <- matrix(nrow=length(peak.parms.vec),ncol=length(peak.parms.vec.summary.names));
  colnames(peak.parms.vec.summary) <- peak.parms.vec.summary.names
  peak.parms.vec.summary <- as.data.frame(peak.parms.vec.summary)
  for (j in 1:length(peak.parms.vec)) {
    peak.parms.vec.summary[j,"Well"]<-peak.parms.vec[[j]]$Well
    peak.parms.vec.summary[j,"baseline"]<-peak.parms.vec[[j]]$baseline
    peak.parms.vec.summary[j,"height"]<-peak.parms.vec[[j]]$height
    peak.parms.vec.summary[j,"amplitude"]<-peak.parms.vec[[j]]$amp
    peak.parms.vec.summary[j,"notch.amp"]<-peak.parms.vec[[j]]$notchpeak
    peak.parms.vec.summary[j,"notch.frac"]<-peak.parms.vec[[j]]$notchpeakfrac
    peak.parms.vec.summary[j,"n.peaks"]<-peak.parms.vec[[j]]$n.peaks
    peak.parms.vec.summary[j,"n.peaks.expected"]<-peak.parms.vec[[j]]$n.peaks.expected
    peak.parms.vec.summary[j,"BPM"]<-peak.parms.vec[[j]]$n.peaks/(max(peak.parms.vec[[j]]$t)/60)
    peak.parms.vec.summary[j,"peak.max.avg"]<-mean(peak.parms.vec[[j]]$peak.max)
    peak.parms.vec.summary[j,"peak.max.sd"]<-sd(peak.parms.vec[[j]]$peak.max)
    peak.parms.vec.summary[j,"peak.max.cv"]<-sd(peak.parms.vec[[j]]$peak.max)/mean(peak.parms.vec[[j]]$peak.max)
    peak.parms.vec.summary[j,"peak.amp.avg"]<-mean(peak.parms.vec[[j]]$peak.max-peak.parms.vec[[j]]$baseline)
    peak.parms.vec.summary[j,"peak.amp.sd"]<-sd(peak.parms.vec[[j]]$peak.max-peak.parms.vec[[j]]$baseline)
    peak.parms.vec.summary[j,"peak.amp.cv"]<-sd(peak.parms.vec[[j]]$peak.max-peak.parms.vec[[j]]$baseline)/
      mean(peak.parms.vec[[j]]$peak.max-peak.parms.vec[[j]]$baseline)
    peak.parms.vec.summary[j,"peak.spacing.avg"]<-mean(peak.parms.vec[[j]]$peak.spaces)
    peak.parms.vec.summary[j,"peak.spacing.sd"]<-sd(peak.parms.vec[[j]]$peak.spaces)
    peak.parms.vec.summary[j,"peak.spacing.cv"]<-sd(peak.parms.vec[[j]]$peak.spaces)/mean(peak.parms.vec[[j]]$peak.spaces)
    peak.parms.vec.summary[j,"peak.freq.avg"]<-mean(peak.parms.vec[[j]]$peak.freqs)
    peak.parms.vec.summary[j,"peak.freq.sd"]<-sd(peak.parms.vec[[j]]$peak.freqs)
    peak.parms.vec.summary[j,"peak.freq.cv"]<-sd(peak.parms.vec[[j]]$peak.freqs)/mean(peak.parms.vec[[j]]$peak.freqs)
    peak.parms.vec.summary[j,"peak.width.avg"]<-mean(peak.parms.vec[[j]]$peak.widths)
    peak.parms.vec.summary[j,"peak.width.sd"]<-sd(peak.parms.vec[[j]]$peak.widths)
    peak.parms.vec.summary[j,"peak.width.cv"]<-sd(peak.parms.vec[[j]]$peak.widths)/mean(peak.parms.vec[[j]]$peak.widths)
    peak.parms.vec.summary[j,"rise.times.avg"]<-mean(peak.parms.vec[[j]]$rise.times)
    peak.parms.vec.summary[j,"rise.times.sd"]<-sd(peak.parms.vec[[j]]$rise.times)
    peak.parms.vec.summary[j,"rise.times.cv"]<-sd(peak.parms.vec[[j]]$rise.times)/mean(peak.parms.vec[[j]]$rise.times)
    peak.parms.vec.summary[j,"decay.times.avg"]<-mean(peak.parms.vec[[j]]$decay.times)
    peak.parms.vec.summary[j,"decay.times.sd"]<-sd(peak.parms.vec[[j]]$decay.times)
    peak.parms.vec.summary[j,"decay.times.cv"]<-sd(peak.parms.vec[[j]]$decay.times)/mean(peak.parms.vec[[j]]$decay.times)
    peak.parms.vec.summary[j,"decay.rise.ratio.avg"]<-mean(peak.parms.vec[[j]]$decay.times/peak.parms.vec[[j]]$rise.times)
    peak.parms.vec.summary[j,"decay.rise.ratio.sd"]<-sd(peak.parms.vec[[j]]$decay.times/peak.parms.vec[[j]]$rise.times)
    peak.parms.vec.summary[j,"decay.rise.ratio.cv"]<-sd(peak.parms.vec[[j]]$decay.times/peak.parms.vec[[j]]$rise.times)/
      mean(peak.parms.vec[[j]]$decay.times/peak.parms.vec[[j]]$rise.times)
  }
  peak.parms.vec.summary
}
