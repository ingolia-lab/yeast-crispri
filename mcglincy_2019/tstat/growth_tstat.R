options(stringsAsFactors=FALSE)
library("RColorBrewer")

datadir <- Sys.getenv("DATADIR")
figuredir <- Sys.getenv("FIGUREDIR")

pal <- brewer.pal(8, "Dark2")

prefix <- "^T"

## timevol [ml / hour] is volume delivered per time of pumping
## ttlvol [ml] is culture volume
## ttlvol / timevol [hr] is time to fill the chamber
## timevol / ttlvol [hr^-1] is dilution rate at constant pumping
leftParam <- list(timevol = 0.255 * 60 * 60, ttlvol = 202)
rightParam <- list(timevol = 0.244 * 60 * 60, ttlvol = 206)

## add transparency to color
setAlpha <- function(color, alpha) {
    x <- col2rgb(color)
    rgb(x[1,], x[2,], x[3,], alpha, maxColorValue=255)
}

samples <- c(166000, 259000, 323000, 385000)
induceColor <- setAlpha(pal[[4]], 128)
sampleColor <- setAlpha(pal[[2]], 128)
sampleColors <- c(induceColor, rep(sampleColor, 3))

## Extract and read turbidostat data
getTlog <- function(screenfile) {
  ## Grep out one turbidostat header followed by all turbidostat data lines
  ## Read turbidostat log into a data.frame
  extractTlog <- function(screenfile, tlogfile) {
    system(sprintf("grep \'%s\' \'%s\' | grep time.s | head -n1 > %s",
                   prefix, screenfile, tlogfile))
    system(sprintf("grep \'%s\' \'%s\' | grep -v time.s >> %s",
                   prefix, screenfile, tlogfile))
    tlog <- read.delim(tlogfile)
    tlog$pumptime.s <- tlog$pumptime.s - min(tlog$pumptime.s)
    tlog
  }
  
  tlogfile <- tempfile(pattern="tlog", tmpdir=dirname(screenfile), fileext=".txt")
  tlog <- tryCatch(extractTlog(screenfile, tlogfile), finally=file.remove(tlogfile))
  tlog
}

## Plot cell density (nephelometry) data from a tlog to the current device
plotNeph <- function(tlog, param) {
  neph <- aggregate(tlog[,c("time.s", "neph", "target", "gain")], 
                    by=list(time.min=round(tlog$time.s / 60)),
                    FUN=median)
  targ <- max(neph$target/neph$gain)
  highest <- max(neph$neph/neph$gain, targ) / targ
  plot(neph$time.s / 3600, (neph$neph / neph$gain) / targ,
       pch=20, cex=0.25, col="black",
       ylim=c(0,1.2*highest), yaxp=c(0,1,2),
       xlim=c(0,108), xaxp=c(0,108,9),
       xlab="Time (hrs)", ylab="Cell Dens (A.U.)")
  abline(v=samples/3600, col=sampleColors, lwd=3)
  abline(h=1, col=setAlpha(pal[[1]], 128), lwd=3)
  text(x=0, y=0.95, labels="Target", adj=c(0,1), col=pal[[1]])
}    

## Plot pumping = dilution = growth data from a tlog to the current device
## Calculates a fit line for the most recent hour (cf recentfit) to estimate rate
plotGrowth <- function(tlog, param) {
  tmin <- aggregate(x=tlog[,c("neph", "pumptime.s")],
                    by=list(time.min = round(tlog$time.s/60)), median, simplify=TRUE)

  ## pumpfit.s is a loess fit of pumptime.s
  range <- 60
  pumpfit <- loess(pumptime.s ~ time.min, data=tmin,
                   span = range / nrow(tmin), degree=1)
  tmin$pumpfit.s <- predict(pumpfit)

  ## pumpduty is duty fraction of pumping in the minute
  ## computed as d(pumping seconds)/d(time)
  tmin[2:(nrow(tmin)-1), "pumpduty"] <- 
    diff(tmin$pumpfit.s,lag=2)/(diff(tmin$time.min,lag=2)*60)
  
  lastmin <- max(tmin$time.min)

  plot(tmin$time.min / 60, tmin$pumpduty*(param$timevol/param$ttlvol),
       type="l", lwd=2, col="black",
       xlab="Time (hrs)", ylab="Dilution Rate (1/hr)",
       xlim=c(0,108), xaxp=c(0,108,9),
       ylim=c(0,0.3), yaxp=c(0,0.3,3))
  
  abline(v=samples/3600, col=sampleColors, lwd=3)
  axis(side=4, at=log(2)/seq(2,8), labels=seq(2,8), cex.axis=0.67)
  axis(side=4, at=log(2)/seq(2.5,7.5), labels=NA, tcl=-0.25)
}

plotDoublings <- function(tlog, param) {
  tmin <- aggregate(x=tlog[,c("neph", "pumptime.s")],
                    by=list(time.min = round(tlog$time.s/60)), median, simplify=TRUE)

  plot(tmin$time.min / 60, (tmin$pumptime.s/3600)*(param$timevol/param$ttlvol),
       type="l", lwd=2,
       xlab="Time (hrs)", ylab="Doublings",
       xlim=c(0,108), xaxp=c(0,108,9),
       ylim=c(0,20), yaxp=c(0,20,4))

  samplePumpTime <- function(stime.s) {
    max(tmin[tmin$time.min * 60 < stime.s,"pumptime.s"])
  }
  
  samp <- data.frame(time.s = samples)
  samp$pump.s <- sapply(samp$time.s, function(t) { samplePumpTime(t) })
  samp$dbl <- (samp$pump.s/3600)*(param$timevol/param$ttlvol)
  samp$color <- brewer.pal(5, "PuBuGn")[2:5]
    
  abline(v=samp$time.s/3600, col=samp$color)
  abline(h=samp$dbl, col=samp$color)
  
  samp
}

handleTstat <- function(subdir, param) {
  tlog <- getTlog(sprintf("%s/tstat_%s.txt", datadir, subdir))

  pdf(sprintf("%s/fig-%s-neph.pdf", figuredir, subdir), useDingbats=FALSE, width=8, height=3)  
  plotNeph(tlog, param)
  dev.off()
  
  pdf(sprintf("%s/fig-%s-growth.pdf", figuredir, subdir), useDingbats=FALSE, width=8, height=4)
  plotGrowth(tlog, param)
  dev.off()
  
  pdf(sprintf("%s/fig-%s-doubling.pdf", figuredir, subdir), useDingbats=FALSE, width=6, height=6)
  dbl <- plotDoublings(tlog, param)
  dev.off()
  write.csv(x=dbl, file=sprintf("%s/doubling-%s.csv", figuredir, subdir), quote=FALSE, row.names=FALSE)
}


handleTstat("left", leftParam)
handleTstat("right", rightParam)
