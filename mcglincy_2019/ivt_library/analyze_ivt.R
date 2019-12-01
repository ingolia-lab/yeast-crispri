options(stringsAsFactors=FALSE)

workdir <- Sys.getenv("WORKDIR")
figuredir <- Sys.getenv("FIGUREDIR")

library(png)
library("RColorBrewer")

l10Decades <- function(xs) {
    log10(c(sapply(xs, function(x) { seq(2,9)*(10**x) })))
}

## Good barcodes
x <- read.delim(sprintf("%s/ivt_good.txt", workdir), row.names=1)

pal <- brewer.pal(4, "PRGn")
ivtcolor <- pal[[1]]
pcrcolor <- pal[[4]]

## Hybrid PDF-PNG scatter plots
## based on https://jonathanchang.org/blog/how-to-partially-rasterize-a-figure-plotted-with-r/
## by jonathan.chang@monash.edu
## We use par(mar=...) and plot(...) rather than plot.window(mar=...) and points(...)

repls <- data.frame(row.names = row.names(x),
                    ivtC = log10(pmax(0.8, x$ivt_c)),
                    ivtD = log10(pmax(0.8, x$ivt_d)),
                    pcrC = log10(pmax(0.8, x$pcr_c)),
                    pcrD = log10(pmax(0.8, x$pcr_d)))

pdf(sprintf("%s/ivt_repl.pdf", figuredir), useDingbats=FALSE, width=5, height=5)
plot(x = repls$ivtC, y = repls$ivtD, type="n", 
     xlim=log10(c(0.7, 50000)), ylim=log10(c(0.7, 50000)),
     axes=FALSE,
     xlab="Read Count, Repl #1", ylab="Read Count, Repl #2")

coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- sprintf("%s/ivt_repl_points.png", workdir)
png(pointsfile, width=width, height=height, units="in", res=300, bg="transparent")
par(mar=c(0,0,0,0))
plot(x = repls$ivtC, y = repls$ivtD,
     pch=20, cex=0.5, col=ivtcolor,
     axes=FALSE, xlim=log10(c(0.7, 50000)), ylim=log10(c(0.7, 50000)))
dev.off()

panel <- readPNG(pointsfile)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
axis(side=1, at=seq(0,4), labels=c("1", "10", "100", "1k", "10k"))
axis(side=1, at=l10Decades(seq(0,3)), labels=F, tcl=-0.25)
axis(side=2, at=seq(0,4), labels=c("1", "10", "100", "1k", "10k"))
axis(side=2, at=l10Decades(seq(0,3)), labels=F, tcl=-0.25)
title(main="IVT-RT", col.main=ivtcolor)
dev.off()

pdf(sprintf("%s/pcr_repl.pdf", figuredir), useDingbats=FALSE, width=5, height=5)
plot(x = repls$pcrC, y = repls$pcrD, type="n", 
     xlim=log10(c(0.7, 50000)), ylim=log10(c(0.7, 50000)),
     axes=FALSE,
     xlab="Read Count, Repl #1", ylab="Read Count, Repl #2")

coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- sprintf("%s/pcr_repl_points.png", workdir)
png(pointsfile, width=width, height=height, units="in", res=300, bg="transparent")
par(mar=c(0,0,0,0))
plot(x = repls$pcrC, y = repls$pcrD,
     pch=20, cex=0.5, col=pcrcolor,
     axes=FALSE, xlim=log10(c(0.7, 50000)), ylim=log10(c(0.7, 50000)))
dev.off()

panel <- readPNG(pointsfile)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
axis(side=1, at=seq(0,4), labels=c("1", "10", "100", "1k", "10k"))
axis(side=1, at=l10Decades(seq(0,3)), labels=F, tcl=-0.25)
axis(side=2, at=seq(0,4), labels=c("1", "10", "100", "1k", "10k"))
axis(side=2, at=l10Decades(seq(0,3)), labels=F, tcl=-0.25)
title(main="PCR", col.main=pcrcolor)
dev.off()

## Averages and replicate ratios for each library type
x$pcrAvg <- (x$pcr_c + x$pcr_d) / 2
x$pcrL2R <- pmin(4, pmax(-4, log2(x$pcr_d / x$pcr_c)))

x$ivtAvg <- (x$ivt_c + x$ivt_d) / 2
x$ivtL2R <- pmin(4, pmax(-4, log2(x$ivt_d / x$ivt_c)))

xhi <- x[x$ivtAvg > 100,]

rbreaks <- seq(-3,3,0.1)
ivt <- hist(pmin(3, pmax(-3, xhi$ivtL2R - median(xhi$ivtL2R, na.rm=TRUE))),
            breaks=rbreaks, plot=FALSE)
pcr <- hist(pmin(3, pmax(-3, xhi$pcrL2R - median(xhi$pcrL2R, na.rm=TRUE))),
            breaks=rbreaks, plot=FALSE)

pdf(sprintf("%s/repl_hist.pdf", figuredir), width=5, height=5, useDingbats=FALSE)
plot(ivt$mids, ivt$density,
     type="l", lwd=2, col=ivtcolor,
     xlim=c(-3,3), ylim=c(0,2),
     axes=FALSE, xlab="Replicate Ratio", ylab=NA)
lines(pcr$mids, pcr$density,
      lwd=2, col=pcrcolor)
axis(side=1, at=seq(-3,3), labels=c("1/8", "1/4", "1/2", "1", "2", "4", "8"))
legend(x="topleft", 
       legend=c("IVT replicates", "PCR replicates"),
       lwd=3, col=c(ivtcolor, pcrcolor), text.col=c(ivtcolor, pcrcolor),
       bty="n")
dev.off()

library(DESeq2)

condIvt <- data.frame(row.names=c("ivt_c", "ivt_d"),
                      libgen=c("ivt", "ivt"))
countsIvt <- x[,colnames(x) %in% row.names(condIvt)]
ddsIvt <- DESeqDataSetFromMatrix(countData = countsIvt,
                                 colData = condIvt,
                                 design = ~ 1)
ddsIvt <- estimateSizeFactors(ddsIvt)
ddsIvt <- estimateDispersions(ddsIvt)

condPcr <- data.frame(row.names=c("pcr_c", "pcr_d"),
                      libgen=c("pcr", "pcr"))
countsPcr <- x[,colnames(x) %in% row.names(condPcr)]
ddsPcr <- DESeqDataSetFromMatrix(countData = countsPcr,
                                 colData = condPcr,
                                 design = ~ 1)
ddsPcr <- estimateSizeFactors(ddsPcr)
ddsPcr <- estimateDispersions(ddsPcr)

pdf(sprintf("%s/compare_cv.pdf", figuredir), useDingbats=FALSE, width=6, height=4)
avgs <- 10^seq(log10(10), log10(5000), length.out=250)
plot(x=log10(avgs),
     y=log10(sqrt(dispersionFunction(ddsIvt)(avgs))),
     type="l", lwd=2, col=ivtcolor,
     xlim=c(1, log10(5000)), ylim=c(-2, 0),
     axes=FALSE,
     xlab="Read Count", ylab="Replicate CV")
lines(x=log10(avgs),
     y=log10(sqrt(dispersionFunction(ddsPcr)(avgs))),
     lwd=2, col=pcrcolor)
axis(side=1, at=seq(1,3), labels=c("10", "100", "1k"))
axis(side=1, at=l10Decades(seq(1,3)), labels=F, tcl=-0.25)
axis(side=2, at=seq(-2, 0), labels=c("1%", "10%", "100%"))
axis(side=2, at=l10Decades(seq(-2, -1)), labels=F, tcl=-0.25)
legend(x="bottomleft", bty="n",
       legend=c("IVT-RT", "PCR"),
       lwd=3, col=c(ivtcolor, pcrcolor), text.col=c(ivtcolor, pcrcolor))
dev.off()
