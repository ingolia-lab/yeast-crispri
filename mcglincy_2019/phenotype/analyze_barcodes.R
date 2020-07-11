options(stringsAsFactors=FALSE)

library(png)
library("RColorBrewer")

barcodingdir <- Sys.getenv("BARCODINGDIR")
workdir <- Sys.getenv("WORKDIR")
figuredir <- Sys.getenv("FIGUREDIR")

l10Decades <- function(xs) {
    log10(c(sapply(xs, function(x) { seq(2,9)*(10**x) })))
}

countFile <- sprintf("%s/nizm005.txt", workdir)

bcToGuideFile <- sprintf("%s/grna-assign-barcode-grna-good.txt", barcodingdir)
bcToGuide <- read.delim(bcToGuideFile)

cts <- read.delim(countFile, row.names=1)
colnames(cts) <- sub("[.]nbhd[.]count[.]txt", "", colnames(cts))
colnames(cts) <- sub(".*[.]", "", colnames(cts))

cts$guide <- bcToGuide[match(row.names(cts), bcToGuide$barcode),"guide"]
cts$guide <- ifelse(is.na(cts$guide), "bad", cts$guide)

write.table(cts, file=sprintf("%s/nizm005-assigned.txt", workdir),
            quote=FALSE, sep="\t")

gatherGuides <- function(x) {
    data.frame(bc1 = row.names(x)[1],
               eff1 = x[1,"eff"],
               bc2 = row.names(x)[2],
               eff2 = x[2,"eff"],
               nbchi = nrow(x))
}

## Selected guide with multiple high-starting-count barcodes
## in both turbidostats
stv1guide <- "YMR054W_08"
sui3guide <- "YPL237W_07"
gensPerSampling <- 3.75

times <- data.frame(row.names=c("d0L1", "d1L1", "d2L1", "d3L1", "d3L2",
                                "d0R1", "d1R1", "d2R1", "d3R1", "d3R2"),
                    times=c(0,1,2,3,3,0,1,2,3,3)*gensPerSampling,
                    jitter=c(-0.03, -0.03, -0.03, -0.03, -0.15,
                             0.09, 0.09, 0.09, 0.09, 0.21),
                    tstat=c("L", "L", "L", "L", "L", "R", "R", "R", "R", "R"))

stv1 <- as.data.frame(t(cts[cts$guide == stv1guide,grep("^d", colnames(cts))]))
stv1$times <- times[match(row.names(stv1), row.names(times)), "times"]
stv1$jitter <- times[match(row.names(stv1), row.names(times)), "jitter"]
stv1$tstat <- times[match(row.names(stv1), row.names(times)), "tstat"]

sui3 <- as.data.frame(t(cts[cts$guide == sui3guide,grep("^d", colnames(cts))]))
sui3$times <- times[match(row.names(sui3), row.names(times)), "times"]
sui3$jitter <- times[match(row.names(sui3), row.names(times)), "jitter"]
sui3$tstat <- times[match(row.names(sui3), row.names(times)), "tstat"]

empty <- as.data.frame(t(cts[cts$guide == "Empty", grep("^d", colnames(cts))]))
empty$times <- times[match(row.names(empty), row.names(times)), "times"]
empty$jitter <- times[match(row.names(empty), row.names(times)), "jitter"]
empty$tstat <- times[match(row.names(empty), row.names(times)), "tstat"]

pal <- brewer.pal(7, "BrBG")
palmore <- brewer.pal(7, "PRGn")

pdf(sprintf("%s/sui3.pdf", figuredir), useDingbats=FALSE, width=3, height=4)
plot(x=sui3$times + sui3$jitter, y=log10(pmax(10, sui3[,1])),
     pch=20, cex=1.5, col=ifelse(sui3$tstat == "L", pal[[1]], pal[[2]]),
     xlim=c(0,3.25*gensPerSampling), ylim=c(1, 4),
     axes=FALSE,
     xlab="Pop Doublings", ylab="Read Count")

points(x=sui3$times + sui3$jitter+0.03, y=log10(pmax(10, sui3[,2])),
       pch=20, cex=1.5, col=ifelse(sui3$tstat == "L", pal[[7]], pal[[6]]))

axis(side=1, at=seq(0,12,4))

axis(side=2, at=seq(1,4), labels=c("10", "100", "1k", "10k"))
axis(side=2, at=l10Decades(seq(1,3)), labels=F, tcl=-0.25)

legend(x="topright", bty="n",
       pch=20, pt.cex=1.5, col=pal[c(1,2,7,6)],
       cex=0.5, legend=c("barcode 1, repl #1", "barcode 1, repl #2", "barcode 2, repl #1", "barcode 2, repl #2"),
       text.col=pal[c(1,2,7,6)])

dev.off()

pdf(sprintf("%s/stv1.pdf", figuredir), useDingbats=FALSE, width=3, height=4)
plot(x=stv1$times + stv1$jitter, y=log10(pmax(10, stv1[,1])),
     pch=20, cex=1.5, col=ifelse(stv1$tstat == "L", pal[[1]], pal[[2]]),
     xlim=c(0,3.25*gensPerSampling), ylim=c(1, 4),
     axes=FALSE,
     xlab="Pop Doublings", ylab="Read Count")

points(x=stv1$times + stv1$jitter+0.03, y=log10(pmax(10, stv1[,2])),
       pch=20, cex=1.5, col=ifelse(stv1$tstat == "L", pal[[7]], pal[[6]]))

axis(side=1, at=seq(0,12,4))

axis(side=2, at=seq(1,4), labels=c("10", "100", "1k", "10k"))
axis(side=2, at=l10Decades(seq(1,3)), labels=F, tcl=-0.25)

dev.off()

pdf(sprintf("%s/empty.pdf", figuredir), useDingbats=FALSE, width=3, height=4)

plot(x=empty$times + empty$jitter, y=log10(pmax(10, empty[,2])),
     pch=20, cex=1.5, col=ifelse(empty$tstat == "L", pal[[1]], pal[[2]]),
     xlim=c(0,3.25*gensPerSampling), ylim=c(1, 4),
     axes=FALSE,
     xlab="Pop Doublings", ylab="Read Count")

points(x=empty$times + empty$jitter+0.03, y=log10(pmax(10, empty[,3])),
       pch=20, cex=1.5, col=ifelse(empty$tstat == "L", pal[[7]], pal[[6]]))

points(x=empty$times + empty$jitter+0.03, y=log10(pmax(10, empty[,4])),
       pch=20, cex=1.5, col=ifelse(empty$tstat == "L", palmore[[1]], palmore[[2]]))

axis(side=1, at=seq(0,12,4))

axis(side=2, at=seq(1,4), labels=c("10", "100", "1k", "10k"))
axis(side=2, at=l10Decades(seq(1,3)), labels=F, tcl=-0.25)

dev.off()
