library("png")
library("RColorBrewer")
options(stringsAsFactors=FALSE)

workdir <- Sys.getenv("WORKDIR")
figuredir <- Sys.getenv("FIGUREDIR")

## All barcodes
counts <- read.delim(sprintf("%s/ivt_all.txt", workdir), row.names=1)
## Discard junk on column names
colnames(counts) <- sub(".nbhd.count.txt", "", colnames(counts))
colnames(counts) <- sub(".*[.]", "", colnames(counts))

## Barcodes present in >1 sample and >32 reads total
good <- counts[rowSums(counts > 0) > 1 & rowSums(counts) > 32,]

## Exclude barcodes with XhoI sites
## Vector can create XhoI sites at the edges
isXho <- (grepl("CTCGAG", row.names(good)) | grepl("^CGAG", row.names(good)) | grepl("CTCGA$", row.names(good)))

xho <- good[isXho,]
noxho <- good[!isXho,]

write.table(noxho, file=sprintf("%s/ivt_good.txt", workdir), row.names=TRUE, sep="\t", quote=FALSE)

## Looking at XhoI drop-outs

l10Decades <- function(xs) {
    log10(c(sapply(xs, function(x) { seq(2,9)*(10**x) })))
}

xhocol <- brewer.pal(3, "PuRd")[[3]]
noxhocol <- brewer.pal(3, "Greys")[[3]]

## Hybrid PDF-PNG scatter plots
## based on https://jonathanchang.org/blog/how-to-partially-rasterize-a-figure-plotted-with-r/
## by jonathan.chang@monash.edu
## We use par(mar=...) and plot(...) rather than plot.window(mar=...) and points(...)

## C replicate
pdf(sprintf("%s/ivt_vs_pcr_c.pdf", figuredir), useDingbats=FALSE, width=5, height=5)
plot(x = pmax(0, log10(noxho$pcr_c)),
     y = pmax(0, log10(noxho$ivt_c)),
     type="n",
     pch=20, cex=0.33, col=noxhocol,
     xlim=c(0, 4.5), ylim=c(0, 4.5),
     axes=FALSE, xlab="PCR", ylab="IVT-RT")

coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- sprintf("%s/ivt_repl_points.png", workdir)
png(pointsfile, width=width, height=height, units="in", res=300, bg="transparent")
par(mar=c(0,0,0,0))
plot(x = pmax(0, log10(noxho$pcr_c)),
     y = pmax(0, log10(noxho$ivt_c)),
     pch=20, cex=0.33, col=noxhocol,
     xlim=c(0, 4.5), ylim=c(0, 4.5),
     axes=FALSE, xlab="PCR", ylab="IVT-RT")
points(x = pmax(0, log10(xho$pcr_c)),
       y = pmax(0, log10(xho$ivt_c)),
       pch=20, cex=0.33, col=xhocol)
dev.off()

panel <- readPNG(pointsfile)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
axis(side=1, at=seq(0,4), labels=c("1", "10", "100", "1k", "10k"))
axis(side=1, at=l10Decades(seq(0,3)), labels=F, tcl=-0.25)
axis(side=2, at=seq(0,4), labels=c("1", "10", "100", "1k", "10k"))
axis(side=2, at=l10Decades(seq(0,3)), labels=F, tcl=-0.25)
legend(x="topleft", bty="n",
       legend=c("No XhoI", "XhoI"),
       pch=20, col=c(noxhocol, xhocol),
       text.col=c(noxhocol, xhocol))
dev.off()

## D replicate

pdf(sprintf("%s/ivt_vs_pcr_d.pdf", figuredir), useDingbats=FALSE, width=5, height=5)
plot(x = pmax(0, log10(noxho$pcr_d)),
     y = pmax(0, log10(noxho$ivt_d)),
     type="n",
     pch=20, cex=0.33, col=noxhocol,
     xlim=c(0, 4.5), ylim=c(0, 4.5),
     axes=FALSE, xlab="PCR", ylab="IVT-RT")

coords <- par("usr")
gx <- grconvertX(coords[1:2], "user", "inches")
gy <- grconvertY(coords[3:4], "user", "inches")
width <- max(gx) - min(gx)
height <- max(gy) - min(gy)

pointsfile <- sprintf("%s/ivt_repl_points.png", workdir)
png(pointsfile, width=width, height=height, units="in", res=300, bg="transparent")
par(mar=c(0,0,0,0))
plot(x = pmax(0, log10(noxho$pcr_d)),
     y = pmax(0, log10(noxho$ivt_d)),
     pch=20, cex=0.33, col=noxhocol,
     xlim=c(0, 4.5), ylim=c(0, 4.5),
     axes=FALSE, xlab="PCR", ylab="IVT-RT")
points(x = pmax(0, log10(xho$pcr_d)),
       y = pmax(0, log10(xho$ivt_d)),
       pch=20, cex=0.33, col=xhocol)
dev.off()

panel <- readPNG(pointsfile)
rasterImage(panel, coords[1], coords[3], coords[2], coords[4])
axis(side=1, at=seq(0,4), labels=c("1", "10", "100", "1k", "10k"))
axis(side=1, at=l10Decades(seq(0,3)), labels=F, tcl=-0.25)
axis(side=2, at=seq(0,4), labels=c("1", "10", "100", "1k", "10k"))
axis(side=2, at=l10Decades(seq(0,3)), labels=F, tcl=-0.25)
legend(x="topleft", bty="n",
       legend=c("No XhoI", "XhoI"),
       pch=20, col=c(noxhocol, xhocol),
       text.col=c(noxhocol, xhocol))
dev.off()
