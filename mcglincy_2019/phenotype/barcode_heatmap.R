options(stringsAsFactors=FALSE)

library(png)
library("RColorBrewer")

workdir <- "work" ## Sys.getenv("WORKDIR")
figuredir <- "figures" ## Sys.getenv("FIGUREDIR")

guide <- read.csv(sprintf("%s/nizm005-guide-deseq.csv", figuredir), row.names=1)

## These are results from DESeq analysis of barcode present >64 in
## both replicates
bcbothl <- read.csv(sprintf("%s/nizm005-barcode-deseq-left.csv", workdir), row.names=1)
bcbothr <- read.csv(sprintf("%s/nizm005-barcode-deseq-right.csv", workdir), row.names=1)

## All barcodes, then filter for >64 in each replicate individually
countFile <- sprintf("%s/nizm005-assigned.txt", workdir)
cts <- read.delim(countFile, row.names=1)
samples <- colnames(cts)

bcboth <- merge(x = bcbothl, y = bcbothr,
              by = "row.names",
              suffixes=c(".l", ".r"))
bcboth <- bcboth[bcboth$guide.l != "bad",]
row.names(bcboth) <- bcboth$Row.names
bcboth$Row.names <- NULL

lfcMean <- function(lfc.l, se.l, lfc.r, se.r) {
    weighted.mean(x=c(lfc.l, lfc.r),
                  w=c(se.l^(-2), se.r^(-2)))
}

bcboth$log2FoldChange.mean <- mapply(FUN=lfcMean,
                                     bcboth$log2FoldChange.l,
                                     bcboth$lfcSE.l,
                                     bcboth$log2FoldChange.r,
                                     bcboth$lfcSE.r)

bc <- merge(x = bcboth, y = cts,
            by = "row.names")
row.names(bc) <- bc$Row.names
bc$Row.names <- NULL

bc$d1L1r <- log2(pmax(1, bc$d1L1)) - log2(bc$d0L1)
bc$d2L1r <- log2(pmax(1, bc$d2L1)) - log2(bc$d0L1)
bc$d3L1r <- log2(pmax(1, bc$d3L1)) - log2(bc$d0L1)
bc$d3L2r <- log2(pmax(1, bc$d3L2)) - log2(bc$d0L1)

bc$d1R1r <- log2(pmax(1, bc$d1R1)) - log2(bc$d0R1)
bc$d2R1r <- log2(pmax(1, bc$d2R1)) - log2(bc$d0R1)
bc$d3R1r <- log2(pmax(1, bc$d3R1)) - log2(bc$d0R1)
bc$d3R2r <- log2(pmax(1, bc$d3R2)) - log2(bc$d0R1)

cor(bc[,c("log2FoldChange.l", "log2FoldChange.r", "log2FoldChange.mean",
          "d1L1r", "d2L1r", "d3L1r",
          "d1R1r", "d2R1r", "d3R1r")])

bc <- bc[order(bc$log2FoldChange.mean),]

relcts <- as.matrix(cbind.data.frame(bc[,c("d1L1r", "d2L1r", "d3L1r")],
                                     zero=0,
                                     bc[,c("d1R1r", "d2R1r", "d3R1r")]))

fitcolor <- colorRamp(brewer.pal(9, "PRGn"))
guide$color <- fitcolor(pmax(0, pmin(1, (guide$lfcMean + 1)/2)))
guide$colorstr <- rgb(guide$color[,1], guide$color[,2], guide$color[,3], maxColorValue=255)

rowcolors <- guide[match(bc$guide, row.names(guide)), "colorstr"]

relcolor <- colorRampPalette(brewer.pal(9, "PRGn"))(100)

zrange <- 3*3.75

pdf(sprintf("%s/barcode_heatmap_both.pdf", figuredir), useDingbats=FALSE)
heatmap(relcts,
        Rowv=NA, Colv=NA, labRow=NA,
        scale="none",
        zlim=c(-zrange, zrange), col=relcolor,
        RowSideColors=rowcolors)
dev.off()

pdf(sprintf("%s/barcode_heatmap_scale.pdf", figuredir), useDingbats=FALSE, width=6, height=2)
scale <- (length(relcolor) - 1) / (2 * zrange)
plot(c(-zrange, zrange), c(0, 10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main='')
axis(1, at=log2(10**seq(-3,3)),
     labels=c("1/1000", "1/100", "1/10", "1", "10", "100", "1000"),
     las=1, lwd=0, lwd.ticks=1)
for (i in 1:(length(relcolor)-1)) {
    x = (i-1)/scale + (-zrange)
    rect(x, 0, x+1/scale, 10, col=relcolor[i], border=NA)
}
dev.off()

## pdf(sprintf("%s/barcode_heatmap_left.pdf", figuredir), useDingbats=FALSE)
## heatmap(as.matrix(bc[,c("d1L1r", "d2L1r", "d3L1r")]),
##         Rowv=NA, Colv=NA, labRow=NA,
##         verbose=TRUE, scale="none", zlim=c(-3*3.75,3*3.75), col=coul)
## dev.off()

## pdf(sprintf("%s/barcode_heatmap_right.pdf", figuredir), useDingbats=FALSE)
## heatmap(as.matrix(bc[,c("d1R1r", "d2R1r", "d3R1r")]),
##         Rowv=NA, Colv=NA, labRow=NA,
##         verbose=TRUE, scale="none", zlim=c(-3*3.75,3*3.75), col=coul)
## dev.off()

