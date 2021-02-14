options(stringsAsFactors=FALSE)

library(png)
library("RColorBrewer")

workdir <- "work" ## Sys.getenv("WORKDIR")
figuredir <- "figures" ## Sys.getenv("FIGUREDIR")

log <- file(sprintf("%s/replicates_out.txt", figuredir), "w")

mergeBarcodes <- function(x) {
    wmean <- weighted.mean(x=x$log2FoldChange, w=1/(x$lfcSE^2))
    data.frame(lfcMean=wmean, 
               lfcSE=sqrt(1/sum(1/(x$lfcSE^2))),
               nbc=nrow(x))
}

## These are results from DESeq analysis of barcode present >64 in
## both replicates
bcbothl<- read.csv(sprintf("%s/nizm005-barcode-deseq-left.csv", workdir), row.names=1)
bcbothr <- read.csv(sprintf("%s/nizm005-barcode-deseq-right.csv", workdir), row.names=1)

byout <- by(data=bcbothl, INDICES=bcbothl$guide, FUN=mergeBarcodes)
guidebothl <- do.call("rbind", byout)

byout <- by(data=bcbothr, INDICES=bcbothr$guide, FUN=mergeBarcodes)
guidebothr <- do.call("rbind", byout)

bcboth <- merge(x = bcbothl, y = bcbothr,
              by = "row.names",
              suffixes=c(".l", ".r"))
bcbothgood <- bcboth[bcboth$guide.l != "bad",]

cat(sprintf("%d barcodes present in both replicates\n", nrow(bcbothgood)), file=log)
cat(sprintf("    Correlation %0.3f between replicates\n", 
            cor(bcbothgood$log2FoldChange.l, bcbothgood$log2FoldChange.r)),
    file=log)

guideboth <- merge(x = guidebothl, y = guidebothr,
                   by = "row.names",
                   suffixes=c(".l", ".r"))

cat(sprintf("%d guides for barcodes present in both replicates\n", nrow(guideboth)), file=log)
cat(sprintf("    Correlation %0.3f between replicates\n", 
            cor(guideboth$lfcMean.l, guideboth$lfcMean.r)),
    file=log)

## All barcodes, then filter for >64 in each replicate individually
countFile <- sprintf("%s/nizm005-assigned.txt", workdir)
cts <- read.delim(countFile, row.names=1)
samples <- colnames(cts)

library("DESeq2")

ctsl <- cts[cts$d0L1 > 64, c("d0L1", "d1L1", "d2L1", "d3L1", "d3L2", "guide")]
condl <- data.frame(gens=c(0,1,2,3,3)*3.75,
                    row.names=c("d0L1", "d1L1", "d2L1", "d3L1", "d3L2"))
ddsl <- DESeqDataSetFromMatrix(countData = ctsl[,row.names(condl)],
                               colData = condl,
                               design = ~ gens)
ddsl <- estimateSizeFactors(ddsl)
ddsl <- estimateDispersions(ddsl)
ddsl <- nbinomWaldTest(ddsl, betaPrior=FALSE)

resl <- results(ddsl, name="gens")
resl <- as.data.frame(resl)
resl$guide <- ctsl$guide

cat(sprintf("%d barcodes in left replicate\n", nrow(resl)), file=log)

byout <- by(data=resl, INDICES=resl$guide, FUN=mergeBarcodes)
guidesl <- do.call("rbind", byout)

cat(sprintf("%d guides in left replicate\n", nrow(guidesl)), file=log)

ctsr <- cts[cts$d0R1 > 64, c("d0R1", "d1R1", "d2R1", "d3R1", "d3R2", "guide")]

condr <- data.frame(gens=c(0,1,2,3,3)*3.75,
                    row.names=c("d0R1", "d1R1", "d2R1", "d3R1", "d3R2"))

ddsr <- DESeqDataSetFromMatrix(countData = ctsr[,row.names(condr)],
                               colData = condr,
                               design = ~ gens)
ddsr <- estimateSizeFactors(ddsr)
ddsr <- estimateDispersions(ddsr)
ddsr <- nbinomWaldTest(ddsr, betaPrior=FALSE)

resr <- results(ddsr, name="gens")
resr <- as.data.frame(resr)
resr$guide <- ctsr$guide

cat(sprintf("%d barcodes in right replicate\n", nrow(resr)), file=log)

byout <- by(data=resr, INDICES=resr$guide, FUN=mergeBarcodes)
guidesr <- do.call("rbind", byout)

cat(sprintf("%d guides in right replicate\n", nrow(guidesr)), file=log)

guides <- merge(x = guidesl, y = guidesr,
                by = "row.names",
                suffixes=c(".l", ".r"))
cat(sprintf("%d guides in left / right intersection\n", nrow(guidesr)), file=log)
cat(sprintf("    Correlation %0.3f between guides\n", 
            cor(guides$lfcMean.l, guides$lfcMean.r)),
    file = log)

guidesmulti <- guides[guides$nbc.l > 1 & guides$nbc.r > 1,]
cat(sprintf("%d guides with >1 barcode in left / right intersection\n", nrow(guidesmulti)),
    file=log)
cat(sprintf("    Correlation %0.3f between guides\n", 
            cor(guidesmulti$lfcMean.l, guidesmulti$lfcMean.r)),
    file = log)
    
write.csv(guides, sprintf("%s/nizm005-guide-deseq-compare.csv", workdir))
