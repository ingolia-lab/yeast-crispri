options(stringsAsFactors=FALSE)

library(png)
library("RColorBrewer")

guidedir <- Sys.getenv("GUIDEDIR")
barcodingdir <- Sys.getenv("BARCODINGDIR")
workdir <- Sys.getenv("WORKDIR")
figuredir <- Sys.getenv("FIGUREDIR")

log <- file(sprintf("%s/compare_out.txt", figuredir), "w")

## SGD annotations for yorf-to-gene mapping
sgdFile <- sprintf("%s/SGD_features.tab", workdir)
if (!file.exists(sgdFile)) {
  sgd <- download.file('https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab', destfile=sgdFile)
}
sgd <- read.delim(sgdFile, header=FALSE, quote="",
                  col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                              "parent", "sgdid2", "chrom", "start", "end",
                              "strand", "genpos", "cver", "sver", "desc"))

## Guide to target yorf mapping
guideToYorfFile <- sprintf("%s/sequence-good-targets.txt", guidedir)
guideToYorfs <- read.delim(guideToYorfFile)

## barcode counts with guides assigned
countFile <- sprintf("%s/nizm005-assigned.txt", workdir)
cts <- read.delim(countFile, row.names=1)

## barcode DESeq results
bcDeseqFile <- sprintf("%s/nizm005-barcode-deseq.csv", workdir)
bcres <- read.csv(bcDeseqFile, row.names=1)
bcresleft <- read.csv(sprintf("%s/nizm005-barcode-deseq-left.csv", workdir), row.names=1)
colnames(bcresleft) <- paste0(colnames(bcresleft), "l")
bcresright <- read.csv(sprintf("%s/nizm005-barcode-deseq-right.csv", workdir), row.names=1)
colnames(bcresright) <- paste0(colnames(bcresright), "r")
if (!all(row.names(bcresleft) == row.names(bcresright))) {
    stop("Mismatch between left and right DESeq")
}
bcrescf <- cbind.data.frame(bcresleft, bcresright)

bcres$d0L <- cts[match(row.names(bcres), row.names(cts)), "d0L1"]
bcres$d0R <- cts[match(row.names(bcres), row.names(cts)), "d0R1"]

bcrescf$d0L <- cts[match(row.names(bcrescf), row.names(cts)), "d0L1"]
bcrescf$d0R <- cts[match(row.names(bcrescf), row.names(cts)), "d0R1"]
bcrescf$d0min <- pmin(bcrescf$d0L, bcrescf$d0R)

compareBarcodes <- function(x) {
    xsorted <- x[order(x$baseMean, decreasing=TRUE),]
    if (nrow(xsorted) < 2) {
        data.frame()
    } else {
        data.frame(baseMean1 = xsorted$baseMean[[1]],
                   baseMean2 = xsorted$baseMean[[2]],
                   lfc1 = xsorted$log2FoldChange[[1]],
                   lfc2 = xsorted$log2FoldChange[[2]],
                   d0L1 = xsorted$d0L[[1]],
                   d0R1 = xsorted$d0R[[1]],
                   d0L2 = xsorted$d0L[[2]],
                   d0R2 = xsorted$d0R[[2]])
    }
}

byout <- by(data=bcres, INDICES=bcres$guide, FUN=compareBarcodes)
barcodeCompare <- do.call("rbind", byout)
barcodeCompare <- barcodeCompare[order(barcodeCompare$baseMean2, decreasing=TRUE),]

barcodeCompare$d0min <- pmin(barcodeCompare$d0L1, barcodeCompare$d0R1, barcodeCompare$d0L2, barcodeCompare$d0R2)

write.csv(barcodeCompare, file=sprintf("%s/nizm005-barcode-compare.csv", figuredir))

for (cutoff in c(64,128,256,512)) {
    cat(sprintf("Comparing barcode pairs at >= %d reads", cutoff), file=log)
    cat(sprintf("  Using %d barcode pairs", sum(barcodeCompare$d0min >= cutoff)), file=log)
    cat(sprintf("  Correlation = %0.3f\n",
                cor(barcodeCompare[barcodeCompare$d0min >= cutoff, "lfc1"],
                    barcodeCompare[barcodeCompare$d0min >= cutoff, "lfc2"])),
        file=log)
}

for (cutoff in c(64,128,256,512)) {
    cat(sprintf("Comparing turbidostats at >= %d reads", cutoff), file=log)
    cat(sprintf("  Using %d barcodes", sum(bcrescf$d0min >= cutoff)), file=log)
    cat(sprintf("  Correlation = %0.3f\n",
                cor(bcrescf[bcrescf$d0min >= cutoff, "log2FoldChangel"],
                    bcrescf[bcrescf$d0min >= cutoff, "log2FoldChanger"])),
        file=log)
}
