options(stringsAsFactors=FALSE)

library(png)
library("RColorBrewer")
library("DESeq2")

guidedir <- Sys.getenv("GUIDEDIR")
barcodingdir <- Sys.getenv("BARCODINGDIR")
workdir <- Sys.getenv("WORKDIR")
figuredir <- Sys.getenv("FIGUREDIR")

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
samples <- colnames(cts)

library("DESeq2")

cond <- data.frame(gens=c(0,1,2,3,3,0,1,2,3,3)*3.75,
                   tstat=c("L", "L", "L", "L", "L", "R", "R", "R", "R", "R"),
                   row.names=c("d0L1", "d1L1", "d2L1", "d3L1", "d3L2",
                               "d0R1", "d1R1", "d2R1", "d3R1", "d3R2"))
ctshi <- cts[cts$d0L1 > 64 | cts$d0R1 > 64,]
dds <- DESeqDataSetFromMatrix(countData = ctshi[,row.names(cond)],
                              colData = cond,
                              design = ~ gens + tstat)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds, betaPrior=FALSE)

gens <- results(dds, name="gens")
gens <- as.data.frame(gens)

gens$guide <- ctshi$guide

write.csv(gens, file=sprintf("%s/nizm005-barcode-deseq.csv", workdir))

mergeBarcodes <- function(x) {
  wmean <- weightedMean(x=x$log2FoldChange, w=1/(x$lfcSE^2))
  data.frame(lfcMean=wmean, 
             lfcMed=median(x$log2FoldChange),
             lfcSE=sqrt(1/sum(1/(x$lfcSE^2))),
             nbc=nrow(x))
}

byout <- by(data=gens, INDICES=gens$guide, FUN=mergeBarcodes)
guideres <- do.call("rbind", byout)

guideres <- cbind.data.frame(guideres,
                             guideToYorfs[match(row.names(guideres), guideToYorfs$Guide),])
guideres$Guide <- NULL
guideres$Gene1 <- sgd[match(guideres$Yorf1, sgd$name),"gene"]
guideres$Gene1 <- ifelse(is.na(guideres$Gene1), guideres$Yorf1, guideres$Gene1)
guideres <- guideres[,c("lfcMean", "lfcSE", "nbc", "Gene1", "Yorf1", "Offset1", "StartType1",
                        "Yorf2", "Offset2", "StartType2", "Yorfs", "TargetLoc", "Oligo")]

write.csv(guideres, file=sprintf("%s/nizm005-guide-deseq.csv", figuredir))

guideTable <- guideToYorfs[,c("Guide", "TargetLoc", "Yorf1", "Offset1", "StartType1", "Yorfs", "Oligo")]
guideTable$Gene1 <- sgd[match(guideTable$Yorf1, sgd$name), "gene"]
guideTable$Gene1 <- ifelse(is.na(guideTable$Gene1), guideTable$Yorf1, guideTable$Gene1)
guideTable$lfcMean <- round(guideres[match(row.names(guideres), guideTable$Guide), "lfcMean"], digits=2)
guideTable$lfcSE <- round(guideres[match(row.names(guideres), guideTable$Guide), "lfcSE"], digits=2)
guideTable$nbc <- round(guideres[match(row.names(guideres), guideTable$Guide), "nbc"], digits=2)

write.csv(guideTable, file=sprintf("%s/guide-table-nizm005-deseq.csv", figuredir))

## DESeq2 estimates are log2 fold-change in 1 (wild-type) generation
## selective coefficient s is 2^(log2 fold-change)
##   or alternately lg(s) is log2 fold-change
## relative fitness w is 1 + s
## s = exp(kg TD0) / exp(k0 TD0)
##   by def'n, k0 TD0 = ln(2) and exp(k0 TD0) = 2
## s = (1/2) exp(kg TD0)
##   by defn'n, kg TDg = ln(2) and kg = ln(2) / TDg
## s = (1/2) exp(ln(2) TD0 / TDg)
##   = (1/2) (exp(ln(2)))^(TD0 / TDg)
##   = (1/2) 2^(TD0 / TDg)
## 2s = 2^(TD0 / TDg)
## lg(2s) = TD0/TDg
## 1 + lg(s) = TD0/TDg
## TDg = TD0 * 1/(1 + lg(s))
## lg(s) = TD0/TDg - 1

