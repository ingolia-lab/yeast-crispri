options(stringsAsFactors=FALSE)

library(png)
library("RColorBrewer")

guidedir <- Sys.getenv("GUIDEDIR")
barcodingdir <- Sys.getenv("BARCODINGDIR")
workdir <- Sys.getenv("WORKDIR")
figuredir <- Sys.getenv("FIGUREDIR")

log <- file(sprintf("%s/analyze_fitness_out.txt", figuredir), "w")

## Empirical doubling time
tdbl <- 3.1

## Essential gene list
essFile <- sprintf("%s/Essential.txt", workdir)
if (!file.exists(essFile)) {
  essential <- download.file('http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt', destfile=essFile)
}
essential <- read.delim(essFile)
essentialYorf <- sub(" *$", "", essential$ORF_name)
cat(sprintf("%d essential genes\n", length(essentialYorf)), file=log)

## Information on guide targeting
guideToYorfFile <- sprintf("%s/sequence-good-targets.txt", guidedir)
guideToYorfs <- read.delim(guideToYorfFile)

## Genes with clean targeting
cleanYorfFile <- sprintf("%s/yorf-clean-target.txt", guidedir)
cleanYorfs <- read.delim(cleanYorfFile, header=FALSE)$V1
cat(sprintf("%d genes with \"clean\" (unique, specific, TSS-based) guide choices\n", length(cleanYorfs)), file=log)

cleanEssYorfs <- essentialYorf[essentialYorf %in% cleanYorfs]
cat(sprintf("  %d essential genes with clean guide choices\n", length(cleanEssYorfs)), file=log)

cleanViaYorfs <- cleanYorfs[!(cleanYorfs %in% essentialYorf)]
cat(sprintf("  %d non-essential genes with clean guide choices\n", length(cleanViaYorfs)), file=log)

## Needed for "Empty" barcodes
bcDeseqFile <- sprintf("%s/nizm005-barcode-deseq.csv", workdir)
bcres <- read.csv(bcDeseqFile, row.names=1)

guideDeseqFile <- sprintf("%s/nizm005-guide-deseq.csv", workdir)
guideres <- read.csv(guideDeseqFile, row.names=1)

breaks <- seq(-1.35, 0.45, 0.05)
guideAllHist <- hist(guideres[guideres$Yorfs %in% cleanYorfs,]$lfcMean, breaks=breaks, plot=FALSE)
guideEssHist <- hist(guideres[guideres$Yorfs %in% cleanEssYorfs,]$lfcMean, breaks=breaks, plot=FALSE)
guideViaHist <- hist(guideres[guideres$Yorfs %in% cleanViaYorfs,]$lfcMean, breaks=breaks, plot=FALSE)
bcEmptyHist <- hist(bcres[grepl("Empty", bcres$guide),]$log2FoldChange, breaks=breaks, plot=FALSE)

pdf(sprintf("%s/guide-fitness.pdf", figuredir), useDingbats=FALSE, width=6, height=4)
pal <- brewer.pal(5, "Dark2")
plot(guideAllHist$mids, guideAllHist$density, 
     type="l", lwd=2, col="black", 
     xlim=c(-1.4, 0.6), ylim=c(0,4), 
     xlab="Fitness (log2(s))", ylab="Frequency of Guides", 
     xaxp=c(-1.5,0.5,4), yaxt="n")
lines(guideEssHist$mids, guideEssHist$density,
      lwd=2, col=pal[[2]])
lines(guideViaHist$mids, guideViaHist$density,
      lwd=2, col=pal[[5]])
lines(bcEmptyHist$mids, bcEmptyHist$density,
      lwd=2, col=pal[[3]])
axis(side=3, at=tdbl/c(12,6,4,3,2)-1, labels=c(12,6,4,3,2))
axis(side=3, at=-1, labels=expression(infinity))
axis(side=3, at=tdbl/seq(2,12)-1, labels=NA)
axis(side=3, at=tdbl/seq(2.5,11.5)-1, labels=NA, tcl=-0.25)
legend(x="topleft", bty="n",
       legend=c("All Genes", "Essential", "Non-essential", "No gRNA"),
       lwd=3, col=c("black", pal[[2]], pal[[5]], pal[[3]]),
       text.col=c("black", pal[[2]], pal[[5]], pal[[3]]))
dev.off()

mergeGuides <- function(x) {
  data.frame(lfcMin = min(x$lfcMean), 
             lfcMax = max(x$lfcMean),
             grnaMin = row.names(x)[which.min(x$lfcMean)],
             grnaMax = row.names(x)[which.max(x$lfcMean)],
             nguide = nrow(x))
}

byout <- by(data=guideres, INDICES=guideres$Yorf1, FUN=mergeGuides)
yorfres <- do.call("rbind", byout)
yorfres$clean <- row.names(yorfres) %in% cleanYorfs
yorfres$essential <- row.names(yorfres) %in% essentialYorf

write.csv(yorfres, file=sprintf("%s/nizm005-yorf-deseq.csv", workdir))

yorflfc <- data.frame(yorf=row.names(yorfres),
                      lfc=yorfres$lfcMin,
                      essential=yorfres$essential)

write.table(row.names(yorfres[yorfres$essential & yorfres$lfcMin > -0.5,]),
            file=sprintf("%s/nizm005-essential-nosick.txt", figuredir),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(row.names(yorfres[yorfres$essential,]),
            file=sprintf("%s/nizm005-essential-all.txt", figuredir),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(row.names(yorfres[!yorfres$essential & yorfres$lfcMin < -0.5,]),
            file=sprintf("%s/nizm005-viable-sick.txt", figuredir),
            sep="\t", quote=FALSE, row.names=FALSE)
write.table(row.names(yorfres[!yorfres$essential,]),
            file=sprintf("%s/nizm005-viable-all.txt", figuredir),
            sep="\t", quote=FALSE, row.names=FALSE)
            
yorfAllC <- ecdf(yorfres[yorfres$clean,]$lfcMin)
yorfEssC <- ecdf(yorfres[yorfres$clean & yorfres$essential,]$lfcMin)
yorfViaC <- ecdf(yorfres[yorfres$clean & !yorfres$essential,]$lfcMin)
bcEmptyC <- ecdf(bcres[grepl("Empty", bcres$guide),]$log2FoldChange)

pdf(sprintf("%s/gene-fitness.pdf", figuredir), useDingbats=FALSE, width=6, height=4)
plot(breaks, yorfAllC(breaks), 
     type="l", lwd=2, col="black", 
     xlim=c(-1.4, 0.6), ylim=c(0,1), yaxt="n", 
     xlab="lg(s)", ylab="Fraction Guides", 
     xaxp=c(-1.5,0.5,4))
lines(breaks, yorfEssC(breaks), lwd=2, col=pal[[2]])
lines(breaks, yorfViaC(breaks), lwd=2, col=pal[[5]])
lines(breaks, bcEmptyC(breaks), lwd=2, col=pal[[3]])
axis(side=3, at=tdbl/c(12,6,4,3,2)-1, labels=c(12,6,4,3,2))
axis(side=3, at=-1, labels=expression(infinity))
axis(side=3, at=tdbl/seq(2,12)-1, labels=NA)
axis(side=3, at=tdbl/seq(2.5,11.5)-1, labels=NA, tcl=-0.25)
legend(x="topleft", bty="n",
       legend=c("All Genes", "Essential", "Non-essential", "No gRNA"),
       lwd=3, col=c("black", pal[[2]], pal[[5]], pal[[3]]),
       text.col=c("black", pal[[2]], pal[[5]]), pal[[3]])
dev.off()

