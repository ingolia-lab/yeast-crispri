options(stringsAsFactors=FALSE)

library(png)
library("RColorBrewer")

guidedir <- Sys.getenv("GUIDEDIR")
barcodingdir <- Sys.getenv("BARCODINGDIR")
workdir <- Sys.getenv("WORKDIR")
figuredir <- Sys.getenv("FIGUREDIR")

log <- file(sprintf("%s/analyze_guides_out.txt", figuredir), "w")

guideToYorfFile <- sprintf("%s/sequence-good-targets.txt", guidedir)
guideToYorfs <- read.delim(guideToYorfFile)

guideChoiceFile <- sprintf("%s/guide_choice/output/target-choices.txt", guidedir)
guideChoice <- read.delim(guideChoiceFile)

guideDeseqFile <- sprintf("%s/nizm005-guide-deseq.csv", figuredir)
guideres <- read.csv(guideDeseqFile, row.names=1)

yorfresFile <- sprintf("%s/nizm005-yorf-deseq.csv", workdir)
yorfres <- read.csv(yorfresFile, row.names=1)

## Essential genes with unambiguous promoters and a strong-effect guide
sickThreshold <- -0.5
sickYorfs <- row.names(yorfres[yorfres$clean & yorfres$essential & yorfres < sickThreshold,])
cat(sprintf("%d essential genes with clean guide choices and strong phenotype\n", length(sickYorfs)), file=log)

## Guides targeting these genes for guide analysis
allGuides <- guideChoice[guideChoice$Yorf %in% sickYorfs,]
row.names(allGuides) <- guideToYorfs[match(allGuides$TargetLoc, guideToYorfs$TargetLoc), "Guide"]
cat(sprintf("%d guides against essential, clean, strong phenotype genes\n", nrow(allGuides)), file=log)

allGuides$fitness <- guideres[match(row.names(allGuides), row.names(guideres)),"lfcMean"]
seenGuides <- allGuides[!is.na(allGuides$fitness),]
cat(sprintf("%d guides analyzed in DESeq\n", nrow(seenGuides)), file=log)

## Features for nucleotide composition correlation
for (nt in seq(1,20)) {
    seenGuides[,sprintf("nt%02d", nt)] <- substr(seenGuides$TargetSeq, nt, nt)
}

ntfit <- lm(seenGuides$fitness ~ seenGuides$nt01 + seenGuides$nt02 + seenGuides$nt03 + seenGuides$nt04 + seenGuides$nt05
            + seenGuides$nt06 + seenGuides$nt07 + seenGuides$nt08 + seenGuides$nt09 + seenGuides$nt10
            + seenGuides$nt11 + seenGuides$nt12 + seenGuides$nt13 + seenGuides$nt14 + seenGuides$nt15
            + seenGuides$nt16 + seenGuides$nt17 + seenGuides$nt18 + seenGuides$nt19 + seenGuides$nt20)

sink(log)
summary(ntfit)
sink()

##
cat(sprintf("fitness vs ATAC = %0.3f", cor(seenGuides$fitness, seenGuides$ATAC)), file=log)

mergeOffsetBins <- function(x) {
  data.frame(n = nrow(x), 
             meanFit = mean(x$fitness), 
             meanAtac = mean(x$ATAC))
}

byout <- by(seenGuides, INDICES=20*round(seenGuides$Offset/20), FUN=mergeOffsetBins)
offsetBins <- do.call("rbind", byout)
offsetBins$bin <- as.integer(row.names(offsetBins))

seenGuides$yorfFwd <- grepl("[+]", seenGuides$YorfLoc)
seenGuides$targetFwd <- grepl("[+]", seenGuides$TargetLoc)

seenGuidesSame <- seenGuides[seenGuides$yorfFwd == seenGuides$targetFwd,]
seenGuidesOpp <- seenGuides[seenGuides$yorfFwd != seenGuides$targetFwd,]

byout <- by(seenGuidesSame, INDICES=20*round(seenGuidesSame$Offset/20), FUN=mergeOffsetBins)
offsetBinsSame <- do.call("rbind", byout)
offsetBinsSame$bin <- as.integer(row.names(offsetBinsSame))

byout <- by(seenGuidesOpp, INDICES=20*round(seenGuidesOpp$Offset/20), FUN=mergeOffsetBins)
offsetBinsOpp <- do.call("rbind", byout)
offsetBinsOpp$bin <- as.integer(row.names(offsetBinsOpp))


offsetBins <- offsetBins[offsetBins$n > 100,]
offsetBinsOpp <- offsetBinsOpp[offsetBinsOpp$n > 50,]
offsetBinsSame <- offsetBinsSame[offsetBinsSame$n > 50,]

pdf(sprintf("%s/guide-offset.pdf", figuredir), useDingbats=FALSE, width=6, height=4)
piyg <- brewer.pal(5, "PiYG")

plot(offsetBins$bin, offsetBins$meanFit,
     type="b", lwd=2, pch=20, cex=0.5, col="#666666",
     xlim=c(-230,30), xaxp=c(-200,0,4),
     ylim=c(-0.65, 0), yaxp=c(-0.6, 0, 3),
     xlab="Offset guide to TSS (bp)", ylab="Average Fitness")
lines(offsetBinsSame$bin, offsetBinsSame$meanFit,
     type="b", lwd=2, pch=20, cex=0.5, col=piyg[[1]])
lines(offsetBinsOpp$bin, offsetBinsOpp$meanFit,
      type="b", lwd=2, pch=20, cex=0.5, col=piyg[[5]])
legend(x="bottomleft", bty="n", lwd=2,
       legend=c("All Guides", "Fwd (NGG)", "Rev (CCN)"),
       col=c("#666666", piyg[[1]], piyg[[5]]),
       text.col=c("#666666", piyg[[1]], piyg[[5]]))
dev.off()
