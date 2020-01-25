options(stringsAsFactors=FALSE)

library(png)
library("RColorBrewer")

guidedir <- Sys.getenv("GUIDEDIR")
barcodingdir <- Sys.getenv("BARCODINGDIR")
workdir <- Sys.getenv("WORKDIR")
figuredir <- Sys.getenv("FIGUREDIR")

log <- file(sprintf("%s/analyze_guides_out.txt", figuredir), "w")

## Essential gene list
essFile <- sprintf("%s/Essential.txt", workdir)
if (!file.exists(essFile)) {
  essential <- download.file('http://www-sequence.stanford.edu/group/yeast_deletion_project/Essential_ORFs.txt', destfile=essFile)
}
essential <- read.delim(essFile)
essentialYorfs <- sub(" *$", "", essential$ORF_name)

## New accessibility data
## Oberbeckmann et al., Genome Res 2019
## Occupancy by DNA Methylation (in vitro) = ODM
odmFile <- sprintf("%s/sequence-good-odms.txt", workdir)
odm <- read.delim(odmFile, header=FALSE, stringsAsFactors=FALSE)

## Information on guide targeting
guideToYorfFile <- sprintf("%s/sequence-good-targets.txt", guidedir)
guideToYorfs <- read.delim(guideToYorfFile)

guideChoiceFile <- sprintf("%s/guide_choice/output/target-choices.txt", guidedir)
guideChoice <- read.delim(guideChoiceFile)

## Information on fitness effects
guideDeseqFile <- sprintf("%s/nizm005-guide-deseq.csv", figuredir)
guideres <- read.csv(guideDeseqFile, row.names=1)

## Genes with clean targeting
cleanYorfFile <- sprintf("%s/yorf-clean-target.txt", guidedir)
cleanYorfs <- read.delim(cleanYorfFile, header=FALSE)$V1

## Empirical cutoff for "sick" guides, 5th percentile of negatives
bcDeseqFile <- sprintf("%s/nizm005-barcode-deseq.csv", workdir)
bcres <- read.csv(bcDeseqFile, row.names=1)
emptyres <- bcres[grepl("Empty", bcres$guide),]
sickLimit <- quantile(emptyres$log2FoldChange, 0.05)
cat(sprintf("Sick limit %0.2f\n", sickLimit), file=log)

## High-confidence training set
## Multiple barcodes, clean targeting, and essential
efficacyres <- guideres[guideres$nbc > 1
                        & guideres$Yorf1 %in% cleanYorfs
                        & guideres$Yorf1 %in% essentialYorfs,]
efficacy <- data.frame(row.names = efficacyres$TargetLoc,
                       lfcMean = efficacyres$lfcMean)

cat(sprintf("%d guides in efficacy analysis\n", nrow(efficacy)), file=log)

cat(sprintf("%d essential genes with clean promoters\n", sum(cleanYorfs %in% essentialYorfs)), file=log) 

write.table(efficacyres,
            file=sprintf("%s/efficacy.training.txt", workdir),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

## Collect data for efficacy modeling
efficacy$active <- efficacy$lfcMean < sickLimit
efficacy$ATAC <- guideChoice[match(row.names(efficacy), guideChoice$TargetLoc), "ATAC"]
efficacy$Offset <- guideChoice[match(row.names(efficacy), guideChoice$TargetLoc), "Offset"]
efficacy$SameStrand <- (grepl("[0-9]W", efficacyres$Yorf1) ==
                        grepl("[+]", efficacyres$TargetLoc))
efficacy$ODM <- odm[match(row.names(efficacy), odm$V1), "V5"]/100.0
efficacy$ODMBin <- cut(efficacy$ODM, breaks=5, labels=FALSE)
print(table(efficacy$ODMBin))

for (opos in seq(21,40)) {
    efficacy[,sprintf("nt%02d", opos - 20)] <- as.factor(substr(efficacyres$Oligo, opos, opos))
}

## LOESS model of fitness effect vs offset from TSS
offsetModel <- loess(lfcMean ~ Offset, data = efficacy, span = 0.25)
efficacy$OffsetPred <- predict(offsetModel)
offsetout <- data.frame(Offset = seq(-220, 20))
offsetout$Score <- predict(offsetModel, newdata=offsetout)
write.table(offsetout, file=sprintf("%s/offset.predictions.txt", figuredir),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

pdf(sprintf("%s/fitness_by_offset.pdf", figuredir),
    useDingbats=FALSE, width=6, height=4)
brbg <- brewer.pal(9, "BrBG")
odmpal <- c(brbg[[8]], brbg[[7]], "#d0d0d0", brbg[[3]], brbg[[2]])
plot(efficacy$Offset, efficacy$lfcMean,
     pch=20, cex=0.5, col=odmpal[efficacy$ODMBin],
     xlim=c(-230,30), xaxp=c(-200,0,4),
     ylim=c(-1.1, 0.3), yaxp=c(-1, 0, 2),
     xlab="Offset guide to TSS (bp)", ylab="Fitness")     
lines(seq(-200,20), predict(offsetModel, seq(-200,20)),
      lwd=2)

xl <- 1.25
yb <- 1
xr <- 1.5
yt <- 2

par(mar=c(5.1,0.5,4.1,0.5))
plot(NA,type="n",ann=FALSE,xlim=c(1,2),ylim=c(1,2),xaxt="n",yaxt="n",bty="n")
rect(xl,head(seq(yb,yt,(yt-yb)/5),-1),
     xr,tail(seq(yb,yt,(yt-yb)/5),-1),
     col=odmpal, lwd=0)
mtext(seq(0, 1, 0.2),side=2,at=seq(yb,yt,(yt-yb)/5),las=2,cex=0.7)
dev.off()

## Strand dependence of fitness effect
offsetModelSame <- loess(lfcMean ~ Offset, data = efficacy[efficacy$SameStrand,], span = 0.25)
offsetModelOppo <- loess(lfcMean ~ Offset, data = efficacy[!efficacy$SameStrand,], span = 0.25)
efficacy$OffsetStrandPred <- ifelse(efficacy$SameStrand,
                                    predict(offsetModelSame, newdata=efficacy),
                                    predict(offsetModelOppo, newdata=efficacy))
pdf(sprintf("%s/fitness_by_offset_strand.pdf", figuredir),
    useDingbats=FALSE, width=6, height=4)
piyg <- brewer.pal(5, "PiYG")
plot(seq(-200,20), predict(offsetModel, seq(-200,20)),
     type="l", lwd=2, lty=2, col="#666666",
     xlim=c(-230,30), xaxp=c(-200,0,4),
     ylim=c(-1.1, 0.3), yaxp=c(-1, 0, 2),
     xlab="Offset guide to TSS (bp)", ylab="Fitness")     
lines(seq(-200,20), predict(offsetModelSame, seq(-200,20)),
      lwd=2, col=piyg[[1]])
lines(seq(-200,20), predict(offsetModelOppo, seq(-200,20)),
      lwd=2, col=piyg[[5]])
legend(x="bottomleft", bty="n", lwd=2,
       legend=c("All Guides", "Fwd (NGG)", "Rev (CCN)"),
       col=c("#666666", piyg[[1]], piyg[[5]]),
       text.col=c("#666666", piyg[[1]], piyg[[5]]))
dev.off()

## Build model of efficacy
efficacyModel <- glm(active ~ OffsetPred + ODM
                     + nt01 + nt02 + nt03 + nt04 + nt05 + nt06 + nt07 + nt08 + nt09 + nt10
                     + nt11 + nt12 + nt13 + nt14 + nt15 + nt16 + nt17 + nt18 + nt19 + nt20,
                     family = binomial(link = 'logit'),
                     data = efficacy)
write.table(coefficients(efficacyModel),
            file=sprintf("%s/model.coefficients.txt", figuredir),
            quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)

## Rescore all guides according to efficacy model
rescore <- guideToYorfs
rescore$Offset <- ifelse(rescore$StartType1 == "CDS",
                         rescore$Offset1 + 30,
                         rescore$Offset1)
rescore$OffsetPred <- predict(offsetModel, newdata = rescore)
rescore$ODM <- odm[match(rescore$TargetLoc, odm$V1), "V5"]/100.0
for (opos in seq(21,40)) {
    rescore[,sprintf("nt%02d", opos - 20)] <- as.factor(substr(rescore$Oligo, opos, opos))
}
rescore$score <- predict(efficacyModel, newdata=rescore)
rescore$response <- predict(efficacyModel, newdata=rescore, type="response")
write.table(rescore[,c("TargetLoc", "score", "response")],
            file=sprintf("%s/rescore.txt", figuredir),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
write.table(rescore, file=sprintf("%s/rescore.verbose.txt", workdir),
            quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

## Testing set of guides against essential genes not included in training
##   Divergent promoters or just one barcode
testres <- guideres[guideres$Yorf1 %in% essentialYorfs
                    & !(guideres$TargetLoc %in% efficacyres$TargetLoc),]
testres$score <- rescore[match(testres$TargetLoc, rescore$TargetLoc),"score"]
testres$scorebin <- pmin(3, pmax(-3, round(testres$score)))

cat(sprintf("%d guides holdout efficacy test\n", nrow(testres)), file=log)

cat(sprintf("  %d guides active\n", sum(testres$lfcMean < sickLimit)), file=log)

testcdf <- ecdf(testres$lfcMean)
zerocdf <- ecdf(testres[testres$scorebin == 0,]$lfcMean)
breaks <- seq(-1.35, 0.45, 0.05)
puor <- brewer.pal(9, "PuOr")

pdf(sprintf("%s/efficacy_score.pdf", figuredir),
    useDingbats=FALSE, width=6, height=4)
plot(breaks, zerocdf(breaks),
     type="l", lwd=2, col="grey",
     xlim=c(-1.4, 0.6), ylim=c(0,1), yaxt="n",
     xlab="lg(s)", ylab="Fraction Guides",
     xaxp=c(-1.5,0.5,4))

for (bin in c(-3,-2,-1,1,2,3)) {
    bincdf <- ecdf(testres[testres$scorebin == bin,]$lfcMean)
    col <- ifelse(bin < 0, puor[[6 - bin]], puor[[4 - bin]])
    lines(breaks, bincdf(breaks), col=col, lwd=2)
}

legend(x="topleft", bty="n",
       legend = seq(3,-3),
       lwd=2, col=c(puor[c(1,2,3)], "grey", puor[c(7,8,9)]),
       text.col=c(puor[c(1,2,3)], "grey", puor[c(7,8,9)]))

dev.off()

## k-folds cross validation of efficacy models
##   Deterministic set assignments
k <- 10
efficacy$kset <- 1 + (seq(1,nrow(efficacy)) %% k)
efficacy$fittedFull <- NA
efficacy$fittedFullStrand <- NA
efficacy$fittedNoODM <- NA
efficacy$fittedNoSeq <- NA
efficacy$fittedOffset <- NA
efficacy$fittedOffsetStrand <- NA
for (i in seq(1, k)) {
    train <- efficacy[efficacy$kset != i,]
    test <- efficacy[efficacy$kset == i,]

    trainOffset <- loess(lfcMean ~ Offset, data=train, span=0.25)
    trainOffsetSame <- loess(lfcMean ~ Offset, data=train[train$SameStrand,], span=0.25)
    trainOffsetOppo <- loess(lfcMean ~ Offset, data=train[!train$SameStrand,], span=0.25)

    train$OffsetPred <- predict(trainOffset, newdata=train)
    train$OffsetStrandPred <- ifelse(train$SameStrand,
                                     predict(trainOffsetSame, newdata=train),
                                     predict(trainOffsetOppo, newdata=train))

    modelFull <- glm(active ~ OffsetPred + ODM
                     + nt01 + nt02 + nt03 + nt04 + nt05 + nt06 + nt07 + nt08 + nt09 + nt10
                     + nt11 + nt12 + nt13 + nt14 + nt15 + nt16 + nt17 + nt18 + nt19 + nt20,
                     family = binomial(link = 'logit'),
                     data = train)
    modelFullStrand <- glm(active ~ OffsetStrandPred + ODM
                           + nt01 + nt02 + nt03 + nt04 + nt05 + nt06 + nt07 + nt08 + nt09 + nt10
                           + nt11 + nt12 + nt13 + nt14 + nt15 + nt16 + nt17 + nt18 + nt19 + nt20,
                           family = binomial(link = 'logit'),
                           data = train)
    modelNoODM <- glm(active ~ OffsetPred
                      + nt01 + nt02 + nt03 + nt04 + nt05 + nt06 + nt07 + nt08 + nt09 + nt10
                      + nt11 + nt12 + nt13 + nt14 + nt15 + nt16 + nt17 + nt18 + nt19 + nt20,
                      family = binomial(link = 'logit'),
                      data = train)
    modelNoSeq <- glm(active ~ OffsetPred + ODM,
                      family = binomial(link = 'logit'),
                      data = train)
    modelOffset <- glm(active ~ OffsetPred, family = binomial(link = 'logit'), data = train)
    modelOffsetStrand <- glm(active ~ OffsetStrandPred, family = binomial(link = 'logit'), data = train)

    test$OffsetPred <- predict(trainOffset, newdata=test)
    test$OffsetStrandPred <- ifelse(test$SameStrand,
                                    predict(trainOffsetSame, newdata=test),
                                    predict(trainOffsetOppo, newdata=test))

    fittedFull <- predict(modelFull, newdata = test, type = 'response')
    fittedFullStrand <- predict(modelFullStrand, newdata = test, type = 'response')
    fittedNoODM <- predict(modelNoODM, newdata = test, type = 'response')
    fittedNoSeq <- predict(modelNoSeq, newdata = test, type = 'response')
    fittedOffset <- predict(modelOffset, newdata = test, type = 'response')
    fittedOffsetStrand <- predict(modelOffsetStrand, newdata = test, type = 'response')

    efficacy$fittedFull <- ifelse(is.na(match(row.names(efficacy), row.names(test))),
                                  efficacy$fittedFull,
                                  fittedFull[match(row.names(efficacy), row.names(test))])
    efficacy$fittedFullStrand <- ifelse(is.na(match(row.names(efficacy), row.names(test))),
                                        efficacy$fittedFullStrand,
                                        fittedFullStrand[match(row.names(efficacy), row.names(test))])
    efficacy$fittedNoODM <- ifelse(is.na(match(row.names(efficacy), row.names(test))),
                                  efficacy$fittedNoODM,
                                  fittedNoODM[match(row.names(efficacy), row.names(test))])
    efficacy$fittedNoSeq <- ifelse(is.na(match(row.names(efficacy), row.names(test))),
                                  efficacy$fittedNoSeq,
                                  fittedNoSeq[match(row.names(efficacy), row.names(test))])
    efficacy$fittedOffset <- ifelse(is.na(match(row.names(efficacy), row.names(test))),
                                    efficacy$fittedOffset,
                                    fittedOffset[match(row.names(efficacy), row.names(test))])
    efficacy$fittedOffsetStrand <- ifelse(is.na(match(row.names(efficacy), row.names(test))),
                                          efficacy$fittedOffsetStrand,
                                          fittedOffsetStrand[match(row.names(efficacy), row.names(test))])
}

roc <- function(score, active) {
    x <- data.frame(score = score, active = active)
    byscore <- x[order(x$score, decreasing=TRUE),]
    data.frame(tpr = cumsum(byscore$active) / sum(byscore$active),
               fpr = cumsum(!byscore$active) / sum(!byscore$active),
               fitted = byscore$score)
}

rocauc <- function(roc) {
    dx <- roc[seq(2,nrow(roc)),]$fpr - roc[seq(1,nrow(roc)-1),]$fpr
    y <- roc[seq(2,nrow(roc)),]$tpr
    sum(dx * y)
}

rocFull <- roc(efficacy$fittedFull, efficacy$active)
rocFullStrand <- roc(efficacy$fittedFullStrand, efficacy$active)
rocNoODM <- roc(efficacy$fittedNoODM, efficacy$active)
rocNoSeq <- roc(efficacy$fittedNoSeq, efficacy$active)
rocOffset <- roc(efficacy$fittedOffset, efficacy$active)
rocOffsetStrand <- roc(efficacy$fittedOffsetStrand, efficacy$active)

cat(sprintf("AUC Full: %0.3f\n", rocauc(rocFull)), file=log)
cat(sprintf("AUC FullStrand: %0.3f\n", rocauc(rocFullStrand)), file=log)
cat(sprintf("AUC NoODM: %0.3f\n", rocauc(rocNoODM)), file=log)
cat(sprintf("AUC NoSeq: %0.3f\n", rocauc(rocNoSeq)), file=log)
cat(sprintf("AUC Offset: %0.3f\n", rocauc(rocOffset)), file=log)
cat(sprintf("AUC OffsetStrand: %0.3f\n", rocauc(rocOffsetStrand)), file=log)

pdf(sprintf("%s/efficacy_roc_all.pdf", figuredir), useDingbats=FALSE, height = 7, width = 7)
pal <- brewer.pal(6, "Dark2")
plot(rocFull$fpr, rocFull$tpr,
     type="s", lwd=2,
     xlab = "False Positive Rate", ylab = "True Positive Rate")
lines(rocFullStrand$fpr, rocFullStrand$tpr, type="s", col=pal[[4]])
lines(rocNoODM$fpr, rocNoODM$tpr, type="s", col=pal[[2]])
lines(rocNoSeq$fpr, rocNoODM$tpr, type="s", col=pal[[3]])
lines(rocOffset$fpr, rocOffset$tpr, type="s", col=pal[[1]])
lines(rocOffsetStrand$fpr, rocOffsetStrand$tpr, type="s", col=pal[[6]])
legend(x="bottomright", bty="n", lwd=3,
       legend=c(sprintf("Full (AUC %0.2f)", rocauc(rocFull)),
                sprintf("Full w/ Strand (AUC %0.2f)", rocauc(rocFullStrand)),
                sprintf("No ODM (AUC %0.2f)", rocauc(rocNoODM)),
                sprintf("No Seq (AUC %0.2f)", rocauc(rocNoSeq)),
                sprintf("Offset only (AUC %0.2f)", rocauc(rocOffset)),
                sprintf("Offset w/ strand (AUC %0.2f)", rocauc(rocOffsetStrand))),
       col=c("black", pal[c(4,2,3,1,6)]),
       text.col=c("black", pal[c(4,2,3,1,6)]))
dev.off()

pdf(sprintf("%s/efficacy_roc.pdf", figuredir), useDingbats=FALSE, height = 4, width = 4)
pal <- brewer.pal(6, "Dark2")
plot(rocFull$fpr, rocFull$tpr,
     type="s", lwd=2,
     xlab = "False Positive Rate", ylab = "True Positive Rate")
lines(rocOffset$fpr, rocOffset$tpr, type="s", col=pal[[1]], lwd=2)
legend(x="bottomright", bty="n", lwd=3,
       legend=c(sprintf("Full (AUC %0.2f)", rocauc(rocFull)),
                sprintf("Offset only (AUC %0.2f)", rocauc(rocOffset))),
       col=c("black", pal[1]),
       text.col=c("black", pal[1]))
dev.off()

## Model selection tests
sink(file = log)
anova(glm(active ~ OffsetPred + OffsetStrandPred,
          family = binomial(link = 'logit'), data = efficacy), test="Chisq")
anova(glm(active ~ OffsetStrandPred + OffsetPred,
          family = binomial(link = 'logit'), data = efficacy), test="Chisq")

anova(glm(active ~ OffsetStrandPred + ATAC + ODM,
          family = binomial(link = 'logit'), data = efficacy), test="Chisq")
anova(glm(active ~ OffsetStrandPred + ODM + ATAC,
          family = binomial(link = 'logit'), data = efficacy), test="Chisq")

anova(glm(active ~ OffsetPred + ODM +
              nt01 + nt02 + nt03 + nt04 + nt05 + nt06 + nt07 + nt08 + nt09 + nt10 +
              nt11 + nt12 + nt13 + nt14 + nt15 + nt16 + nt17 + nt18 + nt19 + nt20 +
              OffsetStrandPred,
          family = binomial(link = 'logit'), data = efficacy), test="Chisq")

anova(glm(active ~ OffsetStrandPred +
              nt01 + nt02 + nt03 + nt04 + nt05 + nt06 + nt07 + nt08 + nt09 + nt10 +
              nt11 + nt12 + nt13 + nt14 + nt15 + nt16 + nt17 + nt18 + nt19 + nt20 +
              ODM,
          family = binomial(link = 'logit'), data = efficacy), test="Chisq")
sink()
