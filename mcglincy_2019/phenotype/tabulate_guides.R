options(stringsAsFactors=FALSE)

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

## Guide phenotype information
guideDeseqFile <- sprintf("%s/nizm005-guide-deseq.csv", figuredir)
guideres <- read.csv(guideDeseqFile, row.names=1)

guideTable <- guideToYorfs[,c("Guide", "TargetLoc", "Yorf1", "Offset1", "StartType1", "Yorfs", "Oligo")]
guideTable$Yorfs <- gsub(",", ";", guideTable$Yorfs)
guideTable$lfcMean <- round(guideres[match(guideTable$Guide, row.names(guideres)), "lfcMean"], digits=2)
guideTable$lfcSE <- round(guideres[match(guideTable$Guide, row.names(guideres)), "lfcSE"], digits=2)
guideTable$nbc <- round(guideres[match(guideTable$Guide, row.names(guideres)), "nbc"], digits=2)
guideTable$Gene1 <- sgd[match(guideTable$Yorf1, sgd$name), "gene"]
guideTable$Gene1 <- ifelse(is.na(guideTable$Gene1), guideTable$Yorf1, guideTable$Gene1)

write.csv(guideTable, file=sprintf("%s/guide-table-nizm005-deseq.csv", figuredir),
          row.names=FALSE, quote=grep("Gene", colnames(guideTable)))

