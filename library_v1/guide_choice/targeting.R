options(stringsAsFactors=FALSE)

## Configurable locations for files

## Directories for this script
targetingDir <- "./"
targetingWorkDir <- sprintf("%s/targeting_work/", targetingDir)
targetingOutputDir <- targetingDir

dir.create(targetingDir)
dir.create(targetingWorkDir)
dir.create(targetingOutputDir)

## Files from yeast_grna
yeastGrnaDir <- "./"
yeastGrnaWorkDir <- sprintf("%s/work/", yeastGrnaDir)
yeastGrnaOutputDir <- sprintf("%s/output/", yeastGrnaDir)

targetChoiceFile <- sprintf("%s/target-choices.txt", yeastGrnaOutputDir)
targetOligoFile <- sprintf("%s/target-oligos.txt", yeastGrnaOutputDir)
allScoreFile <- sprintf("%s/all-scores.txt", yeastGrnaWorkDir)

sgdFile <- sprintf("%s/SGD_features.tab", targetingWorkDir)

log <- file(sprintf("%s/targeting-log.txt", targetingOutputDir), "w")

## Dubious ORF list from SGD
if (!file.exists(sgdFile)) {
  sgd <- download.file('https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab', destfile=sgdFile)
}
sgd <- read.delim(sgdFile, header=FALSE, quote="",
                  col.names=c("sgdid", "type", "qual", "name", "gene", "alias",
                              "parent", "sgdid2", "chrom", "start", "end",
                              "strand", "genpos", "cver", "sver", "desc"))

dubious <- sgd[sgd$qual == "Dubious","name"]
dubious <- dubious[order(dubious)]
write.table(dubious, file=sprintf("%s/dubious.txt", targetingWorkDir), quote=FALSE, row.names=FALSE, col.names=FALSE)
cat(sprintf("Identified %d dubious ORFs\n", length(dubious)), file=log)

## Target choices table output from yeast_grna
targetChoice <- read.delim(targetChoiceFile)
targetChoice$Guide <- sprintf("%s_%02d", targetChoice$Yorf, targetChoice$GuideNo)
cat(sprintf("Read %d target choices\n", nrow(targetChoice)), file=log)
cat(sprintf("  %d unique locations in target choices\n", length(unique(targetChoice$TargetLoc))), file=log)
cat(sprintf("  %d unique target sequences\n", length(unique(targetChoice$TargetSeq))), file=log)

## Target oligos table output from yeast_grna
targetOligos <- read.delim(targetOligoFile)
cat(sprintf("Read %d target oligos\n", nrow(targetOligos)), file=log)
cat(sprintf("  %d unique locations in oligos\n", length(unique(targetOligos$TargetLoc))), file=log)
cat(sprintf("  %d unique target sequences\n", length(unique(targetOligos$Oligo))), file=log)

## Extract genes with "clean" targeting -- unique sequence,
##   location is specific (no other targets), and based on TSS
mergeChoices <- function(choices) {
    data.frame(allUnique = all(as.logical(choices$Unique)),
               allSpecific = all(as.logical(choices$Specific)),
               nGuide = nrow(choices),
               hasTSS = any(grepl("TSS", choices$YorfStart)))
}
byout <- by(data=targetChoice, INDICES=targetChoice$Yorf, FUN=mergeChoices)
yorfChoice <- do.call("rbind", byout)
cat(sprintf("Merged %s genes with guide choices\n", nrow(yorfChoice)), file=log)

cleanYorfs <- row.names(yorfChoice[yorfChoice$allUnique & yorfChoice$allSpecific & yorfChoice$hasTSS,])
cat(sprintf("  %s genes with \"clean\" (unique, specific, TSS-based) guide choices\n", length(cleanYorfs)), file=log)

write.table(cleanYorfs, file="yorf-clean-target.txt",
            quote=FALSE, row.names=FALSE, col.names=FALSE)

## Target scoring table from yeast_grna
allScores <- read.delim(allScoreFile, header=FALSE,
                        col.names=c("Yorf", "TargetSeq", "TargetLoc", "Offset", "Start", "StartLoc"))
targetAndGene <- allScores[,c("TargetLoc", "Yorf", "Offset")]
targetAndGene$StartType <- sub("=.*", "", allScores$Start)
cat(sprintf("Read %d target-to-gene associations\n", nrow(targetAndGene)), file=log)

## Find all (unique) locations used for oligos from target choices table
## Same location can be chosen for two different genes sometimes
chosenTargetLocs <- unique(targetChoice$TargetLoc)
cat(sprintf("%d unique chosen targets\n", length(chosenTargetLocs)), file=log)

## Ensure we have scoring table entries for every oligo location
missing <- chosenTargetLocs[!(chosenTargetLocs %in% targetAndGene$TargetLoc)]
if (length(missing) > 0) {
    stop(sprintf("Missing %d chosen targets in target-and-gene: %s\n",
                 length(missing), toString(missing)))
} else {
    cat(sprintf("All chosen targets in target-and-gene map\n", length(chosenTargetLocs)), file=log)
}

chosenTargetAndGene <- targetAndGene[targetAndGene$TargetLoc %in% chosenTargetLocs,]
cat(sprintf("%d target-and-gene mappings for chosen targets\n", nrow(chosenTargetAndGene)), file=log)

## Discard targets on mitochondrial chromosome
## Discard targeting information against dubious ORFs -- excludes from scoring for "preferred" target
goodChosenTargetAndGene <- chosenTargetAndGene[!grepl("chrM", chosenTargetAndGene$TargetLoc)
                                               & !(chosenTargetAndGene$Yorf %in% dubious),]
cat(sprintf("%d good mappings for %d good locations after filtering (no mito, no dubious)\n",
            nrow(goodChosenTargetAndGene), length(unique(goodChosenTargetAndGene$TargetLoc))),
    file=log)

## Prioritize targeted genes for each target location
## Get all target genes for target location, sort by proximity to optimal upstream distance
mergeGoodTargets <- function(locTargets) {
    locTargets$distance <- abs(locTargets$Offset - ifelse(locTargets$StartType == "TSS", -50, -120))
    targetOrder <- order(locTargets$distance)
    
    targetLoc <- locTargets[targetOrder[[1]], "TargetLoc"]
    locChoice <- targetChoice[targetChoice$TargetLoc == targetLoc,]

    ## cat(sprintf("TargetLoc is %s. Targets = \n", targetLoc))
    ## print(locTargets)
    ## cat("Choice = \n")
    ## print(locChoice)
    
    ## Pick guide name matching highest-priority YORF
    ## Possible that we didn't pick a guide for a YORF even when it's the most likely target
    guide <- NA
    guidesource <- NA
    guide1 <- locChoice[locChoice$Yorf == locTargets[targetOrder[[1]], "Yorf"], "Guide"]
    if (length(guide1) > 0) {
        guide <- guide1[[1]]
        guidesource <- "Yorf1"
        ## cat(sprintf("Picked guide 1: %s\n", guide))
    } else if (length(targetOrder) >= 2) {
        guide2 <- locChoice[locChoice$Yorf == locTargets[targetOrder[[2]], "Yorf"], "Guide"]
        if (length(guide2) > 0) {
            guide <- guide2[[1]]
            guidesource <- "Yorf2"
            ## cat(sprintf("Picked guide 2: %s\n", guide))
        } else if (length(targetOrder) >= 3) {
            guide3 <- locChoice[locChoice$Yorf == locTargets[targetOrder[[3]], "Yorf"], "Guide"]
            if (length(guide3) > 0) {
                guide <- guide3
                guidesource <- "Yorf3"
                ## cat(sprintf("Picked guide 3: %s\n", guide))
            }
        }
    }
    
    if (is.na(guide)) {
        if (length(locChoice) > 0) {
            guide <- locChoice$Guide[[1]]
            guidesource <- "Choice"
            ## cat(sprintf("Picked guide from choice: %s\n", guide))
        } else {
            warning(sprintf("No guide name for %s\n", targetLoc))
        }
    }

    data.frame(TargetLoc = targetLoc,
               Yorf1 = locTargets[targetOrder[[1]], "Yorf"],
               Offset1 = locTargets[targetOrder[[1]], "Offset"],
               StartType1 = locTargets[targetOrder[[1]], "StartType"],
               Yorf2 = if (length(targetOrder) > 1) { 
                           locTargets[targetOrder[[2]], "Yorf"]
                       } else { NA },
               Offset2 = if (length(targetOrder) > 1) {
                             locTargets[targetOrder[[2]], "Offset"]
                         } else { NA },
               StartType2 = if (length(targetOrder) > 1) {
                             locTargets[targetOrder[[2]], "StartType"]
                         } else { NA },
               Yorf3 = if (length(targetOrder) > 2) {
                           locTargets[targetOrder[[3]], "Yorf"]
                       } else { NA },
               Offset3 = if (length(targetOrder) > 2) {
                             locTargets[targetOrder[[3]], "Offset"]
                         } else { NA },
               StartType3 = if (length(targetOrder) > 2) {
                             locTargets[targetOrder[[3]], "StartType"]
                         } else { NA },
               Yorfs = paste(locTargets[targetOrder,"Yorf"], collapse=","),
               PrimaryGuide = guide,
               PrimaryGuideSource = guidesource)
}

byout <- by(goodChosenTargetAndGene,
            INDICES=goodChosenTargetAndGene$TargetLoc,
            FUN=mergeGoodTargets)
locGoodTargets <- do.call("rbind", byout)

write.table(locGoodTargets,
            file=sprintf("%s/loc-targets.txt", targetingWorkDir),
            quote=FALSE, sep="\t", row.names=FALSE)

## Link named guides from choices with targeting information based on location
guides <- targetChoice[,c("Guide", "TargetLoc")]
guides <- cbind.data.frame(guides,
                           locGoodTargets[match(guides$TargetLoc, locGoodTargets$TargetLoc),-1])
guides$PrimaryGuide <- ifelse(is.na(guides$PrimaryGuide), guides$Guide, guides$PrimaryGuide)
guides <- guides[order(guides$Guide),]
write.table(guides,
            file=sprintf("%s/guide-good-targets.txt", targetingOutputDir),
            quote=FALSE, sep="\t", row.names=FALSE)

## Reduce guides down to unique target locations
##   First pick the primary guide name for each entry
##   Remaining cases are all dubious, no defined targets (Yorf1 is NA)
locations <- guides[guides$Guide == guides$PrimaryGuide,]
cat(sprintf("Reduced %d named guides to %d with primary guide name\n",
            nrow(guides), nrow(locations)), file=log)
redundantLocations <- locations[duplicated(locations$TargetLoc) & !is.na(locations$Yorf1),]
stopifnot(nrow(redundantLocations) == 0)
locations <- locations[!duplicated(locations$TargetLoc),]
cat(sprintf("Reduced to %d unique target locations (discarding locations with all dubious targets)\n", nrow(locations)), file=log)
write.table(locations,
            file=sprintf("%s/location-good-targets.txt", targetingOutputDir),
            quote=FALSE, sep="\t", row.names=FALSE)

## Link named guides from oligos
locations$Oligo <- targetOligos[match(locations$Guide, targetOligos$Yorf_GNo),"Oligo"]
redundantSequences <- locations[duplicated(locations$Oligo),]
sequences <- locations[!duplicated(locations$Oligo),]
sequences$Guide <- ifelse(sequences$Oligo %in% redundantSequences$Oligo,
                          paste0(sequences$Guide, "_DEGEN"),
                          sequences$Guide)
redundantSequences$DegenGuide <- sequences[match(redundantSequences$Oligo, sequences$Oligo), "Guide"]
cat(sprintf("Reduced %d target locations to %d unique sequences\n", nrow(locations), nrow(sequences)), file=log)

write.table(sequences,
            file=sprintf("%s/sequence-good-targets.txt", targetingOutputDir),
            quote=FALSE, sep="\t", row.names=FALSE)
write.table(redundantSequences,
            file=sprintf("%s/sequence-redundant-targets.txt", targetingOutputDir),
            quote=FALSE, sep="\t", row.names=FALSE)
reffa <- file(sprintf("%s/sequence-oligos.fa", targetingOutputDir), "w")
cat(do.call("paste0", as.list(sprintf(">%s\n%s\n", sequences$Guide, sequences$Oligo))), file=reffa)
