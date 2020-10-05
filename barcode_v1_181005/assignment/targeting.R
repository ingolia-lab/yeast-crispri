options(stringsAsFactors=FALSE)

barcodeGuides <- read.delim("../grna-assign-barcode-grna-good.txt")
guideTargets <- read.delim("../../library_v1/guide_choice/guide-good-targets.txt")

barcodeTargets <- cbind.data.frame(barcodeGuides,
                                   guideTargets[match(barcodeGuides$guide, guideTargets$Guide),])
barcodeTargets <- barcodeTargets[!is.na(barcodeTargets$Guide),]

barcodeTargets <- barcodeTargets[,c("barcode", "PrimaryGuide", "TargetLoc",
                                    "Yorf1", "Offset1", "StartType1",
                                    "Yorf2", "Offset2", "StartType2",
                                    "Yorfs")]
barcodeTargets <- barcodeTargets[!is.na(barcodeTargets$Yorf1),]
barcodeTargets <- barcodeTargets[order(barcodeTargets$PrimaryGuide),]
write.csv(barcodeTargets, "../barcode-targets-good.txt", row.names=FALSE)
