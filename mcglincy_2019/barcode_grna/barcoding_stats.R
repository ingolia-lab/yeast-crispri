options(stringsAsFactors=FALSE)

barcodingdir <- Sys.getenv("BARCODINGDIR")
figuredir <- Sys.getenv("FIGUREDIR")

log <- file(sprintf("%s/barcoding_stats_out.txt", figuredir), "w")

fatesTable <- read.delim(sprintf("%s/grna-assign-barcode-fates.txt", barcodingdir), row.names=1, header=FALSE)
fates <- fatesTable$V2
names(fates) <- row.names(fatesTable)

fatesTotal <- fates["Good"] + fates["Bad Match"] + fates["No Match"] + fates["Mixed"]
cat(sprintf("%d assigned barcodes\n", fatesTotal), file=log)
cat(sprintf("%d good assignments (%0.1f %%)\n", fates["Good"], 100 * fates["Good"] / fatesTotal), file=log)
cat(sprintf("%d erroneous (%0.1f %%)\n", fates["Bad Match"], 100 * fates["Bad Match"] / fatesTotal), file=log)
cat(sprintf("%d heterogeneous (%0.1f %%)\n", fates["Mixed"], 100 * fates["Mixed"] / fatesTotal), file=log)
cat(sprintf("%d mysterious (%0.1f %%)\n", fates["No Match"], 100 * fates["No Match"] / fatesTotal), file=log)

good <- read.delim(sprintf("%s/grna-assign-barcode-grna-good.txt", barcodingdir))

nonempty <- good[!grepl("Empty", good$guide),]
cat(sprintf("%d guides assigned to barcodes\n", length(unique(nonempty$guide))), file=log)

cat(sprintf("%d empty-guide barcodes\n", nrow(good[grepl("Empty", good$guide),])), file=log)

guideCount <- table(nonempty$guide)
cat(sprintf("%d guides with >1 barcodes\n", length(guideCount[guideCount > 1])), file=log)

pdf(sprintf("%s/barcodes_per_guide.pdf", figuredir), useDingbats=FALSE, width=5, height=4)
bcCount <- data.frame()
for (nbc in seq(1,11)) {
    bcCount <- rbind.data.frame(bcCount,
                                data.frame(nbc = nbc,
                                           nguides = sum(guideCount >= nbc)))
}

sink(log)
print(bcCount)
sink()

plot(x = bcCount$nbc - 0.5,
     y = bcCount$nguides,
     type="s", lwd=2, col="black",
     xlim=c(0.5, 10.5), ylim=c(0, 50000),
     axes=FALSE,
     xlab="Number of barcodes", ylab="Number of guides")
axis(side=1, at=seq(1,10), cex.axis=0.8, lwd=0, lwd.ticks=1)
axis(side=1, at=seq(0.5, 10.5), labels=FALSE, lwd.ticks=0)
axis(side=2, at=seq(0,40000,20000))
axis(side=2, at=seq(10000,50000,20000),
     labels=FALSE, tcl=-0.25)
dev.off()
