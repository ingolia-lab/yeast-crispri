options(stringsAsFactors=FALSE)
library("RColorBrewer")

guidedir <- Sys.getenv("GUIDEDIR")
figuredir <- Sys.getenv("FIGUREDIR")

log <- file(sprintf("%s/guide_stats_out.txt", figuredir), "w")

targets <- read.delim(sprintf("%s/sequence-good-targets.txt", guidedir))

cat(sprintf("targets = %d\n", nrow(targets)), file=log)
cat(sprintf("%d non-specific, unique guides\n", sum(!is.na(targets$Yorf2))), file=log)
cat(sprintf("%d specific, unique guides\n", sum(is.na(targets$Yorf2))), file=log)

uniqueSpecificCount <- table(targets[is.na(targets$Yorf2), "Yorf1"])
uniqueAssignedCount <- table(targets$Yorf1)

sink(log)
table(uniqueSpecificCount)
table(uniqueAssignedCount)
sink()

pdf(sprintf("%s/unique_specific_guides.pdf", figuredir), useDingbats=FALSE, width=4, height=4)
guideCount <- data.frame()
for (nguide in seq(1,11)) {
    guideCount <- rbind.data.frame(guideCount,
                                   data.frame(nguide = nguide,
                                              specificGene = sum(uniqueSpecificCount >= nguide),
                                              assignedGene = sum(uniqueAssignedCount >= nguide)))
}
sink(log)
print(guideCount)
sink()

unambigColors <- brewer.pal(3, "Greys")
assignColors <- brewer.pal(3, "Greens")

plot(x = guideCount$nguide - 0.5,
     y = guideCount$assignedGene,
     type="s", lwd=2, col=assignColors[3],
     xlim=c(0.5, 10.5), ylim=c(0, 6500),
     axes=FALSE,
     xlab="Number of guides", ylab="Number of genes")
lines(x = guideCount$nguide - 0.5,
      y = guideCount$specificGene,
      type="s", lwd=2, col=unambigColors[3])
axis(side=1, at=seq(1,10), cex.axis=0.8, lwd=0, lwd.ticks=1)
axis(side=1, at=seq(0.5, 10.5), labels=FALSE, lwd.ticks=0)
axis(side=2, at=seq(0,6000,2000))
axis(side=2, at=seq(1000,5000,2000),
     labels=FALSE, tcl=-0.25, lwd=0, lwd.ticks=1)
legend(x = "bottomleft", bty="n",
       legend = c("Unambiguous", "Assigned"),
       col = c(unambigColors[3], assignColors[3]),
       pt.bg = c(unambigColors[1], assignColors[1]),
       text.col = c(unambigColors[3], assignColors[3]),
       pch=22, pt.cex=2)
dev.off()
