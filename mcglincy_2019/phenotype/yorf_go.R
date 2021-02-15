options(stringsAsFactors=FALSE)

workdir <- Sys.getenv("WORKDIR")

yorfs <- read.csv(sprintf("%s/nizm005-yorf-deseq.csv", workdir))

write.table(yorfs[grepl("^Y", yorfs$X), c("X", "lfcMin"),],
            file=sprintf("%s/yorf-lfcmin-all.txt", workdir),
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

write.table(yorfs[grepl("^Y", yorfs$X) & yorfs$essential, c("X", "lfcMin"),],
            file=sprintf("%s/yorf-lfcmin-ess.txt", workdir),
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

write.table(yorfs[grepl("^Y", yorfs$X) & yorfs$essential & yorfs$clean, c("X", "lfcMin"),],
            file=sprintf("%s/yorf-lfcmin-ess-clean.txt", workdir),
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

write.table(yorfs[grepl("^Y", yorfs$X) & !yorfs$essential, c("X", "lfcMin"),],
            file=sprintf("%s/yorf-lfcmin-via.txt", workdir),
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

write.table(yorfs[grepl("^Y", yorfs$X) & !yorfs$essential & yorfs$clean, c("X", "lfcMin"),],
            file=sprintf("%s/yorf-lfcmin-via-clean.txt", workdir),
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
