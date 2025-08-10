#!/usr/bin/env Rscript

###################################################
## sequencing statistics
## a) total amount of reads
## b) mapping hits
###################################################

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("No Log.final.out files provided")
}
files <- args
samples <- sapply(files, function(x) gsub(".Log.final.out", "", basename(x)))

cols <- c("TotalReads", "Multihits", "TooManyLoci", "OtherUnmapped", "TooManyMM", "TooShort", "UniqueHits", "TotalSplice")
star <- c("Number of input reads", "% of reads mapped to multiple loci", "% of reads mapped to too many loci", "% of reads unmapped: other", "% of reads unmapped: too many mismatches", "% of reads unmapped: too short", "Uniquely mapped reads %", "Number of splices: Total")

# Read all logs
stats_list <- lapply(files, function(f) {
    df <- read.csv(f, sep="\t", as.is=TRUE)
    colnames(df) <- c("category", "value")
    df$category <- gsub("\\|", "", df$category)
    df$category <- trimws(df$category)
    df <- df[match(star, df$category), ]
    df$value <- gsub("%", "", df$value)
    as.numeric(df$value)
})

stats <- do.call(cbind, stats_list)
rownames(stats) <- cols
colnames(stats) <- samples

# Derived values
stats <- rbind(stats, 
               sumOther = stats[1,]/100*colSums(stats[2:5,]),
               short    = stats[1,]/100*stats[6,],
               unique   = stats[1,]/100*stats[7,])

# Output
dir.create("plots", showWarnings = FALSE)

expname <- "_Project_bc_mut"

pdf(paste0("plots/ReadStats", expname, ".pdf"), height=8, width=5)
par(mar=c(4,8,2,1))
bp <- barplot(stats[11:9,]/10^6, main='Library sizes and mapping', las=1,
              border=NA, horiz=TRUE, xlab='Reads in Mio', ylab='',
              xlim=c(0, max(stats[1,])*1.15/10^6), col=c("yellowgreen", "gray88", "gray48"))
legend("topright", rownames(stats)[11:9], bty='n', col=c("yellowgreen", "gray88", "gray48"), pch=15)
mtext(side=2, las=1, line=-2, at=bp, text=paste(round(stats[11,]/stats[1,]*100), "%", sep=""))
abline(v=30, col="grey30", lty=3)
dev.off()

write.table(stats[11:9,], "star_stat_summary.tsv", sep="\t", quote=FALSE)
