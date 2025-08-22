###################################################
## visualise statistics
###################################################

# Check/install required packages
list.of.packages <- c("tidyr", "dplyr", "data.table", "ggplot2", "gridExtra", "ggpubr", "cowplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  message("Installing missing R packages: ", paste(new.packages, collapse = ", "))
  install.packages(new.packages, repos = "https://cloud.r-project.org")
} else {
  message("All required R packages are already installed.")
}

library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(cowplot)

# Get command-line args: 
# 1 = file2sample.csv, 2 = STATS_FOLDER path
args <- commandArgs(trailingOnly = TRUE)
STATS_FOLDER <- paste0(args[1], "/")
ANNOTATE_VARS_FOLDER <- paste0(args[2], "/")
STAR_STATS_FOLDER <- paste0(args[3], "/")

# load statistics for all mutation and further filter
all = read.csv(paste0(STATS_FOLDER, "stats_hc_all.txt"), sep=":", header=FALSE)
colnames(all) = c("cell_line", "all")
all$cell_line = gsub("\\.hc\\.vcf\\.gz", "", all$cell_line)
pass = read.csv(paste0(STATS_FOLDER, "stats_hc_pass.txt"), sep=":", header=FALSE)
pass2 = read.csv(paste0(STATS_FOLDER, "stats_hc_pass2.txt"), sep=":", header=FALSE)
lcr = read.csv(paste0(STATS_FOLDER, "stats_hc_lcr.txt"), sep=":", header=FALSE)
kG = read.csv(paste0(STATS_FOLDER, "stats_hc_1kG.txt"), sep=":", header=FALSE)
dbsnp = read.csv(paste0(STATS_FOLDER, "stats_hc_dbsnp.txt"), sep=":", header=FALSE)

all$pass = pass$V2
all$edit = pass2$V2
all$lcr = lcr$V2
all$kG = kG$V2
all$dbsnp = dbsnp$V2
all_long <- gather(all, filter, count, all:dbsnp, factor_key = TRUE)

# load mapping statistics
common_mut = fread(paste0(ANNOTATE_VARS_FOLDER, "snv_indel_maf.tsv"))
sample20 = fread(paste0(ANNOTATE_VARS_FOLDER, "snv_indel_maf_less20.tsv"))
tab = common_mut %>% count(sampleID) %>% as.data.frame()
tab = merge(tab, by="sampleID", sample20 %>% count(sampleID) %>% as.data.frame())
colnames(tab) = c("sampleID", "snv_indel", "snv_indel_20")
stats = read.table(paste0(STAR_STATS_FOLDER, "star_stat_summary.tsv"), sep="\t", header=TRUE, row.names=1, check.names=FALSE)
mapped <- as.numeric(stats["unique", ]) / 1e6
names(mapped) <- colnames(stats)  # name the vector by sampleID
# Add 'mapped' to tab for the correct sample
tab$mapped <- mapped[tab$sampleID]  # vector length = nrow(tab)

# safety measure in case no variants in tab
if(nrow(tab) > 0){
    all = cbind(all,tab[,-1])
    all = all[,c(1,10, 2:9)]
    all_long  <- gather(all, filter, count, all:snv_indel, factor_key = TRUE)
    all_long2 <- gather(all, filter, count, all:snv_indel_20, factor_key = TRUE)

    # visualise
    pdf("hc_stats_merged_single_violin.pdf", width = 9, height = 9)
    p <- ggplot(data = all_long[!all_long$filter %in% c("kG", "dbsnp", "snv_indel", "snv_indel_20") ,], aes(x = cell_line, y = count, fill = filter)) +
        geom_bar(stat = "identity", position=position_dodge()) +
        theme_bw(base_size = 12, base_family = "") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position="none") +
        facet_wrap(~filter, ncol = 1) +
        scale_fill_brewer(palette = "Paired") +
        ggtitle("Mutation filter") +
        xlab("") + ylab("Number of mutations")
    p1 <- ggplot(data = all_long2[!all_long$filter %in% "kG",], aes(x = filter, y = count)) +
        geom_violin() + #geom_jitter(width = 0.2) +
        geom_boxplot(width=0.2, color="grey", alpha=0.2) +
        scale_y_log10() +
        theme_bw(base_size = 12, base_family = "") +
        # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggtitle("Filtering steps") +
        xlab("") + ylab("Number of mutations")
    # number of filtered variants correlation to star statistics
    p2 <- ggplot(data = all_long2[all_long2$filter %in% c("all", "pass", "dbsnp", "snv_indel_20"),], aes(x = mapped, y = count)) +
        ggtitle("Reads vs variants") +
        xlab("Mapped reads in million") + ylab("Number of mutations") +
        theme_bw(base_size = 12, base_family = "") +
        facet_wrap( ~ filter, ncol=1, scales="free_y") +
        stat_smooth(method = "lm", color="black", formula = y ~ x) + # linear regression line, by default includes 95% confidence region
        geom_point() +
        stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
                r.accuracy = 0.01)
    pp1 = plot_grid(p, p1, labels = c("A", "B"), ncol=1, rel_heights = c(2, 1))
    final_plot = plot_grid(pp1, p2, ncol=2, labels = c("", "C"), rel_widths = c(5, 2))
    print (final_plot)
    dev.off()

    # for review BMC res notes
    pdf("hc_stats_merged_single_violin2.pdf", width = 9, height = 9)
    p <- ggplot(data = all_long[!all_long$filter %in% c("kG", "dbsnp", "snv_indel", "snv_indel_20") ,], aes(x = cell_line, y = count, fill = filter)) +
        geom_bar(stat = "identity", position=position_dodge()) +
        theme_bw(base_size = 12, base_family = "") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.position="none") +
        facet_wrap(~filter, ncol = 1) +
        scale_fill_brewer(palette = "Paired") +
        ggtitle("Mutation filter") +
        xlab("") + ylab("Number of mutations")
    p1 <- ggplot(data = all_long2[!all_long2$filter %in% "kG",], aes(x = filter, y = count)) +
        geom_violin() + #geom_jitter(width = 0.2) +
        geom_boxplot(width=0.2, color="grey", alpha=0.2) +
        scale_y_log10() +
        expand_limits(y=c(0.1,max(all_long2$count))) +
        theme_bw(base_size = 12, base_family = "") +
        # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ggtitle("Filtering steps") +
        xlab("") + ylab("Number of mutations")
    # number of filtered variants correlation to star statistics
    p2 <- ggplot(data = all_long2[all_long2$filter %in% c("all", "pass", "dbsnp", "snv_indel_20"),], aes(x = mapped, y = count)) +
        ggtitle("Reads vs variants") +
        xlab("Mapped reads in million") + ylab("Number of mutations") +
        theme_bw(base_size = 12, base_family = "") +
        facet_wrap( ~ filter, ncol=1, scales="free_y") +
        stat_smooth(method = "lm", color="black", formula = y ~ x) + # linear regression line, by default includes 95% confidence region
        geom_point() +
        stat_cor(aes(label = paste(..rr.label..)), # adds R^2 value
                r.accuracy = 0.01)
    pp1 = plot_grid(p, p1, labels = c("a", "b"), ncol=1, rel_heights = c(2, 1))
    final_plot = plot_grid(pp1, p2, ncol=2, labels = c("", "c"), rel_widths = c(5, 2))
    print (final_plot)
    dev.off()
}

