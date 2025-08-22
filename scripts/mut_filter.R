###################################################
## Get mutation numbers for each sample
## After filtering common mutants
## Using gnomad, dbsnp, 1000genomes filters
## And after filtering mutations present in >20% of samples
###################################################

# Check/install required packages
list.of.packages <- c("dplyr", "data.table", "ggplot2", "tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
  message("Installing missing R packages: ", paste(new.packages, collapse = ", "))
  install.packages(new.packages, repos = "https://cloud.r-project.org")
} else {
  message("All required R packages are already installed.")
}

# Load libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(tidyr)

# Get command-line args: 
# 1 = file2sample.csv, 2 = STATS_FOLDER path
args <- commandArgs(trailingOnly = TRUE)
file2sample <- args[1]
STATS_FOLDER <- paste0(args[2], "/")

# Expression name for output files (based on current directory)
# expname = paste0("_", basename(getwd()), "_gatk_hc_dbsnp")
expname = ""

# Read samples and remove duplicates by cell_line
samples <- read.csv(file2sample)
samples <- samples[!duplicated(samples$cell_line), ]
Load_vcf_hc <- function(vcf_file, sample){
  VCF <- as.data.frame(fread(cmd = paste("zgrep -v '##' ", vcf_file, sep = ''), sep = '\t'))
  VCF <- dplyr::rename(VCF, CHROM = `#CHROM`)
  VCF <- VCF[VCF$CHROM %in% c(paste('chr',seq(1,22), sep = ''),'chrX','chrY'),]
  
  # keep PASS mutations
  VCF2 <- VCF[VCF$FILTER == 'PASS',]
  
  # Extract maf file
  VCF2 <- tidyr::separate(data = VCF2, col = sample, sep = ':', 
                          into = c('T_GT','T_AD','Depth','T_GQ','T_PL'))
  VCF2$REF_count <- apply(VCF2, 1, function(x) as.numeric(unlist(strsplit(x = x["T_AD"], split = ','))[1]))
  VCF2$ALT_count <- apply(VCF2, 1, function(x) as.numeric(unlist(strsplit(x = x["T_AD"], split = ','))[2]))
  VCF2$Depth <- as.numeric(VCF2$Depth) # We want this value to be considered as a number
  
  # Keep interesting columns
  VCF3 <- VCF2[,c('CHROM','POS','REF','ALT','QUAL','FILTER', 'Depth','REF_count','ALT_count')]
  
  # Compute Variant Allele Frequency (VAF) in the tumour sample
  VCF3$VAF <- VCF3$ALT_count/VCF3$Depth
  return(VCF3)
}

# get variants
Mut_calls <- NULL
count = 1
for (sample in samples$sample_id) {
  bc = Load_vcf_hc(paste(STATS_FOLDER, sample, ".hc.pass2.lcr.dbsnp.vcf.gz", sep=""), sample)
  bc$sampleID = sample
  bc <- bc[bc$Depth >= 5,] # filter depth
  
  # read in maf file
  files_vep = paste(STATS_FOLDER, sample, ".hc.pass2.lcr.dbsnp.vcf.maf.gz", sep="")
  vep = read.csv(gzfile(files_vep), header=T, stringsAsFactors=F, sep = "\t", comment.char="#")
  # filter coding variants
  vep = vep[grep("protein_coding", vep$BIOTYPE), c(1:16, 35:43, 47:75, 77:87, 93:114)]
  vep = vep[vep$HGVSp!="", ]
  # gnomad filter
  vep = vep[vep$FILTER != "common_variant", ]
  # filter 1000 genomes
  vep = vep[is.na(vep$AF) | vep$AF<0.01, ]
  vep$id = paste(vep$Chromosome, vep$Start_Position, sep="_")
  # account for deletion with pos-1
  tmp = which(vep$Variant_Type == "DEL")
  vep$id[tmp] = paste(vep$Chromosome[tmp], vep$Start_Position[tmp]-1, sep="_")
  
  bc$id = paste(bc$CHROM, bc$POS, sep="_")
  bc = bc[bc$id %in% vep$id,]
  bc = merge(bc, vep[,c("id", "Chromosome", "Start_Position", "End_Position", "Strand", "HGVSc", "HGVSp", "HGVSp_Short", "Variant_Classification", "Consequence", "Existing_variation", "Hugo_Symbol", "Gene", "Entrez_Gene_Id", "SIFT", "PolyPhen", "IMPACT")], by="id")
  
  Mut_calls[[count]] <- bc
  count = count + 1
}
# Usually faster than rbind in each iteration of the for loop
Mut_calls <- do.call("rbind", Mut_calls)
# rearrange column order
Mut_calls = Mut_calls[,c(1,12:16, 4:5, 8:11, 17:ncol(Mut_calls))]

# filter mutations occuring in more than 20% of samples
tmp = as.data.frame(table(Mut_calls$id))
tmp = as.vector(tmp$Var1[tmp$Freq >= ceiling(nrow(samples)/5)])
Mut_calls$Freq = "<=20%"
Mut_calls$Freq[Mut_calls$id %in% tmp] = ">20%"
write.table(Mut_calls[,-1], paste('snv_indel_maf',expname,'.tsv',sep=''), row.names=F, quote=F, sep="\t")
write.table(Mut_calls[!Mut_calls$id %in% tmp,-1], paste('snv_indel_maf_less20',expname,'.tsv',sep=''), row.names=F, quote=F, sep="\t")

# store summary numbers in one table
common_mut = data.table::fread(paste('snv_indel_maf',expname,'.tsv',sep=''))
sample20 = data.table::fread(paste('snv_indel_maf_less20',expname,'.tsv',sep=''))
tab = common_mut %>% count(sampleID) %>% as.data.frame()
tab = merge(tab, by="sampleID", sample20 %>% count(sampleID) %>% as.data.frame())
colnames(tab) = c("sampleID", "snv_indel", "snv_indel_20")
write.table(tab, paste('snv_indel_summary',expname,'.tsv',sep=''), row.names=F, quote=F, sep="\t")