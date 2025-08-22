# Bioconductor base image
FROM bioconductor/bioconductor_docker:RELEASE_3_21

# Install R and Bioconductor packages
RUN R -e "BiocManager::install(c( \
  'BSgenome', \
  'BSgenome.Hsapiens.UCSC.hg38', \
  'cowplot', \
  'data.table', \
  'dplyr', \
  'ggplot2', \
  'ggpubr', \
  'gplots', \
  'gridExtra', \
  'gtools', \
  'GenomicRanges', \
  'GenVisR', \
  'MutationalPatterns', \
  'NMF', \
  'plyr', \
  'psych', \
  'RColorBrewer', \
  'tidyr', \
  'readxl' \
))"