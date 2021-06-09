library(ggplot2)
library(dplyr)

setwd("C:/Users/ab/Desktop/hse/bioinf/project")

#1NAME <- 'H3K27me3_A549.ENCFF522WJJ.hg19'
#2NAME <- 'H3K27me3_A549.ENCFF522WJJ.hg38'
#3NAME <- 'H3K27me3_A549.ENCFF684ZZH.hg19'
NAME <- 'H3K27me3_A549.ENCFF684ZZH.hg38'

OUT_DIR <- 'Results/'

bed_df <- read.delim(paste0(NAME, '.bed'), as.is = TRUE, header = FALSE)
colnames(bed_df) <- c('chrom', 'start', 'end', 'name', 'score')
bed_df$len <- bed_df$end - bed_df$start
head(bed_df)

ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('len_hist.', NAME, '.pdf'), path = OUT_DIR)

bed_df <- bed_df %>%
  arrange(-len) %>%
  filter(len < 50000)

head(bed_df)

ggplot(bed_df) +
  aes(x = len) +
  geom_histogram() +
  ggtitle(NAME, subtitle = sprintf('Number of peaks = %s', nrow(bed_df))) +
  theme_bw()
ggsave(paste0('len_hist.', NAME, '.filtered.pdf'), path = OUT_DIR)

bed_df %>%
  select(-len) %>%
  write.table(file='Results/H3K27me3_A549.ENCFF522WJJ.hg19.filtered.bed',
              col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)

