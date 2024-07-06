rm(list=ls())

library(tidyverse)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load data ---------------------------------------------------------------

nucmer.results <- 
  read_tsv('../data/416_Illumina_only_nucmer_results.tsv',
           skip = 1)


# plot nucmer differences between same-donor isolates and different-donor isolates----------

nucmer.plot.data <- 
  nucmer.results %>% 
  filter(donor == 'same' | donor == 'different')
nucmer.plot.data$donor <- as.factor(nucmer.plot.data$donor)

## ANI plot------
require(ggbreak)
ggplot(nucmer.plot.data, mapping = aes(x = donor,
                                       y = AvgIdentity)) + 
  geom_violin(mapping = aes(color = donor),
              scale = 'width') +
  geom_jitter(alpha = 0.1, size = 1) +
  theme_bw() +
  ylim(85, 100) +
  scale_y_cut(breaks = c(95,99.9), which = c(1,2,3), scales = c(3,3,1)) +
  ylab('Average Nucleotide Identity (%)')
ggsave(filename = '../results/FigS2A.ANI.cutoff.pdf', 
       device = 'pdf', width = 6, height = 6)


## TotalGSNPs plot -------
require(ggbreak)
ggplot(nucmer.plot.data, mapping = aes(x = donor,
                                       y = TotalGSNPs)) + 
  geom_violin(mapping = aes(color = donor),
              scale = 'width') +
  geom_jitter(alpha = 0.1, size = 1) +
  theme_bw() +
  # ylim(85, 100) +
  scale_y_cut(breaks = c(1000), which = c(1,2), scales = c(1,1)) +
  ylab('Total Genomic SNPs')
ggsave(filename = '../results/FigS2B.TotalGenomicSNPs.cutoff.pdf', 
       device = 'pdf', width = 6, height = 6)



## TotalGInDels plot ------
require(ggbreak)
ggplot(nucmer.plot.data, mapping = aes(x = donor,
                                       y = TotalGIndels)) + 
  geom_violin(mapping = aes(color = donor),
              scale = 'width') +
  geom_jitter(alpha = 0.1, size = 1) +
  theme_bw() +
  # ylim(85, 100) +
  scale_y_cut(breaks = c(75), which = c(1,2), scales = c(1,1)) +
  ylab('Total Genomic InDels')
ggsave(filename = '../results/FigS2C.TotalGenomicInDels.cutoff.pdf', 
       device = 'pdf', width = 6, height = 6)


