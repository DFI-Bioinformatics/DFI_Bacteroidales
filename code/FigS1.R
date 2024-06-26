# Fig S1, donor 16S metagenomics plot

library(readxl)
library(tidyverse)
library(phyloseq)
library(yingtools2)

t <- read_csv("data/donor.16S.csv")

devtools::source_url("https://github.com/yingeddi2008/DFIutility/blob/master/getRdpPal.R?raw=TRUE")

mpal <- getRdpPal(t)

t %>%
  add_count(new_id,wt = pctseqs, name = "total") %>%
  mutate(pctseqs = pctseqs/total) %>%
  count(new_id,Kingdom,Phylum,Class,Order,Family,Genus, wt = pctseqs, name = "pctseqs") %>%
  arrange(Kingdom, Phylum, Class, Order, Family, Genus) %>%
  mutate(Genus = factor(Genus, levels = unique(Genus))) %>%
  group_by(new_id) %>%
  arrange(Genus) %>%
  ggplot(aes(x=new_id,y=pctseqs)) +
  geom_bar(aes(fill=Genus),stat="identity") +
  scale_fill_manual(values=mpal) +
  labs(x = "", y = "16S % abundance",
       title = "rdp 16s taxonomy for healthy donors") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        legend.position="none") +
  scale_y_continuous(expand = c(0.005,0.005))

ggsave("data/FigS1.donor_barplot.pdf", width = 8.15, height = 5)

