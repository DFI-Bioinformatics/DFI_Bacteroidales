


rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)


# read data --------------
isolate_table <- read_csv("../data/genome.size.csv")

gene_cluster_table <- read_csv("../data/Anvio_geneCluster.matrix.csv") %>% 
  pivot_longer(!msk_id, names_to = "gene_cluster_id", values_to = "gene_count") %>%
  left_join(isolate_table) %>% 
  mutate(cluster_presence = if_else(gene_count > 0, 1, 0))


# count the number of times gene cluster appears at each tax level --------------
gene_cluster_table <- gene_cluster_table %>% 
  group_by(phylum, class, order, family, genus, species) %>% 
  mutate(strain_count = length(unique(strain))) %>% 
  group_by(phylum, class, order, family, genus) %>% 
  mutate(species_count = length(unique(species))) %>% 
  group_by(phylum, class, order, family) %>% 
  mutate(genus_count = length(unique(genus))) %>% 
  group_by(phylum, class, order) %>% 
  mutate(family_count = length(unique(family))) %>% 
  group_by(phylum, class) %>% 
  mutate(order_count = length(unique(order))) %>% 
  ungroup()


# identify tax level for each gene cluster --------------
gene_cluster_table <- gene_cluster_table %>% 
  group_by(gene_cluster_id, phylum, class, order, family, genus, species, strain) %>% 
  mutate(strain_level = if_else(n() == sum(cluster_presence), "strain", NA)) %>% 
  group_by(gene_cluster_id, phylum, class, order, family, genus, species) %>% 
  mutate(species_level = if_else(n() == sum(cluster_presence) & strain_count != 1, "species", NA)) %>% 
  group_by(gene_cluster_id, phylum, class, order, family, genus) %>% 
  mutate(genus_level = if_else(n() == sum(cluster_presence) & species_count != 1, "genus", NA)) %>% 
  group_by(gene_cluster_id, phylum, class, order, family) %>% 
  mutate(family_level = if_else(n() == sum(cluster_presence) & genus_count != 1, "family", NA)) %>% 
  group_by(gene_cluster_id, phylum, class, order) %>% 
  mutate(order_level = if_else(n() == sum(cluster_presence), "order", NA)) %>% 
  ungroup()


# consolidate tax levels and identify the hightest level for each gene cluster --------------
gene_cluster_table <- gene_cluster_table %>% 
  mutate(level = paste(order_level, family_level, genus_level, species_level, strain_level,
                       sep = ",")) %>% 
  mutate(level = str_replace_all(level, "NA,|NA", "")) %>% 
  mutate(level = word(level, sep = ",", start = 1)) %>% 
  filter(cluster_presence > 0) %>% 
  mutate(genus_species = paste0(genus, " ", species)) %>% 
  group_by(genome, genus_species, phylum, class, order, family, genus, species, strain, msk_id, level, type, label) %>% 
  summarize(cluster_count = sum(cluster_presence))


# plot core genome comparison --------------
custom_color <- c("#1EB800", "#4F8FE6", "#A93400", "#2619D1", "#FFAB00", "#00CF91")

gene_cluster_table %>% 
  mutate(level = str_to_title(level)) %>% 
  arrange(type, genome) %>% 
  mutate(label = factor(label, levels = unique(.$label))) %>% 
  mutate(level = factor(level, levels = c("Strain", "Species", "Genus", "Family", "Order"))) %>%
  mutate(type = factor(type, levels = c("Biobank", "Type Strain"))) %>% 
  ggplot(aes(x = label, y = cluster_count)) +
  geom_bar(aes(fill=level), stat="identity") +
  geom_point(aes(x= label, y = 5500, shape = type, color = type), size = 3) +
  scale_y_continuous(breaks=c(0,1000,2000,3000,4000, 5000)) +
  facet_grid(. ~ family+genus_species, scales = "free",space = "free") +
  ylab("# of Gene Clusters") +
  xlab("") +
  scale_fill_manual(values = custom_color) +
  scale_color_manual(values = c("#0E3386","#CC3433")) +
  guides(fill = guide_legend(title = "Taxonomy Level:",
                             keywidth = unit(1, "cm"),
                             keyheight = unit(1, "cm")), 
         color = guide_legend(title = "Cohort:",
                              override.aes = list(size=7),
                              keywidth = unit(1, "cm"),
                              keyheight = unit(1, "cm")),
         shape = guide_legend(title = "Cohort:",
                              keywidth = unit(1, "cm"),
                              keyheight = unit(1, "cm"))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=270, hjust = 0, vjust = 0.5 ),
        strip.text.x = element_text(angle=270, face = "italic"),
        axis.title.y = element_text(size = 12, margin = unit(c(0, 4, 0, 0), "mm")),
        legend.title = element_text(size = 12, face = "bold"),
        panel.background = element_rect(fill = "white"),      
        strip.background = element_rect(fill = "white")) 



ggsave("../results/Fig2.core_genome_comparison.pdf", height = 25, width = 100, units = "cm", limitsize = FALSE)










