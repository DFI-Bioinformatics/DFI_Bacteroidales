# Figure 4. CAZyme and PUL analyses of Bacteroidales isolates.  
# Figure S5. Additional CAZyme and PUL analyses of Bacteroidales isolates.  

### parse hmm results from hmmsearch result with cazyme V12

library(tidyverse)

# script setup ------------------------------------------------------------

## read in hmm output ------------------------------------
new.dbcan <- read_tsv("data/new_output_from_CAZyme_hmmsearch.txt",
                      comment = "#", col_names = F)  %>%
  separate(X1, c("target","tacc","locus_tag", "qacc", "seq_evalue","seq_score",
                 "seq_bias","dom_evalue","dom_score","dom_bias",
                 "exp","reg","clu","ov","env","dom","rep","inc","description"),
           convert = T, extra = "merge",
           sep = "\\s+")

summary(new.dbcan$seq_evalue)
# seq_evalue is been filtered to 1e-18

new.dbcan %>% 
  count(target, sort = T)
# multiple hits with one target

filtered.dbcan <- new.dbcan %>%
  # filter(seq_evalue<1e-10) %>%
  group_by(target) %>%
  top_n(1, wt = -seq_evalue) %>%
  top_n(1, wt = dom_score) %>%
  arrange(desc(seq_evalue)) %>%
  dplyr::slice(1) %>%
  # ungroup() %>%
  # add_count(locus_tag) %>%
  # filter(n > 1)
  select(target, locus_tag, seq_evalue) %>%
  ungroup()

# switch column names from hmmscan
colnames(filtered.dbcan) <- colnames(filtered.dbcan)[c(2,1,3)]

## find old color palette and isolate meta sheet ---------------------------

genome.size <- read_csv("data/genome.size.csv")

anno <- readRDS("data/illumina_typeStrain.prokka.rds")

taxpal <- readRDS("data/sppal.extended.rds")

## cazyme family -----------------------------------------------------------

cazymecat <- tibble(cazyfamily = c("AA","CBM","CE","GH","GT","PL","SLH"),
                    cazyname = c("Auxiliary Activities","Carbohydrate-Binding Modules",
                                 "Carbohydrate Esterases","Glycoside Hydrolases",
                                 "GlycosylTransferases", "Polysaccharide Lyases",
                                 "S-layer homologous"))

clean.dbcan <- filtered.dbcan %>% 
  mutate(cazy=gsub("\\.hmm","",target),
         cazyfamily=str_extract(cazy,pattern="[A-Z]+")) %>% 
  left_join(anno %>% 
              select(locus_tag, msk_id)) %>% 
  left_join(genome.size %>% 
              select(msk_id, species, type))

dbcount <- clean.dbcan %>% 
  filter(!is.na(msk_id)) %>%
  count(msk_id, species, type, cazyfamily, cazy)

dbcount %>% 
  select(-cazy) %>% 
  pivot_wider(names_from = cazyfamily,
              values_from = n,
              values_fn = sum, values_fill = 0) %>% 
  arrange(species) %>% 
  write_csv("data/parsed.CAZyme_hmmsearch.count.csv")

dbcount %>% 
  select(-cazy) %>% 
  mutate(n = 1) %>% 
  pivot_wider(names_from = cazyfamily,
              values_from = n,
              values_fn = sum, values_fill = 0) %>% 
  arrange(species) %>% 
  write_csv("data/parsed.CAZyme_hmmsearch.unique.csv")

genome.size %>% 
  anti_join(dbcount %>% distinct(msk_id))
  
dbcount %>% 
  count(msk_id, species, type, cazyfamily) %>% 
  spread(cazyfamily, n, fill = 0) %>% 
  arrange(species) %>% 
  view

# Figure 4 A CAZyme numbers ---------------------------------------------------

my.cazy <- c("Glycoside Hydrolases", 
             "Polysaccharide Lyases", 
             "Carbohydrate Esterases", 
             "Carbohydrate-Binding Modules", 
             "GlycosylTransferases")

dbpltdf <- dbcount %>%
  left_join(cazymecat) %>%
  # filter(!is.na(cazyname)) %>% 
  mutate(cazyname = factor(cazyname,
                           levels = my.cazy)) 

dbpltdf %>% 
  filter(is.na(cazyname))

dbpltdf %>% 
  filter(cazyname %in% my.cazy) %>% 
  count(msk_id, species, type, cazyname, cazyfamily, wt = n) %>% 
  ggplot(aes(n, reorder(species, n))) +
  geom_jitter(aes(fill = species, shape = type), alpha = 0.65,
              width = 0.15, height = 0.35, size = 3) +
  scale_fill_manual(values = taxpal, guide = "none") +
  scale_shape_manual(values = c(21,24), guide = "none") +
  facet_grid( . ~ cazyname, scales = "free") +
  theme_bw() +
  labs(x = "Number of CAZymes",
       y = "Species") +
  theme(axis.text.y = element_text(face = "italic"),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        strip.text.x = element_text(face = "bold", size = 10.5),
        axis.title = element_text(size = 13, face = "bold"))+
  scale_x_continuous(breaks= scales::pretty_breaks())

ggsave(paste0("results/Fig4_A.cazyme.summary",
              ".pdf"), width = 12.5, height = 6.5)

# Fig S5 B  unique number of cazymes --------------------------------------

# number of unique cazymes
dbcount %>% 
  count(msk_id, species, type, cazyfamily) %>% 
  spread(cazyfamily, n, fill = 0) %>%
  gather("cazyfamily", "n", -c(msk_id,species,type)) %>%
  left_join(cazymecat) %>%
  filter(cazyname %in% my.cazy) %>% 
  mutate(cazyname = factor(cazyname,
                           levels = my.cazy)) %>% 
  ggplot(aes(n, reorder(species, n))) +
  geom_jitter(aes(fill = species, shape = type), alpha = 0.65,
              width = 0.15, height = 0.35, size = 3) +
  scale_fill_manual(values = taxpal, guide = "none") +
  scale_shape_manual(values = c(21,24), guide = "none") +
  facet_grid( . ~ cazyname, scales = "free") +
  theme_bw() +
  labs(x = "Number of Unique CAZymes",
       y = "Species")+
  theme(axis.text.y = element_text(face = "italic"),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        strip.text.x = element_text(face = "bold", size = 10.5),
        axis.title = element_text(size = 13, face = "bold")) +
  scale_x_continuous(breaks= scales::pretty_breaks())

ggsave(paste0("results/FigS5_B.cazyme.unique.summary",
              ".pdf"), width = 12.5, height = 6.5)

# Fig 5S A inverse simpson ------------------------------------------------

## cazyme inverse simpson ---------------------------------------

dbmat <- dbcount %>% 
  select(msk_id, cazy, n ) %>% 
  spread(cazy, n, fill = 0) %>%
  column_to_rownames(var = "msk_id") 

library(vegan)

inv <- diversity(dbmat, index = "inv")

invdf <- tibble(msk_id = names(inv),
                invSimp =  inv)

dbcantop.df <- dbcount %>%
  count(msk_id, species, type, wt = n, name = "nCazyme") %>%
  left_join(invdf) 

## remove linear line and add species label --------------------------------

label.df <- dbcantop.df %>% 
  filter(type == "Type Strain") %>% 
  mutate(label = sub("[a-z]+","\\.", species)) 

dbcantop.df %>%
  ggplot(aes(nCazyme, invSimp)) +
  # geom_smooth(method = "lm", linetype = 2, color = "black") +
  geom_point(aes(fill = species, shape = type), alpha = 0.65, size = 3.5) +
  scale_fill_manual(values = taxpal) +
  scale_shape_manual(values = c(21,24), guide = "none") +
  labs(x = "Number of CAZymes",
       y = "Inverse Simpson") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.y = element_text(face = "italic"),
        strip.background = element_blank(),
        axis.text = element_text(size = 12),
        strip.text.x = element_text(face = "bold", size = 10.5),
        axis.title = element_text(size = 13, face = "bold")) +
  ggrepel::geom_text_repel(data = label.df, 
                           aes(label = label, color = species),
                           bg.color = "black",fontface = "italic",
                           bg.r = 0.05, size = 3.15) +
  scale_color_manual(values = taxpal)

ggsave(paste0("results/FigS5_A.cazyme.invSimp.rmline",
              ".pdf"), width = 6.5, height = 5.8)

# umap of cazyme -------------------------------------------------

library(umap)

custom.config <- umap.defaults
custom.config$n_neighbors <- 85

random.stats <- sample(10000, 10)

# get donor code
donor.df <- read_csv("data/donor.16S.csv") %>% 
  distinct(donor_id, new_id) %>% 
  mutate(community = case_when(
    grepl("^FC", donor_id) ~ paste0("MSK.",gsub("FC0+","",donor_id)),
    grepl("^DFI", donor_id) ~ gsub("0+",".",donor_id),
    TRUE ~ "Type Strain"
  ))

## loop over random state: --------------------------------------------------

for (ri in random.stats){
  
  custom.config$random_state <- ri 
  
  clus1 <- umap(dbmat,config = custom.config)
  
  uplt1 <- clus1$layout %>%
    as.data.frame() %>%
    mutate(msk_id=row.names(.)) %>%
    left_join(genome.size)
  
  ## Fig 4 B: color by species ----------------------------------------------------------
  
  ggumapp1 <- uplt1 %>%
    ggplot(aes(x=V1,y=V2,label = msk_id)) +
    theme_bw() +
    geom_point(aes(fill=species, shape = type, size = type),alpha=0.65) +
    scale_fill_manual(values=taxpal) +
    scale_shape_manual(values=c(21,24), guide=FALSE) +
    scale_size_manual(values=c(3.5,5), guide=FALSE) +
    theme(legend.title=element_blank()) +
    labs(title = "UMAP of Bacteroidales based on presense and absence of cazyme families",
         subtitle = "colored by species",
         x = "", y = "",
         caption = paste0("n_neigh=",custom.config$n_neighbors)) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          axis.text = element_text(size = 12),
          strip.text.x = element_text(face = "bold", size = 10.5),
          axis.title = element_text(size = 13, face = "bold"))  +
    ggrepel::geom_text_repel(data = uplt1 %>%
                               filter(type == "Type Strain") %>%
                               mutate(lab = sub("[a-z]+",".", species)),
                             aes(label = lab, color = species),
                             bg.color = "black", fontface = 'italic',
                             bg.r = 0.05,
                             max.overlaps = 12)+
    scale_color_manual(values=taxpal) 
  
  ggumapp1
  
  ggumapp1 <- ggumapp1 + 
    coord_cartesian(xlim = c(min(uplt1$V1)-1.45, max(uplt1$V1) + 0.25), 
                    ylim = c(min(uplt1$V2)-1.45, max(uplt1$V2) + 0.25))
  
  ggumapp1
  
  uplt1
  
  legendgg <- uplt1 %>%
    ggplot(aes(x=V1,y=V2)) +
    geom_point(aes(color = species),alpha = 0.8, shape = 15, size = 3.8) +
    scale_color_manual(values= taxpal) +
    theme_minimal() +
    guides(color=guide_legend(ncol=2)) +
    theme(legend.text = element_text(face="italic"))
  
  splegend <- ggpubr::get_legend(legendgg)
  # cowplot::get_legend(legendgg)
  ggpubr::ggarrange(ggumapp1, splegend, widths = c(2.15, 1.3))
  
  ggsave(paste0("results/Fig4_B.umap_cazyme.bySpecies",
                ".",ri,".pdf"), height=6.75,width=10.25)
  
  ## Fig S5 C: color by donor ----------------------------------------------------------
  
  uplt2 <- uplt1 %>% 
    mutate(community = gsub("\\.[0-9]+[A-Z]?$","",msk_id)) %>% 
    left_join(donor.df) %>% 
    mutate(new_id = if_else(is.na(new_id),"Type Strain", new_id))
  
  donor_list <- unique(uplt2$new_id)
  
  donorpal <- taxpal[1:length(donor_list)]  
  names(donorpal) <- donor_list
  
  ggumapp2 <- uplt2 %>%
    ggplot(aes(x=V1,y=V2,label = msk_id)) +
    theme_bw() +
    geom_point(aes(fill=new_id, shape = type, size = type),alpha=0.8) +
    scale_fill_manual(values=donorpal) +
    scale_shape_manual(values=c(21,24), guide=FALSE) +
    scale_size_manual(values=c(3.5,5), guide=FALSE) +
    theme(legend.title=element_blank()) +
    labs(title = "UMAP of Bacteroidales based on presense and absence of cazyme families",
         subtitle = "colored by donor",
         x = "", y = "",
         caption = paste0("n_neigh=",custom.config$n_neighbors)) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          axis.text = element_text(size = 12),
          strip.text.x = element_text(face = "bold", size = 10.5),
          axis.title = element_text(size = 13, face = "bold"))
  
  ggumapp2
  
  uplt2
  
  legendgg <- uplt2 %>%
    ggplot(aes(x=V1,y=V2)) +
    geom_point(aes(color = new_id),alpha = 0.8, shape = 15, size = 3.8) +
    scale_color_manual(values= donorpal) +
    theme_minimal() +
    guides(color=guide_legend(ncol=1)) +
    labs(color = "Donor")
  
  splegend <- ggpubr::get_legend(legendgg)
  # cowplot::get_legend(legendgg)
  ggpubr::ggarrange(ggumapp2, splegend, widths = c(2.55, 0.65))
  
  ggsave(paste0("results/FigS5_C.umap_cazyme.byDonor.",ri,".",
                ".pdf"), height=6.75,width=8.75)
}






