# Figure 3. UMAP analysis of the genomic diversity of Bacteroidales isolates
# Figure S4. UMAP analysis of genomic diversity of Bacteroidales isolates within genera or family.  

library(tidyverse)
library(umap)
library(conflicted)

conflict_prefer("n", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("count", "dplyr")
conflicts_prefer(dplyr::filter)

# load data ---------------------------------------------------------------

taxpal <- readRDS("./data/sppal.extended.rds")

tax <- read_csv("data/genome.size.csv") %>% 
  select(msk_id, genus, species, type)

gcl.mat <- read_csv("data/Anvio_geneCluster.matrix.csv") %>% 
  as.matrix()

# umap settings -----------------------------------------------------------

custom.config <- umap.defaults
custom.config$n_neighbors <- 70

random.stats <- sample(10000, 10)

# Fig 3 -------------------------------------------------------------------

## loop over random state: --------------------------------------------------

for (ri in random.stats){
  
  custom.config$random_state = ri 
  
  clus1 <- umap(gcl.mat,config = custom.config)
  
  uplt1 <- clus1$layout %>%
    as.data.frame() %>%
    mutate(msk_id=row.names(.)) %>%
    left_join(tax)
  
  ggumapp1 <- uplt1 %>%
    ggplot(aes(x=V1,y=V2)) +
    theme_minimal() +
    geom_point(aes(fill=species, shape = type, size = type),alpha=0.8) +
    scale_fill_manual(values=taxpal) +
    scale_shape_manual(values=c(21,24), guide=FALSE) +
    scale_size_manual(values=c(3.5,5), guide=FALSE) +
    theme(legend.title=element_blank()) +
    labs(title = "UMAP of 162 Bacteroidales isolates with type strains",
         subtitle = "Anvio Gene Clusters",
         x = "", y = "",
         caption = paste0("n_neighbors=", custom.config$n_neighbors)) +
    theme(legend.position = "none") +
    ggrepel::geom_text_repel(data = uplt1 %>%
                               filter(type == "Type Strain") %>% 
                               mutate(lab = sub("[a-z]+",".", species)),
                             aes(label = lab, color = species),
                             bg.color = "black",fontface = "italic",
                             bg.r = 0.05, size = 3.15
    ) +
    scale_color_manual(values=taxpal)
  
  ggumapp1
  
  ## legend
  legendgg <- uplt1 %>%
    ggplot(aes(x=V1,y=V2)) +
    geom_point(aes(color = species),alpha = 0.8, shape = 15, size = 3.8) +
    scale_color_manual(values= taxpal) +
    theme_minimal() +
    guides(color=guide_legend(ncol=2)) +
    theme(legend.text = element_text(face ="italic"))
  
  splegend <- ggpubr::get_legend(legendgg)
  
  ## assemble graph
  ggpubr::ggarrange(ggumapp1, splegend, widths = c(2.15, 1.25))
  
  ggsave(paste0("results/Fig3.umap.",ri,".pdf"), height=6.75,width=10.65)
  
  # write_csv(uplt1, paste0("./rds/genomic.all.umap.coord.",ri,".csv"))
  
}


## family/genus level sub-UMAP --------------------------------------------------

sub.genus <- list("Alistipes" = "Alistipes",
                  "Bacteroides" = "Bacteroides",
                  #  "Butyricimonas" = "Butyricimonas",
                  "Prevotellaceae" = c("Hoylesella",
                                       "Leyella",
                                       "Segatella"),
                  # "Odoribacter" = "Odoribacter",
                  "Parabacteroides" = "Parabacteroides",
                  "Phocaeicola" = "Phocaeicola")

for (ri in random.stats){
  
  # initiate an empty list
  ri.list <- list()
  
  # UMAP config: random state setting
  custom.config <- umap.defaults
  custom.config$random_state = ri 
  
  for (gi in names(sub.genus)){
    
    tar.genus <- sub.genus[[gi]]
    
    tar.msk <- tax %>% 
      filter(genus %in% tar.genus) %>% 
      pull(msk_id)
    
    gcl.mat.tmp <- gcl.org %>% 
      ungroup() %>%
      filter(msk_id %in% tar.msk) %>% 
      count(msk_id, gene_cluster_id) %>% 
      spread(gene_cluster_id, n, fill = 0) %>% 
      column_to_rownames(var = "msk_id") %>% 
      as.matrix()
    
    dim(gcl.mat.tmp)
    
    # get count of non-ubiquitous genes
    nuc <- gcl.org %>% 
      filter(msk_id %in% tar.msk) %>%
      ungroup() %>% 
      distinct(genome_name, gene_cluster_id) %>% 
      dplyr::count(gene_cluster_id) %>%
      summarise(num = sum(n!=nrow(gcl.mat.tmp))) %>%
      `$`(num)
    
    custom.config$n_neighbors <- ceiling(nrow(gcl.mat.tmp) * 0.45)
    
    clus2 <- umap(gcl.mat.tmp,config = custom.config)
    
    uplt2 <- clus2$layout %>%
      as.data.frame() %>%
      mutate(msk_id=row.names(.)) %>%
      left_join(tax) %>% 
      mutate(random_state = ri,
             n_neighbors = custom.config$n_neighbors,
             # facets = paste0(gi,"\nn_neighbor=", custom.config$n_neighbors),
             facets = gi)
    
    ri.list[[gi]] <- uplt2
    
  }
  
  ri.tmp.df <- bind_rows(ri.list) 
  
  ri.tmp.df %>% 
    ggplot(aes(x=V1,y=V2)) +
    theme_bw() +
    geom_point(aes(fill=species, shape = type, size = type),alpha=0.8) +
    scale_fill_manual(values=taxpal) +
    scale_shape_manual(values=c(21,24), guide=FALSE) +
    scale_size_manual(values=c(3.5,5), guide=FALSE) +
    labs(title = "sub-UMAP of chosen isolates with type strains",
         subtitle = "Anvio Gene Clusters",
         x = "", y = "") +
    theme(legend.position = "none") +
    facet_wrap("facets") +
    ggrepel::geom_text_repel(data = ri.tmp.df %>%
                               filter(type == "Type Strain") %>% 
                               mutate(lab = if_else(facets == "Bacteroides",
                                                    gsub("Bacteroides ", "B. ", species),
                                                    species)),
                             aes(label = lab, color = species),
                             bg.color = "black",fontface = "italic",
                             bg.r = 0.05, size = 3.25)+
    scale_color_manual(values=taxpal) +
    theme(legend.title=element_blank(),
          # panel.grid.major = element_blank(),
          # panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_text(face = "bold.italic", size = 13))
  
  ggsave(paste0("FigS4.sub-umap.",ri,".pdf"), height=8.15,width=11.75)
  
  # write_csv(ri.tmp.df, paste0("./rds/sub-umap.coord.",ri,".csv"))
  
}

