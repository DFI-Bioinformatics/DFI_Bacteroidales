rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(magrittr)
library(ape)

'%!in%' <- function(x,y)!('%in%'(x,y))

# load data ---------

metab.qual.matrix.log2 <- read_csv('../data/metab.qual.matrix.log2.csv')
sppal <- readRDS('../data/sppal.extended.rds')


physub <- read.tree("../data/custom_bac71_phylogenomic-tree.txt")
rootphysub <- root(physub, outgroup = "Escherichia_coli_K_12_MG1655")
  
# clean tree data ------------------------------------------------

treemeta <- tibble(tip.label = rootphysub$tip.label) %>%
  mutate(string1 =lapply(str_split(tip.label, pattern = '_'), '[',1)) %>%
  mutate(string2 =lapply(str_split(tip.label, pattern = '_'), '[',2)) %>% 
  mutate(species = paste(string1, string2)) %>%
  mutate(string3 =lapply(str_split(tip.label, pattern = '_'), '[',3)) %>%
  mutate(string4 =lapply(str_split(tip.label, pattern = '_'), '[',4)) %>%
  mutate(string5 =lapply(str_split(tip.label, pattern = '_'), '[',5)) %>%
  unite(col = msk.id, string3:string5, sep = '.', na.rm = T) %>% 
  mutate(cohort = if_else(grepl('MSK|DFI', msk.id), 'Biobank', 'Type Strain'))

# for custom_tree.txt, need to amend treemeta 
treemeta$species[which(treemeta$msk.id == 'MSK.20.12')] <- 'Bacteroides xylanisolvens'
treemeta$species[which(treemeta$msk.id == 'MSK.9.14')] <- 'Parabacteroides merdae'
treemeta$species[which(treemeta$msk.id == 'DFI.7.71')] <- 'Hoylesella buccalis'
treemeta$species[which(treemeta$msk.id == 'ATCC35310')] <- 'Hoylesella buccalis'
treemeta$species[which(treemeta$msk.id == 'DSM18206')] <- 'Leyella stercorea'
treemeta$species[which(treemeta$msk.id == 'DSM18205')] <- 'Segatella copri'
treemeta$species[which(treemeta$msk.id == 'MSK.20.82')] <- 'Bacteroides thetaiotaomicron'
treemeta$species[which(treemeta$msk.id == 'MSK.18.83')] <- 'Bacteroides ovatus'
treemeta$species[which(treemeta$msk.id == 'MSK.17.76')] <- 'Bacteroides ovatus'
treemeta$species[which(treemeta$msk.id == 'MSK.16.13')] <- 'Bacteroides ovatus'
treemeta$species[which(treemeta$msk.id == 'MSK.21.44')] <- 'Segatella sinica'
treemeta$msk.id[which(treemeta$msk.id == 'DFI.3.60B')] <- 'DFI.3.60'


# list of isolates included in metabolomics

iso.to.remain <- metab.qual.matrix.log2$msk.id
iso.to.drop <- 
  treemeta %>% 
  filter(msk.id %!in% iso.to.remain)

require(tidytree)
tree.reduced <- drop.tip(rootphysub, iso.to.drop$tip.label)
# 111 tips, and isolates with metabolomics have 111 rows
tree.tbl <- tidytree::as_tibble(tree.reduced) %>% 
  left_join(treemeta %>% select(label = tip.label, msk.id), by = 'label') 
# replace label column with  msk.id
# first, fix NA rows in msk.id
tree.tbl$msk.id[is.na(tree.tbl$msk.id)] <- tree.tbl$label[is.na(tree.tbl$msk.id)]
# then replace label column
tree.tbl$label <- tree.tbl$msk.id
tree.parsed <- as.phylo(tree.tbl %>% select(-msk.id))



# Cophylogeny analysis of SCG tree vs. metabolomics tree-------------
require(paco)
require(ape)

# transform tree into pairwise distances
D.SCG <- cophenetic.phylo(tree.parsed)

# get distance matrix from metabolomics matrix
metab.qual.matrix <- metab.qual.matrix.log2 %>% 
  column_to_rownames(var = 'msk.id') %>% 
  as.matrix()
D.metab <- as.matrix(dist(metab.qual.matrix))

# construct SCG-metab (host-parasite) relationship matrix
# row names are hosts, column names are parasites
isolate.name.SCG.metab <- rownames(D.SCG)
link.SCG.metab <- matrix(data = 0, 
                        ncol = length(isolate.name.SCG.metab), 
                        nrow = length(isolate.name.SCG.metab))
isolate.name.SCG.metab -> colnames(link.SCG.metab)
isolate.name.SCG.metab -> rownames(link.SCG.metab)
# link host to parasite, a binary matrix
# Host in rows, parasites in columns
for (i in 1:length(isolate.name.SCG.metab)) {link.SCG.metab[i,i] = 1}


## Procrustes Analysis------------
require(paco)
# construct PACo object
D.SCG.metab <- prepare_paco_data(H = D.SCG, # host
                                P = D.metab, # parasite
                                HP = link.SCG.metab) #host-parasite association
# translate distance matrices into Principal Coordinates
D.SCG.metab <- add_pcoord(D.SCG.metab)
#Procrustes analysis
PACo.SCG.metab <- PACo(D.SCG.metab,
                      nperm = 1000, #permutation time
                      seed = 42, #Answer to the Ultimate Question of Life, the Universe, and Everything
                      method = 'r00', #permutation method
                      shuffled = T) #permutation statistics
print(PACo.SCG.metab$gof)
# p ~ 0.038
# significant co-evolution of SCG vs metabolomics.

residuals_paco(PACo.SCG.metab$proc)
#cophylogeny effect size (ES) = (null_dist - obs_dist) / (null_dist)
ES.PACo <- (mean(PACo.SCG.metab$shuffled) - PACo.SCG.metab$proc[['ss']]) / mean(PACo.SCG.metab$shuffled)
# Effect Size ~ 0.255

## plot Procrustes results------

########  Procrustes Analysis from RStudio: https://rpubs.com/mengxu/procrustes_analysis 

procrustes.diy <- function(A, B){
  # A = PCoA.16S$points %>% as_tibble() %>% select(1:2)
  # B = PCoA.lanI$points%>% as_tibble() %>% select(1:2)
  # center and normalize A 
  A.centered <- t(scale(t(A), center = TRUE, scale = FALSE))
  A.size <- norm(A.centered, type = "F") / (ncol(A) * nrow(A))
  A.normalized <- A.centered / A.size
  
  # center and normalize B
  B.centered <- t(scale(t(B), center = TRUE, scale = FALSE))
  B.size <- norm(B.centered, type = "F") / (ncol(B) * nrow(B))
  B.normalized <- B.centered / B.size
  
  # Rotation matrix T 
  svd.results <- svd(B.normalized %*% t(A.normalized))
  U <- svd.results$u
  V <- svd.results$v
  T <- V %*% t(U)
  
  # B transformed
  B.transformed <- T %*% B.normalized
  
  # Error after superimposition
  RSS <- norm(A.normalized - B.transformed,  type = "F")
  
  # Return
  return(list(A.normalized = A.normalized, B.normalized = B.normalized, rotation.mtx = T, B.transformed = B.transformed, RSS = RSS))
}

require(vegan)
PCoA.SCG.2dim <- wcmdscale(D.SCG) %>% 
  as_tibble(rownames = NA) %>% rownames_to_column('msk_id') %>% 
  arrange(msk_id) %>% select(msk_id, V1, V2) %>% column_to_rownames('msk_id')

PCoA.metab.2dim <- wcmdscale(D.metab) %>% 
  as_tibble(rownames = NA) %>% rownames_to_column('msk_id') %>% 
  arrange(msk_id) %>% select(msk_id, V1, V2) %>% column_to_rownames('msk_id')

# DIY procrustes function
pro <- procrustes.diy(A = t(PCoA.SCG.2dim),B = t(PCoA.metab.2dim)) # two matrices should have same order of rows

procrustes.data <- bind_rows(as.data.frame(t(pro$A.normalized)) %>% rownames_to_column('msk_id')%>% mutate(group = 'SCG'),
                             as.data.frame(t(pro$B.transformed))%>% rownames_to_column('msk_id')%>% mutate(group = 'metab'))

#test significance of similarity
vegan::protest(PCoA.SCG.2dim, PCoA.metab.2dim)
#Significance ~ 0.001

# plot Procrustes PCoA plots with connecting lines
pdf(file = '../results/FigS6C.metab.vs.SCG.procrustes.jitter.pdf', width = 8, height = 6 )  
ggplot(data = procrustes.data,
       aes(x = V1, y = V2, group = group, color = group)) + 
  geom_line(aes(V1, V2, group=msk_id), color='grey', alpha = 0.2) + # plot linking lines
  geom_jitter(width = 0.5, height = 0.5) +
  xlab('PCoA1') + ylab('PCoA2') +
  theme_bw() +
  scale_color_manual(name = "",
                     values = c( "SCG" = "#E69F00", "metab" = "#56B4E9"),
                     labels = c("Single-copy Core Genes", "Metabolomics")
  ) +
  annotate(geom = 'text',label = 'PACo significance = 0.038',
           x= 30, y = 50, size = 4)
dev.off()  


## plot 16S tree with metabolomics tree face-to-face-------------

require(ape)

tree.metab <- njs(D.metab)

require(ggtree)
p.metab <- ggtree(tree.metab,
                  branch.length = 'none') #cladogram

p.SCG <- ggtree(tree.parsed,
                branch.length = 'none') #cladogram

#extract position data of two trees
d.SCG <- p.SCG$data
d.metab <- p.metab$data

# reverse x-axis of metab tree and set offset to make the tree in the right hand side of the SCG tree
d.metab$x <- max(d.metab$x) - d.metab$x + max(d.SCG$x) + 20 #for cladogram

# flip metab tree vertically 
# d.metab$y <- mean(d.metab$y)*2 - d.metab$y

dd.SCG.metab <- bind_rows(d.SCG, d.metab) %>%
  filter(isTip == TRUE)

Figure.SCG.metab <- p.SCG + geom_tree(data=d.metab) + 
  geom_line(aes(x, y, group=label), data=dd.SCG.metab, color='grey') +
  geom_tiplab(align = F, size = 2) +
  geom_tiplab(data = d.metab, hjust=1, align = F, size = 2) +
  ggtitle('SCG cladogram (left) versus metabolomics cladogram (right)') +
  theme(plot.title = element_text(hjust = 0.5))

Figure.SCG.metab
ggsave(plot = Figure.SCG.metab,
       filename = '../results/FigS6B.metab.vs.SCG.cladogram.pdf',
       device = 'pdf', width = 8, height = 8)
dev.off()
