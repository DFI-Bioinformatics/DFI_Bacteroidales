
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(tidyverse)

# draw tree---------------

library(grid)
library(ape)
library(yingtools2)
library(ggtree)

physub <- read.tree("../data/custom_bac71_phylogenomic-tree.txt")
rootphysub <- root(physub, outgroup = "Escherichia_coli_K_12_MG1655")

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
treemeta$species[which(treemeta$msk.id == 'MSK.21.56')] <- 'Segatella brunsvicensis'
treemeta$species[which(treemeta$msk.id == 'MSK.21.64')] <- 'Segatella brasiliensis'


sppal <- readRDS('../data/sppal.extended.rds') 

sspcs <- names(sppal)
hilsize = 3 #species name size
xend = 0.618 #highlight ring size
xadd_var = 0.618 #highlight ring size


tree <- ggtree(rootphysub, 
       layout="circular",
       ladderize=T,
       size= 0.33) %<+% #%>% #weight of tree lines
  # flip(190, 192) %>% rotate(172)  %<+%   #adjust tree orders by flipping two nodes or rotating
  treemeta +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[1]),fill=sppal[1],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[2]),fill=sppal[2],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[3]),fill=sppal[3],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[4]),fill=sppal[4],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[5]),fill=sppal[5],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[6]),fill=sppal[6],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[7]),fill=sppal[7],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[8]),fill=sppal[8],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[9]),fill=sppal[9],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[10]),fill=sppal[10],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[11]),fill=sppal[11],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[12]),fill=sppal[12],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[13]),fill=sppal[13],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[14]),fill=sppal[14],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[15]),fill=sppal[15],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[16]),fill=sppal[16],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[17]),fill=sppal[17],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[18]),fill=sppal[18],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[19]),fill=sppal[19],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[20]),fill=sppal[20],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[21]),fill=sppal[21],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[22]),fill=sppal[22],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[23]),fill=sppal[23],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[24]),fill=sppal[24],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[25]),fill=sppal[25],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[26]),fill=sppal[26],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[27]),fill=sppal[27],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[28]),fill=sppal[28],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[29]),fill=sppal[29],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[30]),fill=sppal[30],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[31]),fill=sppal[31],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[32]),fill=sppal[32],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[33]),fill=sppal[33],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[34]),fill=sppal[34],xend=xend, xadd=xadd_var,size=hilsize) +
  yingtools2::geom_hilight(aes(isTip=isTip,var=species,value=sspcs[35]),fill=sppal[35],xend=xend, xadd=xadd_var,size=hilsize) +
  # geom_text(aes(label = node), check_overlap = F) #+  #find node(s) to flip() or ratate()
  geom_tippoint(aes(shape = cohort, fill = cohort), size= 1, alpha=0.87) + #tip size
  geom_tiplab2(aes(label = msk.id),
               size = 1.68, #tiplabel size
               hjust= -1, #tiplabel position to center
               align= F) +
  # scale_color_manual(values = c(sppal,E.coli = "#822E1C")) +
  # geom_tippoint(size=1,alpha=0.55) +
  xlim(c(-0.01,1.44)) +
  geom_rootpoint() +
  geom_treescale(width=0.2) +
  scale_shape_manual(values = c(21,24)) + #tip shapes
  theme(legend.position = "right") 
tree

pdf(file="../results/Fig1A.PhyloTree.Bacteroidales.pdf",height=12,width=15)
tree
dev.off()
