rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(magrittr)

'%!in%' <- function(x,y)!('%in%'(x,y))


# read data --------------
sppal <- readRDS('../data/sppal.extended.rds')

annotation <- 
  read_csv('../data/genome.size.csv') %>% 
  select(species, msk_id)
  
metab.qual.matrix.log2 <- 
  read_csv('../data/metab.qual.matrix.log2.csv') %>% 
  column_to_rownames(var = 'msk.id') %>% 
  as.matrix()
  

### plot matrix --------

require(pheatmap)
require(RColorBrewer)

# row annotations
row.annot.qual <- annotation %>% 
  column_to_rownames(var = 'msk_id')

pheatmap(metab.qual.matrix.log2,
         cluster_cols = T, 
         # cluster_rows = T,
         cluster_rows = F,
         annotation_row = row.annot.qual, annotation_colors = list(species = c(
           `Alistipes communis` = "#666666",
           `Alistipes finegoldii` = "#FBB4AE",
           `Alistipes onderdonkii` = "#FB8072",
           `Alistipes putredinis` = "#A65628",
           `Alistipes shahii` = "#A6761D",
           `Bacteroides caccae` = "#FDC086",
           `Bacteroides cellulosilyticus` = "#FF7F00",
           `Bacteroides eggerthii` = "#B3CDE3",
           `Bacteroides faecis` = "#D95F02",
           `Bacteroides finegoldii` = "#999999",
           `Bacteroides fragilis` = "#E6AB02",
           `Bacteroides intestinalis` = "#FFD92F",
           `Bacteroides nordii` = "#BEBADA",
           `Bacteroides ovatus` = "#FFFF33",
           `Bacteroides salyersiae` = "#A6D854",
           `Bacteroides stercoris` = "#4DAF4A",
           `Bacteroides thetaiotaomicron` = "#7FC97F",
           `Bacteroides uniformis` = "#1B9E77",
           `Bacteroides xylanisolvens` = "#8DD3C7",
           # `Butyricimonas paravirosa` = "#1F78B4", # no metabolomic data
           `Butyricimonas synergistica` = "#6A3D9A",
           # `Hoylesella buccalis` = "#80B1D3", # no metabolomic data
           `Leyella stercorea` = "#7570B3",
           `Odoribacter splanchnicus` = "#8DA0CB",
           `Parabacteroides distasonis` = "#984EA3",
           `Parabacteroides johnsonii` = "#CBD5E8",
           `Parabacteroides merdae` = "#BEAED4",
           `Phocaeicola coprocola` = "#BC80BD",
           `Phocaeicola dorei` = "#E78AC3",
           `Phocaeicola massiliensis` = "#E7298A",
           `Phocaeicola vulgatus` = "#FCCDE5",
           `Segatella copri` = "#E31A1C",
           `Segatella sinica` = "#FFED6F") #sppal.extended
         ),
         # main = 'Metabolites log2(fold change)', 
         fontsize_row = 5, fontsize_col = 7,
         color = c(colorRampPalette(rev(brewer.pal(n = 7, name ="Blues")))(6), 
                   colorRampPalette(brewer.pal(n = 7, name ="Reds"))(12)), breaks = -6:12, #customize heatmap color ramps
         filename = '../results/Fig5B.semiquant.metabolomics.log2.fold.change.pdf',
         height = 10, width = 8
)

