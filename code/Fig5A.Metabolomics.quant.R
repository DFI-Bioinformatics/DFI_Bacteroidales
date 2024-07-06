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

metab.quant.matrix <- 
  read_csv('../data/metab.quant.matrix.csv') %>% 
  column_to_rownames(var = 'msk.id') %>% 
  as.matrix()


### plot quant matrix --------

# row annotations
row.annot.quant <- annotation %>%
  column_to_rownames(var = 'msk_id')

pheatmap(metab.quant.matrix,
         cluster_cols = F, 
         cluster_rows = F,
         annotation_row = row.annot.quant,
         annotation_colors = list(species = c(
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
         # main = 'SCFA production or consumption',
         fontsize_row = 3, # row names (msk_id) size
         color = c(colorRampPalette(rev(brewer.pal(n = 7, name ="Blues")))(5), colorRampPalette(brewer.pal(n = 7, name ="Reds"))(30)), breaks = -5:30, #customize heatmap color ramps
         display_numbers = T, number_format = '%.1f', number_color = 'black', fontsize_number = 3, #format numbers within each cell
         filename = '../results/Fig5A.quant.metabolomics.pdf',
         cellwidth = 20, cellheight = 5, # adjust size of individual cells
         height = 10, width = 6
)
