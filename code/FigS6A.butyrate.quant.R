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

metab.quant <- 
  read_csv('../data/metab.quant.matrix.csv') 

GMM_pathway_results <- read_csv('../data/GMM_pathway_results.csv') %>% 
  mutate(Module_description = paste(Module, Description)) %>% 
  dplyr::select(-c('Module','Description')) %>% 
  pivot_wider(names_from = 'Module_description',
              values_from = 'pathway_completeness')

GMM_pathway_butyrate <- 
  GMM_pathway_results %>% 
  select(msk_id,contains('butyrate')) %>% 
  select(msk_id, `MF0089 butyrate production II`) %>% 
  na.omit()

quant.butyrate <- 
  metab.quant %>% 
  select( msk_id = msk.id, Butyrate) %>% 
  left_join(annotation) %>% 
  left_join(GMM_pathway_butyrate) %>% 
  mutate(butyrate.production.pathway.type = if_else(`MF0089 butyrate production II` == 1, 'Complete','Absent or Imcomplete')) %>% 
  distinct() %>% 
  na.omit()

# butyrate quant plot-----------

require(ggbreak)
require(egg)
ggplot(quant.butyrate, mapping = aes(x = butyrate.production.pathway.type,
                                     y = Butyrate) ) + 
  geom_jitter(mapping = aes(color = species), 
              width = 0.25, height = 0,
              size = 3) +
  scale_color_manual(values= sppal)+
  guides(color=guide_legend(ncol =1))+
  theme_presentation() + 
  ylim(-0.2, 5) +
  # break points
  # scale_y_break(c(3,6), scales = 'free') + 
  scale_y_break(c(0.5,1), scales = 'free') + 
  theme(axis.text.x = element_text(angle = 0)) +
  ylab('Butyrate production or consumption (mM)') +
  xlab('Butyrate production pathway') +
  theme(legend.text = element_text(face = 'italic')) # italicize legend text
ggsave(filename = '../results/FigS6A.butyrate.quant.pathway.jitterplot.pdf', 
       device = 'pdf', width = 10, height = 10)
dev.off()


