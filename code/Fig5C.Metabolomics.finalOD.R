rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(tidyverse)
library(magrittr)

'%!in%' <- function(x,y)!('%in%'(x,y))


# read data --------------

final.OD.average <- readxl::read_xlsx('../data/Final_ODs.xlsx',
                                      sheet = 'average',
                                      col_types = c('text','text','numeric'))
sppal <- readRDS('../data/sppal.extended.rds')

## plot Final OD600------------------


# arrange Taxonomy
final.OD.plot <- final.OD.average %>% 
  dplyr::rename(Normalized.Average = `Normalized \r\nAverage`) %>% 
  arrange(Taxonomy, `Isolate ID`) %>% # Arrange dataframe by Taxonomy and Isolate ID
  mutate(`Isolate ID` = factor(`Isolate ID`, levels = rev(unique(`Isolate ID`)))) # Convert Isolate ID to factor with levels sorted by Taxonomy

# Summarize data to calculate the average measurement for each Isolate ID and Taxonomy
final.OD.summary <- final.OD.average %>%
  group_by(`Isolate ID`, Taxonomy) %>%
  dplyr::rename(Normalized.Average = `Normalized \r\nAverage`) %>% 
  summarize(average_measurement = mean(Normalized.Average, na.rm = TRUE), .groups = 'drop')%>% 
  arrange(Taxonomy, `Isolate ID`) %>% # Arrange dataframe by Taxonomy and Isolate ID
  mutate(`Isolate ID` = factor(`Isolate ID`, levels = rev(unique(`Isolate ID`))))# Convert Isolate ID to factor with levels sorted by Taxonomy

ggplot()+
  geom_col(data = final.OD.summary, 
           aes(x = `Isolate ID`, y = average_measurement, fill = Taxonomy))+
  scale_fill_manual(values = sppal)+
  geom_point(data = final.OD.plot, 
             aes(x = `Isolate ID`, y = Normalized.Average), 
             position = position_jitter(width = 0.2))+
  scale_y_continuous(breaks = c(0, 0.3, 0.6, 0.9, 1.2)) + # Set y-axis (future x-axis) breaks
  coord_flip() + # Flip coordinates for a vertical plot
  theme_bw() + # Use a white background
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +# Rotate x-axis labels for better readability
  guides(fill = guide_legend(ncol = 1), color = guide_legend(ncol = 1), title = 'Species')+ # Set legend to one column
  labs(fill = 'Species',
       y = expression(OD[600]),
       x = 'Isolate') # Set x-axis title with subscript, set y-axis title

ggsave(filename = '../results/Fig5C.Bacteroidales.final.OD.pdf',
       width = 6, height = 20)
dev.off()


