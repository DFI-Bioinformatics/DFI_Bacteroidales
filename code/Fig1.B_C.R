library(tidyverse)

# loading data ------------------------------------------------------------

sppal <- readRDS("data/sppal.extended.rds")

dat <- read_csv("data/genome.size.csv")

# GC plot -----------------------------------------------

dat %>% 
  ggplot(aes(GC,reorder(species, GC))) +
  geom_jitter(aes(shape = type, fill = species, size = type),
              alpha = 0.78,
              height = 0.25,
              width = 0) +
  scale_shape_manual(values = c(21, 24)) +
  scale_size_manual(values = c(3, 4.5)) +
  scale_fill_manual(values = sppal) +
  theme_bw() +
  labs(x = "GC%", y = "Species") +
  theme(legend.position = "none",
        legend.title=element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.text = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 13))

ggsave("results/Fig1_C.GCcontent.pdf", width = 6.8, height = 6.45)


# Genome size -------------------------------------------------------------

dat %>% 
  ggplot(aes(bp/1000000,reorder(species, bp))) +
  geom_jitter(aes(shape = type, fill = species, size = type),
              alpha = 0.78,
              height = 0.25,
              width = 0) +
  scale_shape_manual(values = c(21, 24)) +
  scale_size_manual(values = c(3, 4.5)) +
  scale_fill_manual(values = sppal) +
  theme_bw() +
  labs(x = "Genome size (Mbp)", y = "Species") +
  theme(legend.position = "none",
        legend.title=element_blank(),
        axis.text.y = element_text(face = "italic"),
        axis.text = element_text(size = 12),
        axis.title = element_text(face = "bold", size = 13))

ggsave("results/Fig1_B.genomeSize.pdf", width = 6.8, height = 6.45)
