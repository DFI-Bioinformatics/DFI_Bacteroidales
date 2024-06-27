# Fig 4 - C number of non-susD-susC PUL numbers

sppal <- readRDS("data/sppal.extended.rds")

tax <- read_csv("data/genome.size.csv") %>% 
  select(msk_id, species, type)

puls <- readRDS("data/pulpy.rds")

tax %>%
  distinct(msk_id) %>% 
  left_join(puls %>%
              filter(!pattern %in% c("susC-susD","susD-susC")) %>%
              mutate(length = end - start + 1)  %>%
              full_join(tax %>% select(species, type, msk_id)) %>% #view
              group_by(species, type, msk_id) %>%
              summarise(min = min(length),
                        aveLen = mean(length),
                        median = median(length),
                        sd = sd(length),
                        max = max(length),
                        nContigs = n()) ) %>%
  mutate(nContigs = ifelse(is.na(max), 0, nContigs)) %>% #view
  ggplot(aes(x=reorder(species,nContigs),y=nContigs)) +
  theme_bw() +
  geom_jitter(aes(fill=species,shape=type), height = 0,
              size=4,alpha=0.8) +
  scale_fill_manual(values=sppal) +
  theme(legend.title=element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.text.y = element_text(face = "italic"),
        axis.title = element_text(face = "bold", size = 13)
  ) +
  coord_flip() +
  scale_shape_manual(values=c(21,24), guide=FALSE) +
  labs(y = "Number of PULs", x = "Species") 

ggsave("results/Fig4_C.PULnoSusD-C.pdf", width = 6.8, height = 6.45)

