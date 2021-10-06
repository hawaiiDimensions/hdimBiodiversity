# packages
#library(tidyverse)
#library(vegan)
#library(readr)

# read in data
beating_abun <- read_csv("data/github/hdimbeat_abund.csv")

### Calculate Chao ###
######### sum by abundances at plot and pool by site #####
# pool by plot abundances
# point to poolPlot$site for pool
poolPlot <- beating_abun %>% 
  select(otu, abundance, site, sitePlot, smplUnit) %>%
  group_by(site, sitePlot, otu) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  pivot_wider(names_from = otu, values_from = abundance)

xPoolsum <- poolPlot %>%
  mutate_if(is.numeric, replace_na, replace = 0)
xPoolsum <- xPoolsum %>%
  select(-c(site, sitePlot))
xPoolsum <- as.matrix(as.data.frame(xPoolsum))

# chao calculated from site abundances summed
chao <- specpool(xPoolsum, pool = poolPlot$site, smallsample = TRUE)
chao <- chao %>%
  rownames_to_column(var = "site") %>%
  select(-c(jack1, jack1.se, jack2, boot, boot.se))

#### compare to no of distinct otus ######
divMetrics <- beating_abun %>%
  group_by(site, flowAge, volcano) %>% #
  # use otu
  summarise(n_zotus = n_distinct(otu)) %>% 
  left_join(chao) 

############# calculate inverse Simpson ##########
# first need matrix at site level
SiteAbun <- beating_abun %>% 
  select(otu, abundance, site, sitePlot, smplUnit) %>%
  group_by(site, otu) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  pivot_wider(names_from = otu, values_from = abundance)

SiteAbunMatrx <- SiteAbun %>%
  mutate_if(is.numeric, replace_na, replace = 0)
SiteAbunMatrx <- SiteAbunMatrx %>%
  select(-site)
SiteAbunMatrx <- as.matrix(as.data.frame(SiteAbunMatrx))

# Calculate inverse simpson based on abun of each otu at each site
invSimpSite <- diversity(SiteAbunMatrx, "inv")
divMetrics <- cbind(divMetrics, invSimpSite)
colnames(divMetrics)[9] <- "invSimp"

#write_csv(divMetrics, "data/github/divMetrics.csv")

divMetrics %>%
  pivot_longer(c(chao, n_zotus, invSimp), names_to = "metric") %>%
  ggplot() +
  geom_point(aes(x = log(flowAge), y = value, group = metric, colour = metric)) +
  geom_smooth(aes(x = log(flowAge), y = value, group = metric, colour = metric),
              method = "lm", se = TRUE) +
  xlab("Log Substrate Age (yrs)") +
  ylab("Richness") 
#save plot
#ggsave("data/github/compare_numzotu_chao_invS.png")

divMetrics %>%
  pivot_longer(c(chao, n_zotus, invSimp), names_to = "metric") %>%
  filter(metric == "invSimp") %>%
  ggplot() +
  geom_point(aes(x = log(flowAge), y = value, group = metric, colour = metric)) +
  geom_smooth(aes(x = log(flowAge), y = value, group = metric, colour = metric),
              method = "lm", se = FALSE) +
  xlab("Log Substrate Age (yrs)") +
  ylab("Richness") 
#save plot
#ggsave("data/github/invS.png")
divMetrics %>%
  pivot_longer(c(chao, n_zotus, invSimp), names_to = "metric") %>%
  filter(metric == "n_zotus") %>%
  ggplot() +
  geom_point(aes(x = log(flowAge), y = value, group = metric, colour = metric)) +
  geom_smooth(aes(x = log(flowAge), y = value, group = metric, colour = metric),
              method = "lm", se = FALSE) +
  xlab("Log Substrate Age (yrs)") +
  ylab("Richness")
#save plot
#ggsave("data/github/n_zotus.png")
