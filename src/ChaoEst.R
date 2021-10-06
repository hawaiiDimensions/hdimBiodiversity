# packages
#library(tidyverse)
#library(vegan)
#library(readr)

# read in data
beating_abun <- read_csv("data/github/hdimbeat_abund.csv")

######### sum by abundances at smplUnit and pool by site #####
# note that abundances are "per otu per sampling unit" (how much of each otu was found on a plant genus within a plot)

# make simple df
# change to wider format and save as pool data 
# for pool data point to pool$site
pool <- beating_abun %>% 
  select(otu, abundance, site, sitePlot, smplUnit) %>%
  pivot_wider(names_from = otu, values_from = abundance) 

# trim pool to matrix format to supply species data
x <- pool %>%
  mutate_if(is.numeric, replace_na, replace = 0)
x <- x %>%
  select(-c(site, sitePlot, smplUnit))
x <- as.matrix(as.data.frame(x))

# chao calculated from sample unit abundances summed
chao_fromSmplsum <- specpool(x, pool = pool$site, smallsample = TRUE)
chao_fromSmplsum

# note when abundacnes are grouped by smplUnit can also calculate by plot
#chaoPlot <- specpool(x, pool = pool$sitePlot, smallsample = TRUE)
#chaoPlot

######### sum by abundances at plot and pool by site #####
# pool by site abundances
# point to poolSite$site for pool
poolSite <- beating_abun %>% 
  select(otu, abundance, site, sitePlot, smplUnit) %>%
  group_by(site, sitePlot, otu) %>%
  summarise(abundance = sum(abundance)) %>%
  ungroup() %>%
  pivot_wider(names_from = otu, values_from = abundance)

xSitesum <- poolSite %>%
  mutate_if(is.numeric, replace_na, replace = 0)
xSitesum <- xSitesum %>%
  select(-c(site, sitePlot))
xSitesum <- as.matrix(as.data.frame(xSitesum))

# chao calculated from site abundances summed
chao_fromSiteSum <- specpool(xSitesum, pool = poolSite$site, smallsample = TRUE)
chao_fromSiteSum

#write_csv(chao_fromSmplsum, "data/github/chao_fromSmplUnitSum.csv")
#write_csv(chao_fromSiteSum, "data/github/chao_fromSiteSum.csv")

#### compare to no of distinct otus #######
chaoSitesum_forbind <- rownames_to_column(chao_fromSiteSum, var = "site") %>%
  rename(chaoSitesum = chao, chaoSitesum.se =  chao.se) %>%
  select(site, Species, chaoSitesum, chaoSitesum.se)

chao_fromSmplsum_forbind <- rownames_to_column(chao_fromSmplsum, var = "site") %>%
  rename(chaoSmplsum = chao, chaoSmplsum.se =  chao.se) %>%
  select(site, Species, chaoSmplsum, chaoSmplsum.se)

chaoPlotdf <- beating_abun %>%
  group_by(site, flowAge, volcano) %>% #
  # use otu
  summarise(n_species = n_distinct(otu)) %>% 
  left_join(chaoSitesum_forbind) %>%
  left_join(chao_fromSmplsum_forbind)

chaoPlotdf %>%
  ggplot() +
  geom_point(aes(x = log(flowAge), y = n_species), colour = "red") + 
  stat_smooth(aes(x = log(flowAge), y = n_species), method = "lm", size = 1, se = FALSE, colour = "red") +
  geom_text(aes(x = 9, y = 500), label = "N zOTUs", color = "red") +
  geom_point(aes(x = log(flowAge), y = chaoSitesum), colour = "green") + 
  stat_smooth(aes(x = log(flowAge), y = chaoSitesum), method = "lm", size = 1, se = FALSE, colour = "green") +
  geom_text(aes(x = 5, y = 1300), label = "Sitesum Chao", color = "green") +
  geom_point(aes(x = log(flowAge), y = chaoSmplsum), colour = "blue") + 
  stat_smooth(aes(x = log(flowAge), y = chaoSmplsum), method = "lm", size = 1, se = FALSE, colour = "blue")+
  geom_text(aes(x = 6, y = 1700), label = "Smplsum Chao", color = "blue") +
  xlab("Log Substrate Age (yrs)") +
  ylab("Richness") 

#write_csv(chaoPlotdf, "data/github/chaoPlotdf.csv")
#ggsave("data/github/choaCompareN.png")