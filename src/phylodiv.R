library(stringr)
library(readr)
library(magrittr)
library(tidyverse)
library(picante)
library(ggtree)
library(tidytree)
library(cowplot)

beating <- read_csv("data/beatingRedone.csv")
predict <- read_csv("data/niclassify-predictions_replacedfalseendemics.csv")

beating <- beating %>%
  left_join(predict) #%>% 
#  left_join(species) 


tree <- read.tree("data/postcutoff.tre")
wo_beating <- beating %>%
  filter(!site == "kohalaYng" & 
           !site == "LSAG" &
           !site == "hippnet") %>%
  # filter out the plots with more than one sd from mean of sampling time
  filter(!sitePlot == "flow1973_110" & 
           !sitePlot == "flow1973_107" &
           !sitePlot == "kohalaOld_05a" &
           !sitePlot == "waikamoi_06" &
           !sitePlot == "kamakou_15" &
           !sitePlot == "kokee_21a") %>%
  arrange(flowAge) %>%
  mutate(volcano = factor(volcano, unique(volcano))) %>%
  mutate(site = factor(site, unique(site))) %>%
  mutate(sitePlot = factor(sitePlot, unique(sitePlot)))

endemic <- wo_beating %>%
  filter(predict == "endemic")

commhdim.all <- endemic %>%
  group_by(volcano, site, flowAge, otu) %>%
  summarise(abundance = sum(abundance)) %>% # basically summing read abundances by site
  unite(siteInfo, 
        c("volcano", "site", "flowAge"), 
        sep = "_", remove = FALSE) %>%
  pivot_wider(id_cols = siteInfo, 
              names_from = otu, 
              values_from = abundance,
              values_fill = 0) # fills in zeros
commhdim.all <- column_to_rownames(commhdim.all, var = "siteInfo") # keep the row information as a row name not a column
# convert to to a matrix
commhdim.all <- as.matrix(as.data.frame(commhdim.all))

commhdim.plot <- endemic %>%
  group_by(sitePlot, otu) %>%
  summarise(abundance = sum(abundance)) %>% # basically summing read abundances by site
  pivot_wider(id_cols = sitePlot, 
              names_from = otu, 
              values_from = abundance,
              values_fill = 0) %>%
  column_to_rownames(var = "sitePlot")

siteAttributes <- data.frame(site = unique(str_split_fixed(endemic$sitePlot, n = 2, pattern = "_")[,1]),
                            flowAge = c(42, 133, 300, 475, 575, 7500, 20500, 365000, 545000, 1400000, 4150000))

# prune the tree based on the subsampled community
length(colnames(commhdim.all)) # 2736 OTUs
length(tree$tip.label) # 3785 tips
sum(colnames(commhdim.all) %in% tree$tip.label) # all 2736 OTUs are represented in tree
prunedphy.all <- prune.sample(commhdim.all, tree)
length(prunedphy.all$tip.label) # 2736 tips

# plot phylogeny
otu_phylo <- ggtree(prunedphy.all, size = 0.2, layout = "circular")
ggsave(otu_phylo, filename = "figures/otu_phylo.pdf")

# count number of OTUs
rowSums(decostand(commhdim.all, method = "pa"))

# plot site composition onto phylogeny
site_phylo <- list()
for(i in rownames(commhdim.all)){
  comm <- commhdim.all[i,]
  otus <- names(comm)[comm > 0]
  otu_phylo$data$test <- ifelse(otu_phylo$data$label %in% otus, "present", "absent") 
  site_phylo[[i]] <- otu_phylo + 
    geom_tippoint(aes(color = test), na.rm = TRUE, size = 0.5) +  # ggtree still plots even if it's NA
    scale_colour_manual(values = c("transparent", "red")) +
    ggtitle(i) +
    theme(legend.position = "none")
}
site_phylo_combined <- plot_grid(plotlist = site_phylo, nrow = 4, ncol= 3)
ggsave(site_phylo_combined, filename = "figures/comm_phylo.pdf", width = 12, height = 16)

# Make a phylogenetic distance matrix
phydist <- cophenetic(prunedphy.all)

# NTI and MTND
nti <- ses.mntd(decostand(commhdim.all, method = "pa"), phydist, null.model="richness", runs = 1000)
nti <- nti %>% rownames_to_column("siteInfo")
nti$flowAge <- as.numeric(str_split_fixed(nti$siteInfo, pattern = "_", n = 3)[,3])
# negative values of mntd.obs.z = species are more closely related than expected by random chance

ggplot(data = nti) + geom_point(aes(y = mntd.obs, x = log(flowAge), size = ntaxa))
ggplot(data = nti) + geom_point(aes(y = mntd.obs.z, x = log(flowAge), size = ntaxa))
ggplot(data = nti) + geom_point(aes(y = mntd.rand.mean, x = ntaxa))

rowSums(decostand(commhdim.all, method = "pa"))
rando_mat <- randomizeMatrix(decostand(commhdim.all, method = "pa"), null.model = "richness")
rowSums(rando_mat)
## QUESTIONS FOR NAT
# first part of code merges with "species" but I don't have access to it
# "endemic" does that include non-endemic natives?


## REDO ANALYSES FOR EACH PLOT SEPARATELY
nti_plot <- ses.mntd(decostand(commhdim.plot, method = "pa"), phydist, null.model = "richness", runs = 1000)
nti_plot <- rownames_to_column(nti_plot, var = "sitePlot")
nti_plot$site <- str_split_fixed(nti_plot$sitePlot, n = 2, pattern = "_")[,1]

nti_plot <- merge(nti_plot, siteAttributes, by = "site")
ggplot(aes(y = mntd.obs, x =log(flowAge)), data = nti_plot) + geom_point(aes(size = ntaxa)) + geom_smooth(method = "lm")

summary(lm(mntd.obs.z ~ log(flowAge), data = nti_plot))
cor.test(nti_plot$mntd.obs.z, nti_plot$flowAge, method = "spearman", exact = FALSE)
