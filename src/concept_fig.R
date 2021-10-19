## concept_fig.R
# To generate some images for Fig 1

## PACKAGES ============
#setwd("~/Dropbox/projects/2017/hawaiiCommunityAssembly/hdimBiodiversity/")
#devtools::install_github("ahasverus/orthomap")
library(orthomap)
library(rgdal)
library(ggplot2)
library(RColorBrewer)
library(ape)
library(ggtree)
library(rmapshaper)
library(sf)

## GENERATE SUBSTRATE AGE MAP ============
# USGS Substrate age map (https://pubs.usgs.gov/of/2007/1089/)
haw_geo <- readOGR("~/Dropbox/projects/2014/Hawaii/Maps/Haw_geo/Haw_St_geo_20070426_region.shp")
haw_geo <- spTransform(haw_geo, CRSobj = CRS("+proj=longlat")) # default datum is now wgs84
haw_geo_small <- ms_simplify(st_as_sf(haw_geo), keep_shapes = TRUE, keep = 0.05)
haw_union <- st_union(haw_geo_small)

## Define numeric ages to broad age categories
# See Table 6 of https://pubs.usgs.gov/of/2007/1089/Hawaii_metadata_for_users.pdf
# age_range <- data.frame(AGE_GROUP = 0:14,
#                         max_age = c(NA, 200, 750, 1500, 3000, 5000, 10E3, 30E3, 50E3, 140E3, 780E3, 1E6, 2E6, 4E6, 6E6),
#                         min_age = c(NA, 0, 200, 750, 1500, 3000, 5000, 10E3, 30E3, 50E3, 140E3, 780E3, 1E6, 2E6, 4E6))
# age_range$midpoint <- (age_range$max_age + age_range$min_age) / 2
# haw_geo@data$id <- rownames(haw_geo@data)
# haw_geo@data <- merge(x = haw_geo@data, y = age_range, all.x = TRUE)

color_ramp  <- colorRampPalette(brewer.pal(9, "YlGnBu"))
haw_geo_plot <- 
  ggplot() + 
    geom_sf(aes(fill = AGE_GROUP), data = haw_geo_small, size = 0.01, colour = NA) +
    geom_sf(data = haw_union, fill = NA, colour = "black", size = 0.2) +
    scale_fill_manual(values = c("grey90", color_ramp(14)),
                      breaks = c(0:14)) +
    # coord_fixed() + # coord_fixed does not work well with geom_sf
    theme(panel.background = element_rect(color = "black", fill = "white"),
          panel.grid = element_blank(),
          axis.text = element_blank(),
        axis.title = element_blank())
ggsave(haw_geo_plot, filename = "figures/haw_geo_plot.pdf")

## GENERATE GLOBE CENTERED ON HAWAII ==========
haw_centroid <- apply(haw_geo@bbox, MARGIN = 1, FUN = mean)

pdf("figures/world_map_hawaii.pdf")
orthomap(centre = c(haw_centroid[2], haw_centroid[1]), grid = TRUE, nx = 18, ny = 10, fill = "grey30", border.size = 3, grid.size = 1)
dev.off()

## GENERATE MOCK PHYLOGENIES ==========
set.seed(2)
tree <- rtree(30)

young_site_tip <- c("t14", "t22", "t17")
young_site <- rep("grey60", 2*30-1)
# let's pick t14, t22, t17
young_site[which(x$data$label %in% young_site_tip)] <- "black"

y <- ggtree(tree, col = young_site, size = 1)
head(y$data)
young_site_phylo <- y + geom_tippoint(size = 3, aes(color = factor(ifelse(label %in% young_site_tip, 1, 0)))) + 
  scale_colour_manual(values = c("grey60", "black")) +
  theme(legend.position = 0)
ggsave(young_site_phylo, filename = "figures/young_site_phylo.pdf")
