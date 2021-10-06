# sem.R
# Structural equation modelling
# Author: Jun Ying Lim

## PACKAGES =========
library(plyr)
library(ggplot2)
library(lavaan)
library(lavaanPlot)
library(DiagrammeRsvg)
#library(semPlot)

## IMPORT DATA =========
metadata <- read.csv("data/divMetrics.csv")
otutable <- read.csv("data/hdimbeat_abund.csv")

## CALCULATE PLANT DIVERSITY =========
# Calculate number of unique OTUs and plant genera sampled for each plot
plot_div <- ddply(.data = otutable, .variables = .(sitePlot), summarize,
      site = site[1],
      n_otu = length(unique(otu)),
      n_plant = length(unique(plant)))

plot_final <- merge(plot_div, metadata[c("site","flowAge", "volcano")])

## STRUCTURAL EQUATION MODELLING =========
# Standardize variables
plot_final$logflowAge  <- log(plot_final$flowAge)
plot_final$logflowAge_scale <- scale(plot_final$logflowAge)
plot_final$n_plant_scale <- scale(plot_final$n_plant)
plot_final$n_otu_scale <- scale(plot_final$n_otu)

# Some exploratory plots
plot(n_otu ~ n_plant, data = plot_final)
plot(n_plant ~ logflowAge, data = plot_final)
plot(n_otu ~ logflowAge, data = plot_final)

# Fit structural equation model
sem_mod <- sem("n_plant_scale ~ logflowAge_scale
            n_otu_scale ~ n_plant_scale + logflowAge_scale", data = plot_final)
summary(sem_mod, stand = T, rsq = T, fit.measures = T)

#semPaths(test, what = "stand", residuals = TRUE, intercepts = FALSE, rotation = 3,nCharNodes = 0)

# Plot 
sem_plot <- lavaanPlot(model = sem_mod, node_options = list(shape = "box", fontname = "Helvetica"),
                edge_options = list(color = "grey"), coefs = T, covs = TRUE, stars = "regress")
embed_plot_pdf(sem_plot, path = "figures/sem.pdf")
