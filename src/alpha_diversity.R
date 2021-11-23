# sem.R
# Structural equation modelling
# Author: Jun Ying Lim

## PACKAGES =========
library(plyr)
library(ggplot2)
library(lavaan)
library(lavaanPlot)
library(DiagrammeRsvg)
library(vegan)
library(MuMIn)
library(lme4)
library(ggeffects) #https://strengejacke.github.io/ggeffects/articles/introduction_randomeffects.html
#library(semPlot)

## IMPORT DATA =========
metadata <- read.csv("data/divMetrics.csv")
otutable <- read.csv("data/hdimbeat_abund.csv")


## GRAPHICAL PARAMETERS =========
berkcol3 = c("#003262", "#3B7EA1", "#9BBEA9", "#00B0DA", "#00A598", brewer.pal(9, "BuGn")[8], "#CFDD45",
             "#859438", "#FDB515", brewer.pal(9, "YlOrRd")[5], "#ED4E33", "#C4820E", "#D9661F",
             "#6C3302")
point_theme <- theme(panel.background = element_rect(color = "black", fill = NA),
                     panel.grid = element_blank(),
                     plot.background = element_blank())

## CALCULATE PLANT DIVERSITY =========
# Calculate number of unique OTUs and plant genera sampled for each plot
plot_div <- ddply(.data = otutable, .variables = .(sitePlot), summarize,
      site = site[1],
      n_otu = length(unique(otu)),
      inv_simp = 1 / sum((abundance / sum(abundance))^2),
      n_plant = length(unique(plant)))

plot_div_final <- merge(plot_div, metadata[c("site","flowAge", "volcano")])

# Calculate site by order 
length(table(otutable$orderGrp))
plot_order_div <- ddply(.data = otutable, .variables = .(sitePlot, orderGrp), summarize,
                        site = site[1],
                        n_otu = length(unique(otu)),
                        inv_simp = 1 / sum((abundance / sum(abundance))^2),
                        n_plant = length(unique(plant)))
plot_order_div_final <- merge(plot_order_div, metadata[c("site","flowAge", "volcano")])
plot_order_div_final$logFlowAge <- log(plot_order_div_final$flowAge)
plot_order_div_final$inv_simp_scale <- scale(plot_order_div_final$inv_simp)
plot_order_div_final$logOTU_scale <- scale(plot_order_div_final$n_otu)

## ALPHA DIVERSITY MODELLING =========
age_bounds <- c(min(plot_order_div_final$logFlowAge), max(plot_order_div_final$logFlowAge))
orderGrp <- unique(plot_order_div_final$orderGrp)

# Inverse Simpson mixed effect model
inv_simp_order_mod <- lmer(inv_simp_scale ~ logFlowAge + (logFlowAge|orderGrp), data = plot_order_div_final) 
hist(residuals(inv_simp_order_mod))
plot(fitted(inv_simp_order_mod) ~residuals(inv_simp_order_mod))
plot(fitted(inv_simp_order_mod) ~ scale(plot_order_div_final$inv_simp))
abline(0,1, col = "red")

r.squaredGLMM(inv_simp_order_mod) # marginal = 0.084, conditional = 0.594
random_slopes <- data.frame(coef(inv_simp_order_mod)$orderGrp)
random_slopes$orderGrp <- rownames(coef(inv_simp_order_mod)$orderGrp)

random_ends <- predict(inv_simp_order_mod, newdata = data.frame(logFlowAge = rep(age_bounds, 14), 
                                                                orderGrp = rep(orderGrp, each = 2)))
random_ends2 <- data.frame("inv_simp_scale" = random_ends,
                           logFlowAge = rep(age_bounds, 14), orderGrp = rep(orderGrp, each = 2))


fixed_effect <- ggeffects::ggpredict(inv_simp_order_mod, "logFlowAge", type = "fixed")

inv_simp_plot <- ggplot() +
  geom_point(aes(y = inv_simp_scale, x = logFlowAge, colour = orderGrp),
             position = "jitter", pch = 16, size = 2, alpha = 0.5, data = plot_order_div_final) +
  #geom_abline(aes(intercept = X.Intercept., slope = logFlowAge, color = orderGrp), data = random_slopes, alpha = 0.7) +
  #geom_path(aes(y = inv_simp_scale, x = logFlowAge, color = orderGrp), data = random_ends2) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, x = x), data = fixed_effect, alpha = 0.3) +
  geom_path(aes(y = predicted, x = x), colour = "black", size = 1, data = fixed_effect) +
  scale_y_continuous(limits = c(-1, 6)) +
  labs(y = "Inverse Simpson (scaled)", x = "Substrate age (Log Ma)") +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  point_theme +
  theme(legend.title = element_blank(),
        legend.key= element_blank()) +
  scale_color_manual(values = berkcol3)
ggsave(inv_simp_plot, filename = "figures/inv_simpson.pdf", width = 6, height = 4)

# OTU richness mixed effect model
otu_div_order_mod <- lmer(logOTU_scale ~ logFlowAge + (logFlowAge|orderGrp), data = plot_order_div_final) 
hist(residuals(otu_div_order_mod))
plot(fitted(otu_div_order_mod) ~residuals(otu_div_order_mod))
plot(fitted(otu_div_order_mod) ~ scale(plot_order_div_final$logOTU_scale))
abline(0,1, col = "red")

r.squaredGLMM(otu_div_order_mod) # marginal = 0.032, conditional = 0.671
random_slopes <- data.frame(coef(otu_div_order_mod)$orderGrp)
random_slopes$orderGrp <- rownames(coef(otu_div_order_mod)$orderGrp)

random_ends <- predict(otu_div_order_mod, newdata = data.frame(logFlowAge = rep(age_bounds, 14), 
                                                                orderGrp = rep(orderGrp, each = 2)))
random_ends2 <- data.frame("logOTU_scale" = random_ends,
                           logFlowAge = rep(age_bounds, 14), orderGrp = rep(orderGrp, each = 2))
fixed_effect <- ggeffects::ggpredict(otu_div_order_mod, "logFlowAge", type = "fixed")

otu_div_plot <- ggplot() +
  geom_point(aes(y = logOTU_scale, x = logFlowAge, colour = orderGrp),
             position = "jitter", pch = 16, size = 2, alpha = 0.5, data = plot_order_div_final) +
  geom_ribbon(aes(ymax = conf.high, ymin = conf.low, x = x), data = fixed_effect, alpha = 0.3) +
  geom_path(aes(y = predicted, x = x), colour = "black", size = 1, data = fixed_effect) +
  scale_y_continuous(limits = c(-1, 6)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(y = "Log OTU richness (scaled)", x = "Substrate age (Log Ma)") +
  scale_color_manual(values = berkcol3) +
  point_theme +
  theme(legend.title = element_blank(),
        legend.key= element_blank())
  
ggsave(otu_div_plot, filename = "figures/otu_div.pdf", width = 6, height = 4)


# Combine plots
library(cowplot)
alpha_div_plot_leg <- get_legend(otu_div_plot+theme(legend.position = "right"))

alpha_div_plot <- plot_grid(plot_grid(otu_div_plot + theme(legend.position = "none"),
                                      inv_simp_plot + theme(legend.position = "none"), 
                                      labels = "auto"),
                            alpha_div_plot_leg,
                            nrow = 1, rel_widths= c(0.9, 0.2))
ggsave(alpha_div_plot, filename = "figures/alpha_div.pdf", height = 5, width = 11)


## STRUCTURAL EQUATION MODELLING =========
# All orders combined
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

# SEM with order as random effect on slope
library(piecewiseSEM)
plot_order_div_final$n_plant_scale <- scale(plot_order_div_final$n_plant)
plot_order_div_final$logFlowAge_scale <- scale(plot_order_div_final$logFlowAge)

sem_mixed  <- psem(PlantDiv = lmer(n_plant_scale ~ logFlowAge_scale + (logFlowAge_scale|orderGrp), data = plot_order_div_final),
                             N_OTU = lmer(logOTU_scale ~ n_plant_scale + logFlowAge_scale + (logFlowAge_scale | orderGrp), data = plot_order_div_final),
                             data = plot_order_div_final)
