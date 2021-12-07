#Correlate acoustic indices with vertebrate diversity - 2021

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(lme4)
library(lmerTest)
library(boot)


# Load indices and richness data, merge -----------------------------------

setwd("C:/Users/jc696551/OneDrive - James Cook University/Projects/acousticindices_vertebratediversity")

load("./outputs/data/2021-12-06_acousticIndices_summary.RData") #load indices

richness <- read_csv("./rawdata/biodiversity/richness.csv") #load richness

#format richness to match indices
richness <- richness %>% rename(Site = site, Sensor = plot)
richness$Sensor <- gsub(" ", "", richness$Sensor)

#merge richness and acoustic indices
acousticIndices_richness <- merge(acousticIndices_summary, richness[richness$day.date == 'all',], by = c("type", "Site", "Sensor", "sampling.period"))


# Mixed-effects models ----------------------------------------------------

Plots_Indices_Richness <- list()
Models_Indices_Richness <- list()
Bootstrap_Indices_Richness <- list()

IndicesRichness_Combinations <- expand.grid(Indices = grep("mean", colnames(acousticIndices_richness), value = TRUE),
                                            type = unique(acousticIndices_richness$type))

for (combination in 1:nrow(IndicesRichness_Combinations)) {
  currentAcousticIndex <- IndicesRichness_Combinations$Indices[combination]
  currentRichness <- IndicesRichness_Combinations$type[combination]
  
  data <- acousticIndices_richness %>% filter(type == currentRichness)
  
  #Mixed-effects model with 'Site' as random effect
  Models_Indices_Richness[[paste0(currentAcousticIndex, "_", currentRichness)]] <- lmer(as.formula(paste0("richness ~ ", currentAcousticIndex, " + (1|Site/Sensor) + (1|sampling.period)")), 
                                                                                        data = data)
  
  #Plot and save visual check of model assumptions
  png(filename = paste0("./outputs/model.checks/", currentAcousticIndex, "_", currentRichness, ".png"), width = 30, height = 20, units = "cm", res = 800)
  performance::check_model(Models_Indices_Richness[[paste0(currentAcousticIndex, "_", currentRichness)]])
  dev.off()
  
  #Model-based (Semi-)Parametric Bootstrap for Mixed Models - Used to calculate confidence intervals
  Bootstrap_Indices_Richness[[paste0(currentAcousticIndex, "_", currentRichness)]] <- bootMer(Models_Indices_Richness[[paste0(currentAcousticIndex, "_", currentRichness)]],
                                                                                              FUN = function(x)predict(x, re.form=NA),
                                                                                              nsim = 1000)
  
  data$lci <- apply(Bootstrap_Indices_Richness[[paste0(currentAcousticIndex, "_", currentRichness)]]$t, 2, quantile, 0.025)
  data$uci <- apply(Bootstrap_Indices_Richness[[paste0(currentAcousticIndex, "_", currentRichness)]]$t, 2, quantile, 0.975)
  data$pred <- predict(Models_Indices_Richness[[paste0(currentAcousticIndex, "_", currentRichness)]], re.form=NA)
  
  ylimits <- if(currentRichness == 'birds') {c(10, 50)} else {c(0, 12)}
  
  #Plot of fixed effects and confidence intervals
  Plots_Indices_Richness[[paste0(currentAcousticIndex, "_", currentRichness)]] <- data %>% 
    ggplot(aes_string(x = paste0(currentAcousticIndex), y = "richness")) + 
    geom_ribbon(aes_string(x = paste0(currentAcousticIndex), ymin = "lci", ymax = "uci"), fill = "black", alpha = 0.1) +
    geom_line(aes_string(x = paste0(currentAcousticIndex), y = "pred"), color = "black", lwd = 1) +
    geom_point(size = 2) +
    scale_y_continuous(limits = ylimits) + 
    annotate(geom = "text", -Inf, Inf, hjust = -0.1, vjust = 1.8, parse = T,
             label = paste0("R^2 == ", format(round(as.numeric(performance::r2(Models_Indices_Richness[[paste0(currentAcousticIndex, "_", currentRichness)]])[[2]]), 2), nsmall = 2))) +
    labs(x = gsub("_mean", "", paste0(currentAcousticIndex)), y = "Species Richness") +
    theme_classic() +
    theme(axis.title.y = element_blank())
}


# Arrange multi-panel plots -----------------------------------------------

#Plot acoustic indices vs bird richness
Plot_Indices_Birds <- plot_grid(plotlist = Plots_Indices_Richness[grep("birds", names(Plots_Indices_Richness))],
                                ncol = 4) %>% annotate_figure(top = "lmer - AcousticIndicesMorning_BirdRichness")

Plot_Indices_Birds <- annotate_figure(Plot_Indices_Birds,
                                      left = "Species Richness",
                                      bottom = "Acoustic Index")

ggsave(filename = paste0("./outputs/figures/", Sys.Date(), "Plots_BirdRichness_MorningIndicesMean_lmer.png"),
       plot = Plot_Indices_Birds,
       width = 30, height = 25, units = "cm", dpi = 1200)

#Plot acoustic indices vs frog richness
Plot_Indices_Frogs <- plot_grid(plotlist = Plots_Indices_Richness[grep("frogs", names(Plots_Indices_Richness))],
                                ncol = 4) %>% annotate_figure(top = "lmer - AcousticIndicesEvening_FrogRichness")

Plot_Indices_Frogs <- annotate_figure(Plot_Indices_Frogs,
                                      left = "Species Richness",
                                      bottom = "Acoustic Index")

ggsave(filename = paste0("./outputs/figures/", Sys.Date(), "Plots_FrogRichness_EveningIndicesMean_lmer.png"),
       plot = Plot_Indices_Frogs,
       width = 30, height = 25, units = "cm", dpi = 1200)