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

IndicesRichness_Combinations <- expand.grid(index = grep("mean", colnames(acousticIndices_richness), value = TRUE),
                                            taxa = unique(acousticIndices_richness$type),
                                            diversity = c("richness", "shannon", "count"))

for (combination in 1:nrow(IndicesRichness_Combinations)) {
  currentAcousticIndex <- IndicesRichness_Combinations$index[combination]
  currentTaxa <- IndicesRichness_Combinations$taxa[combination]
  currentDiversity <- IndicesRichness_Combinations$diversity[combination]
  
  data <- acousticIndices_richness %>% filter(type == currentTaxa)
  
  #Mixed-effects model with 'Site' as random effect
  Models_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]] <- lmer(as.formula(paste0(currentDiversity, " ~ ", currentAcousticIndex, " + (1|Site/Sensor) + (1|sampling.period)")), 
                                                                                                                                  data = data)
  
  #Plot and save visual check of model assumptions
  png(filename = paste0("./outputs/model.checks/", currentDiversity, "_", currentAcousticIndex, "_", currentTaxa, ".png"), width = 30, height = 20, units = "cm", res = 800)
  print(performance::check_model(Models_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]]))
  dev.off()
  
  #Model-based (Semi-)Parametric Bootstrap for Mixed Models - Used to calculate confidence intervals
  Bootstrap_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]] <- bootMer(Models_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]],
                                                                                                                 FUN = function(x)predict(x, re.form=NA),
                                                                                                                 nsim = 1000)
  
  data$lci <- apply(Bootstrap_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]]$t, 2, quantile, 0.025)
  data$uci <- apply(Bootstrap_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]]$t, 2, quantile, 0.975)
  data$pred <- predict(Models_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]], re.form=NA)
  
  ylimits <- if(currentTaxa == 'birds' & currentDiversity == 'richness') {c(10, 50)} else 
    if(currentTaxa == 'frogs' & currentDiversity == 'richness') {c(0, 12)} else
      if(currentTaxa == 'birds' & currentDiversity == 'shannon') {c(1.8, 4.2)} else 
        if(currentTaxa == 'frogs' & currentDiversity == 'shannon') {c(0, 1.2)} else
          if(currentTaxa == 'birds' & currentDiversity == 'count') {c(0, 500)} else 
            if(currentTaxa == 'frogs' & currentDiversity == 'count') {c(0, 60)} 
  
  
  #Plot of fixed effects and confidence intervals
  Plots_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]] <- data %>% 
    ggplot(aes_string(x = paste0(currentAcousticIndex), y = paste0(currentDiversity))) + 
    geom_ribbon(aes_string(x = paste0(currentAcousticIndex), ymin = "lci", ymax = "uci"), fill = "black", alpha = 0.1) +
    geom_line(aes_string(x = paste0(currentAcousticIndex), y = "pred"), color = "black", lwd = 1) +
    geom_point(size = 2) +
    scale_y_continuous(limits = ylimits) + 
    annotate(geom = "text", -Inf, Inf, hjust = -0.1, vjust = 1.8, parse = T,
             label = paste0("R^2 == ", format(round(as.numeric(performance::r2(Models_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]])[[2]]), 2), nsmall = 2))) +
    labs(x = gsub("_mean", "", paste0(currentAcousticIndex)), y = "Diversity Measure") +
    theme_classic() +
    theme(axis.title.y = element_blank())
}


# Arrange multi-panel plots -----------------------------------------------

plot_combinations <- unique(IndicesRichness_Combinations[,c('diversity','taxa')])

for (combination in 1:nrow(plot_combinations)) {
  
  ytitle <- if(plot_combinations$diversity[combination] == 'richness') {c("Species Richness")} else
    if(plot_combinations$diversity[combination] == 'shannon') {c("Shannon Diversity")} else
      if(plot_combinations$diversity[combination] == 'count') {c("Total Abundance")}
  
  Plot <- plot_grid(plotlist = Plots_Indices_Richness[grep(paste0(plot_combinations$diversity[combination], ".*", plot_combinations$taxa[combination]), names(Plots_Indices_Richness))],
                    ncol = 4) %>% 
    annotate_figure(top = paste0("lmer - AcousticIndices_", plot_combinations$taxa[combination], "_", plot_combinations$diversity[combination]),
                    left = ytitle,
                    bottom = "Acoustic Index")
  
  ggsave(filename = paste0("./outputs/figures.local/", Sys.Date(), "Plots_", plot_combinations$taxa[combination], "_", plot_combinations$diversity[combination], ".png"),
         plot = Plot,
         width = 30, height = 25, units = "cm", dpi = 1200)
}

#Save plots of specific indices that we are going to analyse
indicesToUse <- c("ACI", 
                  "ADI", 
                  "AEI", 
                  "BI", 
                  "NDSI", 
                  "EVN")

for (combination in 1:nrow(plot_combinations)) {
  
  ytitle <- if(plot_combinations$diversity[combination] == 'richness') {c("Species Richness")} else
    if(plot_combinations$diversity[combination] == 'shannon') {c("Shannon Diversity")} else
      if(plot_combinations$diversity[combination] == 'count') {c("Total Abundance")}
  
  Plot <- plot_grid(plotlist = Plots_Indices_Richness[grep(paste0(plot_combinations$diversity[combination], "_(", paste0(indicesToUse, collapse = "|"), ")_.*", plot_combinations$taxa[combination]), names(Plots_Indices_Richness), value = T)],
                    ncol = 3) %>% 
    annotate_figure(top = paste0("lmer - AcousticIndices_", plot_combinations$taxa[combination], "_", plot_combinations$diversity[combination]),
                    left = ytitle,
                    bottom = "Acoustic Index")
  
  ggsave(filename = paste0("./outputs/figures/", plot_combinations$taxa[combination], "_", plot_combinations$diversity[combination], ".png"),
         plot = Plot,
         width = 22.5, height = 12.5, units = "cm", dpi = 1200)
}
