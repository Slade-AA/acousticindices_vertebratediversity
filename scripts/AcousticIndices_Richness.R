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
  Models_Indices_Richness[[paste0(currentAcousticIndex, "_", currentRichness)]] <- lmer(as.formula(paste0("richness ~ ", currentAcousticIndex, " + (1|Site)")), 
                                                                                        data = data)
  
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
    theme_classic()
}

#Plot acoustic indices vs bird richness
Plot_Indices_Birds <- plot_grid(plotlist = Plots_Indices_Richness[grep("birds", names(Plots_Indices_Richness))],
                                ncol = 4) %>% annotate_figure(top = "lmer - AcousticIndicesMorning_BirdRichness")

annotate_figure(Plot_Indices_Birds,
                left = "Species Richness",
                bottom = "Acoustic Index")

ggsave(paste0("./outputs/figures/", Sys.Date(), "Plots_BirdRichness_MorningIndicesMean_lmer.png"),
       width = 30, height = 18, units = "cm", dpi = 1200)

#Plot acoustic indices vs frog richness
Plot_Indices_Frogs <- plot_grid(plotlist = Plots_Indices_Richness[grep("frogs", names(Plots_Indices_Richness))],
                                ncol = 4) %>% annotate_figure(top = "lmer - AcousticIndicesEvening_FrogRichness")

ggsave(paste0("./outputs/figures/", Sys.Date(), "Plots_FrogRichness_EveningIndicesMean_lmer.png"),
       width = 30, height = 18, units = "cm", dpi = 1200)








#Mixed-effects model plots
Plots_SensorTotals_MorningIndicesMean_lmer <- list()
model_lmer <- list()
bb.lmer <- list()
for (indice in grep("mean", colnames(IndicesSummary_Morning))) {
  currentAcousticIndex <- gsub("_mean", "", colnames(IndicesSummary_Morning[,indice]))
  
  #Mixed-effects model with 'Site' as random effect
  model_lmer[[currentAcousticIndex]] <- lmer(as.formula(paste0("n ~ ", colnames(IndicesSummary_Morning[,indice]), " + (1|Site)")), 
                                             data = Counts_SiteSensor_IndicesMorning)
  
  #Model-based (Semi-)Parametric Bootstrap for Mixed Models - Used to calculate confidence intervals
  bb.lmer[[currentAcousticIndex]] <- bootMer(model_lmer[[currentAcousticIndex]],
                     FUN = function(x)predict(x, re.form=NA),
                     nsim = 1000)
  
  Counts_SiteSensor_IndicesMorning$lci <- apply(bb.lmer[[currentAcousticIndex]]$t, 2, quantile, 0.025)
  Counts_SiteSensor_IndicesMorning$uci <- apply(bb.lmer[[currentAcousticIndex]]$t, 2, quantile, 0.975)
  Counts_SiteSensor_IndicesMorning$pred <- predict(model_lmer[[currentAcousticIndex]], re.form=NA)
  
  #Plot of fixed effects and confidence intervals
  Plots_SensorTotals_MorningIndicesMean_lmer[[currentAcousticIndex]] <- Counts_SiteSensor_IndicesMorning %>% 
    ggplot(aes_string(x = colnames(IndicesSummary_Morning[,indice]), y = "n")) + 
    geom_ribbon(aes_string(x = colnames(IndicesSummary_Morning[,indice]), ymin = "lci", ymax = "uci"), fill = "black", alpha = 0.1) +
    geom_line(aes_string(x = colnames(IndicesSummary_Morning[,indice]), y = "pred"), color = "black", lwd = 1) +
    geom_point(size = 2) +
    scale_y_continuous(limits = c(10, 50)) + 
    annotate(geom = "text", -Inf, Inf, hjust = -0.1, vjust = 1.8, parse = T,
             label = paste0("R^2 == ", format(round(as.numeric(performance::r2(model_lmer[[currentAcousticIndex]])[[2]]), 2), nsmall = 2))) +
    labs(x = gsub("_mean", "", colnames(IndicesSummary_Morning[,indice])), y = "Species Richness") +
    theme_classic()
}

plot_grid(plotlist = Plots_SensorTotals_MorningIndicesMean_lmer,
          ncol = 4) %>% annotate_figure(top = "lmer - AcousticIndicesMorning_BirdRichness")

ggsave("./Plots/Plots_SensorTotals_MorningIndicesMean_lmer.png",
       width = 30, height = 18, units = "cm", dpi = 1200)