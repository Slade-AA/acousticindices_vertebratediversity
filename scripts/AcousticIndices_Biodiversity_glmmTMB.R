#Correlate acoustic indices with vertebrate diversity - 2021

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(lme4)
library(lmerTest)
library(boot)
library(glmmTMB)


# Load indices and richness data, merge -----------------------------------

setwd("C:/Users/jc696551/OneDrive - James Cook University/Projects/acousticindices_vertebratediversity")

load("./outputs/data/2021-12-06_acousticIndices_summary.RData") #load indices

richness <- read_csv("./rawdata/biodiversity/richness_diversity_updated.csv") #load richness

#format richness to match indices
richness <- richness %>% rename(Site = site, Sensor = original.plot, richness = Richness, shannon = Shannon, type = taxa)
#richness$Sensor <- gsub(" ", "", richness$Sensor) #now redundant?
richness <- richness %>% mutate(sampling.period = paste0(season, ".2021"))
richness <- richness %>% mutate(type = recode(type,
                                              'bird' = 'birds', 'frog' = 'frogs'))

#duplicate indices 'all' and rename to 'not.birds' for comparison to all other species
acousticIndices_summary <- rbind(acousticIndices_summary, 
                                 acousticIndices_summary %>% filter(type == 'all') %>% mutate(type = recode(type, 'all' = 'not.birds')))

#merge richness and acoustic indices
acousticIndices_richness <- merge(acousticIndices_summary, richness[richness$day == 'all',], by = c("type", "Site", "Sensor", "sampling.period"))


#scale all indices to 0-1 and 'censor' min and max values
#range01 <- function(x){(x-min(x))/(max(x)-min(x))}

#acousticIndices_richness <- acousticIndices_richness %>% 
#  mutate_at(vars(matches("_mean")), range01) %>% 
#  mutate_at(vars(matches("_mean")), ~ replace(.x, which(.x == min(.x)), 0.00001)) %>% 
#  mutate_at(vars(matches("_mean")), ~ replace(.x, which(.x == max(.x)), 0.99999))


#divide all indices by their maximum - except NDSI ((NDSI + 1)/2) - (sensu Bradfer-Lawrence et al. 2020)
#what about BGN? - using abs
range02 <- function(x){abs(x)/max(abs(x))}
ndsi_range <- function(x){(x+1)/2}

acousticIndices_richness <- acousticIndices_richness %>% 
  mutate_at(vars(matches("_mean"), -matches("NDSI")), range02) %>% 
  #mutate_at(vars(matches("_mean")), ~ replace(.x, which(.x == min(.x)), 0.00001)) %>% 
  mutate_at(vars(matches("_mean"), -matches("NDSI")), ~ replace(.x, which(.x == max(.x)), 0.99999)) %>% 
  mutate_at(vars(matches("NDSI")), ndsi_range)

# Mixed effects glmmTMB models --------------------------------------------

Plots_Indices_Richness <- list()
Models_Indices_Richness <- list()
Bootstrap_Indices_Richness <- list()

IndicesRichness_Combinations <- expand.grid(index = grep("mean", colnames(acousticIndices_richness), value = TRUE),
                                            taxa = unique(acousticIndices_richness$type),
                                            diversity = c("richness", "shannon", "count"))

pb = txtProgressBar(min = 0, max = nrow(IndicesRichness_Combinations), initial = 0, style = 3) #progress bar for loop

for (combination in 1:nrow(IndicesRichness_Combinations)) {
  currentAcousticIndex <- IndicesRichness_Combinations$index[combination]
  currentTaxa <- IndicesRichness_Combinations$taxa[combination]
  currentDiversity <- IndicesRichness_Combinations$diversity[combination]
  
  data <- acousticIndices_richness %>% filter(type == currentTaxa)
  
  Models_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]] <- glmmTMB(as.formula(paste0(currentAcousticIndex, " ~ ", currentDiversity, " + (1|Site) + (1|sampling.period)")),
                                                                                                              data = data,
                                                                                                              family = beta_family())
  
  Bootstrap_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]] <- bootMer(Models_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]],
                                                                                                                 FUN = function(x)predict(x, re.form=NA, type = "response"),
                                                                                                                 nsim = 100)
  
  data$lci <- apply(Bootstrap_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]]$t, 2, quantile, 0.025)
  data$uci <- apply(Bootstrap_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]]$t, 2, quantile, 0.975)
  data$pred <- predict(Models_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]], re.form=NA, type = "response")
  
  Plots_Indices_Richness[[paste0(currentDiversity, "_", currentAcousticIndex, "_", currentTaxa)]] <- data %>% 
    ggplot(aes_string(x = paste0(currentDiversity), y = paste0(currentAcousticIndex))) +
    geom_ribbon(aes_string(x = paste0(currentDiversity), ymin = "lci", ymax = "uci"), fill = "black", alpha = 0.1) +
    geom_line(aes_string(x = paste0(currentDiversity), y = "pred"), color = "black", lwd = 1) +
    geom_point(size = 2) +
    labs(y = gsub("_mean", "", paste0(currentAcousticIndex)), x = "Diversity Measure") +
    theme_classic() +
    theme(axis.title.x = element_blank())
  
  setTxtProgressBar(pb,combination)
}

# Arrange multi-panel plots -----------------------------------------------

plot_combinations <- unique(IndicesRichness_Combinations[,c('diversity','taxa')])

for (combination in 1:nrow(plot_combinations)) {
  
  xtitle <- if(plot_combinations$diversity[combination] == 'richness') {c("Species Richness")} else
    if(plot_combinations$diversity[combination] == 'shannon') {c("Shannon Diversity")} else
      if(plot_combinations$diversity[combination] == 'count') {c("Total Abundance")}
  
  Plot <- plot_grid(plotlist = Plots_Indices_Richness[grep(paste0(plot_combinations$diversity[combination], ".*mean_", plot_combinations$taxa[combination]), names(Plots_Indices_Richness))],
                    ncol = 4) %>% 
    annotate_figure(top = paste0("glmmTMB_beta - AcousticIndices_", plot_combinations$taxa[combination], "_", plot_combinations$diversity[combination]),
                    bottom = xtitle,
                    left = "Acoustic Index")
  
  ggsave(filename = paste0("./outputs/figures.local/glmmTMB_beta/", Sys.Date(), "Plots_", plot_combinations$taxa[combination], "_", plot_combinations$diversity[combination], ".png"),
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
  
  xtitle <- if(plot_combinations$diversity[combination] == 'richness') {c("Species Richness")} else
    if(plot_combinations$diversity[combination] == 'shannon') {c("Shannon Diversity")} else
      if(plot_combinations$diversity[combination] == 'count') {c("Total Abundance")}
  
  Plot <- plot_grid(plotlist = Plots_Indices_Richness[grep(paste0(plot_combinations$diversity[combination], "_(", paste0(indicesToUse, collapse = "|"), ")_.*mean_", plot_combinations$taxa[combination]), names(Plots_Indices_Richness), value = T)],
                    ncol = 3) %>% 
    annotate_figure(top = paste0("glmmTMB_beta - AcousticIndices_", plot_combinations$taxa[combination], "_", plot_combinations$diversity[combination]),
                    bottom = xtitle,
                    left = "Acoustic Index")
  
  ggsave(filename = paste0("./outputs/figures/glmmTMB_beta/", plot_combinations$taxa[combination], "_", plot_combinations$diversity[combination], ".png"),
         plot = Plot,
         width = 22.5, height = 12.5, units = "cm", dpi = 1200)
}





lme4::bootMer(glmmTMB_model_richness,nsim=10,FUN=function(x) unlist(fixef(x)))


glmmTMB_model_shannon <- glmmTMB(ACI_mean ~ shannon,
                                  data = testData,
                                  family = beta_family())

AICcmodavg::aictab(cand.set = list(richness = glmmTMB_model_richness,
                                   shannon = glmmTMB_model_shannon))