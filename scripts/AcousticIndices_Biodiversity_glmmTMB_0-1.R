#Correlate acoustic indices with vertebrate diversity - 2021

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(lme4)
library(lmerTest)
library(boot)
library(glmmTMB)

# Load indices and biodiversity data --------------------------------------

#summary indices
files <- file.info(list.files("./outputs/data/weekly_summaries/", pattern = ".*_acousticIndices_summary.RData$", full.names = TRUE)) #list files
latestFile <- rownames(files)[which.max(files$mtime)] #determine most recent file to use for loading

load(latestFile)

#richness
richness <- read_csv("./rawdata/biodiversity/richness_diversity_updated.csv") #load richness

#format richness to match indices
richness <- richness %>% rename(Site = site, Sensor = original.plot, richness = Richness, shannon = Shannon, type = taxa)
richness <- richness %>% mutate(sampling.period = paste0(season, ".2021"))
richness <- richness %>% mutate(type = recode(type,
                                              'bird' = 'birds', 'frog' = 'frogs'))

# Merge indices and biodiversity data -------------------------------------

#function to combine indices (all, day, night, morning, afternoon, evening) and biodiversity (all, not.birds, birds, frogs) based on user specified combinations
combineIndicesBiodiversity <- function(indices,
                                       richness,
                                       combinations) {
  frames <- list()
  for (combination in 1:length(combinations)) {
    frames[[combination]] <- merge(indices %>% 
                                     filter(type == combinations[[combination]][1]) %>%
                                     mutate(type = case_when(type == combinations[[combination]][1] ~ names(combinations[combination]))),
                                   richness %>% 
                                     filter(day == 'all') %>% 
                                     filter(type == combinations[[combination]][2]) %>% 
                                     mutate(type = case_when(type == combinations[[combination]][2] ~ names(combinations[combination]))),
                                   by = c("type", "Site", "Sensor", "sampling.period"))
  }
  frames <- do.call(rbind, frames)
  
  return(frames)
}

#combine indices and biodiversity using above function
acousticIndices_richness <- combineIndicesBiodiversity(indices = acousticIndices_summary,
                                                       richness = richness,
                                                       combinations = list(all_all = c('all', 'all'),
                                                                           not.birds_all = c('all', 'not.birds'),
                                                                           birds_morning = c('morning', 'birds'),
                                                                           birds_afternoon = c('afternoon', 'birds'),
                                                                           birds_day = c('day', 'birds'),
                                                                           birds_all = c('all', 'birds'),
                                                                           frogs_evening = c('evening', 'frogs'),
                                                                           frogs_night = c('night', 'frogs')))
acousticIndices_richness <- acousticIndices_richness %>% filter(p > 0.7) #remove datapoints where less than 70% of audio was available

#scale all indices to 0-1 and 'censor' min and max values
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

acousticIndices_richness <- acousticIndices_richness %>% 
  group_by(type) %>% #group_by to scale by taxa - frogs and birds use different indice values than all and not.birds
  mutate_at(vars(matches("_mean")), range01) %>% 
  mutate_at(vars(matches("_mean")), ~ replace(.x, which(.x == min(.x)), 0.00001)) %>% 
  mutate_at(vars(matches("_mean")), ~ replace(.x, which(.x == max(.x)), 0.99999))

# Mixed effects glmmTMB models --------------------------------------------

#indicesToUse <- c('ACI', 'ADI', 'AEI', 'BI', 'NDSI', 'EVN', 'SH', 'LFC', 'MFC', 'HFC') #specify which indices to use in models
indicesToUse <- c('ADI', 'AEI', 'BI', 'NDSI', 'SH', 
                  'Activity', 'EventsPerSecond', 'LowFreqCover', 'MidFreqCover', 'HighFreqCover', 
                  'AcousticComplexity', 'ClusterCount', 'SptDensity')

Plots_Indices_Richness <- list()
Models_Indices_Richness <- list()
Bootstrap_Indices_Richness <- list()

IndicesRichness_Combinations <- expand.grid(index = grep("mean", colnames(acousticIndices_richness), value = TRUE)[which(gsub("_mean","",grep("mean", colnames(acousticIndices_richness), value = TRUE)) %in% indicesToUse)],
                                            comparison = c("all_all", "not.birds_all", "birds_day", "frogs_night"),
                                            measure = c("richness", "shannon", "count"))

pb = txtProgressBar(min = 0, max = nrow(IndicesRichness_Combinations), initial = 0, style = 3) #progress bar for loop

for (combination in 1:nrow(IndicesRichness_Combinations)) {
  currentAcousticIndex <- IndicesRichness_Combinations$index[combination]
  currentComparison <- IndicesRichness_Combinations$comparison[combination]
  currentMeasure <- IndicesRichness_Combinations$measure[combination]
  
  data <- acousticIndices_richness %>% filter(type == currentComparison)
  
  Models_Indices_Richness[[paste0(currentMeasure, "_", currentAcousticIndex, "_", currentComparison)]] <- glmmTMB(as.formula(paste0(currentAcousticIndex, " ~ ", currentMeasure, " + (1|Site) + (1|sampling.period)")),
                                                                                                              data = data,
                                                                                                              family = beta_family())
  
  Bootstrap_Indices_Richness[[paste0(currentMeasure, "_", currentAcousticIndex, "_", currentComparison)]] <- bootMer(Models_Indices_Richness[[paste0(currentMeasure, "_", currentAcousticIndex, "_", currentComparison)]],
                                                                                                                 FUN = function(x)predict(x, re.form=NA, type = "response"),
                                                                                                                 nsim = 100)
  
  data$lci <- apply(Bootstrap_Indices_Richness[[paste0(currentMeasure, "_", currentAcousticIndex, "_", currentComparison)]]$t, 2, quantile, 0.025, na.rm = TRUE)
  data$uci <- apply(Bootstrap_Indices_Richness[[paste0(currentMeasure, "_", currentAcousticIndex, "_", currentComparison)]]$t, 2, quantile, 0.975, na.rm = TRUE)
  data$pred <- predict(Models_Indices_Richness[[paste0(currentMeasure, "_", currentAcousticIndex, "_", currentComparison)]], re.form=NA, type = "response")
  
  Plots_Indices_Richness[[paste0(currentMeasure, "_", currentAcousticIndex, "_", currentComparison)]] <- data %>% 
    ggplot(aes_string(x = paste0(currentMeasure), y = paste0(currentAcousticIndex))) +
    geom_ribbon(aes_string(x = paste0(currentMeasure), ymin = "lci", ymax = "uci"), fill = "black", alpha = 0.1) +
    geom_line(aes_string(x = paste0(currentMeasure), y = "pred"), color = "black", lwd = 1) +
    geom_point(size = 2) +
    scale_y_continuous(limits = c(0, 1), 
                       labels = scales::label_number(accuracy = 0.1),
                       breaks = seq(0, 1, 0.2)) +
    labs(y = gsub("_mean", "", paste0(currentAcousticIndex)), x = "Diversity Measure") +
    theme_classic() +
    theme(axis.title.x = element_blank())
  
  setTxtProgressBar(pb,combination)
}

# Arrange multi-panel plots -----------------------------------------------

plot_combinations <- unique(IndicesRichness_Combinations[,c('measure','comparison')])

for (combination in 1:nrow(plot_combinations)) {
  
  xtitle <- if(plot_combinations$measure[combination] == 'richness') {c("Species Richness")} else
    if(plot_combinations$measure[combination] == 'shannon') {c("Shannon Diversity")} else
      if(plot_combinations$measure[combination] == 'count') {c("Total Abundance")}
  
  Plot <- plot_grid(plotlist = Plots_Indices_Richness[grep(paste0(plot_combinations$measure[combination], ".*mean_", plot_combinations$comparison[combination]), names(Plots_Indices_Richness))],
                    ncol = 4) %>% 
    annotate_figure(top = paste0("glmmTMB_beta - AcousticIndices_", plot_combinations$comparison[combination], "_", plot_combinations$measure[combination]),
                    bottom = xtitle,
                    left = "Acoustic Index")
  
  ggsave(filename = paste0("./outputs/figures.local/glmmTMB_beta_0-1/", Sys.Date(), "Plots_", plot_combinations$comparison[combination], "_", plot_combinations$measure[combination], ".png"),
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
  
  xtitle <- if(plot_combinations$measure[combination] == 'richness') {c("Species Richness")} else
    if(plot_combinations$measure[combination] == 'shannon') {c("Shannon Diversity")} else
      if(plot_combinations$measure[combination] == 'count') {c("Total Abundance")}
  
  Plot <- plot_grid(plotlist = Plots_Indices_Richness[grep(paste0(plot_combinations$measure[combination], "_(", paste0(indicesToUse, collapse = "|"), ")_.*mean_", plot_combinations$comparison[combination]), names(Plots_Indices_Richness), value = T)],
                    ncol = 3) %>% 
    annotate_figure(top = paste0("glmmTMB_beta - AcousticIndices_", plot_combinations$comparison[combination], "_", plot_combinations$measure[combination]),
                    bottom = xtitle,
                    left = "Acoustic Index")
  
  ggsave(filename = paste0("./outputs/figures/glmmTMB_beta_0-1/", plot_combinations$comparison[combination], "_", plot_combinations$measure[combination], ".png"),
         plot = Plot,
         width = 22.5, height = 12.5, units = "cm", dpi = 1200)
}

# Save workspace for later loading ----------------------------------------

save.image(file = "outputs/workspaces/glmmTMB_0-1.RData")