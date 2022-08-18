#Single site linear models of acoustic indices with vertebrate diversity - 2021

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(lme4)
library(lmerTest)

# Load indices and biodiversity data --------------------------------------

#summary indices
files <- file.info(list.files("./outputs/data/weekly_summaries/", pattern = ".*_acousticIndices_summary_sunrise_sunset.RData$", full.names = TRUE)) #list files
latestFile <- rownames(files)[which.max(files$mtime)] #determine most recent file to use for loading

load(latestFile)

#richness
richness <- read_csv("./rawdata/biodiversity/updated.diversity.indices.allonly.csv") #load richness

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
                                                                           #birds_morning = c('morning', 'birds'),
                                                                           #birds_afternoon = c('afternoon', 'birds'),
                                                                           birds_day = c('day', 'birds'),
                                                                           #birds_all = c('all', 'birds'),
                                                                           #frogs_evening = c('evening', 'frogs'),
                                                                           frogs_night = c('night', 'frogs')))
acousticIndices_richness <- acousticIndices_richness %>% filter(p > 0.7) #remove datapoints where less than 70% of audio was available

# Fit single site linear models -------------------------------------------

models_singleSite <- list()
obs_pred_values <- list()
for (site in c("Tarcutta", "Wambiana", "Rinyirru")) {
  for (comparison in c("birds_day")) {
    for (measure in c("richness", "count")) {
      models_singleSite[[paste(site, comparison, measure, sep = "_")]] <- lmer(as.formula(paste0(measure, " ~ ", "ClusterCount_mean + SptDensity_mean + MidFreqCover_mean", " + (1|Sensor)")), 
                                                                               data = acousticIndices_richness[acousticIndices_richness$type == comparison &
                                                                                                                 acousticIndices_richness$Site == site,])
      obs_pred_values[[paste(site, comparison, measure, sep = "_")]] <- data.frame(obs = acousticIndices_richness[acousticIndices_richness$type == comparison &
                                                                                                                             acousticIndices_richness$Site == site,
                                                                                                                  c(paste0(measure))],
                                                                                   pred = predict(models_singleSite[[paste(site, comparison, measure, sep = "_")]]))
      
    }
  }
}


# Observed vs predicted plots ---------------------------------------------

Plots_obs_pred <- list()
for (model in 1:length(obs_pred_values)) {
  #plot obs vs pred
  Plots_obs_pred[[paste0(names(obs_pred_values)[model])]] <- ggplot(data = obs_pred_values[[model]],
                                                                    aes(x = pred, y = obs)) +
    geom_abline(slope = 1, linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point() +
    annotate(geom = "text", 
             x = 0.8*(max(obs_pred_values[[model]]) - min(obs_pred_values[[model]])) + min(obs_pred_values[[model]]), 
             y = min(obs_pred_values[[model]]), 
             vjust = 0, size = 3,
             label = paste0("CCC: ", round(DescTools::CCC(obs_pred_values[[model]]$pred, 
                                                          obs_pred_values[[model]]$obs)$rho.c$est, 2))) +
    scale_x_continuous(limits = c(min(obs_pred_values[[model]]),
                                  max(obs_pred_values[[model]]))) +
    scale_y_continuous(limits = c(min(obs_pred_values[[model]]),
                                  max(obs_pred_values[[model]]))) +
    labs(x = "Predicted", y = "Observed") +
    theme_classic() +
    theme(axis.title = element_blank())
}

#arrange plots
Plot_ObservedPredicted <- egg::ggarrange(as_ggplot(text_grob(label = "richness")),
                                         as_ggplot(text_grob(label = "count")),
                                         ggplot() + theme_void(),
                                         Plots_obs_pred$Tarcutta_birds_day_richness,
                                         Plots_obs_pred$Tarcutta_birds_day_count,
                                         as_ggplot(text_grob(label = "Tarcutta", rot = 270)),
                                         Plots_obs_pred$Wambiana_birds_day_richness,
                                         Plots_obs_pred$Wambiana_birds_day_count,
                                         as_ggplot(text_grob(label = "Wambiana", rot = 270)),
                                         Plots_obs_pred$Rinyirru_birds_day_richness,
                                         Plots_obs_pred$Rinyirru_birds_day_count,
                                         as_ggplot(text_grob(label = "Rinyirru", rot = 270)),
                                         nrow = 4, ncol = 3,
                                         heights = c(0.1,1,1,1),
                                         widths = c(1,1,0.1)) %>% 
  annotate_figure(bottom = "Predicted",
                  left = "Observed")

ggsave(filename = "outputs/figures/singleSite_Birds_obs_pred_sunrise_sunset_fixedbiodiversity.png",
       plot = Plot_ObservedPredicted,
       width = 12, height = 18, units = "cm", dpi = 800)

#New 1000dpi figure for publication (Ecological Indicators)

ggsave(filename = "outputs/figures/publication/Figure06.png",
       plot = Plot_ObservedPredicted,
       width = 90, height = 135, units = "mm", dpi = 1000)

# Save workspace for later loading ----------------------------------------

save.image(file = "outputs/workspaces/SingleSite_lmer_sunrise_sunset_fixedbiodiversity.RData")