# Load packages ----
library(tidyverse)
library(ggcorrplot)
library(ggpubr)
library(cowplot)


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

# Create correlation plots for each vertebrate taxaonomic group ----

Plots_IndicesCorrelations <- list()
for (comparison in unique(acousticIndices_richness$type)) {
  tmp <- cor(acousticIndices_richness[acousticIndices_richness$type == comparison, which(colnames(acousticIndices_richness) %in% paste0(indicesToUse, "_mean"))])
  
  axisLabels_corrplot <- c("Activity_mean" = "ACT", 
                           "EventsPerSecond_mean" = "ENV", 
                           "SpectralCentroid_mean" = "CENT", 
                           "HighFreqCover_mean" = "HFC", 
                           "MidFreqCover_mean" = "MFC", 
                           "LowFreqCover_mean" = "LFC", 
                           "AcousticComplexity_mean" = "ACI", 
                           "TemporalEntropy_mean" = "ENT", 
                           "ClusterCount_mean" = "CLS", 
                           "ThreeGramCount_mean" = "TGC", 
                           "Ndsi_mean" = "NDSI", 
                           "SptDensity_mean" = "SPD",
                           "BI_mean" = "BI",
                           "AEI_mean" = "AEI",
                           "ADI_mean" = "ADI",
                           "NDSI_mean" = "NDSI",
                           "SH_mean" = "SH")
  
  Plots_IndicesCorrelations[[comparison]] <- ggcorrplot(tmp, 
                                                        method = "square", 
                                                        type = "upper",
                                                        show.diag = FALSE,
                                                        colors = c("#6D9EC1", "white", "#E46726"),
                                                        outline.color = "black",
                                                        lab = TRUE) + 
    scale_x_discrete(labels = axisLabels_corrplot) +
    scale_y_discrete(labels = axisLabels_corrplot) +
    scale_size(range = c(1, 6)) +
    theme(legend.position = "none")
}

# ├ Create combined plot ----


plot_grid(plotlist = Plots_IndicesCorrelations,
          labels = names(Plots_IndicesCorrelations))