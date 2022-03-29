# Load packages -----------------------------------------------------------

library(tidyverse)
library(factoextra)

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


# PCA ---------------------------------------------------------------------

indicesToUse <- c('ADI', 'AEI', 'BI', 'NDSI', 'SH', 
                  'Activity', 'EventsPerSecond', 'LowFreqCover', 'MidFreqCover', 'HighFreqCover', 
                  'AcousticComplexity', 'ClusterCount', 'SptDensity')

pca.all <- prcomp(acousticIndices_richness[acousticIndices_richness$type == 'all_all', 
                                           c(which(colnames(acousticIndices_richness) %in% paste0(indicesToUse, "_mean")))], 
                  center = TRUE, scale. = TRUE)

var <- get_pca_var(pca.all)
corrplot::corrplot(var$cos2, is.corr=FALSE)

# Contributions of variables to PC1
fviz_contrib(pca.all, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(pca.all, choice = "var", axes = 2, top = 10)

fviz_eig(pca.all, addlabels = TRUE)

fviz_pca_ind(pca.all,
             col.ind = acousticIndices_richness$richness[acousticIndices_richness$type == 'all_all'],
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

fviz_pca_biplot(pca.all,
                habillage = acousticIndices_richness$Site[acousticIndices_richness$type == 'all_all'],
                addEllipses = TRUE,
                mean.point = FALSE,
                ellipse.level = 0.95,
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                repel = TRUE)