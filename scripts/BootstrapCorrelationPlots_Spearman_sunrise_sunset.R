#Bootstrap correlation analysis of individual indices and vertebrate biodiversity

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(lme4)
library(lmerTest)
library(boot)

# Load indices and biodiversity data --------------------------------------

#summary indices
files <- file.info(list.files("./outputs/data/weekly_summaries/", pattern = ".*_acousticIndices_summary_sunrise_sunset.RData$", full.names = TRUE)) #list files
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

#function to combine indices (all, day, night, morning, evening) and biodiversity (all, not.birds, birds, frogs) based on user specified combinations
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

# Calculate bootstrap correlation between indices and biodiversity --------

bootCor_results <- list()

for (measure in c('count', 'richness', 'shannon')) {
  
  for (comparison in unique(acousticIndices_richness$type)) {
    
    for (index in colnames(select(acousticIndices_richness, ends_with(c("mean"))))) {
      set.seed(1234)#set seed for reproducibility
      bootResults <- boot(acousticIndices_richness[acousticIndices_richness$type == comparison,], 
                          statistic = function(data, i) {
                            cor(data[i, measure], data[i, index], method='spearman')
                          },
                          R = 1000)
      bootResultsCI <- boot.ci(bootResults, 
                               conf = 0.95, type = "bca")
      
      bootCor_results[[paste(measure, comparison, index, sep = "_")]] <- data.frame(Index = gsub("_mean", "", index),
                                                                              Taxa = comparison,
                                                                              Measure = measure,
                                                                              Mean = mean(bootResults$t),
                                                                              Low = bootResultsCI$bca[4],
                                                                              High = bootResultsCI$bca[5])
    }
  }
}

bootCor_results <- do.call(rbind, bootCor_results)

# Plot correlation bootstrap results --------------------------------------

#create acronymns for AP indices
axisLabels <- c("Activity" = "ACT", 
                "EventsPerSecond" = "ENV", 
                "SpectralCentroid" = "CENT", 
                "HighFreqCover" = "HFC", 
                "MidFreqCover" = "MFC", 
                "LowFreqCover" = "LFC", 
                "AcousticComplexity" = "ACI", 
                "TemporalEntropy" = "ENT", 
                "ClusterCount" = "CLS", 
                "ThreeGramCount" = "TGC", 
                "Ndsi" = "NDSI", 
                "SptDensity" = "SPD")

#correlation plot for birds with different ACI frequency bands
indicesToUse <- c("ACI", 
                  "ACI_1000_4000", "ACI_3000_6000", "ACI_5000_8000",
                  "ACI_1000_2000", "ACI_2000_3000", "ACI_3000_4000", "ACI_4000_5000",
                  "ACI_5000_6000", "ACI_6000_7000", "ACI_7000_8000")

correlationPlots_birdsACI <- list()
for (measure in c('richness', 'shannon', 'count')) {
  tmp_data <- bootCor_results[bootCor_results$Measure == measure & 
                                bootCor_results$Index %in% indicesToUse,]
  tmp_data <- tmp_data[grep("^birds*", tmp_data$Taxa),]
  tmp_data$Index <- fct_relevel(tmp_data$Index, indicesToUse)
  
  correlationPlots_birdsACI[[measure]] <- ggplot(data = tmp_data, aes(x = Mean, y = Index, group = Taxa, colour = Taxa)) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4)) +
    geom_hline(yintercept = 1+0.5, linetype = "dotted") +
    geom_hline(yintercept = 4+0.5, linetype = "dotted") +
    scale_x_continuous(limits = c(-0.9, 0.9), breaks = seq(-0.8, 0.8, 0.4)) +
    scale_color_viridis_d() +
    labs(x = "Mean correlation") +
    theme_classic() +
    theme(axis.title = element_blank(),
          legend.position = "none")
}

legend_bottom <- get_legend(
  correlationPlots_birdsACI[['count']] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank())
)

correlationPlot_birdsACI <- plot_grid(plotlist = correlationPlots_birdsACI,
                             ncol = 1,
                             labels = c(paste0("a) ", names(correlationPlots_birdsACI)[1]),
                                        paste0("b) ", names(correlationPlots_birdsACI)[2]),
                                        paste0("c) ", names(correlationPlots_birdsACI)[3])),
                             hjust = 0, label_x = 0.23) %>% 
  annotate_figure(left = "Acoustic index", bottom = "Mean correlation") %>% 
  plot_grid(legend_bottom, ncol = 1, rel_heights = c(1, .1))

ggsave("outputs/figures/bootstrapcorrelations_sunrise_sunset/bootstrap_correlations_birdsACI_spearman.png",
       correlationPlot_birdsACI,
       width = 14, height = 24, units = "cm", dpi = 800)

#correlation plot for birds using indices for morning, afternoon and day
#indicesToUse <- c('ACI', 'ADI', 'AEI', 'BI', 'NDSI', 'EVN', 'SH', 'LFC', 'MFC', 'HFC')
indicesToUse <- c('ADI', 'AEI', 'BI', 'NDSI', 'SH', 
                  'Activity', 'EventsPerSecond', 'LowFreqCover', 'MidFreqCover', 'HighFreqCover', 
                  'AcousticComplexity', 'ClusterCount', 'SptDensity')

correlationPlots_birds <- list()
for (measure in c('richness', 'shannon', 'count')) {
  tmp_data <- bootCor_results[bootCor_results$Measure == measure & 
                                bootCor_results$Index %in% indicesToUse,]
  tmp_data <- tmp_data[grep("^birds_morning|^birds_afternoon|^birds_day", tmp_data$Taxa),]
  tmp_data$Index <- fct_relevel(tmp_data$Index, indicesToUse)
  
  correlationPlots_birds[[measure]] <- ggplot(data = tmp_data, aes(x = Mean, y = Index, group = Taxa, colour = Taxa)) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4)) +
    scale_x_continuous(limits = c(-0.9, 0.9), breaks = seq(-0.8, 0.8, 0.4)) +
    scale_y_discrete(labels = axisLabels) +
    scale_color_viridis_d() +
    labs(x = "Mean correlation") +
    theme_classic() +
    theme(axis.title = element_blank(),
          legend.position = "none")
}

legend_bottom <- get_legend(
  correlationPlots_birds[['count']] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank())
)

correlationPlot_birds <- plot_grid(plotlist = correlationPlots_birds,
                                   ncol = 1,
                                   labels = c(paste0("a) ", names(correlationPlots_birds)[1]),
                                              paste0("b) ", names(correlationPlots_birds)[2]),
                                              paste0("c) ", names(correlationPlots_birds)[3])),
                                   hjust = 0, label_x = 0.12) %>% 
  annotate_figure(left = "Acoustic index", bottom = "Mean correlation") %>% 
  plot_grid(legend_bottom, ncol = 1, rel_heights = c(1, .1))

ggsave("outputs/figures/bootstrapcorrelations_sunrise_sunset/bootstrap_correlations_birds_spearman.png",
       correlationPlot_birds,
       width = 12, height = 24, units = "cm", dpi = 800)

#correlation plot for frogs using indices for evening and night
#indicesToUse <- c('ACI', 'ADI', 'AEI', 'BI', 'NDSI', 'EVN', 'SH', 'LFC', 'MFC', 'HFC')
indicesToUse <- c('ADI', 'AEI', 'BI', 'NDSI', 'SH', 
                  'Activity', 'EventsPerSecond', 'LowFreqCover', 'MidFreqCover', 'HighFreqCover', 
                  'AcousticComplexity', 'ClusterCount', 'SptDensity')

correlationPlots_frogs <- list()
for (measure in c('richness', 'shannon', 'count')) {
  tmp_data <- bootCor_results[bootCor_results$Measure == measure & 
                                bootCor_results$Index %in% indicesToUse,]
  tmp_data <- tmp_data[grep("^frogs*", tmp_data$Taxa),]
  tmp_data$Index <- fct_relevel(tmp_data$Index, indicesToUse)
  
  correlationPlots_frogs[[measure]] <- ggplot(data = tmp_data, aes(x = Mean, y = Index, group = Taxa, colour = Taxa)) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4)) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
    scale_y_discrete(labels = axisLabels) +
    scale_color_viridis_d() +
    labs(x = "Mean correlation") +
    theme_classic() +
    theme(axis.title = element_blank(),
          legend.position = "none")
}

legend_bottom <- get_legend(
  correlationPlots_frogs[['count']] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank())
)

correlationPlot_frogs <- plot_grid(plotlist = correlationPlots_frogs,
                                   ncol = 1,
                                   labels = c(paste0("a) ", names(correlationPlots_frogs)[1]),
                                              paste0("b) ", names(correlationPlots_frogs)[2]),
                                              paste0("c) ", names(correlationPlots_frogs)[3])),
                                   hjust = 0, label_x = 0.12) %>% 
  annotate_figure(left = "Acoustic index", bottom = "Mean correlation") %>% 
  plot_grid(legend_bottom, ncol = 1, rel_heights = c(1, .1))

ggsave("outputs/figures/bootstrapcorrelations_sunrise_sunset/bootstrap_correlations_frogs_spearman.png",
       correlationPlot_frogs,
       width = 12, height = 24, units = "cm", dpi = 800)




#correlation plot facetted by taxa (using birds_day and frogs_night)
#indicesToUse <- c('ACI', 'ADI', 'AEI', 'BI', 'NDSI', 'EVN', 'SH', 'LFC', 'MFC', 'HFC')
indicesToUse <- c('ADI', 'AEI', 'BI', 'NDSI', 'SH', 
                  'Activity', 'EventsPerSecond', 'LowFreqCover', 'MidFreqCover', 'HighFreqCover', 
                  'AcousticComplexity', 'ClusterCount', 'SptDensity')

correlationPlots <- list()
for (taxa in c('all_all', 'not.birds_all', 'birds_day', 'frogs_night')) {
  tmp_data <- bootCor_results[bootCor_results$Taxa == taxa & bootCor_results$Index %in% indicesToUse,]
  tmp_data$Measure <- factor(tmp_data$Measure, levels = c("richness", "shannon", "count"))
  tmp_data$Index <- fct_relevel(tmp_data$Index, indicesToUse)
  
  correlationPlots[[taxa]] <- ggplot(data = tmp_data, aes(x = Mean, y = Index, group = Measure, colour = Measure)) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4)) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
    scale_y_discrete(labels = axisLabels) +
    scale_color_viridis_d(labels = c('Richness', 'Shannon', 'Count')) +
    labs(x = "Mean correlation") +
    theme_classic() +
    theme(axis.title = element_blank(),
          legend.position = "none")
}

legend_bottom <- get_legend(
  correlationPlots[['all_all']] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank())
)

correlationPlot <- plot_grid(plotlist = correlationPlots,
                             ncol = 2,
                             labels = c("a) all vertebrates",
                                        "b) non-avian vertebrates",
                                        "c) birds",
                                        "d) frogs"),
                             hjust = 0, label_x = 0.12, label_y = 1.02) %>% 
  annotate_figure(left = "Acoustic index", bottom = "Mean correlation") %>% 
  plot_grid(legend_bottom, ncol = 1, rel_heights = c(1, .1))

ggsave("outputs/figures/bootstrapcorrelations_sunrise_sunset/bootstrap_correlations_bytaxa_spearman.png",
       correlationPlot,
       width = 24, height = 24, units = "cm", dpi = 800)


# Create table of correlation values --------------------------------------

library(kableExtra)

bootCor_results_table <- bootCor_results %>% 
  filter(Index %in% indicesToUse) %>% 
  filter(Taxa %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'))
bootCor_results_table$Index <- fct_relevel(bootCor_results_table$Index, indicesToUse)


bootCor_results_table$value <- paste0(round(bootCor_results_table$Mean, 2), 
                                      " (", 
                                      round(bootCor_results_table$Low, 2), 
                                      " - ", 
                                      round(bootCor_results_table$High, 2), 
                                      ")")

bootCor_results_table <- bootCor_results_table %>% select(-Mean, -Low, -High) %>% pivot_wider(names_from = Index, values_from = value)

target <- c('richness', 'shannon', 'count')

bootCor_results_table <- bootCor_results_table %>% arrange(factor(Measure, levels = target))
write_csv(bootCor_results_table,
          file = "outputs/figures/bootstrapcorrelations_sunrise_sunset/correlationTable_spearman.csv")



#Calculate spearman rank correlations for best acoustic index per comparison
bootCor_results <- bootCor_results %>% unite("comparison", Taxa:Measure, remove = FALSE)

indicesToUse <- c('ADI', 'AEI', 'BI', 'NDSI', 'SH', 
                  'Activity', 'EventsPerSecond', 'LowFreqCover', 'MidFreqCover', 'HighFreqCover', 
                  'AcousticComplexity', 'ClusterCount', 'SptDensity')

bestindex_spearman <- bootCor_results %>% filter(Index %in% indicesToUse) %>% select(-Low, -High) %>% group_by(comparison) %>% top_n(1, Mean)

save(bestindex_spearman, file = "outputs/figures/bootstrapcorrelations_sunrise_sunset/bestindex_spearman.RData")
  
# Save workspace for later loading ----------------------------------------

save.image(file = "outputs/workspaces/BootstrapCorrelations_spearman_sunrise_sunset.RData")