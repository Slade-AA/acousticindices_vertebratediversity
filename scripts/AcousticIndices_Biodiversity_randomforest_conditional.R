#Random forest model (using party package) of acoustic indices with vertebrate diversity - 2021

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(randomForest)
library(caret)
library(party)
library(permimp)

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

# Fit conditional random forest models ------------------------------------

control <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verbose = FALSE, savePredictions = TRUE)
tunegrid <- expand.grid(.mtry=c(2:10))

RandomForestFits_cforest <- list()
RandomForestPerformance_cforest <- data.frame(comparison = character(),
                                              measure = character(),
                                              ACItype = character(),
                                              mtry = numeric(),
                                              RMSE = numeric(),
                                              Rsquared = numeric(),
                                              MAE = numeric(),
                                              RMSESD = numeric(),
                                              RsquaredSD = numeric(),
                                              MAESD = numeric(),
                                              minResponse = numeric(),
                                              meanResponse = numeric(),
                                              maxResponse = numeric())
RandomForestImportance_cforest <- list()
RandomForestImportance_cforest_conditional <- list()
RandomForestPredictions_cforest <- list()

#indicesToUse <- c('ACI', 'ADI', 'AEI', 'BI', 'NDSI', 'EVN', 'SH', 'LFC', 'MFC', 'HFC') #specify which indices to use in rf models
indicesToUse <- c('ADI', 'AEI', 'BI', 'NDSI', 'SH', 
                  'Activity', 'EventsPerSecond', 'LowFreqCover', 'MidFreqCover', 'HighFreqCover', 
                  'AcousticComplexity', 'ClusterCount', 'SptDensity')

pb = txtProgressBar(min = 0, max = length(unique(acousticIndices_richness$type)) * 3 * 1, initial = 0, style = 3); k <- 0

for (comparison in unique(acousticIndices_richness$type)) {
  tmpdata <- acousticIndices_richness[acousticIndices_richness$type == comparison,]
  
  for (measure in c("richness", "shannon", "count")) {
    
    #forumlas for using just total ACI, 3kHz ACI, and 1kHz ACI values
    #formulas <- list(totalACI = as.formula(paste0(measure, " ~ ", paste(grep(paste(paste0(indicesToUse, "_*"), collapse = "|"), grep("ACI_[0-9].", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), value = TRUE), collapse = " + "))),
    #                 ACI_3kHz = as.formula(paste0(measure, " ~ ", paste(grep(paste(paste0(indicesToUse, "_*"), collapse = "|"), grep("ACI_mean$|ACI_1000_2.|ACI_2000_3.|ACI_3000_4.|ACI_4000_5.|ACI_5000_6.|ACI_6000_7.|ACI_7000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), value = TRUE), collapse = " + "))),
    #                 ACI_1kHz = as.formula(paste0(measure, " ~ ", paste(grep(paste(paste0(indicesToUse, "_*"), collapse = "|"), grep("ACI_mean$|ACI_1000_4.|ACI_3000_6.|ACI_5000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), value = TRUE), collapse = " + "))))
    
    #formula for new indices (AP + Kaleidoscope)
    formulas <- list(totalACI = as.formula(paste0(measure, " ~ ", paste(grep(paste(paste0(indicesToUse, "_*"), collapse = "|"), grep("ACI_[0-9].", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), value = TRUE), collapse = " + "))))
    
    for (formula in 1:length(formulas)) {
      set.seed(1234)#set seed for reproducibility
      RandomForestFits_cforest[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]] <- train(formulas[[formula]],
                                                                                                             data = tmpdata,
                                                                                                             method = "cforest",
                                                                                                             trControl = control,
                                                                                                             tuneGrid = tunegrid,
                                                                                                             controls = cforest_unbiased(ntree = 1000))
      
      RandomForestPerformance_cforest <- rbind(RandomForestPerformance_cforest, data.frame(cbind(comparison = comparison, 
                                                                                                 measure = measure, 
                                                                                                 ACItype = names(formulas)[[formula]],
                                                                                                 RandomForestFits_cforest[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$results[RandomForestFits_cforest[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$results$mtry == RandomForestFits_cforest[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$bestTune[[1]],],
                                                                                                 minResponse = min(tmpdata[measure]),
                                                                                                 meanResponse = mean(tmpdata[[measure]]),
                                                                                                 maxResponse = max(tmpdata[measure]))))
      
      RandomForestImportance_cforest[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]] <- varImp(RandomForestFits_cforest[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]])
      RandomForestImportance_cforest_conditional[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]] <- permimp(RandomForestFits_cforest[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$finalModel, conditional = TRUE, progressBar = FALSE, scaled = TRUE)
      
      RandomForestPredictions_cforest[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]] <- data.frame(predictions = predict(RandomForestFits_cforest[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$finalModel),
                                                                                                                         observations = tmpdata[[measure]],
                                                                                                                         comparison = paste0(comparison, "_", measure),
                                                                                                                         ACItype = names(formulas)[[formula]],
                                                                                                                         stringsAsFactors = FALSE)
      
      
      k <- k+1; setTxtProgressBar(pb, k)
    }
  }
}

#Normalize RMSE and MAE, calculate Scatter Index
RandomForestPerformance_cforest <- RandomForestPerformance_cforest %>% mutate(normRMSE = RMSE/(maxResponse-minResponse),
                                                                              normMAE = MAE/(maxResponse-minResponse),
                                                                              SI = (RMSE/meanResponse)*100,
                                                                              normRMSE_SE = RMSESD/(maxResponse-minResponse)/sqrt(30),
                                                                              normMAE_SE = MAESD/(maxResponse-minResponse)/sqrt(30),
                                                                              SI_SE = (RMSESD/meanResponse)*100/sqrt(30),
                                                                              Rsquared_SE = RsquaredSD/sqrt(30))

RandomForestPredictions_cforest <- do.call(rbind, RandomForestPredictions_cforest)

# Plot random forest performance ------------------------------------------

#compare rf models with different ACI bands used as predictors
ggplot(data = RandomForestPerformance_cforest, aes(x = measure, y = normRMSE, group = ACItype, fill = ACItype)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = normRMSE-normRMSE_SE, ymax = normRMSE+normRMSE_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised RMSE") +
  facet_wrap(~comparison) +
  theme_bw()

#compare bird indices comparisons - using totalACI rf models
Plot_RMSE_birds <- ggplot(data = RandomForestPerformance_cforest[RandomForestPerformance_cforest$ACItype == 'totalACI' & RandomForestPerformance_cforest$comparison %in% grep("^birds*", RandomForestPerformance_cforest$comparison, value = TRUE),], 
                          aes(x = measure, y = normRMSE, group = comparison, fill = comparison)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = normRMSE-normRMSE_SE, ymax = normRMSE+normRMSE_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised RMSE") +
  theme_bw() +
  theme(legend.position = "none")
Plot_SI_birds <- ggplot(data = RandomForestPerformance_cforest[RandomForestPerformance_cforest$ACItype == 'totalACI' & RandomForestPerformance_cforest$comparison %in% grep("^birds*", RandomForestPerformance_cforest$comparison, value = TRUE),], 
                        aes(x = measure, y = SI, group = comparison, fill = comparison)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = SI-SI_SE, ymax = SI+SI_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Taxa", y = "Scatter Index") +
  theme_bw() +
  theme(legend.position = "none")
Plot_R2_birds <- ggplot(data = RandomForestPerformance_cforest[RandomForestPerformance_cforest$ACItype == 'totalACI' & RandomForestPerformance_cforest$comparison %in% grep("^birds*", RandomForestPerformance_cforest$comparison, value = TRUE),], 
                        aes(x = measure, y = Rsquared, group = comparison, fill = comparison)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = Rsquared-Rsquared_SE, ymax = Rsquared+Rsquared_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "R squared") +
  theme_bw() +
  theme(legend.position = "none")

legend_bottom <- get_legend(
  Plot_R2_birds + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank())
)

Plot_Birds_Time <- plot_grid(Plot_RMSE_birds, Plot_SI_birds, Plot_R2_birds, legend_bottom,
                             ncol = 2)
#rel_heights = c(1, 1, 1, 0.1))
ggsave(filename = "outputs/figures/randomforestperformance/Birds_MorningDayAll_cforest.png",
       plot = Plot_Birds_Time,
       width = 24, height = 20, units = "cm", dpi = 800)


#plot performance measures for all comparisons - totalACI
Plot_RMSE <- ggplot(data = RandomForestPerformance_cforest[RandomForestPerformance_cforest$ACItype == 'totalACI' &
                                                     RandomForestPerformance_cforest$comparison %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'),], 
                    aes(x = comparison, y = normRMSE, group = measure, fill = measure)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = normRMSE-normRMSE_SE, ymax = normRMSE+normRMSE_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised RMSE") +
  theme_classic() +
  theme(legend.position = "none")
Plot_MAE <- ggplot(data = RandomForestPerformance_cforest[RandomForestPerformance_cforest$ACItype == 'totalACI' &
                                                    RandomForestPerformance_cforest$comparison %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'),], 
                   aes(x = comparison, y = normMAE, group = measure, fill = measure)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = normMAE-normMAE_SE, ymax = normMAE+normMAE_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised MAE") +
  theme_classic() +
  theme(legend.position = "none")
Plot_SI <- ggplot(data = RandomForestPerformance_cforest[RandomForestPerformance_cforest$ACItype == 'totalACI' &
                                                   RandomForestPerformance_cforest$comparison %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'),], 
                  aes(x = comparison, y = SI, group = measure, fill = measure)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = SI-SI_SE, ymax = SI+SI_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 40)) +
  labs(x = "Taxa", y = "Scatter Index") +
  theme_classic() +
  theme(legend.position = "none")
Plot_R2 <- ggplot(data = RandomForestPerformance_cforest[RandomForestPerformance_cforest$ACItype == 'totalACI' &
                                                   RandomForestPerformance_cforest$comparison %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'),], 
                  aes(x = comparison, y = Rsquared, group = measure, fill = measure)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = Rsquared-Rsquared_SE, ymax = Rsquared+Rsquared_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "R Squared") +
  theme_classic() +
  theme(legend.position = "none")

legend_bottom <- get_legend(
  Plot_R2 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank())
)

Plot_AllRF <- plot_grid(Plot_RMSE, Plot_MAE, Plot_SI, Plot_R2,
                        ncol = 2)
Plot_AllRF <- plot_grid(Plot_AllRF, legend_bottom,
                        ncol = 1, 
                        rel_heights = c(1, 0.1))
ggsave(filename = "outputs/figures/randomforestperformance/AllComparisons_cforest.png",
       plot = Plot_AllRF,
       width = 24, height = 20, units = "cm", dpi = 800)

# Observed vs predicted plots ---------------------------------------------

Plots_ObsPred <- list()

for (predictions in c('all_all_richness', 'all_all_shannon', 'all_all_count',
                      'not.birds_all_richness', 'not.birds_all_shannon', 'not.birds_all_count',
                      'birds_day_richness', 'birds_day_shannon', 'birds_day_count',
                      'frogs_night_richness', 'frogs_night_shannon', 'frogs_night_count')) {
  Plots_ObsPred[[paste0(predictions)]] <- ggplot(RandomForestPredictions_cforest[RandomForestPredictions_cforest$comparison == predictions &
                                                                                   RandomForestPredictions_cforest$ACItype == 'totalACI',], 
                                                 aes(x = .outcome, y = observations)) +
    geom_abline(slope = 1, linetype = 'dashed') +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point() +
    annotate(geom = "text", 
             x = 0.8*(max(RandomForestPredictions_cforest[RandomForestPredictions_cforest$comparison == predictions &
                                                            RandomForestPredictions_cforest$ACItype == 'totalACI',c(1:2)]) -
                        min(RandomForestPredictions_cforest[RandomForestPredictions_cforest$comparison == predictions &
                                                            RandomForestPredictions_cforest$ACItype == 'totalACI',c(1:2)])) + 
               min(RandomForestPredictions_cforest[RandomForestPredictions_cforest$comparison == predictions &
                                                     RandomForestPredictions_cforest$ACItype == 'totalACI',c(1:2)]), 
             y = min(RandomForestPredictions_cforest[RandomForestPredictions_cforest$comparison == predictions &
                                                                    RandomForestPredictions_cforest$ACItype == 'totalACI',c(1:2)]), 
             vjust = 0, size = 3,
             label = paste0("CCC: ", round(DescTools::CCC(RandomForestPredictions_cforest$.outcome[RandomForestPredictions_cforest$comparison == predictions &
                                                                                                     RandomForestPredictions_cforest$ACItype == 'totalACI'], 
                                                          RandomForestPredictions_cforest$observations[RandomForestPredictions_cforest$comparison == predictions &
                                                                                                         RandomForestPredictions_cforest$ACItype == 'totalACI'])$rho.c$est, 2))) +
    scale_x_continuous(limits = c(min(RandomForestPredictions_cforest[RandomForestPredictions_cforest$comparison == predictions &
                                                                        RandomForestPredictions_cforest$ACItype == 'totalACI',c(1:2)]), 
                                  max(RandomForestPredictions_cforest[RandomForestPredictions_cforest$comparison == predictions &
                                                                        RandomForestPredictions_cforest$ACItype == 'totalACI',c(1:2)]))) +
    scale_y_continuous(limits = c(min(RandomForestPredictions_cforest[RandomForestPredictions_cforest$comparison == predictions &
                                                                        RandomForestPredictions_cforest$ACItype == 'totalACI',c(1:2)]), 
                                  max(RandomForestPredictions_cforest[RandomForestPredictions_cforest$comparison == predictions &
                                                                        RandomForestPredictions_cforest$ACItype == 'totalACI',c(1:2)]))) +
    theme_classic() +
    theme(axis.title = element_blank())
}

#plot_grid(plotlist = Plots_ObsPred)

#arrange plot using egg package - even sized plots
Plot_ObservedPredicted <- egg::ggarrange(as_ggplot(text_grob(label = "richness")),
                                         as_ggplot(text_grob(label = "shannon")),
                                         as_ggplot(text_grob(label = "count")),
                                         ggplot() + theme_void(),
                                         Plots_ObsPred[['all_all_richness']],
                                         Plots_ObsPred[['all_all_shannon']],
                                         Plots_ObsPred[['all_all_count']],
                                         as_ggplot(text_grob(label = "all vertebrates", rot = 270)),
                                         Plots_ObsPred[['not.birds_all_richness']],
                                         Plots_ObsPred[['not.birds_all_shannon']],
                                         Plots_ObsPred[['not.birds_all_count']],
                                         as_ggplot(text_grob(label = "non-avian", rot = 270)),
                                         Plots_ObsPred[['birds_day_richness']],
                                         Plots_ObsPred[['birds_day_shannon']],
                                         Plots_ObsPred[['birds_day_count']],
                                         as_ggplot(text_grob(label = "birds", rot = 270)),
                                         Plots_ObsPred[['frogs_night_richness']],
                                         Plots_ObsPred[['frogs_night_shannon']],
                                         Plots_ObsPred[['frogs_night_count']],
                                         as_ggplot(text_grob(label = "frogs", rot = 270)),
                                         ncol = 4,
                                         heights = c(0.1,1,1,1,1),
                                         widths = c(1,1,1,0.1)) %>% 
  annotate_figure(bottom = "Predicted",
                  left = "Observed")

ggsave(filename = "outputs/figures/randomforestobspred/ObservedPredicted_cforest.png",
       plot = Plot_ObservedPredicted,
       width = 20, height = 20, units = "cm", dpi = 800)

# Variable importance -----------------------------------------------------

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

#rename acoustic index names to remove '_mean' & reorder acoustic indices to order specified in 'indicesToUse'
for (df in 1:length(RandomForestImportance_cforest_conditional)) {
  names(RandomForestImportance_cforest_conditional[[df]]$values) <- gsub("_mean", "", names(RandomForestImportance_cforest_conditional[[df]]$values))
  RandomForestImportance_cforest_conditional[[df]]$data <- data.frame(AcousticIndex = names(RandomForestImportance_cforest_conditional[[df]]$values),
                                                                      Value = as.numeric(RandomForestImportance_cforest_conditional[[df]]$values))
  RandomForestImportance_cforest_conditional[[df]]$data$AcousticIndex <- factor(RandomForestImportance_cforest_conditional[[df]]$data$AcousticIndex,
                                                                                levels = indicesToUse)
}

Plots_VariableImportance <- list()
for (randomforest in 1:length(RandomForestImportance_cforest_conditional)) {
  
  Plots_VariableImportance[[names(RandomForestImportance_cforest_conditional)[randomforest]]] <- ggplot(RandomForestImportance_cforest_conditional[[randomforest]]$data, 
                                                                                                        aes(x = Value, 
                                                                                                            y = AcousticIndex)) + 
    geom_vline(xintercept = 0, linetype = "dotted") +
    geom_col() +
    scale_y_discrete(labels = axisLabels) +
    scale_x_continuous(limits = c(-0.03, 1.03)) +
    labs(y = "Acoustic Index") +
    theme_classic() +
    theme(axis.title = element_blank())
}

Plot_VariableImportance <- egg::ggarrange(as_ggplot(text_grob(label = "richness")),
                                          as_ggplot(text_grob(label = "shannon")),
                                          as_ggplot(text_grob(label = "count")),
                                          ggplot() + theme_void(),
                                          Plots_VariableImportance[['all_all_richness_totalACI']],
                                          Plots_VariableImportance[['all_all_shannon_totalACI']] + rremove("y.text"),
                                          Plots_VariableImportance[['all_all_count_totalACI']] + rremove("y.text"),
                                          as_ggplot(text_grob(label = "all vertebrates", rot = 270)),
                                          Plots_VariableImportance[['not.birds_all_richness_totalACI']],
                                          Plots_VariableImportance[['not.birds_all_shannon_totalACI']] + rremove("y.text"),
                                          Plots_VariableImportance[['not.birds_all_count_totalACI']] + rremove("y.text"),
                                          as_ggplot(text_grob(label = "non-avian", rot = 270)),
                                          Plots_VariableImportance[['birds_day_richness_totalACI']],
                                          Plots_VariableImportance[['birds_day_shannon_totalACI']] + rremove("y.text"),
                                          Plots_VariableImportance[['birds_day_count_totalACI']] + rremove("y.text"),
                                          as_ggplot(text_grob(label = "birds", rot = 270)),
                                          Plots_VariableImportance[['frogs_night_richness_totalACI']],
                                          Plots_VariableImportance[['frogs_night_shannon_totalACI']] + rremove("y.text"),
                                          Plots_VariableImportance[['frogs_night_count_totalACI']] + rremove("y.text"),
                                          as_ggplot(text_grob(label = "frogs", rot = 270)),
                                          ncol = 4,
                                          heights = c(0.1,1,1,1,1),
                                          widths = c(1,1,1,0.1)) %>% 
  annotate_figure(bottom = "Overall",
                  left = "Acoustic Index")

ggsave(filename = "outputs/figures/randomforestvariableimportance/VariableImportance_cforest.png",
       plot = Plot_VariableImportance,
       width = 20, height = 20, units = "cm", dpi = 800)

# Save workspace for later loading ----------------------------------------

save.image(file = "outputs/workspaces/RandomForest_cforest.RData")