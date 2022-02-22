#Random forest model of acoustic indices with vertebrate diversity - 2021

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(randomForest)
library(caret)

# Load indices and biodiversity data --------------------------------------

#summary indices
files <- file.info(list.files("./outputs/data/", pattern = ".*_acousticIndices_summary.RData$", full.names = TRUE)) #list files
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
                                                                           birds_morning = c('morning', 'birds'),
                                                                           birds_afternoon = c('afternoon', 'birds'),
                                                                           birds_day = c('day', 'birds'),
                                                                           birds_all = c('all', 'birds'),
                                                                           frogs_evening = c('evening', 'frogs'),
                                                                           frogs_night = c('night', 'frogs')))

# Fit random forest models ------------------------------------------------

control <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verbose = FALSE, savePredictions = TRUE)
tunegrid <- expand.grid(.mtry=c(2:10))

RandomForestFits <- list()
RandomForestPerformance <- data.frame(comparison = character(),
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
RandomForestImportance <- list()

RandomForestPredictions <- list()

indicesToUse <- c('ACI', 'ADI', 'AEI', 'BI', 'NDSI', 'EVN', 'SH', 'LFC', 'MFC', 'HFC')

pb = txtProgressBar(min = 0, max = length(unique(acousticIndices_richness$type)) * 3 * 3, initial = 0, style = 3); k <- 0

for (comparison in unique(acousticIndices_richness$type)) {
  tmpdata <- acousticIndices_richness[acousticIndices_richness$type == comparison,]
  
  for (measure in c("richness", "shannon", "count")) {
    
    #forumlas for using just total ACI, 3kHz ACI, and 1kHz ACI values
    formulas <- list(totalACI = as.formula(paste0(measure, " ~ ", paste(grep(paste(paste0(indicesToUse, "_*"), collapse = "|"), grep("ACI_[0-9].", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), value = TRUE), collapse = " + "))),
                     ACI_3kHz = as.formula(paste0(measure, " ~ ", paste(grep(paste(paste0(indicesToUse, "_*"), collapse = "|"), grep("ACI_mean$|ACI_1000_2.|ACI_2000_3.|ACI_3000_4.|ACI_4000_5.|ACI_5000_6.|ACI_6000_7.|ACI_7000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), value = TRUE), collapse = " + "))),
                     ACI_1kHz = as.formula(paste0(measure, " ~ ", paste(grep(paste(paste0(indicesToUse, "_*"), collapse = "|"), grep("ACI_mean$|ACI_1000_4.|ACI_3000_6.|ACI_5000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), value = TRUE), collapse = " + "))))
    
    for (formula in 1:length(formulas)) {
      set.seed(1234)#set seed for reproducibility
      RandomForestFits[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]] <- train(formulas[[formula]],
                                                                                               data = tmpdata,
                                                                                               method = "rf",
                                                                                               trControl = control,
                                                                                               importance = TRUE,
                                                                                               allowParallel = TRUE,
                                                                                               tuneGrid = tunegrid,
                                                                                               ntree = 1000)
      
      RandomForestPerformance <- rbind(RandomForestPerformance, data.frame(cbind(comparison = comparison, 
                                                                                 measure = measure, 
                                                                                 ACItype = names(formulas)[[formula]],
                                                                                 RandomForestFits[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$results[RandomForestFits[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$results$mtry == RandomForestFits[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$bestTune[[1]],],
                                                                                 minResponse = min(tmpdata[measure]),
                                                                                 meanResponse = mean(tmpdata[[measure]]),
                                                                                 maxResponse = max(tmpdata[measure]))))
      
      RandomForestImportance[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]] <- varImp(RandomForestFits[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]])
      
      RandomForestPredictions[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]] <- data.frame(predictions = predict(RandomForestFits[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$finalModel),
                                                                                                                 observations = tmpdata[[measure]],
                                                                                                                 comparison = paste0(comparison, "_", measure),
                                                                                                                 ACItype = names(formulas)[[formula]],
                                                                                                                 stringsAsFactors = FALSE)
      
      
      k <- k+1; setTxtProgressBar(pb, k)
    }
  }
}

#Normalize RMSE and MAE, calculate Scatter Index
RandomForestPerformance <- RandomForestPerformance %>% mutate(normRMSE = RMSE/(maxResponse-minResponse),
                                                              normMAE = MAE/(maxResponse-minResponse),
                                                              SI = (RMSE/meanResponse)*100,
                                                              normRMSE_SE = RMSESD/(maxResponse-minResponse)/sqrt(30),
                                                              normMAE_SE = MAESD/(maxResponse-minResponse)/sqrt(30),
                                                              SI_SE = (RMSESD/meanResponse)*100/sqrt(30),
                                                              Rsquared_SE = RsquaredSD/sqrt(30))

RandomForestPredictions <- do.call(rbind, RandomForestPredictions)

# Plot random forest performance ------------------------------------------

#compare rf models with different ACI bands used as predictors
ggplot(data = RandomForestPerformance, aes(x = measure, y = normRMSE, group = ACItype, fill = ACItype)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = normRMSE-normRMSE_SE, ymax = normRMSE+normRMSE_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised RMSE") +
  facet_wrap(~comparison) +
  theme_bw()

#compare bird indices comparisons - using totalACI rf models
Plot_RMSE_birds <- ggplot(data = RandomForestPerformance[RandomForestPerformance$ACItype == 'totalACI' & RandomForestPerformance$comparison %in% grep("^birds*", RandomForestPerformance$comparison, value = TRUE),], 
                          aes(x = measure, y = normRMSE, group = comparison, fill = comparison)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = normRMSE-normRMSE_SE, ymax = normRMSE+normRMSE_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised RMSE") +
  theme_bw() +
  theme(legend.position = "none")
Plot_SI_birds <- ggplot(data = RandomForestPerformance[RandomForestPerformance$ACItype == 'totalACI' & RandomForestPerformance$comparison %in% grep("^birds*", RandomForestPerformance$comparison, value = TRUE),], 
                        aes(x = measure, y = SI, group = comparison, fill = comparison)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = SI-SI_SE, ymax = SI+SI_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Taxa", y = "Scatter Index") +
  theme_bw() +
  theme(legend.position = "none")
Plot_R2_birds <- ggplot(data = RandomForestPerformance[RandomForestPerformance$ACItype == 'totalACI' & RandomForestPerformance$comparison %in% grep("^birds*", RandomForestPerformance$comparison, value = TRUE),], 
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
ggsave(filename = "outputs/figures/randomforestperformance/Birds_MorningDayAll.png",
       plot = Plot_Birds_Time,
       width = 24, height = 20, units = "cm", dpi = 800)


#plot performance measures for all comparisons - totalACI
Plot_RMSE <- ggplot(data = RandomForestPerformance[RandomForestPerformance$ACItype == 'totalACI' &
                                                     RandomForestPerformance$comparison %in% c('all_all', 'not.birds_all', 'birds_morning', 'frogs_evening'),], 
                    aes(x = comparison, y = normRMSE, group = measure, fill = measure)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = normRMSE-normRMSE_SE, ymax = normRMSE+normRMSE_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised RMSE") +
  theme_bw() +
  theme(legend.position = "none")
Plot_MAE <- ggplot(data = RandomForestPerformance[RandomForestPerformance$ACItype == 'totalACI' &
                                                    RandomForestPerformance$comparison %in% c('all_all', 'not.birds_all', 'birds_morning', 'frogs_evening'),], 
                   aes(x = comparison, y = normMAE, group = measure, fill = measure)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = normMAE-normMAE_SE, ymax = normMAE+normMAE_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised MAE") +
  theme_bw() +
  theme(legend.position = "none")
Plot_SI <- ggplot(data = RandomForestPerformance[RandomForestPerformance$ACItype == 'totalACI' &
                                                   RandomForestPerformance$comparison %in% c('all_all', 'not.birds_all', 'birds_morning', 'frogs_evening'),], 
                  aes(x = comparison, y = SI, group = measure, fill = measure)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = SI-SI_SE, ymax = SI+SI_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 120)) +
  labs(x = "Taxa", y = "Scatter Index") +
  theme_bw() +
  theme(legend.position = "none")
Plot_R2 <- ggplot(data = RandomForestPerformance[RandomForestPerformance$ACItype == 'totalACI' &
                                                   RandomForestPerformance$comparison %in% c('all_all', 'not.birds_all', 'birds_morning', 'frogs_evening'),], 
                  aes(x = comparison, y = Rsquared, group = measure, fill = measure)) +
  geom_col(position = "dodge", color = "black") +
  geom_errorbar(aes(ymin = Rsquared-Rsquared_SE, ymax = Rsquared+Rsquared_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "R Squared") +
  theme_bw() +
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
ggsave(filename = "outputs/figures/randomforestperformance/AllComparisons.png",
       plot = Plot_AllRF,
       width = 24, height = 20, units = "cm", dpi = 800)


# Observed vs predicted plots ---------------------------------------------

Plots_ObsPred <- list()

for (predictions in c('all_all_richness', 'all_all_shannon', 'all_all_count',
                      'not.birds_all_richness', 'not.birds_all_shannon', 'not.birds_all_count',
                      'birds_morning_richness', 'birds_morning_shannon', 'birds_morning_count',
                      'frogs_night_richness', 'frogs_night_shannon', 'frogs_night_count')) {
  Plots_ObsPred[[paste0(predictions)]] <- ggplot(RandomForestPredictions[RandomForestPredictions$comparison == predictions &
                                                                          RandomForestPredictions$ACItype == 'totalACI',], 
                                                aes(x = predictions, y = observations)) +
    geom_abline(slope = 1, linetype = 'dashed') +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point() +
    scale_x_continuous(limits = c(min(RandomForestPredictions[RandomForestPredictions$comparison == predictions &
                                                                RandomForestPredictions$ACItype == 'totalACI',c(1:2)]), 
                                  max(RandomForestPredictions[RandomForestPredictions$comparison == predictions &
                                                                RandomForestPredictions$ACItype == 'totalACI',c(1:2)]))) +
    scale_y_continuous(limits = c(min(RandomForestPredictions[RandomForestPredictions$comparison == predictions &
                                                                RandomForestPredictions$ACItype == 'totalACI',c(1:2)]), 
                                  max(RandomForestPredictions[RandomForestPredictions$comparison == predictions &
                                                                RandomForestPredictions$ACItype == 'totalACI',c(1:2)]))) +
    theme_classic() +
    theme(axis.title = element_blank())
}

plot_grid(plotlist = Plots_ObsPred)

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
                                         Plots_ObsPred[['birds_morning_richness']],
                                         Plots_ObsPred[['birds_morning_shannon']],
                                         Plots_ObsPred[['birds_morning_count']],
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

ggsave(filename = "outputs/figures/randomforestobspred/ObservedPredicted.png",
       plot = Plot_ObservedPredicted,
       width = 20, height = 20, units = "cm", dpi = 800)

# Variable importance -----------------------------------------------------

#rename acoustic index names to remove '_mean'
for (df in 1:length(RandomForestImportance)) {
  RandomForestImportance[[df]]$importance$AcousticIndex <- gsub("_mean", "", rownames(RandomForestImportance[[df]]$importance))
  RandomForestImportance[[df]]$importance$AcousticIndex <- factor(RandomForestImportance[[df]]$importance$AcousticIndex,
                                                                  levels = c(indicesToUse, setdiff(unique(RandomForestImportance[[df]]$importance$AcousticIndex), indicesToUse)))
}

Plots_VariableImportance <- list()
for (randomforest in 1:length(RandomForestImportance)) {
  
  Plots_VariableImportance[[names(RandomForestImportance)[randomforest]]] <- ggplot(RandomForestImportance[[randomforest]]$importance, 
                                                                                    aes(x = Overall, 
                                                                                        y = AcousticIndex)) + 
    geom_col() + 
    labs(y = "Acoustic Index") +
    theme_bw() +
    theme(axis.title = element_blank())
}

#arrange plot using egg package - even sized plots
Plot_VariableImportance <- egg::ggarrange(as_ggplot(text_grob(label = "richness")),
                                          as_ggplot(text_grob(label = "shannon")),
                                          as_ggplot(text_grob(label = "count")),
                                          ggplot() + theme_void(),
                                          Plots_VariableImportance[['all_all_richness_totalACI']] + rremove("x.text"),
                                          Plots_VariableImportance[['all_all_shannon_totalACI']] + rremove("axis.text"),
                                          Plots_VariableImportance[['all_all_count_totalACI']] + rremove("axis.text"),
                                          as_ggplot(text_grob(label = "all vertebrates", rot = 270)),
                                          Plots_VariableImportance[['not.birds_all_richness_totalACI']] + rremove("x.text"),
                                          Plots_VariableImportance[['not.birds_all_shannon_totalACI']] + rremove("axis.text"),
                                          Plots_VariableImportance[['not.birds_all_count_totalACI']] + rremove("axis.text"),
                                          as_ggplot(text_grob(label = "non-avian", rot = 270)),
                                          Plots_VariableImportance[['birds_morning_richness_totalACI']] + rremove("x.text"),
                                          Plots_VariableImportance[['birds_morning_shannon_totalACI']] + rremove("axis.text"),
                                          Plots_VariableImportance[['birds_morning_count_totalACI']] + rremove("axis.text"),
                                          as_ggplot(text_grob(label = "birds", rot = 270)),
                                          Plots_VariableImportance[['frogs_evening_richness_totalACI']],
                                          Plots_VariableImportance[['frogs_evening_shannon_totalACI']] + rremove("y.text"),
                                          Plots_VariableImportance[['frogs_evening_count_totalACI']] + rremove("y.text"),
                                          as_ggplot(text_grob(label = "frogs", rot = 270)),
                                          ncol = 4,
                                          heights = c(0.1,1,1,1,1),
                                          widths = c(1,1,1,0.1)) %>% 
  annotate_figure(bottom = "Overall",
                  left = "Acoustic Index")

ggsave(filename = "outputs/figures/randomforestvariableimportance/VariableImportance.png",
       plot = Plot_VariableImportance,
       width = 20, height = 20, units = "cm", dpi = 800)

#save everything in environment for later loading
save.image(file = "outputs/workspaces/RandomForest.RData")