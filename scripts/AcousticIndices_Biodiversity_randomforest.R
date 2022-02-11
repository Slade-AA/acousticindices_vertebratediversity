#Random forest model of acoustic indices with vertebrate diversity - 2021

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(randomForest)
library(caret)

#WHAT ABOUT XGBOOST!!!!

# Load indices and biodiversity data --------------------------------------

setwd("C:/Users/jc696551/OneDrive - James Cook University/Projects/acousticindices_vertebratediversity")

#summary indices
files <- file.info(list.files("./outputs/data/", pattern = ".*_acousticIndices_summary.RData$", full.names = TRUE)) #list files
latestFile <- rownames(files)[which.max(files$mtime)] #determine most recent file to use for loading

load(latestFile)

#richness
richness <- read_csv("./rawdata/biodiversity/richness_diversity_updated.csv") #load richness

#format richness to match indices
richness <- richness %>% rename(Site = site, Sensor = original.plot, richness = Richness, shannon = Shannon, type = comparison)
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

for (comparison in unique(acousticIndices_richness$type)) {
  tmpdata <- acousticIndices_richness[acousticIndices_richness$type == comparison,]
  
  for (measure in c("richness", "shannon", "count")) {
    
    formulas <- list(totalACI = as.formula(paste0(measure, " ~ ", paste(grep("ACI_[0-9].", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + "))),
                     ACI_3kHz = as.formula(paste0(measure, " ~ ", paste(grep("ACI_mean$|ACI_1000_2.|ACI_2000_3.|ACI_3000_4.|ACI_4000_5.|ACI_5000_6.|ACI_6000_7.|ACI_7000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + "))))
    
    #formula <- as.formula(paste0(measure, " ~ ", paste(grep("ACI_[0-9].", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + ")))
    #formula <- as.formula(paste0(measure, " ~ ", paste(grep("ACI_mean$|ACI_1000_2.|ACI_2000_3.|ACI_3000_4.|ACI_4000_5.|ACI_5000_6.|ACI_6000_7.|ACI_7000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + ")))
    for (formula in 1:length(formulas)) {
      
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
    }
  }
}

#Normalize RMSE and MAE, calculate Scatter Index
RandomForestPerformance <- RandomForestPerformance %>% mutate(normRMSE = RMSE/(maxResponse-minResponse),
                                                              normMAE = MAE/(maxResponse-minResponse),
                                                              SI = (RMSE/meanResponse)*100)

# Plot random forest performance ------------------------------------------

#compare ACI type
ggplot(data = RandomForestPerformance, aes(x = measure, y = normRMSE, group = ACItype, fill = ACItype)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised RMSE") +
  facet_wrap(~comparison) +
  theme_bw()


Plot_RMSE <- ggplot(data = RandomForestPerformance[RandomForestPerformance$ACItype == 'totalACI',], aes(x = comparison, y = normRMSE, group = measure, fill = measure)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised RMSE") +
  theme_bw() +
  theme(legend.position = "none")
Plot_MAE <- ggplot(data = RandomForestPerformance[RandomForestPerformance$ACItype == 'totalACI',], aes(x = comparison, y = normMAE, group = measure, fill = measure)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised MAE") +
  theme_bw() +
  theme(legend.position = "none")
Plot_R2 <- ggplot(data = RandomForestPerformance[RandomForestPerformance$ACItype == 'totalACI',], aes(x = comparison, y = Rsquared, group = measure, fill = measure)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "R Squared") +
  theme_bw() +
  theme(legend.position = "none")
Plot_SI <- ggplot(data = RandomForestPerformance[RandomForestPerformance$ACItype == 'totalACI',], aes(x = comparison, y = SI, group = measure, fill = measure)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 100)) +
  labs(x = "Taxa", y = "Scatter Index") +
  theme_bw() +
  theme(legend.position = "none")

legend_bottom <- get_legend(
  Plot_R2 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank())
)

plot_grid(Plot_MAE, Plot_R2, legend_bottom,
          ncol = 1,
          rel_heights = c(1, 1, 0.1))

# Variable importance -----------------------------------------------------

#rename acoustic index names to remove '_mean'
for (df in 1:length(RandomForestImportance)) {
  RandomForestImportance[[df]]$importance$AcousticIndex <- gsub("_mean", "", rownames(RandomForestImportance[[df]]$importance))
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
egg::ggarrange(as_ggplot(text_grob(label = "richness")),
               as_ggplot(text_grob(label = "shannon")),
               as_ggplot(text_grob(label = "count")),
               ggplot() + theme_void(),
               Plots_VariableImportance[[1]] + rremove("x.text"),
               Plots_VariableImportance[[3]] + rremove("axis.text"),
               Plots_VariableImportance[[5]] + rremove("axis.text"),
               as_ggplot(text_grob(label = "all", rot = 270)),
               Plots_VariableImportance[[7]] + rremove("x.text"),
               Plots_VariableImportance[[9]] + rremove("axis.text"),
               Plots_VariableImportance[[11]] + rremove("axis.text"),
               as_ggplot(text_grob(label = "not.birds", rot = 270)),
               Plots_VariableImportance[[13]] + rremove("x.text"),
               Plots_VariableImportance[[15]] + rremove("axis.text"),
               Plots_VariableImportance[[17]] + rremove("axis.text"),
               as_ggplot(text_grob(label = "birds", rot = 270)),
               Plots_VariableImportance[[19]],
               Plots_VariableImportance[[21]] + rremove("y.text"),
               Plots_VariableImportance[[23]] + rremove("y.text"),
               as_ggplot(text_grob(label = "frogs", rot = 270)),
               ncol = 4,
               heights = c(0.1,1,1,1,1),
               widths = c(1,1,1,0.1)) %>% 
  annotate_figure(bottom = "Overall",
                  left = "Acoustic Index")




# Random forest using party and conditional importance --------------------

library(party)
library(permimp)

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

for (comparison in unique(acousticIndices_richness$type)) {
  tmpdata <- acousticIndices_richness[acousticIndices_richness$type == comparison,]
  
  for (measure in c("richness", "shannon", "count")) {
    
    formulas <- list(totalACI = as.formula(paste0(measure, " ~ ", paste(grep("ACI_[0-9].", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + "))),
                     ACI_3kHz = as.formula(paste0(measure, " ~ ", paste(grep("ACI_mean$|ACI_1000_2.|ACI_2000_3.|ACI_3000_4.|ACI_4000_5.|ACI_5000_6.|ACI_6000_7.|ACI_7000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + "))))
    
    #formula <- as.formula(paste0(measure, " ~ ", paste(grep("ACI_[0-9].", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + ")))
    #formula <- as.formula(paste0(measure, " ~ ", paste(grep("ACI_mean$|ACI_1000_2.|ACI_2000_3.|ACI_3000_4.|ACI_4000_5.|ACI_5000_6.|ACI_6000_7.|ACI_7000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + ")))
    for (formula in 1:length(formulas)) {
      
      RandomForestFits_cforest[[paste0(comparison, "_", measure)]] <- train(formulas[[formula]],
                                                                      data = tmpdata,
                                                                      method = "cforest",
                                                                      trControl = control,
                                                                      tuneGrid = tunegrid,
                                                                      controls = cforest_unbiased(ntree = 1000))
      
      RandomForestPerformance_cforest <- rbind(RandomForestPerformance_cforest, data.frame(cbind(comparison = comparison, 
                                                                                                 measure = measure, 
                                                                                                 ACItype = names(formulas)[[formula]],
                                                                                                 RandomForestFits_cforest[[paste0(comparison, "_", measure)]]$results[RandomForestFits_cforest[[paste0(comparison, "_", measure)]]$results$mtry == RandomForestFits_cforest[[paste0(comparison, "_", measure)]]$bestTune[[1]],],
                                                                                                 minResponse = min(tmpdata[measure]),
                                                                                                 meanResponse = mean(tmpdata[[measure]]),
                                                                                                 maxResponse = max(tmpdata[measure]))))
      
      RandomForestImportance_cforest[[paste0(comparison, "_", measure)]] <- varImp(RandomForestFits_cforest[[paste0(comparison, "_", measure)]])
      RandomForestImportance_cforest_conditional[[paste0(comparison, "_", measure)]] <- permimp(RandomForestFits_cforest[[paste0(comparison, "_", measure)]]$finalModel, conditional = TRUE)
    }
  }
}



cf_model <- cforest(formula,
                    data = tmpdata,
                    controls = cforest_unbiased(mtry = 4, ntree = 1000, minbucket = 1, minsplit = 3))
CPI <- permimp(cf_model, conditional = TRUE)
#plot(CPI, type = "bar", interval = "quantile")
plot(CPI, type = "box", horizontal = TRUE)

cforestStats(cf_model)