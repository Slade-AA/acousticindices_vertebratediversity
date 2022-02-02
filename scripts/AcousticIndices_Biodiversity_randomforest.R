#Random forest model of acoustic indices with vertebrate diversity - 2021

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(randomForest)
library(caret)

#WHAT ABOUT XGBOOST!!!!

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


# Fit random forest models ------------------------------------------------

control <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verbose = FALSE, savePredictions = TRUE)
tunegrid <- expand.grid(.mtry=c(2:10))

RandomForestFits <- list()
RandomForestPerformance <- data.frame(taxa = character(),
                                      measure = character(),
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

for (taxa in c('all', 'not.birds', 'birds', 'frogs')) {
  tmpdata <- acousticIndices_richness[acousticIndices_richness$type == taxa,]
  
  for (measure in c("richness", "shannon", "count")) {
    
    formula <- as.formula(paste0(measure, " ~ ", paste(grep("*_mean", colnames(acousticIndices_richness), value = TRUE), collapse = " + ")))
    
    RandomForestFits[[paste0(taxa, "_", measure)]] <- train(formula,
                                                            data = tmpdata,
                                                            method = "rf",
                                                            trControl = control,
                                                            importance = TRUE,
                                                            allowParallel = TRUE,
                                                            tuneGrid = tunegrid,
                                                            ntree = 1000)
    
    RandomForestPerformance <- rbind(RandomForestPerformance, data.frame(cbind(taxa = taxa, 
                                                                               measure = measure, 
                                                                               RandomForestFits[[paste0(taxa, "_", measure)]]$results[RandomForestFits[[paste0(taxa, "_", measure)]]$results$mtry == RandomForestFits[[paste0(taxa, "_", measure)]]$bestTune[[1]],],
                                                                               minResponse = min(tmpdata[measure]),
                                                                               meanResponse = mean(tmpdata[[measure]]),
                                                                               maxResponse = max(tmpdata[measure]))))
    
    RandomForestImportance[[paste0(taxa, "_", measure)]] <- varImp(RandomForestFits[[paste0(taxa, "_", measure)]])
  }
}

#Normalize RMSE and MAE, calculate Scatter Index
RandomForestPerformance <- RandomForestPerformance %>% mutate(normRMSE = RMSE/(maxResponse-minResponse),
                                                              normMAE = MAE/(maxResponse-minResponse),
                                                              SI = (RMSE/meanResponse)*100)

# Plot random forest performance ------------------------------------------

Plot_RMSE <- ggplot(data = RandomForestPerformance, aes(x = taxa, y = normRMSE, group = measure, fill = measure)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised RMSE") +
  theme_bw() +
  theme(legend.position = "none")
Plot_MAE <- ggplot(data = RandomForestPerformance, aes(x = taxa, y = normMAE, group = measure, fill = measure)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised MAE") +
  theme_bw() +
  theme(legend.position = "none")
Plot_R2 <- ggplot(data = RandomForestPerformance, aes(x = taxa, y = Rsquared, group = measure, fill = measure)) +
  geom_col(position = "dodge") +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "R Squared") +
  theme_bw() +
  theme(legend.position = "none")
Plot_SI <- ggplot(data = RandomForestPerformance, aes(x = taxa, y = SI, group = measure, fill = measure)) +
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
  rownames(RandomForestImportance[[df]]$importance) <- gsub("_mean", "", rownames(RandomForestImportance[[df]]$importance))
}

Plots_VariableImportance <- list()
for (randomforest in 1:length(RandomForestImportance)) {
  
  Plots_VariableImportance[[names(RandomForestImportance)[randomforest]]] <- ggplot(RandomForestImportance[[randomforest]]$importance, 
                                                                                    aes(x = Overall, 
                                                                                        y = rownames(RandomForestImportance[[randomforest]]$importance))) + 
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
               Plots_VariableImportance[[2]] + rremove("axis.text"),
               Plots_VariableImportance[[3]] + rremove("axis.text"),
               as_ggplot(text_grob(label = "all", rot = 270)),
               Plots_VariableImportance[[4]] + rremove("x.text"),
               Plots_VariableImportance[[5]] + rremove("axis.text"),
               Plots_VariableImportance[[6]] + rremove("axis.text"),
               as_ggplot(text_grob(label = "not.birds", rot = 270)),
               Plots_VariableImportance[[7]] + rremove("x.text"),
               Plots_VariableImportance[[8]] + rremove("axis.text"),
               Plots_VariableImportance[[9]] + rremove("axis.text"),
               as_ggplot(text_grob(label = "birds", rot = 270)),
               Plots_VariableImportance[[10]],
               Plots_VariableImportance[[11]] + rremove("y.text"),
               Plots_VariableImportance[[12]] + rremove("y.text"),
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
RandomForestPerformance_cforest <- data.frame(taxa = character(),
                                              measure = character(),
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

for (taxa in c('all', 'not.birds', 'birds', 'frogs')) {
  tmpdata <- acousticIndices_richness[acousticIndices_richness$type == taxa,]
  
  for (measure in c("richness", "shannon", "count")) {
    
    formula <- as.formula(paste0(measure, " ~ ", paste(grep("*_mean", colnames(acousticIndices_richness), value = TRUE), collapse = " + ")))
    
    RandomForestFits_cforest[[paste0(taxa, "_", measure)]] <- train(formula,
                                                                    data = tmpdata,
                                                                    method = "cforest",
                                                                    trControl = control,
                                                                    tuneGrid = tunegrid,
                                                                    controls = cforest_unbiased(ntree = 1000))
    
    RandomForestPerformance_cforest <- rbind(RandomForestPerformance_cforest, data.frame(cbind(taxa = taxa, 
                                                                                               measure = measure, 
                                                                                               RandomForestFits_cforest[[paste0(taxa, "_", measure)]]$results[RandomForestFits_cforest[[paste0(taxa, "_", measure)]]$results$mtry == RandomForestFits_cforest[[paste0(taxa, "_", measure)]]$bestTune[[1]],],
                                                                                               minResponse = min(tmpdata[measure]),
                                                                                               meanResponse = mean(tmpdata[[measure]]),
                                                                                               maxResponse = max(tmpdata[measure]))))
    
    RandomForestImportance_cforest[[paste0(taxa, "_", measure)]] <- varImp(RandomForestFits_cforest[[paste0(taxa, "_", measure)]])
    RandomForestImportance_cforest_conditional[[paste0(taxa, "_", measure)]] <- permimp(RandomForestFits_cforest[[paste0(taxa, "_", measure)]]$finalModel, conditional = TRUE)
  }
}



cf_model <- cforest(formula,
                    data = tmpdata,
                    controls = cforest_unbiased(mtry = 4, ntree = 1000, minbucket = 1, minsplit = 3))
CPI <- permimp(cf_model, conditional = TRUE)
#plot(CPI, type = "bar", interval = "quantile")
plot(CPI, type = "box", horizontal = TRUE)

cforestStats(cf_model)