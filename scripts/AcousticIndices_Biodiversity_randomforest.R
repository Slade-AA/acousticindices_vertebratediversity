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
                                                                               meanResponse = mean(tmpdata[measure])
                                                                               maxResponse = max(tmpdata[measure]))))
    
    RandomForestImportance[[paste0(taxa, "_", measure)]] <- varImp(RandomForestFits[[paste0(taxa, "_", measure)]])
  }
}


RandomForestPerformance <- RandomForestPerformance %>% mutate(normRMSE = RMSE/(maxResponse-minResponse),
                                                              normMAE = MAE/(maxResponse-minResponse),
                                                              SI = (RMSE/meanResponse)*100)

# tt ----------------------------------------------------------------------


res <- resamples(list(all_richness = RandomForestFits$all_richness, 
                      all_shannon = RandomForestFits$all_shannon, 
                      all_count = RandomForestFits$all_count))

summary(res)



Plots_VariableImportance <- list()
for (randomforest in 1:length(RandomForestImportance)) {
  
  Plots_VariableImportance[[names(RandomForestImportance)[randomforest]]] <- ggplot(RandomForestImportance[[randomforest]]$importance, 
                                                                                    aes(x = Overall, 
                                                                                        y = rownames(RandomForestImportance[[randomforest]]$importance))) + 
    geom_col() + 
    labs(y = "Acoustic Index") +
    theme_bw()
}

plot_grid(plotlist = Plots_VariableImportance)




library(party)
library(permimp)
cf_model <- cforest(formula,
                    data = tmpdata,
                    controls = cforest_unbiased(mtry = 4, ntree = 1000, minbucket = 1, minsplit = 3))
CPI <- permimp(cf_model, conditional = TRUE)
#plot(CPI, type = "bar", interval = "quantile")
plot(CPI, type = "box", horizontal = TRUE)

cforestStats(cf_model)