#Random forest model of acoustic indices with vertebrate diversity - 2021

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(randomForest)
library(caret)




control <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verbose = FALSE, savePredictions = TRUE)
tunegrid <- expand.grid(.mtry=c(2:10))

for (measure in c("richness", "shannon", "count")) {
  
  formula <- as.formula(paste0(measure, " ~ ", paste(grep("*_mean", colnames(acousticIndices_richness), value = TRUE), collapse = " + ")))
  
  fit <- train(formula,
               data = acousticIndices_richness[acousticIndices_richness$type == 'birds',],
               method = "rf",
               trControl = control,
               importance = TRUE,
               allowParallel = TRUE,
               tuneGrid = tunegrid,
               ntree = 1000
  )
  fit
  plot(fit)
  varImp(fit)
}

