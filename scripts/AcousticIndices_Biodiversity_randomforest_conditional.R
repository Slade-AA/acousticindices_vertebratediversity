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

pb = txtProgressBar(min = 0, max = length(unique(acousticIndices_richness$type)) * 3 * 3, initial = 0, style = 3); k <- 0

for (comparison in unique(acousticIndices_richness$type)) {
  tmpdata <- acousticIndices_richness[acousticIndices_richness$type == comparison,]
  
  for (measure in c("richness", "shannon", "count")) {
    
    #forumlas for using just total ACI, 3kHz ACI, and 1kHz ACI values
    formulas <- list(totalACI = as.formula(paste0(measure, " ~ ", paste(grep("ACI_[0-9].", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + "))),
                     ACI_3kHz = as.formula(paste0(measure, " ~ ", paste(grep("ACI_mean$|ACI_1000_2.|ACI_2000_3.|ACI_3000_4.|ACI_4000_5.|ACI_5000_6.|ACI_6000_7.|ACI_7000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + "))),
                     ACI_1kHz = as.formula(paste0(measure, " ~ ", paste(grep("ACI_mean$|ACI_1000_4.|ACI_3000_6.|ACI_5000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + "))))
    
    for (formula in 1:length(formulas)) {
      
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
      RandomForestImportance_cforest_conditional[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]] <- permimp(RandomForestFits_cforest[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$finalModel, conditional = TRUE, progressBar = FALSE)
      
      k <- k+1; setTxtProgressBar(pb, k)
      
    }
  }
}

#Normalize RMSE and MAE, calculate Scatter Index
RandomForestPerformance_cforest <- RandomForestPerformance_cforest %>% mutate(normRMSE = RMSE/(maxResponse-minResponse),
                                                                              normMAE = MAE/(maxResponse-minResponse),
                                                                              SI = (RMSE/meanResponse)*100)