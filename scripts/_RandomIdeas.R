#plot biodiversity by site - how can I visualise this with index values?
ggplot(data = acousticIndices_richness[acousticIndices_richness$type == 'birds',], 
       aes(x = fct_relevel(Site, "Tarcutta", "Duval", "Mourachan", "Wambiana", "Undara", "Rinyirru"), 
           y = count)) +
  geom_point(position = position_dodge(width = 0.2)) +
  facet_wrap(~sampling.period) +
  theme_bw()


ggplot(data = acousticIndices_richness[acousticIndices_richness$type == 'birds',], 
       aes(x = fct_relevel(Site, "Tarcutta", "Duval", "Mourachan", "Wambiana", "Undara", "Rinyirru"), 
           y = richness)) +
  geom_point(position = position_dodge(width = 0.2)) +
  facet_wrap(~sampling.period) +
  theme_bw()

ggplot(data = acousticIndices_richness[acousticIndices_richness$type == 'birds',], 
       aes(x = fct_relevel(Site, "Tarcutta", "Duval", "Mourachan", "Wambiana", "Undara", "Rinyirru"), 
           y = richness)) +
  geom_boxplot() +
  facet_wrap(~sampling.period) +
  theme_bw()

ggplot(data = acousticIndices_richness[acousticIndices_richness$type == 'birds',], 
       aes(x = fct_relevel(Site, "Tarcutta", "Duval", "Mourachan", "Wambiana", "Undara", "Rinyirru"), 
           y = BI_mean)) +
  geom_boxplot() +
  facet_wrap(~sampling.period) +
  theme_bw()



#ale plots
predfun <- function(X.model, newdata) predict(X.model, as.matrix(newdata))

ALEPlot::ALEPlot(X = acousticIndices_richness[acousticIndices_richness$type == 'birds',], 
                 X.model = RandomForestFits$birds_richness$finalModel, 
                 J = c(7,10),
                 pred.fun = predfun)


#multiple correlation
test_model <- lm(richness ~ ACI_mean + BI_mean + MFC_mean + NDSI_mean, data = acousticIndices_richness[acousticIndices_richness$type == 'birds',])
cor(test_model$model$richness, test_model$fitted.values)

#heatmap of two continous variables using interpolation

library(akima)
data <- data.frame(BI = acousticIndices_richness$BI_mean[acousticIndices_richness$type == 'birds'],
                   ACI = acousticIndices_richness$ACI_mean[acousticIndices_richness$type == 'birds'],
                   richness = acousticIndices_richness$richness[acousticIndices_richness$type == 'birds'])
resolution <- 0.5 # you can increase the resolution by decreasing this number (warning: the resulting dataframe size increase very quickly)
a <- interp(x=data$BI, y=data$ACI, z=data$richness, 
            xo=seq(min(data$BI),max(data$BI),by=resolution), 
            yo=seq(min(data$ACI),max(data$ACI),by=resolution), duplicate="mean")


res <- a$z %>% 
  magrittr::set_colnames(a$y) %>% 
  as_tibble() %>% 
  mutate(x=a$x) %>% 
  gather(y, z, -x, convert=TRUE)

res %>% 
  ggplot(aes(x, y)) +
  geom_tile(aes(fill=z)) +
  geom_point(data = data, aes(x = BI, y = ACI)) +
  scale_fill_viridis_c() +
  labs(x = "Bioacoustic Index", y = "Acoustic Complexity Index") +
  theme_classic()




#ggplot version - would work well if we had more data points
ggplot(data, aes(BI, ACI, z = richness)) +
  stat_summary_2d() +
  geom_point(shape = 1, col = 'white') +
  viridis::scale_fill_viridis() +
  theme_classic()


library(factoextra)

pca.all <- prcomp(acousticIndices_richness[acousticIndices_richness$type == 'birds',c(5:17)], center = TRUE, scale. = TRUE)

fviz_eig(pca.all)

fviz_pca_ind(pca.all,
             col.ind = acousticIndices_richness$richness[acousticIndices_richness$type == 'birds'],
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

fviz_pca_var(pca.all,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             title = NULL,
             ggtheme = theme_classic()) + ggtitle(NULL) + coord_cartesian() + theme(legend.position = c(0.3, 0.3), legend.background = element_blank(), legend.key.width = unit(0.5, "cm"), legend.key.size = unit(0.5, "cm"))








# quantile regression random forest for prediction intervals

library(quantregForest)

control <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verbose = FALSE, savePredictions = TRUE)
tunegrid <- expand.grid(.mtry=c(2:10))


RandomForestQuantRegFits <- list()
RandomForestQuantRegPerformance <- data.frame(taxa = character(),
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
RandomForestQuantRegImportance <- list()

for (taxa in c('all', 'not.birds', 'birds', 'frogs')) {
  tmpdata <- acousticIndices_richness[acousticIndices_richness$type == taxa,]
  
  for (measure in c("richness", "shannon", "count")) {
    
    formulas <- list(totalACI = as.formula(paste0(measure, " ~ ", paste(grep("ACI_[0-9].", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + "))),
                     ACI_3kHz = as.formula(paste0(measure, " ~ ", paste(grep("ACI_mean$|ACI_1000_2.|ACI_2000_3.|ACI_3000_4.|ACI_4000_5.|ACI_5000_6.|ACI_6000_7.|ACI_7000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + "))))
    
    #formula <- as.formula(paste0(measure, " ~ ", paste(grep("ACI_[0-9].", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + ")))
    #formula <- as.formula(paste0(measure, " ~ ", paste(grep("ACI_mean$|ACI_1000_2.|ACI_2000_3.|ACI_3000_4.|ACI_4000_5.|ACI_5000_6.|ACI_6000_7.|ACI_7000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + ")))
    for (formula in 1:length(formulas)) {
      
      RandomForestQuantRegFits[[paste0(taxa, "_", measure, "_", names(formulas)[[formula]])]] <- train(formulas[[formula]],
                                                                                                       data = tmpdata,
                                                                                                       method = "qrf",
                                                                                                       trControl = control,
                                                                                                       importance = TRUE,
                                                                                                       allowParallel = TRUE,
                                                                                                       tuneGrid = tunegrid,
                                                                                                       ntree = 1000,
                                                                                                       keep.inbag = TRUE)
      
      RandomForestQuantRegPerformance <- rbind(RandomForestQuantRegPerformance, data.frame(cbind(taxa = taxa, 
                                                                                                 measure = measure, 
                                                                                                 ACItype = names(formulas)[[formula]],
                                                                                                 RandomForestQuantRegFits[[paste0(taxa, "_", measure, "_", names(formulas)[[formula]])]]$results[RandomForestQuantRegFits[[paste0(taxa, "_", measure, "_", names(formulas)[[formula]])]]$results$mtry == RandomForestQuantRegFits[[paste0(taxa, "_", measure, "_", names(formulas)[[formula]])]]$bestTune[[1]],],
                                                                                                 minResponse = min(tmpdata[measure]),
                                                                                                 meanResponse = mean(tmpdata[[measure]]),
                                                                                                 maxResponse = max(tmpdata[measure]))))
      
      RandomForestQuantRegImportance[[paste0(taxa, "_", measure, "_", names(formulas)[[formula]])]] <- varImp(RandomForestQuantRegFits[[paste0(taxa, "_", measure, "_", names(formulas)[[formula]])]])
    }
  }
}

#Normalize RMSE and MAE, calculate Scatter Index
RandomForestQuantRegPerformance <- RandomForestQuantRegPerformance %>% mutate(normRMSE = RMSE/(maxResponse-minResponse),
                                                                              normMAE = MAE/(maxResponse-minResponse),
                                                                              SI = (RMSE/meanResponse)*100)




#observed vs predicted plot
test <- data.frame(predictions = predict(RandomForestQuantRegFits$birds_count_totalACI$finalModel),
                   observations = acousticIndices_richness$count[acousticIndices_richness$type == 'birds'])

ggplot(test, aes(x = observations, y = predictions)) +
  geom_point() +
  scale_x_continuous(limits = c(min(test), max(test))) +
  scale_y_continuous(limits = c(min(test), max(test))) +
  geom_abline(slope = 1) +
  theme_bw()




test2 <- data.frame(cbind(predict(RandomForestQuantRegFits$birds_count_totalACI$finalModel),
                          observations = acousticIndices_richness$count[acousticIndices_richness$type == 'birds']))

ggplot(test2, aes(x = observations, y = quantile..0.5)) +
  geom_ribbon(aes(x = observations, ymin = quantile..0.1, ymax = quantile..0.9)) +
  geom_point() +
  scale_x_continuous(limits = c(min(test), max(test))) +
  scale_y_continuous(limits = c(min(test), max(test))) +
  geom_abline(slope = 1) +
  theme_bw()