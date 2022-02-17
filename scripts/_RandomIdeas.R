#Look at correlation between mean and median indices summarised for different time periods
library(ggcorrplot)
library(cowplot)
library(ggpubr)

Plot_MeanMedian_Cor <- list()
IndiceCorrelations <- list()
for (type in unique(acousticIndices_summary$type)) {
  tmp <- acousticIndices_summary[acousticIndices_summary$type == type,]
  
  res <- cor(tmp[c(5:17, 28:40)])
  IndiceCorrelations[[paste0(type)]] <- res[1:13,14:26] #only keep comparisons between mean and median - i.e. no self comparisons
  Plot_MeanMedian_Cor[[paste0(type)]] <- ggcorrplot(IndiceCorrelations[[paste0(type)]], 
                                                    method = "circle", 
                                                    colors = c("#6D9EC1", "white", "#E46726"),
                                                    outline.color = "black") + 
    scale_size(range = c(1, 6)) +
    theme(legend.position = "none")
}

legend_bottom <- get_legend(
  Plot_MeanMedian_Cor[[paste0(type)]] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.direction = "horizontal")
)

Plot_MeanMedian_Cor_Arranged <- egg::ggarrange(as_ggplot(text_grob(label = "all")),
                                               as_ggplot(text_grob(label = "day")),
                                               as_ggplot(text_grob(label = "night")),
                                               Plot_MeanMedian_Cor$all + rremove("x.text"),
                                               Plot_MeanMedian_Cor$day + rremove("xy.text"),
                                               Plot_MeanMedian_Cor$night + rremove("xy.text"),
                                               as_ggplot(text_grob(label = "morning")),
                                               as_ggplot(text_grob(label = "afternoon")),
                                               as_ggplot(text_grob(label = "evening")),
                                               Plot_MeanMedian_Cor$morning,
                                               Plot_MeanMedian_Cor$afternoon + rremove("y.text"),
                                               Plot_MeanMedian_Cor$evening + rremove("y.text"),
                                               ncol = 3,
                                               heights = c(0.07,1,0.07,1))

Plot_MeanMedian_Cor_Arranged <- plot_grid(Plot_MeanMedian_Cor_Arranged, legend_bottom,
                                          ncol = 1, rel_heights = c(1, 0.1))

ggsave(filename = "outputs/figures/IndicesCorrelation_MeanVsMedian.png",
       plot = Plot_MeanMedian_Cor_Arranged,
       width = 30, height = 20, units = "cm", dpi = 800)




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




# Observed vs predicted plots - standard rf models ------------------------


##### SHOULD THESE PREDICTIONS BE USING '$finalModel' ??????

#observed vs predicted plot - birds_morning_richness
obs_pred_birds_morning_richness <- data.frame(predictions = predict(RandomForestFits$birds_morning_richness_totalACI),
                                              observations = acousticIndices_richness$richness[acousticIndices_richness$type == 'birds_morning'])

Plot_obs_pred_Birds_Richness <- ggplot(obs_pred_birds_morning_richness, aes(x = observations, y = predictions)) +
  geom_point() +
  scale_x_continuous(limits = c(min(obs_pred_birds_morning_richness), max(obs_pred_birds_morning_richness))) +
  scale_y_continuous(limits = c(min(obs_pred_birds_morning_richness), max(obs_pred_birds_morning_richness))) +
  geom_abline(slope = 1, linetype = 'dashed') +
  theme_classic()

#observed vs predicted plot - birds_morning_count
obs_pred_birds_morning_count <- data.frame(predictions = predict(RandomForestFits$birds_morning_count_totalACI),
                                           observations = acousticIndices_richness$count[acousticIndices_richness$type == 'birds_morning'])

Plot_obs_pred_Birds_Count <- ggplot(obs_pred_birds_morning_count, aes(x = observations, y = predictions)) +
  geom_point() +
  scale_x_continuous(limits = c(min(obs_pred_birds_morning_count), max(obs_pred_birds_morning_count))) +
  scale_y_continuous(limits = c(min(obs_pred_birds_morning_count), max(obs_pred_birds_morning_count))) +
  geom_abline(slope = 1, linetype = 'dashed') +
  theme_classic()

Plot_obs_pred_Birds <- plot_grid(Plot_obs_pred_Birds_Richness,
                                 Plot_obs_pred_Birds_Count,
                                 ncol = 1,
                                 labels = c("a) richness",
                                            "b) count"),
                                 hjust = 0, label_x = 0.15)

ggsave(filename = "outputs/figures/randomforestobspred/birds.png",
       plot = Plot_obs_pred_Birds,
       width = 10, height = 20, units = "cm", dpi = 800)


#observed vs predicted plot - all_all_richness
obs_pred_all_all_richness <- data.frame(predictions = predict(RandomForestFits$all_all_richness_totalACI),
                                        observations = acousticIndices_richness$richness[acousticIndices_richness$type == 'all_all'])

Plot_obs_pred_All_Richness <- ggplot(obs_pred_all_all_richness, aes(x = observations, y = predictions)) +
  geom_point() +
  scale_x_continuous(limits = c(min(obs_pred_all_all_richness), max(obs_pred_all_all_richness))) +
  scale_y_continuous(limits = c(min(obs_pred_all_all_richness), max(obs_pred_all_all_richness))) +
  geom_abline(slope = 1, linetype = 'dashed') +
  theme_classic()

#observed vs predicted plot - all_all_count
obs_pred_all_all_count <- data.frame(predictions = predict(RandomForestFits$all_all_count_totalACI),
                                     observations = acousticIndices_richness$count[acousticIndices_richness$type == 'all_all'])

Plot_obs_pred_All_Count <- ggplot(obs_pred_all_all_count, aes(x = observations, y = predictions)) +
  geom_point() +
  scale_x_continuous(limits = c(min(obs_pred_all_all_count), max(obs_pred_all_all_count))) +
  scale_y_continuous(limits = c(min(obs_pred_all_all_count), max(obs_pred_all_all_count))) +
  geom_abline(slope = 1, linetype = 'dashed') +
  theme_classic()

Plot_obs_pred_All <- plot_grid(Plot_obs_pred_All_Richness,
                               Plot_obs_pred_All_Count,
                               ncol = 1,
                               labels = c("a) richness",
                                          "b) count"),
                               hjust = 0, label_x = 0.15)

ggsave(filename = "outputs/figures/randomforestobspred/all.png",
       plot = Plot_obs_pred_All,
       width = 10, height = 20, units = "cm", dpi = 800)






# Quantile regression random forest for prediction intervals --------------

library(quantregForest)

control <- trainControl(method = "repeatedcv", number = 10, repeats = 3, verbose = FALSE, savePredictions = TRUE)
tunegrid <- expand.grid(.mtry=c(2:10))

RandomForestQuantRegFits <- list()
RandomForestQuantRegPerformance <- data.frame(comparison = character(),
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

pb = txtProgressBar(min = 0, max = length(unique(acousticIndices_richness$type)) * 3 * 3, initial = 0, style = 3); k <- 0

for (comparison in unique(acousticIndices_richness$type)) {
  tmpdata <- acousticIndices_richness[acousticIndices_richness$type == comparison,]
  
  for (measure in c("richness", "shannon", "count")) {
    
    #forumlas for using just total ACI, 3kHz ACI, and 1kHz ACI values
    formulas <- list(totalACI = as.formula(paste0(measure, " ~ ", paste(grep("ACI_[0-9].", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + "))),
                     ACI_3kHz = as.formula(paste0(measure, " ~ ", paste(grep("ACI_mean$|ACI_1000_2.|ACI_2000_3.|ACI_3000_4.|ACI_4000_5.|ACI_5000_6.|ACI_6000_7.|ACI_7000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + "))),
                     ACI_1kHz = as.formula(paste0(measure, " ~ ", paste(grep("ACI_mean$|ACI_1000_4.|ACI_3000_6.|ACI_5000_8.", grep("*_mean", colnames(acousticIndices_richness), value = TRUE), value = TRUE, invert = TRUE), collapse = " + "))))
    
    for (formula in 1:length(formulas)) {
      
      RandomForestQuantRegFits[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]] <- train(formulas[[formula]],
                                                                                                             data = tmpdata,
                                                                                                             method = "qrf",
                                                                                                             trControl = control,
                                                                                                             importance = TRUE,
                                                                                                             allowParallel = TRUE,
                                                                                                             tuneGrid = tunegrid,
                                                                                                             ntree = 1000,
                                                                                                             keep.inbag = TRUE)
      
      RandomForestQuantRegPerformance <- rbind(RandomForestQuantRegPerformance, data.frame(cbind(comparison = comparison, 
                                                                                                 measure = measure, 
                                                                                                 ACItype = names(formulas)[[formula]],
                                                                                                 RandomForestQuantRegFits[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$results[RandomForestQuantRegFits[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$results$mtry == RandomForestQuantRegFits[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]]$bestTune[[1]],],
                                                                                                 minResponse = min(tmpdata[measure]),
                                                                                                 meanResponse = mean(tmpdata[[measure]]),
                                                                                                 maxResponse = max(tmpdata[measure]))))
      
      RandomForestQuantRegImportance[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]] <- varImp(RandomForestQuantRegFits[[paste0(comparison, "_", measure, "_", names(formulas)[[formula]])]])
      
      k <- k+1; setTxtProgressBar(pb, k)
    }
  }
}

#Normalize RMSE and MAE, calculate Scatter Index
RandomForestQuantRegPerformance <- RandomForestQuantRegPerformance %>% mutate(normRMSE = RMSE/(maxResponse-minResponse),
                                                                              normMAE = MAE/(maxResponse-minResponse),
                                                                              SI = (RMSE/meanResponse)*100,
                                                                              normRMSE_SE = RMSESD/(maxResponse-minResponse)/sqrt(30),
                                                                              normMAE_SE = MAESD/(maxResponse-minResponse)/sqrt(30),
                                                                              SI_SE = (RMSESD/meanResponse)*100/sqrt(30),
                                                                              Rsquared_SE = RsquaredSD/sqrt(30))





test2 <- data.frame(cbind(predict(RandomForestQuantRegFits$birds_morning_count_totalACI$finalModel)),
                          observations = acousticIndices_richness$count[acousticIndices_richness$type == 'birds_morning'])

ggplot(test2, aes(x = observations, y = quantile..0.5)) +
  geom_ribbon(aes(x = observations, ymin = quantile..0.1, ymax = quantile..0.9)) +
  geom_point() +
  scale_x_continuous(limits = c(min(test2), max(test2))) +
  scale_y_continuous(limits = c(min(test2), max(test2))) +
  geom_abline(slope = 1) +
  theme_bw()



#rfinterval

model.rfinterval <- rfinterval::rfinterval(formula = formulas$totalACI,
                                           train_data = acousticIndices_richness[acousticIndices_richness$type == 'birds_morning',])