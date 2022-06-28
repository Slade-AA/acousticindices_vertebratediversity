#Get individual performance values for boxplot
RMSE_values <- list()
for (comparison in 1:length(RandomForestFits_cforest)) {
  RMSE_values[[names(RandomForestFits_cforest[comparison])]] <- data.frame(comparison = gsub("(_count|_richness|_shannon).*", "", names(RandomForestFits_cforest[comparison])),
                                                                           measure = gsub(".*(day_|night_|all_|morning_|afternoon_)([a-z]+)_.*", "\\2", names(RandomForestFits_cforest[comparison])),
                                                                           ACItype = gsub(".*(count_|richness_|shannon_)", "", names(RandomForestFits_cforest[comparison])),
                                                                           RMSE = RandomForestFits_cforest[[comparison]]$resample$RMSE/(max(RandomForestFits_cforest[[comparison]]$trainingData$.outcome)-min(RandomForestFits_cforest[[comparison]]$trainingData$.outcome)))
}
RMSE_values <- do.call(rbind, RMSE_values)

ggplot(data = RMSE_values[RMSE_values$ACItype == 'totalACI' &
                            RMSE_values$comparison %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'),], 
       aes(x = comparison, y = RMSE, fill = measure)) +
  geom_boxplot(position = "dodge", color = "black") +
  #geom_point(position = position_jitterdodge(0.75), color = "black") + #show points as well
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "RMSE") +
  theme_classic() +
  theme(legend.position = "none")



MAE_values <- list()
for (comparison in 1:length(RandomForestFits_cforest)) {
  MAE_values[[names(RandomForestFits_cforest[comparison])]] <- data.frame(comparison = gsub("(_count|_richness|_shannon).*", "", names(RandomForestFits_cforest[comparison])),
                                                                          measure = gsub(".*(day_|night_|all_|morning_|afternoon_)([a-z]+)_.*", "\\2", names(RandomForestFits_cforest[comparison])),
                                                                          ACItype = gsub(".*(count_|richness_|shannon_)", "", names(RandomForestFits_cforest[comparison])),
                                                                          MAE = RandomForestFits_cforest[[comparison]]$resample$MAE,
                                                                          MAE_norm = RandomForestFits_cforest[[comparison]]$resample$MAE/(max(RandomForestFits_cforest[[comparison]]$trainingData$.outcome)-min(RandomForestFits_cforest[[comparison]]$trainingData$.outcome)))
}
MAE_values <- do.call(rbind, MAE_values)

Plot_MAE_richness <- ggplot(data = MAE_values[MAE_values$ACItype == 'totalACI' &
                                                MAE_values$measure == 'richness' &
                                                MAE_values$comparison %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'),], 
                            aes(x = comparison, y = MAE, fill = measure)) +
  geom_boxplot(position = "dodge", color = "black", fill = "#440154FF") +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  labs(x = "Taxa", y = "Mean Absolute Error") +
  theme_classic() +
  theme(legend.position = "none")

Plot_MAE_shannon <- ggplot(data = MAE_values[MAE_values$ACItype == 'totalACI' &
                                                MAE_values$measure == 'shannon' &
                                                MAE_values$comparison %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'),], 
                            aes(x = comparison, y = MAE, fill = measure)) +
  geom_boxplot(position = "dodge", color = "black", fill = "#21908CFF") +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  labs(x = "Taxa", y = "Mean Absolute Error") +
  theme_classic() +
  theme(legend.position = "none")

Plot_MAE_count <- ggplot(data = MAE_values[MAE_values$ACItype == 'totalACI' &
                                                MAE_values$measure == 'count' &
                                                MAE_values$comparison %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'),], 
                            aes(x = comparison, y = MAE, fill = measure)) +
  geom_boxplot(position = "dodge", color = "black", fill = "#FDE725FF") +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  labs(x = "Taxa", y = "Mean Absolute Error") +
  theme_classic() +
  theme(legend.position = "none")

Plot_MAE_norm <- ggplot(data = MAE_values[MAE_values$ACItype == 'totalACI' &
                           MAE_values$comparison %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'),], 
       aes(x = comparison, y = MAE_norm, fill = measure)) +
  geom_boxplot(position = "dodge", color = "black") +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "Normalised Mean Absolute Error") +
  theme_classic() +
  theme(legend.position = "none")


plot_grid(Plot_MAE_richness, Plot_MAE_shannon, Plot_MAE_count, Plot_MAE_norm,
          nrow = 1)

ggsave(filename = "outputs/figures/randomforestperformance/MAE_norm_cforest.png",
       plot = Plot_MAE_norm + theme(legend.position = "bottom",
                                    legend.title = element_blank()),
       width = 10, height = 10, units = "cm", dpi = 800)





SI_values <- list()
for (comparison in 1:length(RandomForestFits_cforest)) {
  SI_values[[names(RandomForestFits_cforest[comparison])]] <- data.frame(comparison = gsub("(_count|_richness|_shannon).*", "", names(RandomForestFits_cforest[comparison])),
                                                                           measure = gsub(".*(day_|night_|all_|morning_|afternoon_)([a-z]+)_.*", "\\2", names(RandomForestFits_cforest[comparison])),
                                                                           ACItype = gsub(".*(count_|richness_|shannon_)", "", names(RandomForestFits_cforest[comparison])),
                                                                           SI = (RandomForestFits_cforest[[comparison]]$resample$RMSE/mean(RandomForestFits_cforest[[comparison]]$trainingData$.outcome))*100)
}
SI_values <- do.call(rbind, SI_values)

ggplot(data = SI_values[SI_values$ACItype == 'totalACI' &
                          SI_values$comparison %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'),], 
       aes(x = comparison, y = SI, fill = measure)) +
  geom_boxplot(position = "dodge", color = "black") +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, 40)) +
  labs(x = "Taxa", y = "Scatter Index") +
  theme_classic() +
  theme(legend.position = "none")


Rsquared_values <- list()
for (comparison in 1:length(RandomForestFits_cforest)) {
  Rsquared_values[[names(RandomForestFits_cforest[comparison])]] <- data.frame(comparison = gsub("(_count|_richness|_shannon).*", "", names(RandomForestFits_cforest[comparison])),
                                                                               measure = gsub(".*(day_|night_|all_|morning_|afternoon_)([a-z]+)_.*", "\\2", names(RandomForestFits_cforest[comparison])),
                                                                               ACItype = gsub(".*(count_|richness_|shannon_)", "", names(RandomForestFits_cforest[comparison])),
                                                                               Rsquared = RandomForestFits_cforest[[comparison]]$resample$Rsquared)
}
Rsquared_values <- do.call(rbind, Rsquared_values)






ggplot(data = RandomForestPerformance_cforest[RandomForestPerformance_cforest$ACItype == 'totalACI' &
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

ggplot(data = RandomForestPerformance_cforest[RandomForestPerformance_cforest$ACItype == 'totalACI' &
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

ggplot(data = Rsquared_values[Rsquared_values$ACItype == 'totalACI' &
                                Rsquared_values$comparison %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'),], 
                  aes(x = comparison, y = Rsquared, fill = measure)) +
  geom_boxplot(position = "dodge", color = "black") +
  #geom_errorbar(aes(ymin = Rsquared-Rsquared_SE, ymax = Rsquared+Rsquared_SE), width = 0.2, size = 1, position = position_dodge(0.9)) +
  scale_fill_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Taxa", y = "R Squared") +
  theme_classic() +
  theme(legend.position = "none")








#Dot and whisker plot version
ggplot(data = RandomForestPerformance_cforest[RandomForestPerformance_cforest$ACItype == 'totalACI' &
                                                RandomForestPerformance_cforest$comparison %in% c('all_all', 'not.birds_all', 'birds_day', 'frogs_night'),], 
       aes(x = comparison, y = SI, fill = measure, colour = measure)) +
  geom_errorbar(aes(ymin = SI-SI_SE, ymax = SI+SI_SE), width = 0.3, size = 1, position = position_dodge(0.9)) +
  geom_point(position = position_dodge(0.9), size = 2, shape = 21) +
  scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  scale_x_discrete(labels = c('All vertebrates', 'Non-avian', 'Birds', 'Frogs')) +
  scale_y_continuous(limits = c(0, 200), breaks = seq(0, 200, 40)) +
  labs(x = "Taxa", y = "Scatter Index") +
  theme_classic() +
  theme(legend.position = "none")