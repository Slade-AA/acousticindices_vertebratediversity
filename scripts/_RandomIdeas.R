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
