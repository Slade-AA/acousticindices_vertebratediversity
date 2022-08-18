library(vegan)
library(tidyverse)
library(cowplot)

data <- list(all = read.csv("rawdata/biodiversity/sp.all.csv"),
             anuran = read.csv("rawdata/biodiversity/sp.anura.csv"),
             avian = read.csv("rawdata/biodiversity/sp.aves.survey.csv"),
             not.avian = read.csv("rawdata/biodiversity/sp.not.aves.csv"))

#split combined information column into separate columns
all <- data$all %>% mutate(sampling.period = paste0(gsub(" .*", "", season.site.plot), ".2021"),
                      Site = gsub("[A-Za-z]{4,6} ([A-Za-z]{5,9}) .*", "\\1", season.site.plot),
                      Sensor = gsub("[A-Za-z]{4,6} [A-Za-z]{5,9} ([A-Za-z]{4}) .*", "\\1", season.site.plot),
                      type = gsub("[A-Za-z]{4,6} [A-Za-z]{5,9} [A-Za-z]{4} ([A-Za-z]*)", "\\1", season.site.plot))

all %>% select(-season.site.plot, -sampling.period, -Site, -Sensor, -type) %>% rowSums()


#filter out plots that had less than 70% of recordings available so that data matches other analyses

#Summary statistics
speciesPerPlot <- rowSums(all[,2:328])
mean(speciesPerPlot)
sd(speciesPerPlot)

#PCA
# Calculating ordination
ordi <- rda(all[,2:328])

# Ploting ordination diagram:
ordiplot(ordi, type = 'n')
points(ordi, display = 'species', pch = '+', col = 'red')
points(ordi, display = 'sites')
ordihull(ordi, groups = all$Site, draw = 'polygon', alpha = 50)

#NMDS
# Calculating ordination
ordi <- metaMDS(all[,2:328], distance = "bray")

# Ploting ordination diagram:
ordiplot(ordi, type = 'n')
points(ordi, display = 'species', pch = '+', col = 'red')
points(ordi, display = 'sites')
ordihull(ordi, groups = all$Site, draw = 'polygon', alpha = 50)




datasets <- list(all = read.csv("rawdata/biodiversity/sp.all.csv"),
                 not.avian = read.csv("rawdata/biodiversity/sp.not.aves.csv"),
                 avian = read.csv("rawdata/biodiversity/sp.aves.survey.csv"),
                 anuran = read.csv("rawdata/biodiversity/sp.anura.csv"))

#filter out plots that had less than 70% of recordings available so that data matches other analyses

NMDSPlots <- list()
ANOSIM <- list()
SummaryStats <- list()
for (data in 1:length(datasets)) {
  #split combined information column into separate columns
  tmp <- datasets[[data]] %>% mutate(sampling.period = paste0(gsub(" .*", "", season.site.plot), ".2021"),
                          Site = gsub("[A-Za-z]{4,6} ([A-Za-z]{5,9}) .*", "\\1", season.site.plot),
                          Sensor = gsub("[A-Za-z]{4,6} [A-Za-z]{5,9} ([A-Za-z]{4}) .*", "\\1", season.site.plot),
                          type = gsub("[A-Za-z]{4,6} [A-Za-z]{5,9} [A-Za-z]{4} ([A-Za-z]*)", "\\1", season.site.plot))
  
  ordi <- metaMDS(tmp[2:(length(tmp)-4)], distance = "bray")
  
  sites <- data.frame(ordi$points)
  sites$sites <- tmp$Site
  

  
  NMDSPlots[[names(datasets)[data]]] <- ggplot() + 
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
    xlab("NMDS1") +
    ylab("NMDS2") +  
    scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
    scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  
    ggforce::geom_mark_ellipse(data = sites, 
                               aes(x = MDS1, y = MDS2, colour = sites, fill = after_scale(alpha(colour, 0.2))),
                               expand = 0, size = 0.2, show.legend = FALSE) +
    geom_point(data = sites, 
               aes(x = MDS1, y = MDS2, colour = sites, shape = sites), 
               size = 5) +
    #coord_fixed(ratio = 1) +
    theme_classic() +
    theme(legend.position = "none")
  
  #ANOSIM
  ANOSIM[[names(datasets)[data]]] <- anosim(tmp[2:(length(tmp)-4)], tmp$Site, distance = "jaccard")
  
  #Summary statistics
  speciesPerPlot <- rowSums(tmp[2:(length(tmp)-4)])
  SummaryStats[[names(datasets)[data]]] <- data.frame(taxa = names(datasets)[data],
                                                      total = ncol(tmp[2:(length(tmp)-4)]),
                                                      mean = mean(speciesPerPlot),
                                                      sd = sd(speciesPerPlot))
}

legend_bottom <- get_legend(
  NMDSPlots$all + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank()))

BasePlot <- plot_grid(plotlist = NMDSPlots,
                      labels = c("all", "non-avian", "avian", "anuran"), 
                      hjust = 0, label_x = 0.15, label_y = 0.95)

plot_grid(BasePlot, legend_bottom,
          ncol = 1,
          rel_heights = c(1, 0.1))



#NMDS
# Calculating ordination
ordi <- metaMDS(all[2:(length(all)-4)], distance = "bray")

sites <- data.frame(ordi$points)
sites$sites <- all$Site

Plot_NMDS <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab("NMDS1") +
  ylab("NMDS2") +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +  
  ggforce::geom_mark_ellipse(data = sites, 
                             aes(x = MDS1, y = MDS2, colour = sites, fill = after_scale(alpha(colour, 0.2))),
                             expand = 0, size = 0.2, show.legend = FALSE) +
  geom_point(data = sites, 
             aes(x = MDS1, y = MDS2, colour = sites, shape = sites), 
             size = 5) +
  coord_fixed(ratio = 1) +
  theme_bw()