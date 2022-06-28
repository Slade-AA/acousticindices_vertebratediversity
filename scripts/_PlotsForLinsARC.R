# Figure 01 - Map of TERN sites & TimeSeries ACI --------------------------

#Map of TERN Supersites

# Load packages -----------------------------------------------------------

library(sf)
library(ozmaps)
library(tidyverse)
library(cowplot)

# Load in site location data and map data ---------------------------------

oz_states <- ozmaps::ozmap_states

tern_supersites <- read.csv("rawdata/geographic/TERNSupersites.csv")

tern_supersites <- st_as_sf(tern_supersites, coords = c("lon", "lat"), crs = "EPSG:4283")

# Plot map ----------------------------------------------------------------

Plot_SupersiteLocations <- ggplot() +
  geom_sf(data = oz_states) +
  geom_sf(data = tern_supersites, colour = "red", size = 1) +
  scale_y_continuous(breaks = c(-45,-35,-25,-15)) +
  scale_x_continuous(breaks = c(110,130,150)) +
  coord_sf(xlim = c(110, 160)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic()


# Acoustic Index Plots for Lin's ARC grant

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(viridis)
library(RColorBrewer)

# Load data ---------------------------------------------------------------

files <- file.info(list.files("./outputs/data/fullsurveyperiod", pattern = ".*_acousticIndices.RData$", full.names = TRUE)) #list files
latestFile <- rownames(files)[which.max(files$mtime)] #determine most recent file to use for loading

load(latestFile)

files <- file.info(list.files("./outputs/data/fullsurveyperiod", pattern = ".*_acousticIndices_AP.RData$", full.names = TRUE)) #list files
latestFile <- rownames(files)[which.max(files$mtime)] #determine most recent file to use for loading

load(latestFile)

#join AP and Kaleidoscope indices
AP_Kaleidoscope <- full_join(acousticIndices_surveys_AP, acousticIndices_surveys,
                             by = c('sampling.period',
                                    'Site',
                                    'Sensor',
                                    'DATETIME',
                                    'TIME_NEW'))
AP_Kaleidoscope <- AP_Kaleidoscope[complete.cases(AP_Kaleidoscope), ]


# Time Groupings ----------------------------------------------------------

AP_Kaleidoscope$hour <- lubridate::hour(AP_Kaleidoscope$DATETIME)
AP_Kaleidoscope$DATETIME_floored <- lubridate::floor_date(AP_Kaleidoscope$DATETIME, unit = "hour")
AP_Kaleidoscope$DATETIME_hour3 <- lubridate::round_date(AP_Kaleidoscope$DATETIME, "3 hours")
AP_Kaleidoscope$hour3 <- lubridate::hour(AP_Kaleidoscope$DATETIME_hour3)


# Plots -------------------------------------------------------------------

# Same plot, different season, different richness (Tarcutta WetA fall vs spring)

Plot_Tarcutta_WetA_ACI <- ggplot() +
  #26 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Tarcutta' & 
                                       AP_Kaleidoscope$sampling.period == 'fall.2021' &
                                       AP_Kaleidoscope$Sensor == 'WetA',],
              aes(x = TIME_NEW, y = ACI, group = Sensor),
              linetype = "dashed",
              colour = brewer.pal(n = 8, name = 'RdBu')[8], se = FALSE) +
  #54 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Tarcutta' & 
                                       AP_Kaleidoscope$sampling.period == 'spring.2021' &
                                       AP_Kaleidoscope$Sensor == 'WetA',],
              aes(x = TIME_NEW, y = ACI, group = Sensor),
              linetype = "solid",
              colour = brewer.pal(n = 8, name = 'RdBu')[1], se = FALSE) +
  scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                              "15:00:00", "18:00:00", "21:00:00", "24:00:00"),
                   labels = c("00:00", "03:00", "06:00", "09:00", "12:00", 
                              "15:00", "18:00", "21:00", "24:00")) +
  annotate(geom = "text", label = "Bird richness = 54", x = "19:00:00", y = 174, size = 3) + #\nMean ACI = 163
  annotate(geom = "text", label = "Bird richness = 26", x = "19:00:00", y = 149, size = 3) + #\nMean ACI = 158
  xlab("Time") +
  ylab("Acoustic Complexity (ACI)") +
  coord_cartesian(ylim = c(145, 175)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Same season, different site, different richness (Wambiana WetB vs Undara DryA in spring)

Plot_Wambiana_WetB_Undara_DryA_ACI <- ggplot() +
  #59 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Wambiana' & 
                                       AP_Kaleidoscope$sampling.period == 'spring.2021' &
                                       AP_Kaleidoscope$Sensor == 'WetB',],
              aes(x = TIME_NEW, y = ACI, group = Sensor),
              linetype = "twodash",
              colour = brewer.pal(n = 9, name = "Oranges")[5], se = FALSE) +
  #24 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Undara' & 
                                       AP_Kaleidoscope$sampling.period == 'spring.2021' &
                                       AP_Kaleidoscope$Sensor == 'DryA',],
              aes(x = TIME_NEW, y = ACI, group = Sensor),
              linetype = "dotted",
              colour = brewer.pal(n = 9, name = "YlGnBu")[5], se = FALSE) +
  scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                              "15:00:00", "18:00:00", "21:00:00", "24:00:00"),
                   labels = c("00:00", "03:00", "06:00", "09:00", "12:00", 
                              "15:00", "18:00", "21:00", "24:00")) +
  annotate(geom = "text", label = "Bird richness = 59", x = "19:00:00", y = 165, size = 3) + #\nMean ACI = 159
  annotate(geom = "text", label = "Bird richness = 24", x = "08:30:00", y = 147, size = 3) + #\nMean ACI = 152
  xlab("Time") +
  ylab("Acoustic Complexity (ACI)") +
  coord_cartesian(ylim = c(145, 175)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Join map and index plot -------------------------------------------------

Figure01_Vertical <- egg::ggarrange(Plot_SupersiteLocations, Plot_Tarcutta_WetA_ACI,
                                    nrow = 2,
                                    labels = c("A", "B"))

ggsave(filename = "outputs/figures.local/LinsARC/Figure01_Vertical.png",
       plot = Figure01_Vertical,
       width = 8, height = 16, units = "cm", dpi = 800)

Figure01_Horizontal <- egg::ggarrange(Plot_SupersiteLocations, Plot_Tarcutta_WetA_ACI,
                                      nrow = 1, widths = c(1.2,1),
                                      labels = c("A", "B"))

ggsave(filename = "outputs/figures.local/LinsARC/Figure01_Horizontal.png",
       plot = Figure01_Horizontal,
       width = 16, height = 8, units = "cm", dpi = 800)



ggsave(filename = "outputs/figures.local/LinsARC/FigureX_TERNSupersiteLocations.png",
       plot = Plot_SupersiteLocations,
       width = 8, height = 8, units = "cm", dpi = 800)
ggsave(filename = "outputs/figures.local/LinsARC/FigureX_SameSiteDifferentSeason.png",
       plot = Plot_Tarcutta_WetA_ACI,
       width = 8, height = 8, units = "cm", dpi = 800)
ggsave(filename = "outputs/figures.local/LinsARC/FigureX_SameSeasonDifferentSite.png",
       plot = Plot_Wambiana_WetB_Undara_DryA_ACI,
       width = 8, height = 8, units = "cm", dpi = 800)


FigureX_ACIHorizontal <- egg::ggarrange(Plot_Tarcutta_WetA_ACI, Plot_Wambiana_WetB_Undara_DryA_ACI,
                                        nrow = 1, widths = c(1,1),
                                        labels = c("A", "B"))
ggsave(filename = "outputs/figures.local/LinsARC/FigureX_ACIHorizontal_diffLineTypes.png",
       plot = FigureX_ACIHorizontal,
       width = 16, height = 8, units = "cm", dpi = 800)

FigureX_ACIHorizontal_Caption <- gridExtra::grid.arrange(FigureX_ACIHorizontal,
                                                         nrow = 1,
                                                         bottom = grid::textGrob(
                                                           "Fig 3. Index of Acoustic Complexity showing real data. A. Same site in 2 seasons, showing the response of animal noise to \nphotoperiod [a physical variable], the longer photoperiod in the spring, (red), is evident by the difference in the sound peaks \nfrom fall (blue). B. At different sites in the same season, note how the site with lower diversity (turquiose) shows the same \npattern but always lower than orange.",
                                                           gp = grid::gpar(fontface = 3, fontsize = 8),
                                                           just = "left", x = 0, check.overlap = TRUE
                                                         ))
ggsave(filename = "outputs/figures.local/LinsARC/FigureX_ACIHorizontal_Caption.png",
       plot = FigureX_ACIHorizontal_Caption,
       width = 16, height = 10, units = "cm", dpi = 800)


gridExtra::grid.arrange(FigureX_ACIHorizontal,
                        nrow = 1,
                        bottom = grid::textGrob(
                          "Fig 3. Index of Acoustic Complexity showing real data. A. Same site in 2 seasons, showing the response of animal noise to \nphotoperiod [a physical variable], the longer photoperiod in the spring, (red), is evident by the difference in the sound peaks \nfrom fall (blue). B. At different sites in the same season, note how the site with lower diversity (turquiose) shows the same \npattern but always lower than orange.",
                          gp = grid::gpar(fontface = 3, fontsize = 8),
                          just = "left", x = 0, check.overlap = TRUE
                        ))


# Figure 02 - Single site observed vs predicted plots ---------------------


# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(lme4)
library(lmerTest)

# Load indices and biodiversity data --------------------------------------

#summary indices
files <- file.info(list.files("./outputs/data/weekly_summaries/", pattern = ".*_acousticIndices_summary.RData$", full.names = TRUE)) #list files
latestFile <- rownames(files)[which.max(files$mtime)] #determine most recent file to use for loading

load(latestFile)

#richness
richness <- read_csv("./rawdata/biodiversity/richness_diversity_updated.csv") #load richness

#format richness to match indices
richness <- richness %>% rename(Site = site, Sensor = original.plot, richness = Richness, shannon = Shannon, type = taxa)
richness <- richness %>% mutate(sampling.period = paste0(season, ".2021"))
richness <- richness %>% mutate(type = recode(type,
                                              'bird' = 'birds', 'frog' = 'frogs'))

# Merge indices and biodiversity data -------------------------------------

#function to combine indices (all, day, night, morning, afternoon, evening) and biodiversity (all, not.birds, birds, frogs) based on user specified combinations
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
                                                                           birds_afternoon = c('afternoon', 'birds'),
                                                                           birds_day = c('day', 'birds'),
                                                                           birds_all = c('all', 'birds'),
                                                                           frogs_evening = c('evening', 'frogs'),
                                                                           frogs_night = c('night', 'frogs')))
acousticIndices_richness <- acousticIndices_richness %>% filter(p > 0.7) #remove datapoints where less than 70% of audio was available

# Fit single site linear models -------------------------------------------

models_singleSite <- list()
obs_pred_values <- list()
for (site in c("Tarcutta", "Wambiana", "Rinyirru")) {
  for (comparison in c("birds_day")) {
    for (measure in c("richness", "count")) {
      models_singleSite[[paste(site, comparison, measure, sep = "_")]] <- lmer(as.formula(paste0(measure, " ~ ", "ClusterCount_mean + SptDensity_mean + MidFreqCover_mean", " + (1|Sensor)")), 
                                                                               data = acousticIndices_richness[acousticIndices_richness$type == comparison &
                                                                                                                 acousticIndices_richness$Site == site,])
      obs_pred_values[[paste(site, comparison, measure, sep = "_")]] <- data.frame(obs = acousticIndices_richness[acousticIndices_richness$type == comparison &
                                                                                                                    acousticIndices_richness$Site == site,
                                                                                                                  c(paste0(measure))],
                                                                                   pred = predict(models_singleSite[[paste(site, comparison, measure, sep = "_")]]))
      
    }
  }
}


# Observed vs predicted plots ---------------------------------------------

Plots_obs_pred <- list()
for (model in 1:length(obs_pred_values)) {
  #plot obs vs pred
  Plots_obs_pred[[paste0(names(obs_pred_values)[model])]] <- ggplot(data = obs_pred_values[[model]],
                                                                    aes(x = pred, y = obs)) +
    geom_abline(slope = 1, linetype = "dashed") +
    geom_smooth(method = "lm", se = FALSE) +
    geom_point() +
    annotate(geom = "text", 
             x = 0.8*(max(obs_pred_values[[model]]) - min(obs_pred_values[[model]])) + min(obs_pred_values[[model]]), 
             y = min(obs_pred_values[[model]]), 
             vjust = 0, size = 3,
             label = paste0("CCC: ", round(DescTools::CCC(obs_pred_values[[model]]$pred, 
                                                          obs_pred_values[[model]]$obs)$rho.c$est, 2))) +
    scale_x_continuous(limits = c(min(obs_pred_values[[model]]),
                                  max(obs_pred_values[[model]]))) +
    scale_y_continuous(limits = c(min(obs_pred_values[[model]]),
                                  max(obs_pred_values[[model]]))) +
    labs(x = "Predicted Richness", y = "Observed Richness") +
    theme_classic()
}

#arrange plots
Figure02_Vertical <- egg::ggarrange(Plots_obs_pred$Tarcutta_birds_day_richness, Plots_obs_pred$Rinyirru_birds_day_richness,
                                    nrow = 2,
                                    labels = c("A", "B"))

ggsave(filename = "outputs/figures.local/LinsARC/Figure02_Vertical.png",
       plot = Figure02_Vertical,
       width = 8, height = 16, units = "cm", dpi = 800)

Figure02_Horizontal <- egg::ggarrange(Plots_obs_pred$Tarcutta_birds_day_richness, Plots_obs_pred$Rinyirru_birds_day_richness,
                                      nrow = 1,
                                      labels = c("A", "B"))

ggsave(filename = "outputs/figures.local/LinsARC/Figure02_Horizontal.png",
       plot = Figure02_Horizontal,
       width = 16, height = 8, units = "cm", dpi = 800)



ggsave(filename = "outputs/figures.local/LinsARC/FigureX_TarcuttaBirdRichness.png",
       plot = Plots_obs_pred$Tarcutta_birds_day_richness,
       width = 8, height = 8, units = "cm", dpi = 800)
ggsave(filename = "outputs/figures.local/LinsARC/FigureX_RinyirruBirdRichness.png",
       plot = Plots_obs_pred$Rinyirru_birds_day_richness,
       width = 8, height = 8, units = "cm", dpi = 800)




# Figure 03 - Manual vs Audio ---------------------------------------------

library(tidyverse)






frogs <- read.csv("C:/Users/jc696551/OneDrive - James Cook University/Papers/A2O_BiodiversitySurveys/Mammal&GeckoRecognition/ZooniverseProject/ClassificationReports/frogs.csv")
mammals <- read.csv("C:/Users/jc696551/OneDrive - James Cook University/Papers/A2O_BiodiversitySurveys/Mammal&GeckoRecognition/ZooniverseProject/ClassificationReports/mammals.csv")

frogs_mammals <- rbind(frogs, mammals)

subjects <- read.csv("C:/Users/jc696551/OneDrive - James Cook University/Papers/A2O_BiodiversitySurveys/Mammal&GeckoRecognition/ZooniverseProject/ClassificationReports/mammal-and-frog-sounds-subjects_20220427.csv")

subjects <- subjects[,1:4]

joined <- left_join(frogs_mammals, subjects, by = c("workflow_id", "subject_id"))

joined <- joined %>% 
  mutate(survey = case_when(
    subject_set_id %in% c(101806, 101807, 101808, 101809) ~ "Tarcutta_April_May_2021",
    subject_set_id %in% c(101810, 101811, 101812, 101813) ~ "Undara",
    subject_set_id %in% c(101814, 101815, 101816, 101817) ~ "Wambiana_Jul_2021",
    subject_set_id %in% c(101818, 101819, 101820, 101821) ~ "Wambiana_Nov_2021"
  ))

joined %>% filter(species != 'OTHERUNSURE') %>% group_by(survey, taxa) %>% summarise(n = n_distinct(species))



frogs_surveys <- read.csv("C:/Users/jc696551/OneDrive - James Cook University/Papers/A2O_BiodiversitySurveys/Mammal&GeckoRecognition/ZooniverseProject/ClassificationReports/frogs_surveys.csv")
mammals_surveys <- read.csv("C:/Users/jc696551/OneDrive - James Cook University/Papers/A2O_BiodiversitySurveys/Mammal&GeckoRecognition/ZooniverseProject/ClassificationReports/mammals_surveys.csv")

manual_surveys <- rbind(frogs_surveys, mammals_surveys)

manual_surveys <- manual_surveys %>% filter(site %in% c('Tarcutta', 'Undara', 'Wambiana'))

manual_surveys <- manual_surveys %>% 
  mutate(date = as.Date(date, format = "%d/%m/%Y")) %>% 
  mutate(survey = case_when(
    site == 'Tarcutta' & date < "2021-06-01" ~ "Tarcutta_April_May_2021",
    site == 'Undara' & date < "2021-07-01" ~ "Undara_Jun_2021", #unsure if this is what has started to be labelled
    site == 'Wambiana' & date < "2021-09-01" ~ "Wambiana_Jul_2021",
    site == 'Wambiana' & date > "2021-09-01" ~ "Wambiana_Nov_2021"
  ))

manual_surveys %>% group_by(survey, class) %>% summarise(n = n_distinct(scientific.name))


library(ggsci)
library(png)
library(patchwork)

bird_image <- readPNG("outputs/figures.local/LinsARC/bird.png", native = TRUE)
frog_image <- readPNG("outputs/figures.local/LinsARC/frog.png", native = TRUE)
mammal_image <- readPNG("outputs/figures.local/LinsARC/fox.png", native = TRUE)

# manual bird counts vs Hayden's cluster labels (Tarcutta) and Naomi's cluster labels (Undara)
data_birds <- data.frame(Method = c("Manual survey", "Audio survey", "Manual survey", "Audio survey"),
                         Site = c("Tarcutta", "Tarcutta", "Undara", "Undara"),
                         SpeciesRichness = c(54, 46, 51, 37))

Plot_Birds <- ggplot(data = data_birds, aes(x = Site, y = SpeciesRichness, group = Method, fill = Method)) +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_jco() +
  scale_y_continuous(limits = c(0, max(data_birds$SpeciesRichness)*1.3)) +
  labs(y = "Species richness") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        legend.position = "none") + inset_element(bird_image, 
                                                  left = 0.05, right = 0.25, 
                                                  bottom = 0.8, top = 1)

# manual frog counts vs Seb's volunteers labelling of night clusters
data_frogs <- data.frame(Method = c("Manual survey", "Audio survey", "Manual survey", "Audio survey", "Manual survey", "Audio survey"),
                         Site = c("Tarcutta_April_May_2021", "Tarcutta_April_May_2021", "Wambiana_Jul_2021", "Wambiana_Jul_2021", "Wambiana_Nov_2021", "Wambiana_Nov_2021"),
                         SpeciesRichness = c(2, 2, 3, 3, 7, 6))
data_frogs <- data.frame(Method = c("Manual survey", "Audio survey", "Manual survey", "Audio survey"),
                         Site = c("Wambiana_Jul_2021", "Wambiana_Jul_2021", "Wambiana_Nov_2021", "Wambiana_Nov_2021"),
                         SpeciesRichness = c(3, 3, 7, 6))

Plot_Frogs <- ggplot(data = data_frogs, aes(x = Site, y = SpeciesRichness, group = Method, fill = Method)) +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_jco() +
  scale_y_continuous(limits = c(0, max(data_frogs$SpeciesRichness)*1.3),
                     breaks = seq(0,8,2)) +
  labs(y = "Species richness") +
  theme_classic() +
  theme(axis.text.x = element_blank(),legend.position = "none") + inset_element(frog_image, 
                                                                                left = 0.05, right = 0.25, 
                                                                                bottom = 0.8, top = 1)

# manual mammal counts vs Seb's volunteers labelling of night clusters
data_mammals <- data.frame(Method = c("Manual survey", "Audio survey", "Manual survey", "Audio survey", "Manual survey", "Audio survey"),
                           Site = c("Tarcutta_April_May_2021", "Tarcutta_April_May_2021", "Wambiana_Jul_2021", "Wambiana_Jul_2021", "Undara", "Undara"),
                           SpeciesRichness = c(8, 4, 6, 7, 3, 2))
data_mammals <- data.frame(Method = c("Manual survey", "Audio survey", "Manual survey", "Audio survey"),
                           Site = c("Wambiana_Jul_2021", "Wambiana_Jul_2021", "Undara", "Undara"),
                           SpeciesRichness = c(6, 7, 3, 2))

Plot_Mammals <- ggplot(data = data_mammals, aes(x = Site, y = SpeciesRichness, group = Method, fill = Method)) +
  geom_col(position = "dodge", colour = "black") +
  scale_fill_jco() +
  scale_y_continuous(limits = c(0, max(data_mammals$SpeciesRichness)*1.3),
                     breaks = seq(0,8,2)) +
  labs(y = "Species richness") +
  theme_classic() +
  theme(axis.text.x = element_blank(),legend.position = c(0.75, 0.9),
        legend.title = element_blank()) + inset_element(mammal_image, 
                                                        left = 0.05, right = 0.25, 
                                                        bottom = 0.8, top = 1)

#Old non-patchwork plots
FigureX_ManualSurvey_Vs_Audio <- egg::ggarrange(Plot_Birds, Plot_Frogs, Plot_Mammals,
                                                nrow = 1, widths = c(1,1,1),
                                                labels = c("A", "B", "C"))
#Patchwork plots
FigureX_ManualSurvey_Vs_Audio <- Plot_Birds + Plot_Frogs + Plot_Mammals + plot_annotation(tag_levels = list(c('A','','B','','C','')))

ggsave(filename = "outputs/figures.local/LinsARC/FigureX_ManualSurvey_Vs_Audio.png",
       plot = FigureX_ManualSurvey_Vs_Audio,
       width = 24, height = 8, units = "cm", dpi = 800)