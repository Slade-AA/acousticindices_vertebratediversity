# Plots for Lin's ARC grant

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(viridis)

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
              colour = "brown", se = FALSE) +
  #54 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Tarcutta' & 
                                       AP_Kaleidoscope$sampling.period == 'spring.2021' &
                                       AP_Kaleidoscope$Sensor == 'WetA',],
              aes(x = TIME_NEW, y = ACI, group = Sensor),
              colour = "blue", se = FALSE) +
  scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                              "15:00:00", "18:00:00", "21:00:00", "24:00:00")) +
  annotate(geom = "text", label = "Bird richness = 54\nMean ACI = 163", x = "19:00:00", y = 174, size = 3) +
  annotate(geom = "text", label = "Bird richness = 26\nMean ACI = 158", x = "19:00:00", y = 149, size = 3) +
  xlab("Time") +
  ylab("Acoustic Complexity (ACI)") +
  coord_cartesian(ylim = c(145, 175)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Plot_Tarcutta_WetA_CLS <- ggplot() +
  #26 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Tarcutta' & 
                                       AP_Kaleidoscope$sampling.period == 'fall.2021' &
                                       AP_Kaleidoscope$Sensor == 'WetA',],
              aes(x = TIME_NEW, y = ClusterCount, group = Sensor),
              colour = "brown", se = FALSE) +
  #54 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Tarcutta' & 
                                       AP_Kaleidoscope$sampling.period == 'spring.2021' &
                                       AP_Kaleidoscope$Sensor == 'WetA',],
              aes(x = TIME_NEW, y = ClusterCount, group = Sensor),
              colour = "blue", se = FALSE) +
  scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                              "15:00:00", "18:00:00", "21:00:00", "24:00:00")) +
  annotate(geom = "text", label = "Bird richness = 54\nMean CLS = 19.60", x = "19:00:00", y = 26, size = 3) +
  annotate(geom = "text", label = "Bird richness = 26\nMean CLS = 13.10", x = "14:30:00", y = 6, size = 3) +
  xlab("Time") +
  ylab("Cluster Count (CLS)") +
  coord_cartesian(ylim = c(0, 30)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Same plot, different season, same richness (Undara DryA fall vs spring)

Plot_Undara_DryA <- ggplot() +
  #23 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Undara' & 
                                       AP_Kaleidoscope$sampling.period == 'fall.2021' &
                                       AP_Kaleidoscope$Sensor == 'DryA',],
              aes(x = TIME_NEW, y = ClusterCount, group = Sensor),
              colour = "blue") +
  #24 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Undara' & 
                                       AP_Kaleidoscope$sampling.period == 'spring.2021' &
                                       AP_Kaleidoscope$Sensor == 'DryA',],
              aes(x = TIME_NEW, y = ClusterCount, group = Sensor),
              colour = "brown") +
  scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                              "15:00:00", "18:00:00", "21:00:00", "24:00:00")) +
  annotate(geom = "text", label = "Bird richness = 23\nMean CLS = 13.63", x = "13:00:00", y = 20) +
  annotate(geom = "text", label = "Bird richness = 24\nMean CLS = 8.86", x = "08:00:00", y = 6) +
  xlab("Time") +
  ylab("Cluster Count (CLS)") +
  coord_cartesian(ylim = c(0, 30)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Sample site & season, different plot, different richness (Wambiana WetB vs DryB in spring)

Plot_Wambiana_WetBDryB <- ggplot() +
  #59 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Wambiana' & 
                                       AP_Kaleidoscope$sampling.period == 'spring.2021' &
                                       AP_Kaleidoscope$Sensor == 'WetB',],
              aes(x = TIME_NEW, y = ClusterCount, group = Sensor),
              colour = "blue") +
  #26 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Wambiana' & 
                                       AP_Kaleidoscope$sampling.period == 'spring.2021' &
                                       AP_Kaleidoscope$Sensor == 'DryB',],
              aes(x = TIME_NEW, y = ClusterCount, group = Sensor),
              colour = "brown") +
  scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                              "15:00:00", "18:00:00", "21:00:00", "24:00:00")) +
  annotate(geom = "text", label = "Bird richness = 59\nMean CLS = 20.84", x = "20:00:00", y = 26) +
  annotate(geom = "text", label = "Bird richness = 26\nMean CLS = 17.00", x = "20:00:00", y = 8) +
  xlab("Time") +
  ylab("Cluster Count (CLS)") +
  coord_cartesian(ylim = c(0, 30)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plot_grid(Plot_Tarcutta_WetA,
          Plot_Undara_DryA,
          Plot_Wambiana_WetBDryB,
          labels = c("a) Tarcutta WetA fall vs spring",
                     "b) Undara DryA fall vs spring",
                     "c) Wambiana WetB vs DryB (spring)"),
          nrow = 1)





Plot_Wambiana_WetB_Undara_DryA_ACI <- ggplot() +
  #59 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Wambiana' & 
                                       AP_Kaleidoscope$sampling.period == 'spring.2021' &
                                       AP_Kaleidoscope$Sensor == 'WetB',],
              aes(x = TIME_NEW, y = ACI, group = Sensor),
              colour = "blue", se = FALSE) +
  #24 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Undara' & 
                                       AP_Kaleidoscope$sampling.period == 'spring.2021' &
                                       AP_Kaleidoscope$Sensor == 'DryA',],
              aes(x = TIME_NEW, y = ACI, group = Sensor),
              colour = "brown", se = FALSE) +
  scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                              "15:00:00", "18:00:00", "21:00:00", "24:00:00")) +
  annotate(geom = "text", label = "Bird richness = 59\nMean ACI = 159", x = "20:00:00", y = 165, size = 3) +
  annotate(geom = "text", label = "Bird richness = 24\nMean ACI = 152", x = "08:30:00", y = 147, size = 3) +
  xlab("Time") +
  ylab("Acoustic Complexity (ACI)") +
  coord_cartesian(ylim = c(145, 175)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

Plot_Wambiana_WetB_Undara_DryA_CLS <- ggplot() +
  #59 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Wambiana' & 
                                       AP_Kaleidoscope$sampling.period == 'spring.2021' &
                                       AP_Kaleidoscope$Sensor == 'WetB',],
              aes(x = TIME_NEW, y = ClusterCount, group = Sensor),
              colour = "blue", se = FALSE) +
  #24 bird species
  geom_smooth(data = AP_Kaleidoscope[AP_Kaleidoscope$Site == 'Undara' & 
                                       AP_Kaleidoscope$sampling.period == 'spring.2021' &
                                       AP_Kaleidoscope$Sensor == 'DryA',],
              aes(x = TIME_NEW, y = ClusterCount, group = Sensor),
              colour = "brown", se = FALSE) +
  scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                              "15:00:00", "18:00:00", "21:00:00", "24:00:00")) +
  annotate(geom = "text", label = "Bird richness = 59\nMean CLS = 20.84", x = "20:00:00", y = 28, size = 3) +
  annotate(geom = "text", label = "Bird richness = 24\nMean CLS = 8.86", x = "08:30:00", y = 4, size = 3) +
  xlab("Time") +
  ylab("Cluster Count (CLS)") +
  coord_cartesian(ylim = c(0, 30)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


CombinedPlot <- plot_grid(Plot_Tarcutta_WetA_ACI,
                          Plot_Wambiana_WetB_Undara_DryA_ACI,
                          Plot_Tarcutta_WetA_CLS,
                          Plot_Wambiana_WetB_Undara_DryA_CLS,
                          labels = c("a) Same site, different season",
                                     "b) Different site, same season",
                                     "c) Same site, different season",
                                     "d) Different site, same season"),
                          hjust = 0, label_x = 0.13,
                          nrow = 2)

CombinedPlot <- egg::ggarrange(as_ggplot(text_grob(label = "Same site, different season")),
                               as_ggplot(text_grob(label = "Different site, same season")),
                               Plot_Tarcutta_WetA_ACI,
                               Plot_Wambiana_WetB_Undara_DryA_ACI,
                               Plot_Tarcutta_WetA_CLS,
                               Plot_Wambiana_WetB_Undara_DryA_CLS,
                               heights = c(0.1, 1, 1))

ggsave(filename = "C:/Users/jc696551/Documents/TimeSeries_ACI&CLS.png",
       plot = CombinedPlot,
       width = 20, height = 20, units = "cm", dpi = 800)





#load rf workspace

Plot_Birds_rf <- egg::ggarrange(as_ggplot(text_grob(label = "richness")),
                                as_ggplot(text_grob(label = "count")),
                                Plots_ObsPred$birds_day_richness,
                                Plots_ObsPred$birds_day_count,
                                ncol = 2, heights = c(0.1, 1)) %>% 
  annotate_figure(bottom = "Predicted",
                  left = "Observed")

ggsave(filename = "C:/Users/jc696551/Documents/Plot_Birds_rf.png",
       plot = Plot_Birds_rf,
       width = 14, height = 7, units = "cm", dpi = 800)