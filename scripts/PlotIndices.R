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

# Create plots ------------------------------------------------------------

SamplingTrips <- unique(AP_Kaleidoscope[c("Site", "sampling.period")])

for (trip in 1:nrow(SamplingTrips)) {
  data <- AP_Kaleidoscope[AP_Kaleidoscope$Site == SamplingTrips$Site[trip] & AP_Kaleidoscope$sampling.period == SamplingTrips$sampling.period[trip],]
  
  
  Plot_ACI <- ggplot(data = data,
                     aes(x = TIME_NEW, y = ACI, group = Sensor, colour = Sensor, linetype = Sensor)) +
    geom_smooth() +
    xlab("Time") +
    scale_color_manual(values = c("DryA" = viridis(10)[7], "DryB" = viridis(10)[7], "WetA" = viridis(10)[3], "WetB" = viridis(10)[3])) +
    scale_linetype_manual(values = c("DryA" = "solid", "DryB" = "dashed", "WetA" = "solid", "WetB" = "dashed")) + 
    scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                                "15:00:00", "18:00:00", "21:00:00", "24:00:00")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  Plot_BI <- ggplot(data = data,
                    aes(x = TIME_NEW, y = BI, group = Sensor, colour = Sensor, linetype = Sensor)) +
    geom_smooth() +
    xlab("Time") +
    scale_color_manual(values = c("DryA" = viridis(10)[7], "DryB" = viridis(10)[7], "WetA" = viridis(10)[3], "WetB" = viridis(10)[3])) +
    scale_linetype_manual(values = c("DryA" = "solid", "DryB" = "dashed", "WetA" = "solid", "WetB" = "dashed")) + 
    scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                                "15:00:00", "18:00:00", "21:00:00", "24:00:00")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  Plot_ADI <- ggplot(data = data,
                     aes(x = TIME_NEW, y = ADI, group = Sensor, colour = Sensor, linetype = Sensor)) +
    geom_smooth() +
    xlab("Time") +
    scale_color_manual(values = c("DryA" = viridis(10)[7], "DryB" = viridis(10)[7], "WetA" = viridis(10)[3], "WetB" = viridis(10)[3])) +
    scale_linetype_manual(values = c("DryA" = "solid", "DryB" = "dashed", "WetA" = "solid", "WetB" = "dashed")) + 
    scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                                "15:00:00", "18:00:00", "21:00:00", "24:00:00")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  Plot_AEI <- ggplot(data = data,
                     aes(x = TIME_NEW, y = AEI, group = Sensor, colour = Sensor, linetype = Sensor)) +
    geom_smooth() +
    xlab("Time") +
    scale_color_manual(values = c("DryA" = viridis(10)[7], "DryB" = viridis(10)[7], "WetA" = viridis(10)[3], "WetB" = viridis(10)[3])) +
    scale_linetype_manual(values = c("DryA" = "solid", "DryB" = "dashed", "WetA" = "solid", "WetB" = "dashed")) + 
    scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                                "15:00:00", "18:00:00", "21:00:00", "24:00:00")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  Plot_ACT <- ggplot(data = data,
                     aes(x = TIME_NEW, y = ACT, group = Sensor, colour = Sensor, linetype = Sensor)) +
    geom_smooth() +
    xlab("Time") +
    scale_color_manual(values = c("DryA" = viridis(10)[7], "DryB" = viridis(10)[7], "WetA" = viridis(10)[3], "WetB" = viridis(10)[3])) +
    scale_linetype_manual(values = c("DryA" = "solid", "DryB" = "dashed", "WetA" = "solid", "WetB" = "dashed")) + 
    scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                                "15:00:00", "18:00:00", "21:00:00", "24:00:00")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  Plot_EVN <- ggplot(data = data,
                     aes(x = TIME_NEW, y = EVN, group = Sensor, colour = Sensor, linetype = Sensor)) +
    geom_smooth() +
    xlab("Time") +
    scale_color_manual(values = c("DryA" = viridis(10)[7], "DryB" = viridis(10)[7], "WetA" = viridis(10)[3], "WetB" = viridis(10)[3])) +
    scale_linetype_manual(values = c("DryA" = "solid", "DryB" = "dashed", "WetA" = "solid", "WetB" = "dashed")) + 
    scale_x_discrete(breaks = c("00:00:00", "03:00:00", "06:00:00", "09:00:00", "12:00:00", 
                                "15:00:00", "18:00:00", "21:00:00", "24:00:00")) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  #Arrange plots
  Plots_Combined <- plot_grid(Plot_ACI + theme(legend.position="none"), 
                              Plot_BI + theme(legend.position="none"), 
                              Plot_ADI + theme(legend.position="none"), 
                              Plot_AEI + theme(legend.position="none"), 
                              Plot_ACT + theme(legend.position="none"), 
                              Plot_EVN + theme(legend.position="none"),
                              ncol = 2)
  
  
  legend_b <- get_legend(
    Plot_ACI + 
      guides(color = guide_legend(nrow = 1)) +
      theme(legend.position = "bottom")
  )
  
  # add the legend underneath the row we made earlier. Give it 10%
  # of the height of one plot (via rel_heights).
  Plots_Combined <- plot_grid(Plots_Combined, legend_b, ncol = 1, rel_heights = c(1, .1))
  
  ggsave(filename = paste0("./outputs/figures/timeseries/", Sys.Date(), "_", SamplingTrips$Site[trip], "_", SamplingTrips$sampling.period[trip], ".png"),
         plot = Plots_Combined,
         width = 15, height = 20, units = "cm", dpi = 1000)
}