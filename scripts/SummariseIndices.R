#Calculate mean and median values for each acoustic index for whole 7 days, morning(6am-9am), and evening(6pm-9pm)

# Load packages -----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(ggpubr)
library(cowplot)


# Read in pre-cleaned indices ---------------------------------------------

setwd("C:/Users/jc696551/OneDrive - James Cook University/Projects/acousticindices_vertebratediversity")

load("./outputs/data/2021-12-06_acousticIndices.RData")


# Summarise acoustic indices ----------------------------------------------

#Summarise over all 7 days
acousticIndices_7days <- acousticIndices %>% 
  group_by(Site, Sensor, sampling.period) %>% 
  summarise_at(c("SH", "NDSI", "ACI", "ADI", "AEI", "BI", "BGN", "SNR", "ACT", "EVN", "LFC", "MFC", "HFC"), list(mean = mean, median = median))

#Summarise morning 6am-9am
acousticIndices_morning <- acousticIndices %>% 
  filter(hms(TIME_NEW) >= hms("06:00:00") & hms(TIME_NEW) < hms("09:00:00")) %>% 
  group_by(Site, Sensor, sampling.period) %>% 
  summarise_at(c("SH", "NDSI", "ACI", "ADI", "AEI", "BI", "BGN", "SNR", "ACT", "EVN", "LFC", "MFC", "HFC"), list(mean = mean, median = median))

#Summarise evening 6pm-9pm
acousticIndices_evening <- acousticIndices %>% 
  filter(hms(TIME_NEW) >= hms("18:00:00") & hms(TIME_NEW) < hms("21:00:00")) %>% 
  group_by(Site, Sensor, sampling.period) %>% 
  summarise_at(c("SH", "NDSI", "ACI", "ADI", "AEI", "BI", "BGN", "SNR", "ACT", "EVN", "LFC", "MFC", "HFC"), list(mean = mean, median = median))


acousticIndices_summary <- bind_rows(list(all = acousticIndices_7days, 
                                          birds = acousticIndices_morning, 
                                          frogs = acousticIndices_evening),
                                     .id = "type")

save(acousticIndices_summary, file = paste0("./outputs/data/", Sys.Date(), "_acousticIndices_summary.RData"))