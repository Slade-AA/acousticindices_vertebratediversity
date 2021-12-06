#Calculate mean and median values for each indice for whole 7 days, 6am-9am, and 6pm-9pm

library(tidyverse)
library(lubridate)
library(ggpubr)
library(cowplot)

#Read in pre-cleaned indices
setwd("C:/Users/jc696551/OneDrive - James Cook University/Papers/A2O_BiodiversitySurveys/AcousticIndices")

acousticIndices_Trip1 <- read_csv("./_AcousticIndices_Kaleidoscope/2021-11-15/acousticIndices_Trip1.csv")

#Summarise over all 7 days
IndicesSummary_7days <- acousticIndices_Trip1 %>% 
  group_by(Site, Sensor) %>% 
  summarise_at(c("SH", "NDSI", "ACI", "ADI", "AEI", "BI", "ACT", "EVN"), list(mean = mean, median = median))

#Summarise morning 6am-9am
IndicesSummary_Morning <- acousticIndices_Trip1 %>% 
  filter(TIME_NEW >= hms("06:00:00") & TIME_NEW < hms("09:00:00")) %>% 
  group_by(Site, Sensor) %>% 
  summarise_at(c("SH", "NDSI", "ACI", "ADI", "AEI", "BI", "ACT", "EVN"), list(mean = mean, median = median))

#Summarise evening 6pm-9pm
IndicesSummary_Evening <- acousticIndices_Trip1 %>% 
  filter(TIME_NEW >= hms("18:00:00") & TIME_NEW < hms("21:00:00")) %>% 
  group_by(Site, Sensor) %>% 
  summarise_at(c("SH", "NDSI", "ACI", "ADI", "AEI", "BI", "ACT", "EVN"), list(mean = mean, median = median))