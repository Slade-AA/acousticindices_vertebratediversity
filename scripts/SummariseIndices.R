#Calculate mean and median values for each acoustic index for whole 7 days, morning(6am-9am), and evening(6pm-9pm)

# Load packages -----------------------------------------------------------

library(tidyverse)
library(lubridate)

# Read in pre-cleaned indices ---------------------------------------------

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

# Summarise acoustic indices ----------------------------------------------

APIndices <- c('Activity',
               'EventsPerSecond', 
               'SpectralCentroid', 
               'HighFreqCover', 
               'MidFreqCover', 
               'LowFreqCover', 
               'AcousticComplexity', 
               'TemporalEntropy', 
               'EntropyOfAverageSpectrum',
               'EntropyOfVarianceSpectrum',
               'EntropyOfPeaksSpectrum',
               'EntropyOfCoVSpectrum',
               'ClusterCount',
               'ThreeGramCount',
               'Ndsi',
               'SptDensity')
KaleidoscopeIndices <- c("SH", "NDSI", "ACI", "ADI", "AEI", "BI", "BGN", "SNR", "ACT", "EVN", "LFC", "MFC", "HFC")

#Summarise over all 7 days
acousticIndices_7days <- AP_Kaleidoscope %>% 
  group_by(Site, Sensor, sampling.period) %>% 
  mutate(n = n(), p = n()/10080) %>% 
  group_by(Site, Sensor, sampling.period, n, p) %>% 
  summarise_at(vars(all_of(APIndices), all_of(KaleidoscopeIndices),
                    starts_with("ACI")), 
               list(mean = mean, median = median))

#Summarise morning 6am-9am
acousticIndices_morning <- AP_Kaleidoscope %>% 
  filter(hms(TIME_NEW) >= hms("06:00:00") & hms(TIME_NEW) < hms("09:00:00")) %>% 
  group_by(Site, Sensor, sampling.period) %>% 
  mutate(n = n(), p = n()/1260) %>% 
  group_by(Site, Sensor, sampling.period, n, p) %>% 
  summarise_at(vars(all_of(APIndices), all_of(KaleidoscopeIndices),
                    starts_with("ACI")), 
               list(mean = mean, median = median))

#Summarise afternoon 3pm-6pm
acousticIndices_afternoon <- AP_Kaleidoscope %>% 
  filter(hms(TIME_NEW) >= hms("15:00:00") & hms(TIME_NEW) < hms("18:00:00")) %>% 
  group_by(Site, Sensor, sampling.period) %>% 
  mutate(n = n(), p = n()/1260) %>% 
  group_by(Site, Sensor, sampling.period, n, p) %>% 
  summarise_at(vars(all_of(APIndices), all_of(KaleidoscopeIndices),
                    starts_with("ACI")), 
               list(mean = mean, median = median))

#Summarise evening 6pm-9pm
acousticIndices_evening <- AP_Kaleidoscope %>% 
  filter(hms(TIME_NEW) >= hms("18:00:00") & hms(TIME_NEW) < hms("21:00:00")) %>% 
  group_by(Site, Sensor, sampling.period) %>% 
  mutate(n = n(), p = n()/1260) %>% 
  group_by(Site, Sensor, sampling.period, n, p) %>% 
  summarise_at(vars(all_of(APIndices), all_of(KaleidoscopeIndices),
                    starts_with("ACI")), 
               list(mean = mean, median = median))

#Summarise day 6am-6pm
acousticIndices_day <- AP_Kaleidoscope %>% 
  filter(hms(TIME_NEW) >= hms("06:00:00") & hms(TIME_NEW) < hms("18:00:00")) %>% 
  group_by(Site, Sensor, sampling.period) %>% 
  mutate(n = n(), p = n()/5040) %>% 
  group_by(Site, Sensor, sampling.period, n, p) %>% 
  summarise_at(vars(all_of(APIndices), all_of(KaleidoscopeIndices),
                    starts_with("ACI")), 
               list(mean = mean, median = median))

#Summarise night 6pm-6am
acousticIndices_night <- AP_Kaleidoscope %>% 
  filter(hms(TIME_NEW) >= hms("18:00:00") | hms(TIME_NEW) < hms("06:00:00")) %>% 
  group_by(Site, Sensor, sampling.period) %>% 
  mutate(n = n(), p = n()/5040) %>% 
  group_by(Site, Sensor, sampling.period, n, p) %>% 
  summarise_at(vars(all_of(APIndices), all_of(KaleidoscopeIndices),
                    starts_with("ACI")), 
               list(mean = mean, median = median))

#bind in single data frame
acousticIndices_summary <- bind_rows(list(all = acousticIndices_7days, 
                                          morning = acousticIndices_morning,
                                          afternoon = acousticIndices_afternoon,
                                          evening = acousticIndices_evening,
                                          day = acousticIndices_day,
                                          night = acousticIndices_night),
                                     .id = "type") %>% ungroup()

save(acousticIndices_summary, file = paste0("./outputs/data/weekly_summaries/", Sys.Date(), "_acousticIndices_summary.RData"))