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



#Join each sampling trips average sunrise and sunset time to indices
files <- file.info(list.files("./outputs/data", pattern = ".*_AverageSunriseSunsetPerTrip.RData$", full.names = TRUE)) #list files
latestFile <- rownames(files)[which.max(files$mtime)] #determine most recent file to use for loading

load(latestFile)

AP_Kaleidoscope_sunrise_sunset <- left_join(AP_Kaleidoscope, AverageSunriseSunsetPerTrip, 
                                            by = c("Site" = "Site", "sampling.period" = "samplingTrip"))

#Calculate summaries using sunrise and sunset times
#Summarise over all 7 days
acousticIndices_7days <- AP_Kaleidoscope %>% 
  group_by(Site, Sensor, sampling.period) %>% 
  mutate(n = n(), p = n()/10080) %>% 
  group_by(Site, Sensor, sampling.period, n, p) %>% 
  summarise_at(vars(all_of(APIndices), all_of(KaleidoscopeIndices),
                    starts_with("ACI")), 
               list(mean = mean, median = median, iqr = IQR, sd = sd))

#Summarise day sunrise-sunset
acousticIndices_day <- AP_Kaleidoscope_sunrise_sunset %>% 
  group_by(Site, sampling.period) %>% 
  filter(hms(TIME_NEW) >= hms(sunrise) & hms(TIME_NEW) < hms(sunset)) %>% 
  group_by(Site, Sensor, sampling.period) %>% 
  mutate(n = n(), p = n()/(period_to_seconds(hm(gsub("\\:[0-9]{2}$","",sunset[1])) - hm(gsub("\\:[0-9]{2}$","",sunrise[1])))/60*7)) %>% 
  group_by(Site, Sensor, sampling.period, n, p) %>% 
  summarise_at(vars(all_of(APIndices), all_of(KaleidoscopeIndices),
                    starts_with("ACI")), 
               list(mean = mean, median = median, iqr = IQR, sd = sd))

#Summarise night sunset-sunrise
acousticIndices_night <- AP_Kaleidoscope_sunrise_sunset %>% 
  group_by(Site, sampling.period) %>% 
  filter(hms(TIME_NEW) >= hms(sunset) | hms(TIME_NEW) < hms(sunrise)) %>% 
  group_by(Site, Sensor, sampling.period) %>% 
  mutate(n = n(), p = n()/(period_to_seconds(hm("24:00") - hm(gsub("\\:[0-9]{2}$","",sunset[1])) + hm(gsub("\\:[0-9]{2}$","",sunrise[1])))/60*7)) %>% 
  group_by(Site, Sensor, sampling.period, n, p) %>% 
  summarise_at(vars(all_of(APIndices), all_of(KaleidoscopeIndices),
                    starts_with("ACI")), 
               list(mean = mean, median = median, iqr = IQR, sd = sd))

#bind in single data frame
acousticIndices_summary <- bind_rows(list(all = acousticIndices_7days, 
                                          day = acousticIndices_day,
                                          night = acousticIndices_night),
                                     .id = "type") %>% ungroup()

save(acousticIndices_summary, file = paste0("./outputs/data/weekly_summaries/", Sys.Date(), "_acousticIndices_summary_sunrise_sunset.RData"))