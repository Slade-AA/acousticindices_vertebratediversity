library(tidyverse)
library(lubridate)
library(suncalc)

audio_sensors <- read.csv("rawdata/geographic/audio_sensors.csv")

#Extract mean Site locations from individual sensor locations and make column for timezone
audio_sensors %>% mutate(Site = str_to_title(gsub("\\..*", "", name))) %>% 
  group_by(Site) %>% summarise(lat = mean(lat),
                               lon = mean(lon)) %>% 
  mutate(tz = case_when(
    Site %in% c('Tarcutta', 'Duval') ~ 'Australia/Sydney',
    Site %in% c('Mourachan', 'Wambiana', 'Undara', 'Rinyirru') ~ 'Australia/Brisbane'
  )) -> SiteLocations

#Data frame of sampling trips at each site and dates sampled
SamplingTrips <- data.frame(Site = c(rep('Duval', 8),
                                     rep('Mourachan', 8),
                                     rep('Rinyirru', 8),
                                     rep('Tarcutta', 8),
                                     rep('Undara', 8),
                                     rep('Wambiana', 8),
                                     rep('Rinyirru', 8),
                                     rep('Undara', 8),
                                     rep('Wambiana', 8),
                                     rep('Tarcutta', 8)),
                            date = c(seq(ymd('2021-04-18'), ymd('2021-04-25'), by='days'),
                                     seq(ymd('2021-05-09'), ymd('2021-05-16'), by='days'),
                                     seq(ymd('2021-06-14'), ymd('2021-06-21'), by='days'),
                                     seq(ymd('2021-04-29'), ymd('2021-05-06'), by='days'),
                                     seq(ymd('2021-06-03'), ymd('2021-06-10'), by='days'),
                                     seq(ymd('2021-07-05'), ymd('2021-07-12'), by='days'),
                                     seq(ymd('2021-10-09'), ymd('2021-10-16'), by='days'),
                                     seq(ymd('2021-09-29'), ymd('2021-10-06'), by='days'),
                                     seq(ymd('2021-11-09'), ymd('2021-11-16'), by='days'),
                                     seq(ymd('2021-10-18'), ymd('2021-10-25'), by='days')),
                            samplingTrip = c(rep('fall.2021', 8),
                                             rep('fall.2021', 8),
                                             rep('fall.2021', 8),
                                             rep('fall.2021', 8),
                                             rep('fall.2021', 8),
                                             rep('fall.2021', 8),
                                             rep('spring.2021', 8),
                                             rep('spring.2021', 8),
                                             rep('spring.2021', 8),
                                             rep('spring.2021', 8)))

#Merge sensor locations and sampling dates
merged <- left_join(SamplingTrips, SiteLocations) %>% mutate(date = as.Date(date))

#Calculate sunrise and senset times based on lat, lon, date
merged <- left_join(merged, rbind(getSunlightTimes(data = merged[merged$tz == 'Australia/Sydney',],
                                                   keep = c('sunrise', 'sunset'),
                                                   tz = "Australia/Sydney") %>% mutate(sunrise = as.character(sunrise),
                                                                                       sunset = as.character(sunset)),
                                  getSunlightTimes(data = merged[merged$tz == 'Australia/Brisbane',],
                                                   keep = c('sunrise', 'sunset'),
                                                   tz = "Australia/Brisbane") %>% mutate(sunrise = as.character(sunrise),
                                                                                         sunset = as.character(sunset))))


#Calculate average sunrise and sunset time for each sampling trip
merged %>% group_by(Site, samplingTrip) %>% summarise(sunrise = format(mean(strptime(gsub(".* ([0-9:])", "\\1", sunrise), "%H:%M:%S")), "%H:%M:%S"),
                                                      sunset = format(mean(strptime(gsub(".* ([0-9:])", "\\1", sunset), "%H:%M:%S")), "%H:%M:%S"))
