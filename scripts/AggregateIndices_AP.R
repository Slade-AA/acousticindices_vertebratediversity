#Extract summary indices from QUT outputs

# Load packages -----------------------------------------------------------

library(tidyverse)
library(lubridate)


# Read in indices, add columns for site, sensor etc. ----------------------

summaryIndicesFiles <- list.files(path = "G:/_AcousticIndices_AP",
                                  pattern = ".*Towsey.Acoustic.Indices.csv$",
                                  recursive = TRUE,
                                  full.names = TRUE)

summaryIndices <- list()
for (file in summaryIndicesFiles) {
  tmp <- read.csv(file = file)
  
  #extract site name from folder structure
  tmp$Site <- gsub(".*_AP\\/(.+?)_.*", "\\1", dirname(file))
  
  #extract sensor name from folder structure
  tmp$Sensor <- gsub(".*_[0-9]{4}\\/(.{4})\\/.*", "\\1", dirname(file))
  
  #extract date and time from filename and column in csv
  tmp$DATETIME <- as.POSIXlt(paste0(gsub("T.*", "", basename(file)), 
                                    gsub(".*T(.+)\\+.*", "\\1", basename(file))), 
                             format = "%Y%m%d%H%M%S") + tmp$ResultStartSeconds
  
  #remove times that don't line up with the minute
  tmp <- tmp[!second(tmp$DATETIME) <50 | !second(tmp$DATETIME) > 10,]
  
  #round times to nearest 1 minute interval
  tmp$DATETIME <- round_date(tmp$DATETIME, "1 minute")
  
  #extract hour and 1 min interval from datetime
  tmp$TIME_NEW <- strftime(tmp$DATETIME, format = "%H:%M:%S")
  
  summaryIndices[[paste0(tmp$Site[1], "_", tmp$Sensor[1], "_", gsub("\\+.*", "", basename(file)))]] <- tmp
}

summaryIndices <- do.call(rbind, summaryIndices)


# Subset indices to survey periods ----------------------------------------

acousticIndices_surveys <- rbind(summaryIndices %>% filter(Site == "Duval" & DATETIME >= "2021-04-18 12:00:00" & DATETIME < "2021-04-25 12:00:00"), #fall.2021
                                 summaryIndices %>% filter(Site == "Mourachan" & DATETIME >= "2021-05-09 12:00:00" & DATETIME < "2021-05-16 12:00:00"),
                                 summaryIndices %>% filter(Site == "Rinyirru" & DATETIME >= "2021-06-14 12:00:00" & DATETIME < "2021-06-21 12:00:00"),
                                 summaryIndices %>% filter(Site == "Tarcutta" & DATETIME >= "2021-04-29 12:00:00" & DATETIME < "2021-05-06 12:00:00"),
                                 summaryIndices %>% filter(Site == "Undara" & DATETIME >= "2021-06-03 12:00:00" & DATETIME < "2021-06-10 12:00:00"),
                                 summaryIndices %>% filter(Site == "Wambiana" & DATETIME >= "2021-07-05 12:00:00" & DATETIME < "2021-07-12 12:00:00"))
                                 
acousticIndices_surveys %>% group_by(Site) %>% summarise(n = n()) #each site/sampling.period should have ~40320 rows


# Export data set ---------------------------------------------------------

save(acousticIndices_surveys, file = paste0("./outputs/data/", Sys.Date(), "_acousticIndices_AP.RData"))