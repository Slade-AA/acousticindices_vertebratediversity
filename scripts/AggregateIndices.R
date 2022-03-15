# Load packages -----------------------------------------------------------

library(tidyverse)
library(lubridate)

# Read in and clean indices -----------------------------------------------

#all indices - general settings
files_all <- list.files(path = "./rawdata/acousticindices/allindices",
                        pattern = "acousticindex.csv$",
                        full.names = TRUE, 
                        recursive = TRUE)

acousticIndices_all <- lapply(files_all, read_csv)
names(acousticIndices_all) <- c("fall.2021", "spring.2021", "spring.2021") #Tarcutta Spring.2021 indices are in a seperate folder - read in as 3rd file
acousticIndices_all <- bind_rows(acousticIndices_all, .id = "sampling.period")
acousticIndices_all$sampling.period <- factor(acousticIndices_all$sampling.period)

#aci index - different frequency bands
files_aci <- list.files(path = "./rawdata/acousticindices/acionly",
                        pattern = "acousticindex.csv$",
                        full.names = TRUE, 
                        recursive = TRUE)
acousticIndices_aci <- lapply(files_aci, read_csv)
names(acousticIndices_aci) <- paste0(gsub(".*acionly\\/(.*)\\/2.*", "\\1", dirname(files_aci)), "_", 
                                     gsub(".*ACI", "", dirname(files_aci)))
acousticIndices_aci <- lapply(seq_along(acousticIndices_aci), function(i) acousticIndices_aci[[i]] %>% rename(!!paste0("ACI_",gsub(".*2021_(.*)$", "\\1", names(acousticIndices_aci)[i])) := ACI))
names(acousticIndices_aci) <- paste0(gsub(".*acionly\\/(.*)\\/2.*", "\\1", dirname(files_aci)), "_", 
                                     gsub(".*ACI", "", dirname(files_aci)))

acousticIndices_aci <- bind_rows(list(fall.2021 = acousticIndices_aci[grep("fall", names(acousticIndices_aci))] %>% reduce(full_join, by = grep("ACI*", colnames(acousticIndices_aci[[1]]), value = TRUE, invert = TRUE)),
                                      spring.2021 = bind_rows(list(acousticIndices_aci[grep("spring", names(acousticIndices_aci))][1:10] %>% reduce(full_join, by = grep("ACI*", colnames(acousticIndices_aci[[1]]), value = TRUE, invert = TRUE)),
                                                                   acousticIndices_aci[grep("spring", names(acousticIndices_aci))][11:20] %>% reduce(full_join, by = grep("ACI*", colnames(acousticIndices_aci[[1]]), value = TRUE, invert = TRUE))))),
                                 .id = "sampling.period")
acousticIndices_aci$sampling.period <- factor(acousticIndices_aci$sampling.period)


# Merge all acoustic indices with aci indices -----------------------------

acousticIndices_all$Site <- factor(gsub("G:\\\\.*\\\\(.+?)_.*", "\\1", acousticIndices_all$FOLDER))
acousticIndices_all$Sensor <- factor(gsub(".*[0-9|ING]\\\\(.+?)\\\\.*", "\\1", acousticIndices_all$FOLDER))

acousticIndices_aci$Site <- factor(gsub("G:\\\\.*\\\\(.+?)_.*", "\\1", acousticIndices_aci$FOLDER))
acousticIndices_aci$Sensor <- factor(gsub(".*[0-9|ING]\\\\(.+?)\\\\.*", "\\1", acousticIndices_aci$FOLDER))

acousticIndices <- full_join(acousticIndices_all, acousticIndices_aci, by = c('sampling.period', 'IN FILE', 'CHANNEL', 'OFFSET', 'DURATION', 'DATE', 'TIME', 'HOUR', 'Site', 'Sensor'))

#NOTE - as of 2022/03/15 frequency specific ACI values haven't been calculated for all Tarcutta recordings available but they are for the 1 week survey period

# Data cleaning and new columns -------------------------------------------

#remove any rows with durations less than 55 seconds? - from short recordings?
acousticIndices <- acousticIndices[!acousticIndices$DURATION < 55,]


#create proper date time column
acousticIndices$DATETIME <- as.POSIXlt(paste0(gsub("*T.*", "", acousticIndices$`IN FILE`), 
                                              gsub(".*T([0-9]{6})\\+.*", "\\1", acousticIndices$`IN FILE`)), 
                                       format = "%Y%m%d%H%M%S") + acousticIndices$TIME

#remove times that don't line up with the minute?
acousticIndices <- acousticIndices[!second(acousticIndices$DATETIME) <50 | !second(acousticIndices$DATETIME) > 10,]

#round times to nearest 1 minute interval
acousticIndices$DATETIME <- round_date(acousticIndices$DATETIME, "1 minute")

#extract hour and 1 min interval from datetime
acousticIndices$TIME_NEW <- strftime(acousticIndices$DATETIME, format = "%H:%M:%S")



# Subset indices to survey periods ----------------------------------------

acousticIndices_surveys <- rbind(acousticIndices %>% filter(Site == "Duval" & DATETIME >= "2021-04-18 12:00:00" & DATETIME < "2021-04-25 12:00:00"), #fall.2021
                                 acousticIndices %>% filter(Site == "Mourachan" & DATETIME >= "2021-05-09 12:00:00" & DATETIME < "2021-05-16 12:00:00"),
                                 acousticIndices %>% filter(Site == "Rinyirru" & DATETIME >= "2021-06-14 12:00:00" & DATETIME < "2021-06-21 12:00:00"),
                                 acousticIndices %>% filter(Site == "Tarcutta" & DATETIME >= "2021-04-29 12:00:00" & DATETIME < "2021-05-06 12:00:00"),
                                 acousticIndices %>% filter(Site == "Undara" & DATETIME >= "2021-06-03 12:00:00" & DATETIME < "2021-06-10 12:00:00"),
                                 acousticIndices %>% filter(Site == "Wambiana" & DATETIME >= "2021-07-05 12:00:00" & DATETIME < "2021-07-12 12:00:00"),
                                 acousticIndices %>% filter(Site == "Rinyirru" & DATETIME >= "2021-10-09 12:00:00" & DATETIME < "2021-10-16 12:00:00"), #spring.2021
                                 acousticIndices %>% filter(Site == "Undara" & DATETIME >= "2021-09-29 12:00:00" & DATETIME < "2021-10-06 12:00:00"),
                                 acousticIndices %>% filter(Site == "Wambiana" & DATETIME >= "2021-11-09 12:00:00" & DATETIME < "2021-11-16 12:00:00"),
                                 acousticIndices %>% filter(Site == "Tarcutta" & DATETIME >= "2021-10-18 12:00:00" & DATETIME < "2021-10-25 12:00:00"))

acousticIndices_surveys %>% group_by(Site, sampling.period) %>% summarise(n = n()) #each site/sampling.period should have ~40320 rows

SiteSensorDataAvailability <- acousticIndices_surveys %>% group_by(Site, Sensor, sampling.period) %>% summarise(n = n())
write.csv(SiteSensorDataAvailability, file = "outputs/SiteSensorDataAvailability.csv", quote = F, row.names = F)

# Export data set ---------------------------------------------------------

save(acousticIndices_surveys, file = paste0("./outputs/data/fullsurveyperiod/", Sys.Date(), "_acousticIndices.RData"))