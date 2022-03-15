# Set up ------------------------------------------------------------------

library(tidyverse)
library(tuneR) #used to read in audio
library(seewave) #used to generate spectrogram

indiceDuration <- 60

#Determine audio files to use
baseInputDirectory <- "G:/audio_recorder_downloads_trip2_wavs/Tarcutta_Oct_2021"
audioDirs <- list.dirs(baseInputDirectory, recursive = FALSE)
audioDirs <- audioDirs[grep("^G:/audio_recorder_downloads_trip2_wavs/Tarcutta_Oct_2021/_.*", audioDirs, invert = TRUE)]

audioFiles <- list.files(audioDirs, full.names = TRUE, recursive = TRUE)

audioFiles <- audioFiles[which(as.POSIXlt(paste0(gsub("T.*", "", basename(audioFiles)), 
                                                 gsub(".*T([0-9]{6})\\+.*", "\\1", basename(audioFiles))), 
                                          format = "%Y%m%d%H%M%S") >= "2021-10-18 12:00:00" & 
                                 as.POSIXlt(paste0(gsub("T.*", "", basename(audioFiles)), 
                                                   gsub(".*T([0-9]{6})\\+.*", "\\1", basename(audioFiles))), 
                                            format = "%Y%m%d%H%M%S") < "2021-10-25 12:00:00")]

#Set the frequency limits to use
#aciLimits <- list(default = c(0, 0)) #default is c(0, 0)
#aciLimits <- list(low = c(1000, 4000), mid = c(3000, 6000), high = c(5000, 8000))
aciLimits <- list(f01 = c(1000, 2000),
                  f02 = c(2000, 3000),
                  f03 = c(3000, 4000),
                  f04 = c(4000, 5000),
                  f05 = c(5000, 6000),
                  f06 = c(6000, 7000),
                  f07 = c(7000, 8000))


# Generate indices ----------------------------------------------------------

for (limit in aciLimits) {
  #Create output directory and write file list and settings.ini
  baseOutputDirectory <- "G:/_AcousticIndices_Kaleidoscope/spring.2021"
  acousticIndicesOutputDirectory <- paste0(baseOutputDirectory, "/", Sys.Date(), 
                                           "_ACI", limit[1], "_", limit[2])
  dir.create(acousticIndicesOutputDirectory)
  
  #ACIonly settings file
  settingsFile <- "C:/Users/jc696551/OneDrive - James Cook University/Papers/A2O_BiodiversitySurveys/AcousticIndices/KaleidoscopeSettingsFiles/settings_ACIonly.ini"
  
  write.table(audioFiles,
              file = paste0(acousticIndicesOutputDirectory, "/filelist.txt"),
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
  
  #Create settings file with correct output directory
  data <- readLines(settingsFile)
  
  data[22] <- paste0("directory=", gsub("/", "\\\\\\\\", acousticIndicesOutputDirectory))
  cat(data[22])
  
  data[26] <- paste0("max=", indiceDuration)
  
  #ACI limits
  data[132] <- paste0("fmin=", limit[1])
  data[133] <- paste0("fmax=", limit[2])
  
  writeLines(data, paste0(acousticIndicesOutputDirectory, "/", "customsettings.ini"))
  
  #Run Kaleidoscope CLI
  #Need to get full command to have each file path have single backslashes rather than double backslashes
  #Need to get full command to have each file path be enclosed in quotes
  baseCommand <- '"C:\\Program Files\\Wildlife Acoustics\\Kaleidoscope\\kaleidoscope-cli.exe" --batch "'
  
  fullCommand <- paste0(baseCommand,
                        list.files(path = acousticIndicesOutputDirectory,
                                   pattern = "*customsettings.ini$",
                                   full.names = TRUE),
                        '" --file-list "',
                        list.files(path = acousticIndicesOutputDirectory,
                                   pattern = "*filelist.txt$",
                                   full.names = TRUE),
                        '"')
  
  fullCommand <- gsub("/", "\\\\", fullCommand)
  
  system(command = fullCommand, invisible = FALSE)
}

