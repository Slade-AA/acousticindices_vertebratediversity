# Set up ------------------------------------------------------------------

library(tidyverse)
library(tuneR) #used to read in audio
library(seewave) #used to generate spectrogram

indiceDuration <- 60

#Determine audio files to use
baseInputDirectory <- "G:/audio_recorder_downloads_wavs"
audioDirs <- list.dirs(baseInputDirectory, recursive = FALSE)
audioDirs <- audioDirs[grep("^G:/audio_recorder_downloads_wavs/_.*", audioDirs, invert = TRUE)]

audioFiles <- list.files(audioDirs, full.names = TRUE, recursive = TRUE)

#Create output directory and write file list and settings.ini
baseOutputDirectory <- "G:/_AcousticIndices_Kaleidoscope"
acousticIndicesOutputDirectory <- paste0(baseOutputDirectory, "/", Sys.Date())
dir.create(acousticIndicesOutputDirectory, recursive = TRUE)

settingsFile <- "C:/Users/jc696551/OneDrive - James Cook University/Papers/A2O_BiodiversitySurveys/AcousticIndices/KaleidoscopeSettingsFiles/settings.ini"

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

writeLines(data, paste0(acousticIndicesOutputDirectory, "/", "customsettings.ini"))


# Generate indices ----------------------------------------------------------

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