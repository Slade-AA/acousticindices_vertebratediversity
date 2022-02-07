# Load packages -----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(ggpubr)
library(cowplot)
library(viridis)

# Read in and clean indices -----------------------------------------------

setwd("C:/Users/jc696551/OneDrive - James Cook University/Projects/acousticindices_vertebratediversity")

#all indices - general settings
files_all <- list.files(path = "./rawdata/acousticindices/allindices",
                        pattern = "acousticindex.csv$",
                        full.names = TRUE, 
                        recursive = TRUE)

acousticIndices_all <- lapply(files_all, read_csv)
names(acousticIndices_all) <- c("fall.2021", "spring.2021")
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
                                      spring.2021 = acousticIndices_aci[grep("spring", names(acousticIndices_aci))] %>% reduce(full_join, by = grep("ACI*", colnames(acousticIndices_aci[[1]]), value = TRUE, invert = TRUE))),
                                 .id = "sampling.period")
acousticIndices_aci$sampling.period <- factor(acousticIndices_aci$sampling.period)

#merge all acoustic indices with aci indices
acousticIndices_all$Site <- factor(gsub("G:\\\\.*\\\\(.+?)_.*", "\\1", acousticIndices_all$FOLDER))
acousticIndices_all$Sensor <- factor(gsub(".*[0-9|ING]\\\\(.+?)\\\\.*", "\\1", acousticIndices_all$FOLDER))

acousticIndices_aci$Site <- factor(gsub("G:\\\\.*\\\\(.+?)_.*", "\\1", acousticIndices_aci$FOLDER))
acousticIndices_aci$Sensor <- factor(gsub(".*[0-9|ING]\\\\(.+?)\\\\.*", "\\1", acousticIndices_aci$FOLDER))

acousticIndices <- full_join(acousticIndices_all, acousticIndices_aci, by = c('sampling.period', 'IN FILE', 'CHANNEL', 'OFFSET', 'DURATION', 'DATE', 'TIME', 'HOUR', 'Site', 'Sensor'))



#remove any rows with durations less than 58 seconds? - from short recordings?
acousticIndices <- acousticIndices[!acousticIndices$DURATION < 58,]


#create proper date time column
acousticIndices$DATETIME <- as.POSIXlt(paste0(gsub("*T.*", "", acousticIndices$`IN FILE`), 
                                              gsub(".*T(.+)\\+.*", "\\1", acousticIndices$`IN FILE`)), 
                                       format = "%Y%m%d%H%M%S") + acousticIndices$TIME

#remove times that don't line up with the minute?
acousticIndices <- acousticIndices[!second(acousticIndices$DATETIME) <50 | !second(acousticIndices$DATETIME) > 10,]

#round times to nearest 1 minute interval
acousticIndices$DATETIME <- round_date(acousticIndices$DATETIME, "1 minute")

#extract hour and 1 min interval from datetime
acousticIndices$TIME_NEW <- strftime(acousticIndices$DATETIME, format = "%H:%M:%S")

#Only keep data corresponding to 7-day period of biodiversity surveys (i.e. 12pm first survey day to 12pm last survey day)

acousticIndices <- rbind(acousticIndices %>% filter(Site == "Duval" & DATETIME >= "2021-04-18 12:00:00" & DATETIME < "2021-04-25 12:00:00"), #fall.2021
                         acousticIndices %>% filter(Site == "Mourachan" & DATETIME >= "2021-05-09 12:00:00" & DATETIME < "2021-05-16 12:00:00"),
                         acousticIndices %>% filter(Site == "Rinyirru" & DATETIME >= "2021-06-14 12:00:00" & DATETIME < "2021-06-21 12:00:00"),
                         acousticIndices %>% filter(Site == "Tarcutta" & DATETIME >= "2021-04-29 12:00:00" & DATETIME < "2021-05-06 12:00:00"),
                         acousticIndices %>% filter(Site == "Undara" & DATETIME >= "2021-06-03 12:00:00" & DATETIME < "2021-06-10 12:00:00"),
                         acousticIndices %>% filter(Site == "Wambiana" & DATETIME >= "2021-07-05 12:00:00" & DATETIME < "2021-07-12 12:00:00"),
                         acousticIndices %>% filter(Site == "Rinyirru" & DATETIME >= "2021-10-09 12:00:00" & DATETIME < "2021-10-16 12:00:00"), #spring.2021
                         acousticIndices %>% filter(Site == "Undara" & DATETIME >= "2021-09-29 12:00:00" & DATETIME < "2021-10-06 12:00:00"),
                         acousticIndices %>% filter(Site == "Wambiana" & DATETIME >= "2021-11-09 12:00:00" & DATETIME < "2021-11-16 12:00:00"))

acousticIndices %>% group_by(sampling.period, Site) %>% summarise(n = n()) #each site/sampling.period should have ~40320 rows


# Create plots ------------------------------------------------------------

SamplingTrips <- unique(acousticIndices[c("Site", "sampling.period")])

for (trip in 1:nrow(SamplingTrips)) {
  data <- acousticIndices[acousticIndices$Site == SamplingTrips$Site[trip] & acousticIndices$sampling.period == SamplingTrips$sampling.period[trip],]
  
  
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
  
  ggsave(filename = paste0("./outputs/figures/", Sys.Date(), "_", SamplingTrips$Site[trip], "_", SamplingTrips$sampling.period[trip], ".png"),
         plot = Plots_Combined,
         width = 15, height = 20, units = "cm", dpi = 1200)
}

save(acousticIndices, file = paste0("./outputs/data/", Sys.Date(), "_acousticIndices.RData"))