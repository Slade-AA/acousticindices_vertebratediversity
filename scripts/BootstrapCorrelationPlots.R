
# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(lme4)
library(lmerTest)
library(boot)


# Load indices and richness data, merge -----------------------------------

setwd("C:/Users/jc696551/OneDrive - James Cook University/Projects/acousticindices_vertebratediversity")

load("./outputs/data/2021-12-06_acousticIndices_summary.RData") #load indices

richness <- read_csv("./rawdata/biodiversity/richness_diversity_updated.csv") #load richness

#format richness to match indices
richness <- richness %>% rename(Site = site, Sensor = original.plot, richness = Richness, shannon = Shannon, type = taxa)
#richness$Sensor <- gsub(" ", "", richness$Sensor) #now redundant?
richness <- richness %>% mutate(sampling.period = paste0(season, ".2021"))
richness <- richness %>% mutate(type = recode(type,
                                              'bird' = 'birds', 'frog' = 'frogs'))

#duplicate indices 'all' and rename to 'not.birds' for comparison to all other species
acousticIndices_summary <- rbind(acousticIndices_summary, 
                                 acousticIndices_summary %>% filter(type == 'all') %>% mutate(type = recode(type, 'all' = 'not.birds')))

#merge richness and acoustic indices
acousticIndices_richness <- merge(acousticIndices_summary, richness[richness$day == 'all',], by = c("type", "Site", "Sensor", "sampling.period"))


# Calculate bootstrap correlation between indices and biodiversity --------

bootCor_results <- list()

for (measure in c('count', 'richness', 'shannon')) {
  
  for (taxa in c('all', 'birds', 'frogs', 'not.birds')) {
    
    for (index in colnames(select(acousticIndices_richness, ends_with(c("mean"))))) {
      bootResults <- boot(acousticIndices_richness[acousticIndices_richness$type == taxa,], 
                          statistic = function(data, i) {
                            cor(data[i, measure], data[i, index], method='pearson')
                          },
                          R = 1000)
      bootResultsCI <- boot.ci(bootResults, 
                               conf = 0.95, type = "bca")
      
      bootCor_results[[paste(measure, taxa, index, sep = "_")]] <- data.frame(Index = gsub("_mean", "", index),
                                                                              Taxa = taxa,
                                                                              Measure = measure,
                                                                              Mean = mean(bootResults$t),
                                                                              Low = bootResultsCI$bca[4],
                                                                              High = bootResultsCI$bca[5])
    }
  }
}

bootCor_results <- do.call(rbind, bootCor_results)


# Plot correlation bootstrap results --------------------------------------

indicesToUse <- c("ACI", 
                  "ADI", 
                  "AEI", 
                  "BI", 
                  "NDSI", 
                  "EVN")

#plot facetted by biodiversity measure
correlationPlots <- list()
for (measure in c('richness', 'shannon', 'count')) {
  tmp_data <- bootCor_results[bootCor_results$Measure == measure & bootCor_results$Index %in% indicesToUse,]
  
  correlationPlots[[measure]] <- ggplot(data = tmp_data, aes(x = Mean, y = Index, group = Taxa, colour = Taxa)) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4)) +
    scale_x_continuous(limits = c(-0.9, 0.9), breaks = seq(-0.8, 0.8, 0.4)) +
    scale_color_viridis_d(labels = c('All vertebrates', 'Birds', 'Frogs', 'Not birds')) +
    labs(x = "Mean correlation") +
    theme_classic() +
    theme(axis.title = element_blank(),
          legend.position = "none")
}

legend_bottom <- get_legend(
  correlationPlots[['count']] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank())
)

correlationPlot <- plot_grid(plotlist = correlationPlots,
                             ncol = 1,
                             labels = c(paste0("a) ", names(correlationPlots)[1]),
                                        paste0("b) ", names(correlationPlots)[2]),
                                        paste0("c) ", names(correlationPlots)[3])),
                             hjust = 0, label_x = 0.12) %>% 
  annotate_figure(left = "Acoustic index", bottom = "Mean correlation") %>% 
  plot_grid(legend_bottom, ncol = 1, rel_heights = c(1, .1))

ggsave("outputs/figures/bootstrap_correlations.png",
       correlationPlot,
       width = 12, height = 24, units = "cm", dpi = 800)


#plot facetted by taxa
correlationPlots <- list()
for (taxa in c('all', 'not.birds', 'birds', 'frogs')) {
  tmp_data <- bootCor_results[bootCor_results$Taxa == taxa & bootCor_results$Index %in% indicesToUse,]
  tmp_data$Measure <- factor(tmp_data$Measure, levels = c("richness", "shannon", "count"))
  
  correlationPlots[[taxa]] <- ggplot(data = tmp_data, aes(x = Mean, y = Index, group = Measure, colour = Measure)) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    geom_pointrange(aes(xmin = Low, xmax = High), position = position_dodge(width = 0.4)) +
    scale_x_continuous(limits = c(-0.9, 0.9), breaks = seq(-0.8, 0.8, 0.4)) +
    scale_color_viridis_d(labels = c('Richness', 'Shannon', 'Count')) +
    labs(x = "Mean correlation") +
    theme_classic() +
    theme(axis.title = element_blank(),
          legend.position = "none")
}

legend_bottom <- get_legend(
  correlationPlots[['birds']] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.title = element_blank())
)

correlationPlot <- plot_grid(plotlist = correlationPlots,
                             ncol = 2,
                             labels = c(paste0("a) ", names(correlationPlots)[1]),
                                        paste0("b) ", names(correlationPlots)[2]),
                                        paste0("c) ", names(correlationPlots)[3]),
                                        paste0("d) ", names(correlationPlots)[4])),
                             hjust = 0, label_x = 0.12) %>% 
  annotate_figure(left = "Acoustic index", bottom = "Mean correlation") %>% 
  plot_grid(legend_bottom, ncol = 1, rel_heights = c(1, .1))

ggsave("outputs/figures/bootstrap_correlations_bytaxa.png",
       correlationPlot,
       width = 24, height = 24, units = "cm", dpi = 800)




# Create table of correlation values --------------------------------------

library(kableExtra)

bootCor_results_table <- bootCor_results

bootCor_results_table$value <- paste0(round(bootCor_results_table$Mean, 2), 
                                      " (", 
                                      round(bootCor_results_table$Low, 2), 
                                      " - ", 
                                      round(bootCor_results_table$High, 2), 
                                      ")")

bootCor_results_table <- bootCor_results_table %>% select(-Mean, -Low, -High) %>% pivot_wider(names_from = Index, values_from = value)

target <- c('richness', 'shannon', 'count')

bootCor_results_table <- bootCor_results_table %>% arrange(factor(Measure, levels = target))


bootCor_results_table %>%
  kbl() %>%
  kable_classic_2(full_width = F)