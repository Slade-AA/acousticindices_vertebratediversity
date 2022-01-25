#Correlate acoustic indices with vertebrate diversity - 2021

# Load packages -----------------------------------------------------------

library(tidyverse)
library(ggpubr)
library(cowplot)
library(lme4)
library(lmerTest)
library(boot)
library(glmmTMB)


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


#scale all indices to 0-1 and 'censor' min and max values
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

acousticIndices_richness <- acousticIndices_richness %>% 
  mutate_at(vars(matches("_mean")), range01) %>% 
  mutate_at(vars(matches("_mean")), ~ replace(.x, which(.x == min(.x)), 0.00001)) %>% 
  mutate_at(vars(matches("_mean")), ~ replace(.x, which(.x == max(.x)), 0.99999))



# Mixed effects glmmTMB models --------------------------------------------

Plots_Indices_Richness <- list()
Models_Indices_Richness <- list()
Bootstrap_Indices_Richness <- list()

IndicesRichness_Combinations <- expand.grid(index = grep("mean", colnames(acousticIndices_richness), value = TRUE),
                                            taxa = unique(acousticIndices_richness$type),
                                            diversity = c("richness", "shannon", "count"))

############UP TO HERE##############



glmmTMB_model_richness <- glmmTMB(ACI_mean ~ richness + (1|Site) + (1|sampling.period),
                                  data = testData,
                                  family = beta_family())

bootstrapPredictions <- bootMer(glmmTMB_model_richness,
                                FUN = function(x)predict(x, re.form=NA, type = "response"),
                                nsim = 100)

testData$lci <- apply(bootstrapPredictions$t, 2, quantile, 0.025)
testData$uci <- apply(bootstrapPredictions$t, 2, quantile, 0.975)
testData$pred <- predict(glmmTMB_model_richness, re.form=NA, type = "response")

ggplot(data = testData, aes(x = richness, y = ACI_mean)) +
  geom_ribbon(aes_string(x = "richness", ymin = "lci", ymax = "uci"), fill = "black", alpha = 0.1) +
  geom_line(aes_string(x = "richness", y = "pred"), color = "black", lwd = 1) +
  geom_point() +
  theme_classic()


lme4::bootMer(glmmTMB_model_richness,nsim=10,FUN=function(x) unlist(fixef(x)))


glmmTMB_model_shannon <- glmmTMB(ACI_mean ~ shannon,
                                  data = testData,
                                  family = beta_family())

AICcmodavg::aictab(cand.set = list(richness = glmmTMB_model_richness,
                                   shannon = glmmTMB_model_shannon))