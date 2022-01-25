#glmmTMB

library(glmmTMB)

testData <- acousticIndices_richness[acousticIndices_richness$type == 'birds',]

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

testData$ACI_mean <- range01(testData$ACI_mean)

#'censor' min and max values to not be 0 or 1
testData$ACI_mean[which.min(testData$ACI_mean)] <- 0.00001
testData$ACI_mean[which.max(testData$ACI_mean)] <- 0.99999

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