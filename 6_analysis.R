# 01/11/2014
# Team PREDICTS-PD
# How does PD respond to human impacts?

# START
cat (paste0 ('\nStage 6 started at [', Sys.time (), ']'))

# LIBS
library (dplyr)
library(lme4)
library (ggplot2)
library (multcomp)
source (file.path ('tools', 'plotting_tools.R'))

# DIRS
input.dir <- '5_metrics'
output.dir <- '6_analysis'
if (!file.exists (output.dir)) {
  dir.create (output.dir)
}

# INPUT
cat ('\nReading in data ....')
# read in RDS
predicts.data <- readRDS (file.path(input.dir, 'predictsdata_wpd.rds'))
# convert to by site
site_ids <- paste0 (predicts.data$SSID, predicts.data$Site_number)
predicts.data <- predicts.data[!duplicated (site_ids), ]

# PROCESS

# how many species made it into the trees?
cat ('Mean [', (1 - mean (predicts.data$PD_pdropped)) * 100,
     '%] species in tree distritbution', sep = '')
sum (predicts.data$Use_intensity == 'Cannot decide')/ nrow (predicts.data)
sum (predicts.data$Use_intensity == 'Minimal use')/ nrow (predicts.data)
sum (predicts.data$Use_intensity == 'Light use')/ nrow (predicts.data)
sum (predicts.data$Use_intensity == 'Intense use')/ nrow (predicts.data)

# CRUDE MODELLING OF PD
plot (predicts.data$Est_mean_PD~predicts.data$Use_intensity)
hist (predicts.data$Est_mean_PD)
m1.pd <- lmer(Est_mean_PD ~ 1 + (1|SSID), data = predicts.data, REML = FALSE)
m2.pd <- lmer(Est_mean_PD ~ Use_intensity + (1|SSID),
           data = predicts.data, REML = FALSE)
m3.pd <- lmer(Est_mean_PD ~ Predominant_habitat + (1|SSID),
              data = predicts.data, REML = FALSE)
m4.pd <- lmer(Est_mean_PD ~ Use_intensity + (1|Predominant_habitat) + (1|SSID),
              data = predicts.data, REML = FALSE)
summary(m1.pd)
summary(m2.pd)
summary(m3.pd)
summary(m4.pd)
qqnorm(resid(m2.pd))
qqline(resid(m2.pd))
plot(resid(m2.pd))
plot(resid(m2.pd) ~ fitted(m2.pd))
abline(h=0)
plot(predicts.data$Est_mean_PD ~ fitted(m2.pd))
abline(a=0,b=1)
dotplot(ranef(m2,pd))
anova(m1.pd, m2.pd, m3.pd, m4.pd)  # m4 is best model?
# plot -- http://stackoverflow.com/questions/9447329/how-to-plot-the-results-of-a-mixed-model
confints <- as.data.frame (confint (glht (m4.pd))$confint)
# make relative to [1, ]
confints[2, ] <- confints[1, ] + confints[2, ]
confints[3, ] <- confints[1, ] + confints[3, ]
confints[4, ] <- confints[1, ] + confints[4, ]
confints$use <- c ('Minimal use', 'Light use', 'Intense use', 'Cannot decide')
ggplot(confints, aes(x = use, y = Estimate, ymin = lwr, ymax = upr)) +
  geom_errorbar() + geom_point() + xlab ('Use Intensity') + ylab ('PD')

# CRUDE MODELLING OF PSV
plot (predicts.data$Est_mean_PSV~predicts.data$Use_intensity)
hist (predicts.data$Est_mean_PSV)
m1.psv <- lmer(Est_mean_PSV ~ 1 + (1|SSID), data = predicts.data, REML = FALSE)
m2.psv <- lmer(Est_mean_PSV~ Use_intensity + (1|SSID),
           data = predicts.data, REML = FALSE)
m3.psv <- lmer(Est_mean_PSV~ Predominant_habitat + (1|SSID),
               data = predicts.data, REML = FALSE)
m4.psv <- lmer(Est_mean_PSV~ Use_intensity + (1|SSID) + (1|Predominant_habitat),
               data = predicts.data, REML = FALSE)
summary(m1.psv)
summary(m2.psv)
anova(m1.psv,m2.psv, m3.psv, m4.psv)  # m4 is best again?
# plot
confints <- as.data.frame (confint (glht (m4.psv))$confint)
# make relative to [1, ]
confints[2, ] <- confints[1, ] + confints[2, ]
confints[3, ] <- confints[1, ] + confints[3, ]
confints[4, ] <- confints[1, ] + confints[4, ]
confints$use <- c ('Minimal use', 'Light use', 'Intense use', 'Cannot decide')
ggplot(confints, aes(x = use, y = Estimate, ymin = lwr, ymax = upr)) +
  geom_errorbar() + geom_point() + xlab ('Use Intensity') + ylab ('PSV')

# CRUDE MODELLING OF PSE
m1.pse <- lmer(as.factor(Est_mean_PSE) ~ 1 + (1|SSID), data = predicts.data, REML = FALSE)
m2.pse <- lmer(as.factor(Est_mean_PSE)~ as.factor(Use_intensity) + (1|SSID),
           data = predicts.data, REML = FALSE)
summary(m1.pse)
summary(m2.pse)
anova(m1.pse,m2.pse)

# OUTPUT

# FINISH
cat (paste0 ('\nStage 6 finished at [', Sys.time (), ']'))