rm(list=ls(all=TRUE))

library(season)
library(quantreg)
#library(boot)

setwd("C:/Users/IOtte/Desktop/") 
prep <- read.csv("rainfall_enso_iod.csv", header = TRUE, sep = ";")

prep$P_SQRT2 <- prep$P_MOSHI_NEW**0.25
prep$mnth <- substr(prep$YEAR, 6, 7)
prep$mnth <- as.numeric(as.character(prep$mnth))
prep$mnth_precent <- prep$mnth / 12
prep$ann <- substr(prep$YEAR, 1, 4)
prep$ann <- as.numeric(as.character(prep$ann))
prep$date_precent <- prep$ann + prep$mnth_precent

pacf(prep$P_MOSHI_NEW)

qqnorm(diff(prep$P_MOSHI_NEW**0.25))


### test without sin/cos model
## enso/iod in phase wetdry in wetdry
# 1 = dry (Jan/Feb), 2 = long rains, 3 = dry(Jun/Jul/Aug/Sep), 4 = wet(short rains)
# glm
mod_glm <- glm(prep$P_MEAN_NEW ~ prep$mnth_precent + as.factor(prep$phase) 
               + as.factor(prep$wd_numeric))
summary(mod_glm)

#lm
mod_lm <- lm(prep$P_MEAN_NEW ~ prep$mnth_precent + as.factor(prep$phase) 
             + as.factor(prep$wd_numeric))
summary(mod_lm)

## enso/iod und wetdry in phase_wd
# glm
mod_glm <- glm(prep$P_MEAN_NEW ~ prep$mnth_precent + as.factor(prep$phase_wd))
summary(mod_glm)

mod_lm <- lm(prep$P_MEAN_NEW ~ prep$mnth_precent + as.factor(prep$phase_wd))
summary(mod_lm)
