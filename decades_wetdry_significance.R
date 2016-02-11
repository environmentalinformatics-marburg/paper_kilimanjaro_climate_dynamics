rm(list=ls(all=TRUE))

library(season)
library(quantreg)
#library(boot)

setwd("C:/Users/IOtte/Desktop/") 
prep <- read.csv("metoffice_1973-2013_season.csv", header = TRUE, sep = ";")

#prep$P_SQRT2 <- prep$P_MOSHI_NEW**0.25
prep$mnth <- substr(prep$YEAR, 6, 7)
prep$mnth <- as.numeric(as.character(prep$mnth))
prep$mnth_precent <- prep$mnth / 12
prep$ann <- substr(prep$YEAR, 1, 4)
prep$ann <- as.numeric(as.character(prep$ann))
prep$date_precent <- prep$ann + prep$mnth_precent

#pacf(prep$P_MOSHI_NEW)

#qqnorm(diff(prep$P_MOSHI_NEW**0.25))


### test without sin/cos model
## decade in season and wetdry in wetdry
# dry = 1, wet= 2
# glm
mod_glm <- glm(prep$P_MEAN_NEW ~ prep$mnth_precent + as.factor(prep$season) 
               + as.factor(prep$wd_all))
summary(mod_glm)

#lm
mod_lm <- lm(prep$P_MEAN_NEW ~ prep$mnth_precent + as.factor(prep$season) 
             + as.factor(prep$wd_all))
summary(mod_lm)

## enso/iod und wetdry in phase_wd
# decade and wetdry in season_wd combined
# glm
mod_glm <- glm(prep$P_MEAN_NEW ~ prep$mnth_precent + as.factor(prep$season_wd))
summary(mod_glm)

mod_lm <- lm(prep$P_MEAN_NEW ~ prep$mnth_precent + as.factor(prep$season_wd))
summary(mod_lm)
