rm(list=ls(all=TRUE))

library(season)
library(quantreg)
#library(boot)

setwd("C:/Users/IOtte/documents/Desktop/kilimanjaro/data/metoffice/data/") 
prep.enid <- read.csv("metoffice_1973-2013_ensoiod_manipulated.csv", header = TRUE, sep = " ")


# data manipulation
prep.enid.all <- prep.enid
prep.enid.all$type <- "meanPre"
prep.enid.all$phase <- 1

# all El Nino - 12 years
prep.enid.all.nino <- subset(prep.enid, 
                             prep.enid$ENSO == "ME" |
                               prep.enid$ENSO == "SE" |
                               prep.enid$ENSO == "WE")
prep.enid.all.nino$type <- "all"
prep.enid.all.nino$phase <- 2
prep.enid.alles <- rbind(prep.enid.all, prep.enid.all.nino)

# pure El Nino - 7 years
prep.enid.pure.enso <- subset(prep.enid, 
                              prep.enid$ENSO == "ME" &
                                prep.enid$IOD != "P" &
                                prep.enid$IOD != "M" |  
                                prep.enid$ENSO == "SE" &
                                prep.enid$IOD != "P" &
                                prep.enid$IOD != "M" |  
                                prep.enid$ENSO == "WE" &
                                prep.enid$IOD != "P" &
                                prep.enid$IOD != "M")
prep.enid.pure.enso$type <- "pureEnso"
prep.enid.pure.enso$phase <- 3
prep.enid.alles <- rbind(prep.enid.alles, prep.enid.pure.enso)

# El Nino IOD+ - 5 years
prep.enid.enso.iod <- subset(prep.enid, 
                             prep.enid$ENSO == "ME" &
                               prep.enid$IOD == "P" |  
                               prep.enid$ENSO == "SE" &
                               prep.enid$IOD == "P" |
                               prep.enid$ENSO == "WE" &
                               prep.enid$IOD == "P")
prep.enid.enso.iod$type <- "EnsoIOD"
prep.enid.enso.iod$phase <- 4
prep.enid.alles <- rbind(prep.enid.alles, prep.enid.enso.iod)

# purest IOD+ (incl. weak El Nino & weak La Nina) - 3 years
prep.enid.pure.iod <- subset(prep.enid, 
                             prep.enid$IOD == "P" &
                               prep.enid$ENSO == "WE"|
                               prep.enid$IOD == "P" &
                               prep.enid$ENSO == "WL")
prep.enid.pure.iod$type <- "purestIOD"
prep.enid.pure.iod$phase <- 5
prep.enid.alles <- rbind(prep.enid.alles, prep.enid.pure.iod)

prep <- prep.enid.alles

prep$P_SQRT2 <- prep$P_MOSHI_NEW**0.25
prep$mnth <- substr(prep$YEAR, 6, 7)
prep$mnth <- as.numeric(as.character(prep$mnth))
prep$mnth_precent <- prep$mnth / 12
prep$ann <- substr(prep$YEAR, 1, 4)
prep$ann <- as.numeric(as.character(prep$ann))
prep$date_precent <- prep$ann + prep$mnth_precent

pacf(prep$P_MOSHI_NEW)

qqnorm(diff(prep$P_MOSHI_NEW**0.25))

################################
##
## Version Ausprobieren
##
################################
sin_month_percent <- sin(2 * pi * prep$mnth_precent)
cos_month_percent <- cos(2 * pi * prep$mnth_precent)

### GLM
# significant, aber nicht normalverteilt
# season as.factor() testet alle Dekaden in einem Modell
mod_glm <- glm(prep$P_MOSHI_NEW ~ sin_month_percent + 
                 cos_month_percent + as.factor(prep$phase))
summary(mod_glm)

# normalverteilt, aber keine significance
mod_glm_sqrt2 <- glm(prep$P_SQRT2 ~ sin_month_percent + 
                       cos_month_percent + as.factor(prep$phase))
summary(mod_glm_sqrt2)


### LM
# significant, aber nicht normalverteilt
# season as.factor() testet alle Dekaden in einem Modell
# gleich wie glm
mod_lm <- lm(prep$P_MOSHI_NEW ~ sin_month_percent + 
               cos_month_percent + as.factor(prep$phase))
summary(mod_lm)

# normalverteilt, aber keine significance
# gleich wie glm
mod_lm_sqrt2 <- lm(prep$P_SQRT2 ~ sin_month_percent + 
                     cos_month_percent + as.factor(prep$phase))
summary(mod_lm_sqrt2)

### Quantile regression
# Test auf Ver??nderung der maximalen (0.95) und minimalen (0.15) Niederschl??ge
reg <- rq(prep$P_MOSHI_NEW ~ prep$date_precent, tau = c(0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99))
summary(reg, se = "nid")

reg15 <- rq(prep$P_MOSHI_NEW ~ prep$date_precent, tau = c(0.15))
summary(reg15, se = "nid")

reg95 <- rq(prep$P_MOSHI_NEW ~ prep$date_precent, tau = c(0.95))
summary(reg95, se = "nid")

# kein (kaum) Trend in den min Niederschlagen, negativer Trend in den max Niederschl??gen (significant)
plot(prep$P_MOSHI_NEW ~ prep$date_precent)
abline(reg95, col = "blue")
abline(reg15, col = "green")

summary(reg15, se = "nid")
summary(reg95, se = "nid")


# Bootstrap (z. B. f|r El Nino Zeugsel)
# Bootstrap (1000 oder 2000 mal zufaellig verteilen)
# jeweils lm und estimate des in Frage kommenden Faktors
# Vergleich des tatsdchlichen estimates mit 95% (einseitig) bzw. 
# 97.5 (zweiseitig) des randominiserten
mod_lm_sqrt2_org <- lm(prep$P_SQRT2 ~ sin_month_percent + 
                         cos_month_percent + as.factor(prep$phase))
summary(mod_lm_sqrt2_org, se = "nid")
mod_lm_sqrt2_org <- as.data.frame(mod_lm_sqrt2_org$coefficients)



mod_bs <- lapply(seq(2000), function(x){
  # Jetzt sampeln, dann rechnen
  set.seed <- x
  boot <- sample(nrow(prep), replace = FALSE) 
  phase <- as.factor(prep$phase[boot])
  act_mod <- lm(prep$P_SQRT2 ~ sin_month_percent + 
                  cos_month_percent + phase)
  return(act_mod$coefficients)
})
mod_bs <- as.data.frame(do.call("rbind", mod_bs))

# all El Nino - 12 years
quantile(mod_bs$phase2, probs = 0.975) 
quantile(mod_bs$phase2, probs = 0.95)
quantile(mod_bs$phase2, probs = 0.06)
quantile(mod_bs$phase2, probs = 0.05)
quantile(mod_bs$phase2, probs = 0.025)
mod_lm_sqrt2_org
sub <- sort(mod_bs$phase2)

# pure El Nino - 7 years
quantile(mod_bs$phase3, probs = 0.975) 
quantile(mod_bs$phase3, probs = 0.95)
quantile(mod_bs$phase3, probs = 0.06)
quantile(mod_bs$phase3, probs = 0.05)
quantile(mod_bs$phase3, probs = 0.025)
mod_lm_sqrt2_org
sort(mod_bs$phase3)

# El Nino IOD+ - 5 years
quantile(mod_bs$phase4, probs = 0.975)
quantile(mod_bs$phase4, probs = 0.95)
quantile(mod_bs$phase4, probs = 0.08)
quantile(mod_bs$phase4, probs = 0.05)
quantile(mod_bs$phase4, probs = 0.025)
mod_lm_sqrt2_org
sort(mod_bs$phase4)

# purest IOD+ (incl. weak El Nino & weak La Nina) - 3 years
quantile(mod_bs$phase5, probs = 0.975)
quantile(mod_bs$phase5, probs = 0.95)
quantile(mod_bs$phase5, probs = 0.06)
quantile(mod_bs$phase5, probs = 0.05)
quantile(mod_bs$phase5, probs = 0.025)
mod_lm_sqrt2_org
sort(mod_bs$phase5)


