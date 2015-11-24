rm(list=ls(all=TRUE))

library(season)
library(quantreg)
#library(boot)

setwd("C:/Users/IOtte/documents/Desktop/kilimanjaro/data/cloud_cover_monthly_mod06/cluster/") 
cld.enid <- read.csv("clust_mean_series_terra_ensoiod.csv", header = TRUE, sep = ",")


# data manipulation
cld.enid.all <- cld.enid
cld.enid.all$type <- "meanPre"
cld.enid.all$phase <- 1

# all La Nina - 5 years
cld.enid.all.nina <- subset(cld.enid, 
                            cld.enid$ENSO == "ML" |
                              cld.enid$ENSO == "SL" |
                              cld.enid$ENSO == "WL")
cld.enid.all.nina$type <- "all"
cld.enid.all.nina$phase <- 2
cld.enid.alles <- rbind(cld.enid.all, cld.enid.all.nina)

# pure El Nino - 4 years
cld.enid.pure.enso <- subset(cld.enid, 
                              cld.enid$ENSO == "ML" &
                                cld.enid$IOD != "P" &
                                cld.enid$IOD != "M" |  
                                cld.enid$ENSO == "SL" &
                                cld.enid$IOD != "P" &
                                cld.enid$IOD != "M" |  
                                cld.enid$ENSO == "WL" &
                                cld.enid$IOD != "P" &
                                cld.enid$IOD != "M")
cld.enid.pure.enso$type <- "pureEnso"
cld.enid.pure.enso$phase <- 3
cld.enid.alles <- rbind(cld.enid.alles, cld.enid.pure.enso)

# El Nino IOD+ - 1 year
cld.enid.enso.iod <- subset(cld.enid, 
                             cld.enid$ENSO == "ML" &
                               cld.enid$IOD == "P" |  
                               cld.enid$ENSO == "SL" &
                               cld.enid$IOD == "P" |
                               cld.enid$ENSO == "WL" &
                               cld.enid$IOD == "P")
cld.enid.enso.iod$type <- "EnsoIOD"
cld.enid.enso.iod$phase <- 4
cld.enid.alles <- rbind(cld.enid.alles, cld.enid.enso.iod)

# purest IOD+ (incl. weak El Nino & weak La Nina) - 1 year
#cld.enid.pure.iod <- subset(cld.enid, 
#                             cld.enid$IOD == "P" &
#                               cld.enid$ENSO == "WE"|
#                               cld.enid$IOD == "P" &
#                               cld.enid$ENSO == "WL")
#cld.enid.pure.iod$type <- "purestIOD"
#cld.enid.pure.iod$phase <- 5
#cld.enid.alles <- rbind(cld.enid.alles, cld.enid.pure.iod)

cld <- cld.enid.alles

cld$P_SQRT2 <- cld$clust1**0.25
cld$mnth <- substr(cld$YEAR, 6, 7)
cld$mnth <- as.numeric(as.character(cld$mnth))
cld$mnth_precent <- cld$mnth / 12
cld$ann <- substr(cld$YEAR, 1, 4)
cld$ann <- as.numeric(as.character(cld$ann))
cld$date_precent <- cld$ann + cld$mnth_precent

pacf(cld$clust1)

qqnorm(diff(cld$clust1**0.25))

################################
##
## Version Ausprobieren
##
################################
sin_month_percent <- sin(2 * pi * cld$mnth_precent)
cos_month_percent <- cos(2 * pi * cld$mnth_precent)

### GLM
# significant, aber nicht normalverteilt
# season as.factor() testet alle Dekaden in einem Modell
mod_glm <- glm(cld$clust1 ~ sin_month_percent + 
                 cos_month_percent + as.factor(cld$phase))
summary(mod_glm)

# normalverteilt, aber keine significance
mod_glm_sqrt2 <- glm(cld$P_SQRT2 ~ sin_month_percent + 
                       cos_month_percent + as.factor(cld$phase))
summary(mod_glm_sqrt2)


### LM
# significant, aber nicht normalverteilt
# season as.factor() testet alle Dekaden in einem Modell
# gleich wie glm
mod_lm <- lm(cld$clust1 ~ sin_month_percent + 
               cos_month_percent + as.factor(cld$phase))
summary(mod_lm)

# normalverteilt, aber keine significance
# gleich wie glm
mod_lm_sqrt2 <- lm(cld$P_SQRT2 ~ sin_month_percent + 
                     cos_month_percent + as.factor(cld$phase))
summary(mod_lm_sqrt2)

### Quantile regression
# Test auf Ver??nderung der maximalen (0.95) und minimalen (0.15) Niederschl??ge
reg <- rq(cld$clust1 ~ cld$date_precent, tau = c(0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99))
summary(reg, se = "nid")

reg15 <- rq(cld$clust1 ~ cld$date_precent, tau = c(0.15))
summary(reg15, se = "nid")

reg95 <- rq(cld$clust1 ~ cld$date_precent, tau = c(0.95))
summary(reg95, se = "nid")

# kein (kaum) Trend in den min Niederschlagen, negativer Trend in den max Niederschl??gen (significant)
plot(cld$clust1 ~ cld$date_precent)
abline(reg95, col = "blue")
abline(reg15, col = "green")

summary(reg15, se = "nid")
summary(reg95, se = "nid")


# Bootstrap (z. B. f|r El Nino Zeugsel)
# Bootstrap (1000 oder 2000 mal zufaellig verteilen)
# jeweils lm und estimate des in Frage kommenden Faktors
# Vergleich des tatsdchlichen estimates mit 95% (einseitig) bzw. 
# 97.5 (zweiseitig) des randominiserten
mod_lm_sqrt2_org <- lm(cld$P_SQRT2 ~ sin_month_percent + 
                         cos_month_percent + as.factor(cld$phase))
summary(mod_lm_sqrt2_org, se = "nid")
mod_lm_sqrt2_org <- as.data.frame(mod_lm_sqrt2_org$coefficients)



mod_bs <- lapply(seq(2000), function(x){
  # Jetzt sampeln, dann rechnen
  set.seed <- x
  boot <- sample(nrow(cld), replace = FALSE) 
  phase <- as.factor(cld$phase[boot])
  act_mod <- lm(cld$P_SQRT2 ~ sin_month_percent + 
                  cos_month_percent + phase)
  return(act_mod$coefficients)
})
mod_bs <- as.data.frame(do.call("rbind", mod_bs))

# all La Nina - 5 years
quantile(mod_bs$phase2, probs = 0.975) 
quantile(mod_bs$phase2, probs = 0.95)
quantile(mod_bs$phase2, probs = 0.05)
quantile(mod_bs$phase2, probs = 0.025)
mod_lm_sqrt2_org
sort(mod_bs$phase2)

# pure La Nina - 4 years
quantile(mod_bs$phase3, probs = 0.975) 
quantile(mod_bs$phase3, probs = 0.95)
quantile(mod_bs$phase3, probs = 0.05)
quantile(mod_bs$phase3, probs = 0.025)
mod_lm_sqrt2_org
sort(mod_bs$phase3)

# La Nina IOD+ - 1 year
quantile(mod_bs$phase4, probs = 0.975)
quantile(mod_bs$phase4, probs = 0.95)
quantile(mod_bs$phase4, probs = 0.05)
quantile(mod_bs$phase4, probs = 0.025)
mod_lm_sqrt2_org
sort(mod_bs$phase4)

