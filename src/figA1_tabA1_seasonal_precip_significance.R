rm(list=ls(all=TRUE))

library(season)
library(quantreg)
library(ggplot2)
#library()

setwd("C:/Users/IOtte/documents/Desktop/kilimanjaro/data/metoffice/data/")
prep <- read.csv2("precip.ct.all.csv", header = TRUE)

prep$P_SQRT2 <- prep$P_MEAN_NEW**0.25
prep$mnth <- as.numeric(as.character(prep$mnth))
prep$mnth_precent <- prep$mnth / 12
prep$ann <- as.numeric(as.character(prep$ann))
prep$date_precent <- prep$ann + prep$mnth_precent

pacf(prep$P_MEAN_NEW)

qqnorm(diff(prep$P_MEAN_NEW**0.25))

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
mod_glm <- glm(prep$P_MEAN_NEW ~ sin_month_percent + 
                 cos_month_percent + as.factor(prep$season))
summary(mod_glm)

# normalverteilt, aber keine significance
mod_glm_sqrt2 <- glm(prep$P_SQRT2 ~ sin_month_percent + 
                       cos_month_percent + as.factor(prep$season))
summary(mod_glm_sqrt2)


### LM
# significant, aber nicht normalverteilt
# season as.factor() testet alle Dekaden in einem Modell
# gleich wie glm
mod_lm <- lm(prep$P_MEAN_NEW ~ sin_month_percent + 
               cos_month_percent + as.factor(prep$season))
summary(mod_lm)

# normalverteilt, aber keine significance
# gleich wie glm
mod_lm_sqrt2 <- lm(prep$P_SQRT2 ~ sin_month_percent + 
                     cos_month_percent + as.factor(prep$season))
summary(mod_lm_sqrt2)

### Quantile regression
# Test auf Ver??nderung der maximalen (0.95) und minimalen (0.15) Niederschl??ge
reg <- rq(prep$P_MEAN_NEW ~ prep$date_precent, tau = c(0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99))
summary(reg, se = "nid")

reg15 <- rq(prep$P_MEAN_NEW ~ prep$date_precent, tau = c(0.15))
summary(reg15, se = "nid")

reg95 <- rq(prep$P_MEAN_NEW ~ prep$date_precent, tau = c(0.95))
summary(reg95, se = "nid")

# kein (kaum) Trend in den min Niederschlagen, negativer Trend in den max Niederschl??gen (significant)
plot(prep$P_MEAN_NEW ~ prep$date_precent)
abline(reg95, col = "blue")
abline(reg15, col = "green")

prep.reg <- ggplot(prep, aes(x = date_precent, y = P_MEAN_NEW)) +
         geom_point() + 
  geom_abline(intercept = 24.04765, slope = -0.01059, color = "green") +
  geom_abline(intercept = 7828.56524, slope = -3.81220, color = "blue") +
  xlab("") +
  ylab("Precipitation") +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    legend.background = element_rect(color = "black", fill = "white")) +
  labs(fill = "")

png("prep.reg.png", width = 30, height = 10, 
    units = "cm", res = 300, pointsize = 15)
print(prep.reg)
dev.off()

summary(reg15, se = "nid")
summary(reg95, se = "nid")


# Bootstrap (z. B. f|r El Nino Zeugsel)
# Bootstrap (1000 oder 2000 mal zufaellig verteilen)
# jeweils lm und estimate des in Frage kommenden Faktors
# Vergleich des tatsdchlichen estimates mit 95% (einseitig) bzw. 
# 97.5 (zweiseitig) des randominiserten
mod_lm_sqrt2_org <- lm(prep$P_SQRT2 ~ sin_month_percent + 
                     cos_month_percent + as.factor(prep$season))
summary(mod_lm_sqrt2_org, se = "nid")
mod_lm_sqrt2_org <- as.data.frame(mod_lm_sqrt2_org$coefficients)



mod_bs <- lapply(seq(2000), function(x){
  # Jetzt sampeln, dann rechnen
  set.seed <- x
  boot <- sample(nrow(prep), replace = FALSE) 
  season <- as.factor(prep$season[boot])
  act_mod <- lm(prep$P_SQRT2 ~ sin_month_percent + 
                  cos_month_percent + season)
  return(act_mod$coefficients)
})
mod_bs <- as.data.frame(do.call("rbind", mod_bs))

quantile(mod_bs$season2003, probs = 0.975)
quantile(mod_bs$season2003, probs = 0.95)
#quantile(mod_bs$season2003, probs = 0.06)
quantile(mod_bs$season2003, probs = 0.05)
quantile(mod_bs$season2003, probs = 0.025)
mod_lm_sqrt2_org
sort(mod_bs$season2003)

quantile(mod_bs$season1983, probs = 0.975)
quantile(mod_bs$season1983, probs = 0.95)
quantile(mod_bs$season1983, probs = 0.05)
quantile(mod_bs$season1983, probs = 0.025)
mod_lm_sqrt2_org
sort(mod_bs$season1983)

quantile(mod_bs$season1993, probs = 0.975)
quantile(mod_bs$season1993, probs = 0.95)
quantile(mod_bs$season1993, probs = 0.05)
quantile(mod_bs$season1993, probs = 0.025)
mod_lm_sqrt2_org
sort(mod_bs$season1993)

quantile(mod_bs$season2003, probs = 0.975)
quantile(mod_bs$season2003, probs = 0.95)
quantile(mod_bs$season2003, probs = 0.05)
quantile(mod_bs$season2003, probs = 0.025)
mod_lm_sqrt2_org
sort(mod_bs$season2003)

quantile(mod_bs$season2013, probs = 0.975)
quantile(mod_bs$season2013, probs = 0.95)
quantile(mod_bs$season2013, probs = 0.05)
quantile(mod_bs$season2013, probs = 0.025)
mod_lm_sqrt2_org
sort(mod_bs$season2013)

###### bootstrap Insa
season <- as.factor(prep$season)

bs <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample 
  fit <- lm(formula, data=d)
  return(coef(fit)) 
}

mod_bs_in <- boot(data = prep, statistic = bs, 
                  R = 2000, formula = prep$P_SQRT2 ~ sin_month_percent + 
                    cos_month_percent + season)

mod_bs_in <- as.data.frame(do.call("rbind", mod_bs_in)) # sample?

mod_bs_in
plot(mod_bs_in, index=1) # intercept 
plot(mod_bs_in, index=2) # sin_month_percent
plot(mod_bs_in, index=3) # cos_month_precent
plot(mod_bs_in, index=4) # season


boot.ci(mod_bs_in, index=1) 
boot.ci(mod_bs_in, type="bca", index=2) 
boot.ci(mod_bs_in, type="bca", index=3) 
boot.ci(mod_bs_in, index=4)

mod_lm_sqrt2_org
########################################


################################
##
## Version Alberich
##
################################
test<-cosinor(P_MEAN_NEW ~ 1 + as.factor(season), date = "YEAR", data = prep)
test
summary(test)
# Achtung:
# Stationaritdt 
# Symmetrisches 

#andere Variante
t <-  as.numeric(prep$mnth)/12
sint <- sin(2*pi*t)
cost <- cos(2*pi*t)
test_f<-lm(prep$P_MEAN_NEW ~ sint + cost + as.factor(prep$season))
summary(test_f) 
#Vergleiche estimate von test und test_f f|r Faktor as.factor(seasonal_precip$SEASON)
#Vergleiche intercept



################################
##
## Version Mausi
##
################################

# Signifikanz der Dekanden als erklaerende Variable ----------------------------
# Lineares Model mit wurzeltransformiertem Niederschlag
p <- prep$P_MEAN_NEW**0.25
qqnorm(diff(p))

t <-  prep$mnth/12
f <- prep$season
sint <- sin(2*pi*t)
cost <- cos(2*pi*t)
sinf <- sin(2*pi*f)
cosf <- cos(2*pi*f)
l <- lm(p ~ sint + cost + sinf + cosf + f)
summary(l)


# GLM mit nicht-transformiertem Niederschlag
g <- glm(prep$P_MEAN_NEW ~ t + f)
summary(g)
gf <- glm(prep$P_MEAN_NEW ~ sint + cost + as.factor(f))
summary(gf)

# Library season ---------------------------------------------------------------
# Nur 1973-1982
p73 <- prep$P_MEAN_NEW[prep$season == "1973"]
m73 <- prep$mnth[prep$season == "1973"]
s73 <- prep$season[prep$season == "1973"]
df73 <- data.frame(P = p73, mnth = m73, season = s73)


# Nur 2003-2012
p03 <- prep$P_MEAN_NEW[prep$season == "2003"]
m03 <- prep$mnth[prep$season == "2003"]
s03 <- prep$season[prep$season == "2003"]
df03 <- data.frame(P = p03, mnth = m03, season = s03)

# Vergleich Amplitude und Phase (= Peak) fuer beide Dekaden
cnr73 <- cosinor(P ~ season, date = "mnth", data = df73, type = 'monthly', family = gaussian())
cnr03 <- cosinor(P ~ season, date = "mnth", data = df03, type = 'monthly', family = gaussian())

summary(cnr73)
summary(cnr03)

# Quantile regression ----------------------------------------------------------
# Quantile
m <- prep$YEAR + prep$mnth_precent

reg15 <- rq(prep$P_MEAN_NEW ~ m, tau = c(0.15))
summary(reg15, se = "nid")

reg95 <- rq(prep$P_MEAN_NEW ~ m, tau = c(0.95))
summary(reg95, se = "nid")


plot(prep$P_MEAN_NEW ~ m)
abline(reg95, col = "blue")
abline(reg15, col = "green")

