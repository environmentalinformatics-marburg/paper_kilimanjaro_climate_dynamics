#### Set working environment ###################################################

library(dplR)
library(ggplot2)
library(multitaper)
library(quantreg)
library(tseries)

setwd("D:/active/paper_kilimanjaro_climate_dynamics")
atk <- read.table("noaa_gsod_kia_1973_2013_reformat_sd_li_jl_ssa.csv", 
                  header = TRUE, sep = ",", dec = ".")
atn <- read.table("noaa_gsod_nairobi_1973_2013_reformat_sd_li_jl_ssa.csv", 
                  header = TRUE, sep = ",", dec = ".")
atkgaps <- read.table("noaa_gsod_kia_1973_2013_reformat.csv", 
                      header = TRUE, sep = ",", dec = ".")
atngaps <- read.table("noaa_gsod_nairobi_1973_2013_reformat.csv", 
                      header = TRUE, sep = ",", dec = ".")
rf <- read.table("metoffice_1973-2013.csv", header = TRUE, sep = ";", dec = ".")
rfs <- read.table("metoffice_1973-2013_season.csv", header = TRUE, sep = ";", dec = ".")
rfei <- read.table("rainfall_enso_iod.csv", header = TRUE, sep = ";")
ensoiod <- read.table("enso_and_iod.csv", header = TRUE, sep = ",")


atk$Datetime <- as.Date(atk$Datetime)
atn$Datetime <- as.Date(atn$Datetime)
atkgaps$Datetime <- as.Date(atkgaps$Datetime)
atngaps$Datetime <- as.Date(atngaps$Datetime)
rf$YEAR <- as.Date(as.character(rf$YEAR))
colnames(ensoiod)[2] <- "ENSO"


atka <- lapply(seq(9, length(colnames(atk))), function(x){
  ret <- aggregate(atk[, x], by = list(substr(atk$Datetime, 1, 4)), FUN = "mean")
  ret$Group.1 <- as.Date(ret$Group.1, format = "%Y")
  colnames(ret) <- c("Year", "Temperature")
  ret$Parameter <- colnames(atk)[x]
  return(ret)
})
atka <- do.call("rbind", atka)

atna <- lapply(seq(9, length(colnames(atn))), function(x){
  ret <- aggregate(atn[, x], by = list(substr(atn$Datetime, 1, 4)), FUN = "mean")
  ret$Group.1 <- as.Date(ret$Group.1, format = "%Y")
  colnames(ret) <- c("Year", "Temperature")
  ret$Parameter <- colnames(atn)[x]
  return(ret)
})
atna <- do.call("rbind", atna)

atkagaps <- lapply(seq(9, length(colnames(atkgaps))), function(x){
  ret <- aggregate(atkgaps[, x], by = list(substr(atkgaps$Datetime, 1, 4)), 
                   FUN = "mean", na.rm=TRUE, na.action=NULL)
  ret$Group.1 <- as.Date(ret$Group.1, format = "%Y")
  colnames(ret) <- c("Year", "Temperature")
  ret$Parameter <- colnames(atkgaps)[x]
  return(ret)
})
atkagaps <- do.call("rbind", atkagaps)

atnagaps <- lapply(seq(9, length(colnames(atngaps))), function(x){
  ret <- aggregate(atngaps[, x], by = list(substr(atngaps$Datetime, 1, 4)), 
                   FUN = "mean", na.rm=TRUE, na.action=NULL)
  ret$Group.1 <- as.Date(ret$Group.1, format = "%Y")
  colnames(ret) <- c("Year", "Temperature")
  ret$Parameter <- colnames(atngaps)[x]
  return(ret)
})
atnagaps <- do.call("rbind", atnagaps)

rfa <- lapply(seq(2, length(colnames(rf))), function(x){
  ret <- aggregate(rf[, x], by = list(substr(rf$YEAR, 1, 4)), FUN = "sum")
  ret$Group.1 <- as.Date(ret$Group.1, format = "%Y")
  colnames(ret) <- c("Year", "Rainfall")
  ret$Parameter <- colnames(rf)[x]
  return(ret)
})
rfa <- do.call("rbind", rfa)

# Air temperature KIA
for (x in unique(atka$Parameter)){
  print(x)
  print(summary(lm(Temperature ~ time(Year), data = atka[atka$Parameter == x,])))
}

# Air temperature KIA gaps
for (x in unique(atkagaps$Parameter)){
  print(x)
  print(summary(lm(Temperature ~ time(Year), data = atkagaps[atkagaps$Parameter == x,])))
}

# Air temperature Nairobi
for (x in unique(atna$Parameter)){
  print(x)
  print(summary(lm(Temperature ~ time(Year), data = atna[atna$Parameter == x,])))
}

# Air temperature Nairobi gaps
for (x in unique(atnagaps$Parameter)){
  print(x)
  print(summary(lm(Temperature ~ time(Year), data = atnagaps[atnagaps$Parameter == x,])))
}


# Rainfall pre-processing trends
rf$mnth <- as.numeric(substr(rf$YEAR, 6, 7))
rf$mnth_precent <- rf$mnth / 12
rf$ann <- as.numeric(substr(rf$YEAR, 1, 4))
rf$date_precent <- rf$ann + rf$mnth_precent
rf$P_MEAN_NEW_LOG <- log(rf$P_MEAN_NEW + 1)
rf$P_MEAN_NEW_SQRT2 <- rf$P_MEAN_NEW**0.25

P_MEAN_NEW_mmean <- aggregate(rf$P_MEAN_NEW, by = list(rf$mnth), FUN = mean)
rf$P_MEAN_NEW_ds <- rf$P_MEAN_NEW - rep(sapply(1:12, function(x){
  P_MEAN_NEW_mmean[x, "x"]}), nrow(rf)/12)

df <- lapply(seq(12), function(x){
  m <- rep(0, nrow(rf))
  m[rf$mnth == x] <- 1
  df <- data.frame(M = m)
  colnames(df) <- paste0("M", as.character(x))
  return(df)
})
rf <- cbind(rf, do.call("cbind", df))

rf$short_rains <- 0
rf$short_rains[rf$mnth >= 10 & rf$mnth <= 12 | rf$mnth == 1] <- 1

rf$long_rains <- 0
rf$long_rains[rf$mnth >= 3 & rf$mnth <= 5] <- 1

rf$drywet <- "dry"
rf$drywet[rf$long_rains == 1] <- "long_rains"
rf$drywet[rf$short_rains == 1] <- "shor_rains"
rf$drywet[rf$mnth >= 1 & rf$mnth <= 2] <- "dry_mid"

rf <- merge(rf, ensoiod[, 1:3], by.x = "ann", by.y = "Season")

rf$Season <- (c(0,0,0,0,0,0, rep(seq(1:40),each=12), 0,0,0,0,0,0))

rf$ElNino <- 0
rf$ElNino[substr(rf$ENSO, 2, 2) == "E"] <- 1

rf$LaNina <- 0
rf$LaNina[substr(rf$ENSO, 2, 2) == "L"] <- 1


rfts <- ts(rf$P_MEAN_NEW, start = 1973, end = 2013, frequency = 12)
plot(stl(rfts, 10, t.window = 120))


# Rainfall annual trends
for (x in unique(rfa$Parameter)){
  print(x)
  print(summary(lm(Rainfall ~ time(Year), data = rfa[rfa$Parameter == x,])))
}
print(summary(rq(Rainfall ~ time(Year), tau = seq(0.8, 1, 0.05),
                 data = rfa[rfa$Parameter == "P_MEAN_NEW",])))

# Time series stationarity
kpss.test(rf$P_MEAN_NEW)
kpss.test(rf$P_MEAN_NEW, null = "Trend")


# Rainfall frequency
k = kernel("daniell", c(3,3))
rfspec <- spec.pgram(rf$P_MEAN_NEW, kernel = k, taper = 0, plot = FALSE, detrend = TRUE)

rfspecdf <- data.frame(freq = rfspec$freq, spec = rfspec$spec)
ticks <- c(30, 20, 10, 5, 3, 1, 0.5, 0.25)
ticks <- c(30, 20, seq(12, 1, -1), seq(1, 0, -0.25))
labels <- as.character(ticks)

breaks <- 1/(ticks * 12)

ggplot(rfspecdf, aes(x = freq, y = spec)) + 
  geom_line() + 
  scale_x_log10("Period (years)", breaks = breaks, labels = labels) + 
  scale_y_log10()


rfwave <- morlet(y1 = rf$P_MEAN_NEW, x1 = time(rf$YEAR), p2 = 8, dj = 0.1, siglvl = 0.95)
rfwave$period <- rfwave$period/12
levels <- quantile(rfwave$Power, c(0, 0.25, 0.5, 0.75, 0.95, 1))
wavelet.plot(rfwave, wavelet.levels = levels, crn.ylim = c(22.5, 30))

rfwave_avg <- data.frame(power = apply(wave.out$Power, 2, mean), period = (wave.out$period))
ggplot(rfwave_avg, aes(x = period, y = power)) + 
  geom_line() + 
  scale_x_continuous(breaks = seq(25)) + 
  scale_y_continuous()

# Rainfall frequency de-seasoned
k = kernel("daniell", c(3,3))
rfspec <- spec.pgram(rf$P_MEAN_NEW_ds, kernel = k, taper = 0, plot = FALSE, detrend = TRUE)

rfspecdf <- data.frame(freq = rfspec$freq, spec = rfspec$spec)
ticks <- c(30, 20, 10, 5, 3, 1, 0.5, 0.25)
ticks <- c(30, 20, seq(12, 1, -1), seq(1, 0, -0.25))
labels <- as.character(ticks)

breaks <- 1/(ticks * 12)

ggplot(rfspecdf, aes(x = freq, y = spec)) + 
  geom_line() + 
  scale_x_log10("Period (years)", breaks = breaks, labels = labels) + 
  scale_y_log10()


rfwave <- morlet(y1 = rf$P_MEAN_NEW_ds, x1 = time(rf$YEAR), p2 = 8, dj = 0.1, siglvl = 0.95)
rfwave$period <- rfwave$period/12
levels <- quantile(rfwave$Power, c(0, 0.25, 0.5, 0.75, 0.95, 1))
wavelet.plot(rfwave, wavelet.levels = levels, crn.ylim = c(22.5, 30))

rfwave_avg <- data.frame(power = apply(wave.out$Power, 2, mean), period = (wave.out$period))
ggplot(rfwave_avg, aes(x = period, y = power)) + 
  geom_line() + 
  scale_x_continuous(breaks = seq(25)) + 
  scale_y_continuous()



# Yearly rainfall sums
ggplot(rf, aes(x = mnth, y = P_MEAN_NEW, colour = as.factor(Season))) +
  geom_bar(stat = "identity", position = "dodge") +
  xlab("Year") +
  ylab("Precipitation")

# Rainfall trends
rflm <- lm(P_MEAN_NEW_ds ~  time(YEAR), data = rf[rf$mnth == 4,])
summary(rflm)
anova(rflm)

plot(rflm)

rflm <- lm(P_MEAN_NEW_SQRT2 ~  time(YEAR) + short_rains * time(YEAR), data = rf)
summary(rflm)
anova(rflm)

rflm <- lm(P_MEAN_NEW_SQRT2 ~  time(YEAR) + long_rains * time(YEAR), data = rf)
summary(rflm)
anova(rflm)

rflm <- lm(P_MEAN_NEW_SQRT2 ~  time(YEAR) + long_rains * time(YEAR) + short_rains * time(YEAR), data = rf)
summary(rflm)
anova(rflm)



rflm <- lm(P_MEAN_NEW_ds ~  time(YEAR) + as.factor(mnth) * time(YEAR), data = rf)
summary(rflm)
anova(rflm)

rflm <- lm(P_MEAN_NEW_LOG ~  time(YEAR) + drywet * time(YEAR), data = rf)
summary(rflm)
anova(rflm)








rflm <- lm(P_MEAN_NEW_ds ~  time(YEAR) + long_rains * time(YEAR) + short_rains * time(YEAR) +
             M1 * time(YEAR) +
             M3 * time(YEAR) + M4 * time(YEAR) +  M5 * time(YEAR) + 
             M6 * time(YEAR) + M7 * time(YEAR) + M8 * time(YEAR) + 
             M9 * time(YEAR) + M10 * time(YEAR) + M11 * time(YEAR) +
             M12 * time(YEAR), data = rf)
summary(rflm)
anova(rflm)

rflm <- lm(P_MEAN_NEW_ds ~  time(YEAR) + long_rains * time(YEAR) + short_rains * time(YEAR) +
             M1 * time(YEAR) + M2 * time(YEAR) +
             M3 * time(YEAR) + M4 * time(YEAR) +  M5 * time(YEAR) + 
             M6 * time(YEAR) + M7 * time(YEAR) + M8 * time(YEAR) + 
             M9 * time(YEAR) + M10 * time(YEAR) + M11 * time(YEAR) +
             M12 * time(YEAR) + 
             IOD * time(YEAR) + ENSO * time(YEAR), data = rf)
summary(rflm)
anova(rflm)

rflm <- lm(P_MEAN_NEW_LOG ~  time(YEAR) + long_rains * time(YEAR) + short_rains * time(YEAR) +
             M1 * time(YEAR) + M2 * time(YEAR) +
             M3 * time(YEAR) + M4 * time(YEAR) +  M5 * time(YEAR) + 
             M6 * time(YEAR) + M7 * time(YEAR) + M8 * time(YEAR) + 
             M9 * time(YEAR) + M10 * time(YEAR) + M11 * time(YEAR) +
             M12 * time(YEAR) + 
             IOD * time(YEAR) + ENSO * time(YEAR), data = rf)
summary(rflm)
anova(rflm)

rflm <- lm(P_MEAN_NEW ~  time(YEAR) + long_rains * time(YEAR) + short_rains * time(YEAR) +
             M1 * time(YEAR) + M2 * time(YEAR) +
             M3 * time(YEAR) + M4 * time(YEAR) +  M5 * time(YEAR) + 
             M6 * time(YEAR) + M7 * time(YEAR) + M8 * time(YEAR) + 
             M9 * time(YEAR) + M10 * time(YEAR) + M11 * time(YEAR) +
             M12 * time(YEAR) + 
             IOD * time(YEAR) + ElNino * time(LaNina), data = rf)
summary(rflm)
anova(rflm)
