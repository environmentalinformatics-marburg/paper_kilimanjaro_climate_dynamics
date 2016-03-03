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
colnames(rf)[1] <- "Year"
rf$Year <- as.Date(as.character(rf$Year))
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

# Rainfall pre-processing trends
rfa <- lapply(seq(2, length(colnames(rf))), function(x){
  ret <- aggregate(rf[, x], by = list(substr(rf$Year, 1, 4)), FUN = "sum")
  ret$Group.1 <- as.Date(ret$Group.1, format = "%Y")
  colnames(ret) <- c("Year", "Rainfall")
  ret$Parameter <- colnames(rf)[x]
  return(ret)
})
rfa <- do.call("rbind", rfa)
rfa$ann <- as.numeric(substr(rfa$Year, 1, 4))
rfa <- merge(rfa, ensoiod[, 1:3], by.x = "ann", by.y = "Season")
rfa$ElNino <- 0
rfa$ElNino[substr(rfa$ENSO, 2, 2) == "E"] <- 1
rfa$LaNina <- 0
rfa$LaNina[substr(rfa$ENSO, 2, 2) == "L"] <- 1


rf$mnth <- as.numeric(substr(rf$Year, 6, 7))
rf$mnth_precent <- rf$mnth / 12
rf$ann <- as.numeric(substr(rf$Year, 1, 4))
rf$date_precent <- rf$ann + rf$mnth_precent
rf$P_MEAN_NEW_LOG <- log(rf$P_MEAN_NEW + 1)
rf$P_MEAN_NEW_SQRT2 <- rf$P_MEAN_NEW**0.25

P_MEAN_NEW_mmean <- aggregate(rf$P_MEAN_NEW, by = list(rf$mnth), FUN = mean)
rf$P_MEAN_NEW_ds <- rf$P_MEAN_NEW - rep(sapply(1:12, function(x){
  P_MEAN_NEW_mmean[x, "x"]}), nrow(rf)/12)

rf$P_MEAN_NEW_ds_SQRT2 <- (rf$P_MEAN_NEW_ds + 500)**0.25

rf$P_MEAN_NEW_ds_LOG <- (rf$P_MEAN_NEW_LOG + 500)**0.25

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
rf$drywet[rf$mnth >= 6 & rf$mnth <= 9] <- "dry_mid"

rf <- merge(rf, ensoiod[, 1:3], by.x = "ann", by.y = "Season")

rf$Season <- (c(0,0,0,0,0,0, rep(seq(1:40),each=12), 0,0,0,0,0,0))

rf$ElNino <- 0
rf$ElNino[substr(rf$ENSO, 2, 2) == "E"] <- 1

rf$LaNina <- 0
rf$LaNina[substr(rf$ENSO, 2, 2) == "L"] <- 1

save(atk, atka, atkagaps, atkgaps, atn, atna, atnagaps, atngaps, ensoiod, 
     P_MEAN_NEW_mmean, rf, rfa, rfei, rfs, file = "trend_analysis.RData")
