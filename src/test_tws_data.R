#### Set working environment ###################################################

library(ggplot2)
library(quantreg)

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


atk$Datetime <- as.Date(atk$Datetime)
atn$Datetime <- as.Date(atn$Datetime)
atkgaps$Datetime <- as.Date(atkgaps$Datetime)
atngaps$Datetime <- as.Date(atngaps$Datetime)
rf$YEAR <- as.Date(as.character(rf$YEAR))



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

# Rainfall
for (x in unique(rfa$Parameter)){
  print(x)
  print(summary(lm(Rainfall ~ time(Year), data = rfa[rfa$Parameter == x,])))
}


print(summary(rq(Rainfall ~ time(Year), tau = seq(0.8, 1, 0.05),
                 data = rfa[rfa$Parameter == "P_MEAN_NEW",])))

print(summary(rq(P_MEAN_NEW ~ time(YEAR), tau = seq(0.5, 1, 0.05), data = rf)))
