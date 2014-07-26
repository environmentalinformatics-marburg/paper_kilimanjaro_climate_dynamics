# Working directory
switch(Sys.info()[["sysname"]], 
       "Linux" = setwd("D:/active/kilimanjaro"), 
       "Windows" = setwd("D:/active/kilimanjaro"))

library(kza)
library(latticeExtra)
library(Kendall)
library(foreach)
library(ade4)
library(reshape2)
library(ggplot2)
source("scripts/paper_kilimanjaro_climate_dynamics/src/vectorHarmonics.R")
source("scripts/paper_kilimanjaro_climate_dynamics/src/combineAOPI.R")
source("scripts/paper_kilimanjaro_climate_dynamics/src/seasonPlotByAOI.R")
source("scripts/paper_kilimanjaro_climate_dynamics/src/corPlotByAOI.R")

outLayer <- function(x, y) {
  x + as.layer(y)
}

#### Long-term precipitation analysis ##########################################
# Compute long-term anomalies and trends for KIA precipitation and create
# publication quality figure.

# Read precipitation data and create a continous time series.
precip <- 
  read.csv("kilimanjaro_gsod_dynamics/gsod_precip/kia_prcp_1975_2014_mnthly.csv", 
           stringsAsFactors = FALSE)
precip[, 1] <- as.Date(precip[, 1])
st <- precip[1, 1]
nd <- precip[nrow(precip), 1]
ts <- seq(st, nd, "month")
precip <- merge(data.frame(ts), precip, by = 1, all.x = TRUE)

# Compute 3 month running mean of original precipitation values using a
# Kolmogorov-Zurbenko filter with one iteration
precip$kz03k01 <- kz(precip$P_RT_NRT, m = 3, k = 1)

# Compute deseasoned precipitation time series and corresponding 3 month running
# mean using a Kolmogorov-Zurbenko filter with two iterations to close gaps
precip$ssn <- precip$P_RT_NRT - rep(sapply(1:12, function(i) {
  mean(precip$P_RT_NRT[seq(i, nrow(precip), 12)], na.rm = TRUE)
}), nrow(precip) / 12)

precip$ssn_kz03k01 <- kz(precip$ssn, m = 3, k = 1)
precip$ssn_kz03k02 <- kz(precip$ssn, m = 3, k = 2)

clr <- as.character(ifelse(precip$ssn_kz03k02 > 0, "blue", "red"))
clr.blue <- brewer.pal(9, "Blues")[2]
clr.red <- brewer.pal(9, "Reds")[1]

# Create publication quality figure of long term precipitation trends
plot.precip.ssn_kz03k02 <- 
  xyplot(ssn_kz03k02 ~ ts, data = precip, origin = 0, type = "h",
         border = "transparent", col = clr, asp = 0.25,
         xlab = "", ylab = "Precipitation [mm]", 
         lwd = 1.5, ylim = c(-250, 350), as.table = TRUE,
         scales = list(x = list(axs = "i")),
         xscale.components = xscale.components.subticks,
         yscale.components = yscale.components.subticks,
         panel = function(x, y, ...) {
           panel.xblocks(x, y > 0, 
                         col = clr.blue)
           panel.xblocks(x, y < 0, 
                         col = clr.red)
           panel.xyplot(x, y, ...)
           panel.smoother(x, y, method = "lm", 
                          col = "black", 
                          col.se = "black",
                          alpha.se = 0.3, lty = 2)
         })
plot(plot.precip.ssn_kz03k02)
MannKendall(precip$ssn_kz03k02)



#### Seasonal precipitation analysis ###########################################
# Compute seasonal precipitation dynamics for each year of the time series
# centered arround Dezember and create publication quality figure. 

# Compute a three year running mean of the seasonal precipitation dynamics
precip.seasonal03yr <- 
  foreach(st = 1975:2003, nd = 1985:2013, .combine = "rbind") %do% {
    st.process <- grep(st, precip[, 1])[7]
    nd.process <- grep(nd, precip[, 1])[length(grep(nd, precip[, 1]))-6]
    precip.process <- precip[st.process:nd.process, ]
    data.frame(season = paste(st, nd, sep = "-"), 
               p_dyn = vectorHarmonics(precip.process$kz03k01, 
                                       st = c(st, 1), 
                                       nd = c(nd, 12), 
                                       m = 5))
  }

# Split the running mean data set by year and create publication quality figure;
# to smooth the cycle while not reducing the rainfall amounts significantly, 
# use a spline prediction.
precip.seasonal03yr.split <- 
  split(precip.seasonal03yr, precip.seasonal03yr$season)

colors <- colorRampPalette(
  brewer.pal(11, "RdYlGn"))(length(precip.seasonal03yr.split))

plot.precip.seasonal03yr.all <- seasonPlotByAOI(
  lapply(precip.seasonal03yr.split, function(x){x$p_dyn}), colors)
plot(plot.precip.seasonal03yr.all)

# Non-spline version
# plot.precip.seasonal03yr.split <- 
#   lapply(seq(precip.seasonal03yr.split), function(i) {
#     xyplot(p_dyn ~ factor(c(7:12, 1:6), levels = c(7:12, 1:6)), 
#            data = precip.seasonal03yr.split[[i]], type=c("l"), 
#            col = colors(length(precip.seasonal03yr.split))[i], 
#            ylim = c(0,170),
#            xlab = "Month", ylab = "Precipitation")
#   })
# plot.precip.seasonal03yr.all <- 
#   Reduce("outLayer", plot.precip.seasonal03yr.split)



#### Precipitation analysis vs ENSO ############################################
# Prepare aoi record and classifiy years as La Nina (L), El Nino (E) or 
# normal (N); weak ENSO cycles are classified as normal
aoi <- read.csv("enso/enso_and_iod.csv", skip = 1, header = TRUE)
aoi$Season <- paste0(aoi$Season, aoi$X, aoi$X.1)
aoi <- aoi[, -grep("X", names(aoi))]

aoi$TypeClass <- "N"
aoi$TypeClass[grep("L", aoi$Type)] <- "LS"
aoi$TypeClass[grep("WL", aoi$Type)] <- "LW"
aoi$TypeClass[grep("E", aoi$Type)] <- "ES"
aoi$TypeClass[grep("WE", aoi$Type)] <- "EW"
# aoi$TypeClass[grep("W", aoi$Type)] <- "N"
# aoi$TypeClass[grep("M", aoi$Type)] <- "X"

precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.aoi.split.median <- combineAOPI(aoi, precip.shift06m)

if(length(precip.shift06m.aoi.split.median) == 5){
  colors <- c("blue", "lightblue", "red", "bisque", "black")
} else {
  colors <- c("blue", "red", "black")  
}

plot.precip.shift06m.aoi.split.median.all <- 
  seasonPlotByAOI(precip.shift06m.aoi.split.median, colors)

plot(plot.precip.shift06m.aoi.split.median.all)

# Non-spline version
# plot.precip.shift06m.aoi.split.median <- 
#   lapply(seq(precip.shift06m.aoi.split.median), function(x){
#     xyplot(precip.shift06m.aoi.split.median[[x]] ~ 
#              factor(c(7:12, 1:6), levels = c(7:12, 1:6)), type = "l", 
#            ylim = c(0, 200), col = colors[x])
#   })
# 
# plot.precip.shift06m.aoi.split.median.all <- 
#   Reduce("outLayer", plot.precip.shift06m.aoi.split.median)
# plot.precip.shift06m.aoi.split.median.all <- 
#   Reduce("outLayer", plot.precip.shift06m.aoi.split.median)

# Just for background information: the individual seasons which form the above
# combined El Nino, La Nina and normal situation spline.

precip.shift06m.aoi.split <- combineAOPI(aoi, precip.shift06m, rt = "split")

precip.shift06m.aoi.split.season <- lapply(precip.shift06m.aoi.split, function(x){
  split(x, x$Season)
})

plot.precip.shift06m.aoi.split.season.all <- 
  seasonPlotByAOI(precip.shift06m.aoi.split.season, colors, individual = TRUE)
plot(plot.precip.shift06m.aoi.split.season.all)

# Alternative boxplot visualization
## bwplots
precip.shift06m.aoi <- combineAOPI(aoi, precip.shift06m, rt = "org")
precip.shift06m.aoi$month <- substr(precip.shift06m.aoi$ts, 6, 7)
precip.shift06m.aoi.melt <- melt(precip.shift06m.aoi, 
                                 id.vars = c("month", "TypeClass"),
                                 measure.vars = "P_RT_NRT")
bw.plot <- ggplot(data = precip.shift06m.aoi.melt, 
                  aes(y = value, x = month, fill = TypeClass, 
                      dodge = TypeClass))
bw.plot + geom_boxplot()

# Create publication quality figures of correlations between aoi and 
# precipitation
years <- unique(substr(precip$ts, 1, 4))
precip.shift06m.mat.ssn_kz03k01 <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.aoi[grep(i, precip.shift06m.aoi$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$ssn_kz03k01[j]
  }
}

aoi.reshape.shift06m.mat <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.aoi[grep(i, precip.shift06m.aoi$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$oni[j]
  }
}

plot.precip.shift06m.oni.ssn_kz03k01 <- 
  corPlotByAOI(aoi.reshape.shift06m.mat, precip.shift06m.mat.ssn_kz03k01)
plot(plot.precip.shift06m.oni.ssn_kz03k01)



#### Precipitation analysis vs IOD #############################################
# Prepare aoi record and classifiy years as IOD plus (P), IOD minus (M) or
# normal (N)
aoi$TypeClass <- "N"
aoi$TypeClass[grep("E", aoi$Type)] <- "X"
aoi$TypeClass[grep("L", aoi$Type)] <- "X"
aoi$TypeClass[grep("P", aoi$IOD)] <- "IP"
aoi$TypeClass[grep("M", aoi$IOD)] <- "IM"

precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.iod.split.median <- combineAOPI(aoi, precip.shift06m)

if(length(precip.shift06m.iod.split.median) == 5){
  colors <- c("blue", "lightblue", "red", "bisque", "black")
} else {
  colors <- c("blue", "red", "black")  
}

plot.precip.shift06m.iod.split.median <- 
  seasonPlotByAOI(precip.shift06m.iod.split.median, colors)
plot(plot.precip.shift06m.iod.split.median)

# Alternative boxplot visualization
## bwplots
precip.shift06m.aoi <- combineAOPI(aoi, precip.shift06m, rt = "org")
precip.shift06m.aoi$month <- substr(precip.shift06m.aoi$ts, 6, 7)
precip.shift06m.aoi.melt <- melt(precip.shift06m.aoi, 
                                 id.vars = c("month", "TypeClass"),
                                 measure.vars = "P_RT_NRT")
bw.plot <- ggplot(data = precip.shift06m.aoi.melt, 
                  aes(y = value, x = month, fill = TypeClass, 
                      dodge = TypeClass))
bw.plot + geom_boxplot()

# Create publication quality figures of correlations between aoi and 
# precipitation
years <- unique(substr(precip$ts, 1, 4))
precip.shift06m.mat.ssn_kz03k01 <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.aoi[grep(i, precip.shift06m.aoi$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$ssn_kz03k01[j]
  }
}

aoi.reshape.shift06m.mat <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.aoi[grep(i, precip.shift06m.aoi$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$IOD[j]
  }
}

plot.precip.shift06m.oni.ssn_kz03k01 <- 
  corPlotByAOI(aoi.reshape.shift06m.mat, precip.shift06m.mat.ssn_kz03k01)
plot(plot.precip.shift06m.oni.ssn_kz03k01)


# test <- subset(precip.shift06m.aoi, ts >= "1997-01-01" & ts <= "2000-12-31")
# findMaxCCF(as.numeric(scale(test$ssn_kz03k02)), test$aoi, na.action = na.exclude)
# 
# foreach(i = list(c("-03-", "-04-", "-05-"), 
#                  c("-06-", "-07-", "-08-"), 
#                  c("-09-", "-10-", "-11-"), 
#                  c("-12-", "-01-", "-02-"))) %do% {
#   sub <- subset(precip.shift06m.aoi, ts >= "1997-01-01" & ts <= "2000-12-31")                 
#   index <- substr(sub$ts, 5, 8) %in% i
#   test2 <- sub[index, ]
#   findMaxCCF(as.numeric(scale(test2$ssn_kz03k02)), test2$aoi, na.action = na.exclude)
# }
# 
# foreach(i = c("-04-", "-08-", "-12-")) %do% {
#   sub <- subset(precip.shift06m.aoi, ts >= "1975-01-01" & ts <= "2013-12-31")                 
#   index <- substr(sub$ts, 5, 8) %in% i
#   test2 <- sub[index, ]
#   findMaxCCF(as.numeric(scale(test2$ssn_kz03k02)), test2$aoi, na.action = na.exclude)
# }

# aoi$IODAgg <- "N"
# aoi$IODAgg[grep("P", aoi$IOD)] <- "P"
# aoi$IODAgg[grep("M", aoi$IOD)] <- "M"
# 
# aoi.reshape <- do.call("rbind", lapply(seq_len(nrow(aoi)), function(i) {
#   st <- as.Date(paste0(substr(aoi[i, 3], 1, 4), "-07-01"))
#   nd <- as.Date(paste0(substr(aoi[i, 3], 6, 9), "-06-01"))
#   data.frame(aoi[i, c(1:3, (ncol(aoi)-1), ncol(aoi))], 
#              aoi = as.numeric(aoi[i, 4:(ncol(aoi)-2)]), 
#              month = seq(st, nd, "month"))
# }))
# 
# precip.shift06m.aoi <- 
#   merge(precip.shift06m, aoi.reshape, all.x = TRUE, by.x = "ts", by.y = "month")
# 
# precip.shift06m.aoi.split <- 
#   split(precip.shift06m.aoi, precip.shift06m.aoi$IODAgg)
# 
# precip.shift06m.aoi.split.median <- 
#   foreach(i = precip.shift06m.aoi.split) %do% {
#     sapply(1:12, function(j) {
#       median(i$P_RT_NRT[seq(j, nrow(i), 12)], na.rm = TRUE)
#     })
#   }
# 
# plot(precip.shift06m.aoi.split[[3]][1:12, "P_RT_NRT"], type = "l", ylim = c(0, 350), 
#      col = "blue")
# lines(precip.shift06m.aoi.split[[3]][13:36, "P_RT_NRT"], col = "red")
# lines(precip.shift06m.aoi.split[[3]][37:48, "P_RT_NRT"], col = "orange")
# lines(precip.shift06m.aoi.split.median[[3]], col = "black")
# 
# vectorHarmaoics(precip.shift06m.aoi.split[[3]]$P_RT_NRT, st = c(1990, 1), 
#                 nd = c(2002, 12))

