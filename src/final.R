#### Set working environment ###################################################
switch(Sys.info()[["sysname"]], 
       "Linux" = setwd("/media/permanent/active/kilimanjaro"), 
       "Windows" = setwd("D:/active/kilimanjaro"))

sourcePath <- "scripts/paper_kilimanjaro_climate_dynamics/src/"
dataPath <- "data/"
graphicsPath <- "graphics/"

printToFile = TRUE

library(kza)
library(latticeExtra)
library(Kendall)
library(foreach)
library(ade4)
library(reshape2)
library(ggplot2)
source(paste0(sourcePath, "vectorHarmonics.R"))
source(paste0(sourcePath, "combineAOPI.R"))
source(paste0(sourcePath, "seasonPlotByAOI.R"))
source(paste0(sourcePath, "corPlotByAOI.R"))



#### Functions #################################################################
outLayer <- function(x, y) {
  x + as.layer(y)
}



#### Read data sets ############################################################
# Read precipitation data and create a continous time series.
precip <- 
  read.csv(
    paste0(dataPath, 
           "kilimanjaro_gsod_dynamics/gsod_precip/kia_prcp_1975_2014_mnthly.csv"), 
           stringsAsFactors = FALSE)
precip[, 1] <- as.Date(precip[, 1])
st <- precip[1, 1]
nd <- precip[nrow(precip), 1]
ts <- seq(st, nd, "month")
precip <- merge(data.frame(ts), precip, by = 1, all.x = TRUE)

# Read ENSO ONI index
oni <- read.csv(paste0(dataPath, "aoi/enso_and_iod.csv"),
                skip = 1, header = TRUE)
oni$Season <- paste0(oni$Season, oni$X, oni$X.1)
oni <- oni[, -grep("X", names(oni))]

# Read ENSO MEI index
mei <- read.csv(paste0(dataPath, "aoi/mei.txt"), skip = 9, header = TRUE, 
                sep = "\t", nrows = 65)
mei <- mei[,c(1,8:13,2:7)]
mei <- do.call(rbind, lapply(seq(nrow(mei)-1), function(x){
  modrow <- mei[x,]
  modrow[,8:13] <- mei[x+1,8:13]
  return(modrow)
}))
mei$Season <- paste0(mei$YEAR,"-",mei$YEAR+1)
mei <- cbind(oni[,1:3],mei[,2:13])

# Read DMI index
dmi <- read.csv(paste0(dataPath, "aoi/dmi.txt"), header = FALSE, sep = ":")
colnames(dmi) <- c("Year", "Month", "Day", "Hour", "DMI")
dmi$YearMonth <- paste0(dmi$Year,"-",sprintf("%02d", dmi$Month))
dmi <- aggregate(dmi$DMI, by=list(dmi$YearMonth), FUN = "mean")
colnames(dmi) <- c("YearMonth", "DMI")
years <- unique(substr(dmi$YearMonth, 1, 4))
dmi <- as.data.frame(foreach(i = years, .combine = "rbind") %do% {
  dat <- dmi[grep(i, dmi$YearMonth), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$DMI[j]
  }
})
dmi$Year <- years
dmi <- dmi[,c(13,7:12,1:6)]
dmi <- do.call(rbind, lapply(seq(nrow(dmi)-1), function(x){
  modrow <- dmi[x,]
  modrow[,8:13] <- dmi[x+1,8:13]
  return(modrow)
}))
dmi$Season <- paste0(as.numeric(dmi$Year),"-",as.numeric(dmi$Year)+1)
fr <- grep(dmi$Season[1], oni$Season)
lr <- grep(dmi$Season[nrow(dmi)], oni$Season)
dmi <- cbind(oni[fr:lr,1:3],dmi[,2:13])
colnames(dmi)[4:15] <- colnames(oni)[4:15]

# Read DMI-sstHAD index
dmiHAD <- read.csv(paste0(dataPath, "aoi/dmi_hadisst.txt"), 
                   header = FALSE, sep = " ")
colnames(dmiHAD) <- c("Year", "Month", "DMIHAT")
dmiHAD$YearMonth <- paste0(dmiHAD$Year,"-",sprintf("%02d", dmiHAD$Month))
dmiHAD <- dmiHAD[,c(4,3)]
years <- unique(substr(dmiHAD$YearMonth, 1, 4))
dmiHAD <- as.data.frame(foreach(i = years, .combine = "rbind") %do% {
  dat <- dmiHAD[grep(i, dmiHAD$YearMonth), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$DMIHAT[j]
  }
})
dmiHAD$Year <- years
dmiHAD <- dmiHAD[,c(13,7:12,1:6)]
dmiHAD <- do.call(rbind, lapply(seq(nrow(dmiHAD)-1), function(x){
  modrow <- dmiHAD[x,]
  modrow[,8:13] <- dmiHAD[x+1,8:13]
  return(modrow)
}))
dmiHAD$Season <- paste0(as.numeric(dmiHAD$Year),"-",as.numeric(dmiHAD$Year)+1)
fr <- grep(dmiHAD$Season[1], oni$Season)
lr <- grep(dmiHAD$Season[nrow(dmiHAD)], oni$Season)
dmiHAD <- cbind(oni[fr:lr,1:3],dmiHAD[,2:13])
colnames(dmiHAD)[4:15] <- colnames(oni)[4:15]

# Combine aoi data sets
aoi.list <- list(ONI = oni, MEI = mei, DMI = dmi, DMIHAD = dmiHAD)



#### Long-term precipitation analysis ##########################################
# Compute long-term anomalies and trends for KIA precipitation and create
# publication quality figure.
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
MannKendall(precip$ssn_kz03k02)
if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, "plot.precip.ssn_kz03k02.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.ssn_kz03k02)
  dev.off()
  
} else {
  plot(plot.precip.ssn_kz03k02)  
}



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

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, "plot.precip.seasonal03yr.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.seasonal03yr.all)
  dev.off()
  
} else {
  plot(plot.precip.seasonal03yr.all)
}

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
oni <- aoi.list$ONI
oni$TypeClass <- "N"
oni$TypeClass[grep("L", oni$Type)] <- "LS"
oni$TypeClass[grep("WL", oni$Type)] <- "LW"
oni$TypeClass[grep("E", oni$Type)] <- "ES"
oni$TypeClass[grep("WE", oni$Type)] <- "EW"
# oni$TypeClass[grep("W", oni$Type)] <- "N"
# oni$TypeClass[grep("M", oni$Type)] <- "X"

precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.oni.split.median <- combineAOPI(oni, precip.shift06m)

if(length(precip.shift06m.oni.split.median) == 5){
  colors <- c("blue", "lightblue", "red", "bisque", "black")
} else {
  colors <- c("blue", "red", "black")  
}

plot.precip.shift06m.oni.split.median.all <- 
  seasonPlotByAOI(precip.shift06m.oni.split.median, colors)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.oni.split.median.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.oni.split.median.all)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.oni.split.median.all)
}

# Non-spline version
# plot.precip.shift06m.oni.split.median <- 
#   lapply(seq(precip.shift06m.oni.split.median), function(x){
#     xyplot(precip.shift06m.oni.split.median[[x]] ~ 
#              factor(c(7:12, 1:6), levels = c(7:12, 1:6)), type = "l", 
#            ylim = c(0, 200), col = colors[x])
#   })
# 
# plot.precip.shift06m.oni.split.median.all <- 
#   Reduce("outLayer", plot.precip.shift06m.oni.split.median)
# plot.precip.shift06m.oni.split.median.all <- 
#   Reduce("outLayer", plot.precip.shift06m.oni.split.median)

# Just for background information: the individual seasons which form the above
# combined El Nino, La Nina and normal situation spline.

precip.shift06m.oni.split <- combineAOPI(oni, precip.shift06m, rt = "split")

precip.shift06m.oni.split.season <- lapply(precip.shift06m.oni.split, function(x){
  split(x, x$Season)
})

plot.precip.shift06m.oni.split.season.all <- 
  seasonPlotByAOI(precip.shift06m.oni.split.season, colors, individual = TRUE)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.oni.split.season.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.oni.split.season.all)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.oni.split.season.all)
}



# Alternative boxplot visualization
## bwplots
precip.shift06m.oni <- combineAOPI(oni, precip.shift06m, rt = "org")
precip.shift06m.oni$month <- substr(precip.shift06m.oni$ts, 6, 7)
precip.shift06m.oni.melt <- melt(precip.shift06m.oni, 
                                 id.vars = c("month", "TypeClass"),
                                 measure.vars = "P_RT_NRT")
bw.plot <- ggplot(data = precip.shift06m.oni.melt, 
                  aes(y = value, x = month, fill = TypeClass, 
                      dodge = TypeClass))
plot.precip.shift06m.oni.split.median.all.bw <- bw.plot + geom_boxplot() +
  xlab("Month") + ylab("Preciptiation")

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.oni.split.median.all.bw.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.oni.split.median.all.bw)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.oni.split.median.all.bw)
}

# Create publication quality figures of correlations between oni and 
# precipitation
years <- unique(substr(precip$ts, 1, 4))
precip.shift06m.mat.ssn_kz03k01 <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.oni[grep(i, precip.shift06m.oni$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$ssn_kz03k01[j]
  }
}

oni.reshape.shift06m.mat <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.oni[grep(i, precip.shift06m.oni$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$oni[j]
  }
}

plot.precip.shift06m.oni.ssn_kz03k01 <- 
  corPlotByAOI(oni.reshape.shift06m.mat, precip.shift06m.mat.ssn_kz03k01)


if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.oni.ssn_kz03k01.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.oni.ssn_kz03k01)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.oni.ssn_kz03k01)
}




#### Precipitation analysis vs IOD #############################################
# Prepare aoi record and classifiy years as IOD plus (P), IOD minus (M) or
# normal (N)
dmi <- aoi.list$DMIHAD
dmi$TypeClass <- "N"
dmi$TypeClass[grep("P", dmi$IOD)] <- "I1P"
dmi$TypeClass[grep("M", dmi$IOD)] <- "I2M"
#dmi$TypeClass[grep("E", dmi$Type)] <- "X"
#dmi$TypeClass[grep("L", dmi$Type)] <- "X"

precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.iod.split.median <- combineAOPI(dmi, precip.shift06m)

if(length(precip.shift06m.iod.split.median) == 5){
  colors <- c("blue", "lightblue", "red", "bisque", "black")
} else if (length(precip.shift06m.iod.split.median) == 4){
  colors <- c("blue", "red", "black", "bisque")
} else {
  colors <- c("blue", "red", "black")  
}

plot.precip.shift06m.dmi.split.median <- 
  seasonPlotByAOI(precip.shift06m.iod.split.median, colors)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.dmi.split.median.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.dmi.split.median)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.dmi.split.median)
}

# Alternative boxplot visualization
## bwplots
precip.shift06m.dmi <- combineAOPI(dmi, precip.shift06m, rt = "org")
precip.shift06m.dmi$month <- substr(precip.shift06m.dmi$ts, 6, 7)
precip.shift06m.dmi.melt <- melt(precip.shift06m.dmi, 
                                 id.vars = c("month", "TypeClass"),
                                 measure.vars = "P_RT_NRT")
bw.plot <- ggplot(data = precip.shift06m.dmi.melt, 
                  aes(y = value, x = month, fill = TypeClass, 
                      dodge = TypeClass))
plot.precip.shift06m.dmi.split.median.all.bw <- bw.plot + geom_boxplot() +
  xlab("Month") + ylab("Preciptiation")

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.dmi.split.median.all.bw.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.dmi.split.median.all.bw)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.dmi.split.median.all.bw)
}


# Create publication quality figures of correlations between dmi and 
# precipitation
years <- unique(substr(precip.shift06m.dmi$ts, 1, 4))
precip.shift06m.mat.ssn_kz03k01 <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.dmi[grep(i, precip.shift06m.dmi$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$ssn_kz03k01[j]
  }
}

dmi.reshape.shift06m.mat <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.dmi[grep(i, precip.shift06m.dmi$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$oni[j]
  }
}

plot.precip.shift06m.dmi.ssn_kz03k01 <- 
  corPlotBydmi(dmi.reshape.shift06m.mat, precip.shift06m.mat.ssn_kz03k01)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.oni.ssn_kz03k01.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.oni.ssn_kz03k01)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.dmi.ssn_kz03k01)
}


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

