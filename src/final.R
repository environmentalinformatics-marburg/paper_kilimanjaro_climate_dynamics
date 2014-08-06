#### Set working environment ###################################################
switch(Sys.info()[["sysname"]], 
       "Linux" = setwd("/media/permanent/active/kilimanjaro"), 
       "Windows" = setwd("D:/active/kilimanjaro"))
#setwd("C:/Users/Dogbert/Desktop/kilimanjaro")) 

sourcePath <- "scripts/paper_kilimanjaro_climate_dynamics/src/"
dataPath <- "data/"
graphicsPath <- "graphics/"

printToFile <- TRUE
plot2file <- printToFile

library(kza)
library(latticeExtra)
library(Kendall)
library(foreach)
library(ade4)
library(reshape2)
library(ggplot2)
source(paste0(sourcePath, "readData.R"))
source(paste0(sourcePath, "longTermDynamicsPlot.R"))
source(paste0(sourcePath, "seasonalMean.R"))
source(paste0(sourcePath, "vectorHarmonics.R"))
source(paste0(sourcePath, "combineAOPI.R"))
source(paste0(sourcePath, "seasonPlotByAOI.R"))
source(paste0(sourcePath, "visCorPlotTimeSeries.R"))
source(paste0(sourcePath, "mergeAOIwTS.R"))
source(paste0(sourcePath, "visSeasonPlotByAOI.R"))

#### Functions #################################################################
outLayer <- function(x, y) {
  x + as.layer(y)
}


#### Read data sets ############################################################
# Read temperature data and aggregate it to monthly values
ta.list <- readData("temperature")

# Read precipitation data and create a continous time series.
precip.list <- readData("precipitation")

# Read EOT
cloudEOT.list <- readData("cloudEOT")

# Read AOI
aoi.list <- readData("aoi")


#### Long-term temperature analysis ############################################
# Compute long-term anomalies and trends for KIA temperature and create
# publication quality figure.
# Compute 3 month running mean of original temperature values using a
# Kolmogorov-Zurbenko filter with one iteration
ta <- ta.list$KIA
longTermDynamicsPlot(parameter = "temperature", printToFile = printToFile)


#### Long-term precipitation analysis ##########################################
# Compute long-term anomalies and trends for KIA precipitation and create
# publication quality figure.
# Compute 3 month running mean of original precipitation values using a
# Kolmogorov-Zurbenko filter with one iteration
precip <- precip.list$KIA
precip[1:6, 2] <- NA
precip$kz03k01 <- kz(precip$P_RT_NRT, m = 3, k = 1)

# Compute deseasoned precipitation time series and corresponding 3 month running
# mean using a Kolmogorov-Zurbenko filter with two iterations to close gaps
precip$ssn <- precip$P_RT_NRT - rep(sapply(1:12, function(i) {
    mean(precip$P_RT_NRT[seq(i, nrow(precip), 12)], na.rm = TRUE)
  }), nrow(precip) / 12)

precip$ssnmed <- precip$P_RT_NRT - rep(sapply(1:12, function(i) {
  median(precip$P_RT_NRT[seq(i, nrow(precip), 12)], na.rm = TRUE)
}), nrow(precip) / 12)

precip$ssnpmed <- (precip$P_RT_NRT / rep(sapply(1:12, function(i) {
  median(precip$P_RT_NRT[seq(i, nrow(precip), 12)], na.rm = TRUE)
}), nrow(precip) / 12)) - 1.0
precip$ssnpmed[!is.finite(precip$ssnpmed)] <- 0.0

precip$ssnpmean <- (precip$P_RT_NRT / rep(sapply(1:12, function(i) {
  mean(precip$P_RT_NRT[seq(i, nrow(precip), 12)], na.rm = TRUE)
}), nrow(precip) / 12)) - 1.0

precip$ssn_kz03k01 <- kz(precip$ssn, m = 3, k = 1)
precip$ssn_kz03k02 <- kz(precip$ssn, m = 3, k = 2)

precip$ssnmed_kz03k01 <- kz(precip$ssnmed, m = 3, k = 1)
precip$ssnmed_kz03k02 <- kz(precip$ssnmed, m = 3, k = 2)

precip$ssnpmed_kz03k01 <- kz(precip$ssnpmed, m = 3, k = 1)
precip$ssnpmed_kz03k02 <- kz(precip$ssnpmed, m = 3, k = 2)

precip$ssnpmean_kz03k01 <- kz(precip$ssnpmean, m = 3, k = 1)
precip$ssnpmean_kz03k02 <- kz(precip$ssnpmean, m = 3, k = 2)

longTermDynamicsPlot(parameter = "precipitation", printToFile = printToFile)


#### Seasonal precipitation analysis ###########################################
# Split the original precipitation data set by the three main wet/dry phases
# identified in the long-term trend figure (see above)and create publication
# quality figure; to smooth the cycle while not reducing the rainfall amounts 
# significantly, use a spline prediction.
# Split the running mean data set by year and create publication quality figure;
# to smooth the cycle while not reducing the rainfall amounts significantly, 
# use a spline prediction.
colors <- c("black", "blue3", "red", "cornflowerblue")
yminmax<- c(0, 250)
#yminmax<- c(-150, 150)

# 12 month season
precip.seasonalwetdry <- seasonalMean(
  st = c(1975,1976,1992,2001), nd = c(2013,1992,2000,2013),
  st.shift = 7, nd.shift = 0, timespan = 12, fun = "median", prm = "P_RT_NRT")

precip.seasonalwetdry.split <- 
  split(precip.seasonalwetdry, precip.seasonalwetdry$season)
precip.seasonal.normal <- 
  list(lapply(precip.seasonalwetdry.split, function(x){x$p_dyn})$"1975-2013")

plot.precip.seasonalwetdry.all <- seasonPlotByAOI(
  lapply(precip.seasonalwetdry.split, function(x){x$p_dyn}), colors,
  linetype = c(2,1,1,1), ymin = yminmax[1], ymax = yminmax[2])

# 18 month season
precip.18m.seasonalwetdry <- seasonalMean(
  st = c(1975,1976,1992,2001), nd = c(2013,1992,2000,2013), 
  st.shift = 7, nd.shift = 0, timespan = 18, fun = "median", prm = "P_RT_NRT")

precip.18m.seasonalwetdry.split <- 
  split(precip.18m.seasonalwetdry, precip.18m.seasonalwetdry$season)
precip.18m.seasonal.normal <- 
  list(lapply(precip.18m.seasonalwetdry.split, function(x){x$p_dyn})$"1975-2013")

plot.precip18m.seasonalwetdry.all <- 
  visSeasonPlotByAOI(
    lapply(precip.18m.seasonalwetdry.split, function(x){x$p_dyn}), colors,
    linetype = c(2,1,1,1), ymin = yminmax[1], ymax = yminmax[2], timespan = 18)

# 24 month season
precip.24m.seasonalwetdry <- seasonalMean(
  st = c(1975,1976,1992,2001), nd = c(2013,1992,2000,2013), 
  st.shift = 7, nd.shift = 0, timespan = 24, fun = "median", prm = "P_RT_NRT")

precip.24m.seasonalwetdry.split <- 
  split(precip.24m.seasonalwetdry, precip.24m.seasonalwetdry$season)
precip.24m.seasonal.normal <- 
  list(lapply(precip.24m.seasonalwetdry.split, function(x){x$p_dyn})$"1975-2013")

plot.precip24m.seasonalwetdry.all <- 
  visSeasonPlotByAOI(
  lapply(precip.24m.seasonalwetdry.split, function(x){x$p_dyn}), colors,
  linetype = c(2,1,1,1), ymin = yminmax[1], ymax = yminmax[2], timespan = 24)


if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, "plot.precip.seasonalwetdry.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.seasonalwetdry.all)
  dev.off()
  tiff(filename = paste0(graphicsPath, "plot.precip18m.seasonalwetdry.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip18m.seasonalwetdry.all)
  dev.off()
  tiff(filename = paste0(graphicsPath, "plot.precip24m.seasonalwetdry.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip24m.seasonalwetdry.all)
  dev.off()
} else {
  plot(plot.precip.seasonalwetdry.all)
  plot(plot.precip18m.seasonalwetdry.all)
  plot(plot.precip24m.seasonalwetdry.all)
}


#### Precipitation analysis vs ENSO ############################################
# Prepare aoi record and classifiy years as La Nina (L), El Nino (E) or 
# normal (N); weak ENSO cycles are classified as normal
enso <- aoi.list$MEI
enso$TypeClass <- "N"
enso$TypeClass[grep("E", enso$Type)] <- "El Nino"
enso$TypeClass[grep("WE", enso$Type)] <- "El Nino W"
enso$TypeClass[grep("L", enso$Type)] <- "La Nina"
enso$TypeClass[grep("WL", enso$Type)] <- "La Nina W"
#enso$TypeClass[grep("W", enso$Type)] <- "N"
enso$TypeClass[grep("P", enso$IOD)] <- "N"
enso$TypeClass[grep("M", enso$IOD)] <- "N"

# Compute plot for long-term normal distribution
yminmax = c(0, 350)
#yminmax = c(-200,200)
colors <- c("black")

precip.shift06m <- precip[7:(nrow(precip)-6), ]
precip.shift06m.enso.split.median <- combineAOPI(enso, precip.shift06m)

plot.precip.shift06m.enso.split.median.normal <- 
  seasonPlotByAOI(precip.seasonal.normal, colors,
                  linetype = c(2), ymin = yminmax[1], ymax = yminmax[2])


precip.18m <- precip[7:(nrow(precip)-0), ]
precip.18m.enso.split.median <- mergeAOIwTS(enso, precip.18m, 
                                            timespan = 18,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")

plot.precip.18m.enso.split.median.normal <- 
  visSeasonPlotByAOI(precip.18m.seasonal.normal, colors,
                     linetype = c(2), ymin = yminmax[1], ymax = yminmax[2])


precip.24m <- precip[7:(nrow(precip)-0), ]
precip.24m.enso.split.median <- mergeAOIwTS(enso, precip.24m, 
                                            timespan = 24,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")

plot.precip.24m.enso.split.median.normal <- 
  visSeasonPlotByAOI(precip.24m.seasonal.normal, colors,
                     linetype = c(2), ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 24)


# Compute seasonal distribution by major aoi situation
red <- brewer.pal(4, "Reds")
blue <- brewer.pal(4, "Blues")
colors <- c(blue[4], blue[2], red[4], red[2], "black")

plot.precip.shift06m.enso.split.median.all <- 
  seasonPlotByAOI(precip.shift06m.enso.split.median, colors,
                  linetype = c(1,2,1,2,1),
                  normal = plot.precip.shift06m.enso.split.median.normal,
                  ymin = yminmax[1], ymax = yminmax[2])

plot.precip.18m.enso.split.median.all <- 
  visSeasonPlotByAOI(precip.18m.enso.split.median, colors,
                  linetype = c(1,2,1,2,1),
                  normal = plot.precip.18m.enso.split.median.normal,
                  ymin = yminmax[1], ymax = yminmax[2],
                  vline.pos = 501)

plot.precip.24m.enso.split.median.all <- 
  visSeasonPlotByAOI(precip.24m.enso.split.median, colors,
                     linetype = c(1,2,1,2,1),
                     normal = plot.precip.24m.enso.split.median.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 24,
                     vline.pos = 501)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.enso.split.median.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.enso.split.median.all)
  dev.off()
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.18m.enso.split.median.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.18m.enso.split.median.all)
  dev.off()
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.24m.enso.split.median.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.24m.enso.split.median.all)
  dev.off()
} else {
  plot(plot.precip.shift06m.enso.split.median.all)
  plot(plot.precip.18m.enso.split.median.all)
  plot(plot.precip.24m.enso.split.median.all)
}

# Alternative boxplot visualization
precip.shift06m.enso <- combineAOPI(enso, precip.shift06m, rt = "org")
precip.shift06m.enso$month <- substr(precip.shift06m.enso$ts, 6, 7)
precip.shift06m.enso.melt <- melt(precip.shift06m.enso, 
                                  id.vars = c("month", "TypeClass"),
                                  measure.vars = "P_RT_NRT")
bw.plot <- ggplot(data = precip.shift06m.enso.melt, 
                  aes(y = value, x = month, fill = TypeClass, 
                      dodge = TypeClass))
plot.precip.shift06m.enso.split.median.all.bw <- bw.plot + geom_boxplot() +
  xlab("Month") + ylab("Preciptiation")

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.enso.split.median.all.bw.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.enso.split.median.all.bw)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.enso.split.median.all.bw)
}

# Create publication quality figures of correlations between enso and 
# precipitation
#test <- precip.shift06m.enso[precip.shift06m.enso$TypeClass == "N",]
precip.shift06m.enso <- combineAOPI(enso, precip.shift06m, rt = "org")
precip.shift06m.enso$StartSeason <- substr(precip.shift06m.enso$Season,1,4)
test <- precip.shift06m.enso
m12 <- visCorPlotTimeSeries(df = test, 
                     x.prm = "aoi",
                     y.prm = "P_RT_NRT",
                     t.prm = "StartSeason",
                     lable.nbrs = c(c(7:12),c(1:6)),
                     p.thv = 0.08,
                     plot.filepath = "plot.precip.shift06m.enso.cor.tif",
                     plot2file = plot2file,
                     rt = TRUE)

precip.18m.enso <- mergeAOIwTS(enso, precip.18m, 
                               timespan = 18, rt = "org")
test <- precip.18m.enso # [precip.18m.enso$TypeClass == "El Nino", ]
m18 <- visCorPlotTimeSeries(df = test, 
                     x.prm = "aoi",
                     y.prm = "P_RT_NRT",
                     t.prm = "StartSeason",
                     lable.nbrs = c(c(7:12),c(1:12)),
                     p.thv = 0.08,
                     plot.filepath = "plot.precip18m.enso.cor.tif",
                     plot2file = plot2file,
                     rt = TRUE)

precip.24m.enso <- mergeAOIwTS(enso, precip.24m, 
                               timespan = 24, rt = "org")
test <- precip.24m.enso # [precip.24m.enso$TypeClass == "El Nino",]
m24 <- visCorPlotTimeSeries(df = test, 
                     x.prm = "aoi",
                     y.prm = "P_RT_NRT",
                     t.prm = "StartSeason",
                     lable.nbrs = c(c(7:12),c(1:12),c(1:6)),
                     p.thv = 0.08,
                     plot.filepath = "plot.precip24m.enso.cor.tif",
                     plot2file = plot2file,
                     rt = TRUE)

testm12m18 <- m12 - m18[1:12,1:12]
testm18m24 <- m18 - m24[1:18,1:18]
range(testm12m18, na.rm = TRUE)
range(testm18m24, na.rm = TRUE)

#### Precipitation analysis vs IOD #############################################
# Prepare aoi record and classifiy years as IOD plus (P), IOD minus (M) or
# normal (N)
iod <- aoi.list$DMIHAD
iod$TypeClass <- "X"
iod$TypeClass[grep("P", iod$IOD)] <- "C1 P"
iod$TypeClass[grep("M", iod$IOD)] <- "C2 M"

# Compute plot for long-term normal distribution
precip.shift06m <- precip[7:(nrow(precip)-6), ]
precip.shift06m.iod.split.median <- combineAOPI(iod, precip.shift06m)
yminmax = c(0, 250)
colors <- c("black")
plot.precip.shift06m.iod.split.median.normal <- 
  seasonPlotByAOI(precip.seasonal.normal, colors,
                  linetype = c(2), ymin = yminmax[1], ymax = yminmax[2])

# Compute seasonal distribution by major aoi situation
if(length(precip.shift06m.iod.split.median) == 5){
  colors <- c("blue", "lightblue", "red", "bisque", "black")
} else if(length(precip.shift06m.iod.split.median) == 1){
  colors <- c("black")
} else {
  colors <- c("blue", "red", "black")  
}
colors <- c("blue", "red")  
plot.precip.shift06m.iod.split.median.all <- 
  seasonPlotByAOI(precip.shift06m.iod.split.median, colors,
                  normal = plot.precip.shift06m.iod.split.median.normal,
                  ymin = yminmax[1], ymax = yminmax[2])
if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.iod.split.median.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.iod.split.median.all)
  dev.off()
} else {
  plot(plot.precip.shift06m.iod.split.median.all)
}

# Alternative boxplot visualization
precip.shift06m.iod <- combineAOPI(iod, precip.shift06m, rt = "org")
precip.shift06m.iod$month <- substr(precip.shift06m.iod$ts, 6, 7)
precip.shift06m.iod.melt <- melt(precip.shift06m.iod, 
                                 id.vars = c("month", "TypeClass"),
                                 measure.vars = "P_RT_NRT")
bw.plot <- ggplot(data = precip.shift06m.iod.melt, 
                  aes(y = value, x = month, fill = TypeClass, 
                      dodge = TypeClass))
plot.precip.shift06m.iod.split.median.all.bw <- bw.plot + geom_boxplot() +
  xlab("Month") + ylab("Preciptiation")

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.iod.split.median.all.bw.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.iod.split.median.all.bw)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.iod.split.median.all.bw)
}

# Create publication quality figures of correlations between iod and 
# precipitation
test <- precip.shift06m.iod
visCorPlotTimeSeries(df = test, 
                     x.prm = "aoi",
                     y.prm = "ssn_kz03k01",
                     t.prm = "ts",
                     x.lable = "DMI",
                     lable.nbrs = c(c(7:12),c(1:6)),
                     p.thv = 0.05,
                     plot.filepath = "plot.precip.shift06m.iod.ssn_kz03k01.tif",
                     plot2file = plot2file)



#### Precipitation analysis vs enid #############################################
# Prepare aoi record and classifiy years as enid plus (P), enid minus (M) or
# normal (N)
enidENSO <- aoi.list$MEI
enidDMI <- aoi.list$DMIHAD
st <- enidDMI$Season[1]
nd <- enidDMI$Season[nrow(enidDMI)]
fr <- grep(st, enidENSO$Season)
lr <- grep(nd, enidENSO$Season)
enidENSO <- enidENSO[fr:lr,]
enidENSO[4:15] <- enidENSO[4:15] / max(abs(enidENSO[4:15]))
enidDMI[4:15] <- enidDMI[4:15]/max(abs(enidDMI[4:15]))
enid <- cbind(enidENSO[,1:3], enidDMI[,4:15] + enidENSO[,4:15])
colnames(dmi)[4:15] <- colnames(oni)[4:15]

enid$TypeClass <- "N"
enid$TypeClass[grep("P", enid$enid)] <- "I1P"
enid$TypeClass[grep("M", enid$enid)] <- "I2M"
enid$TypeClass[grep("L", enid$Type)] <- "LS"
enid$TypeClass[grep("WL", enid$Type)] <- "LW"
enid$TypeClass[grep("E", enid$Type)] <- "ES"
enid$TypeClass[grep("WE", enid$Type)] <- "EW"

precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.enid.split.median <- combineAOPI(enid, precip.shift06m)

precip.shift06m.enid <- combineAOPI(enid, precip.shift06m, rt = "org")

# Create publication quality figures of correlations between enid and 
# precipitation
test <- precip.shift06m.enid
visCorPlotTimeSeries(df = test, 
                     x.prm = "aoi",
                     y.prm = "ssn_kz03k01",
                     t.prm = "ts",
                     x.lable = "ENID",
                     lable.nbrs = c(c(7:12),c(1:6)),
                     p.thv = 0.05,
                     plot.filepath = "plot.precip.shift06m.enid.ssn_kz03k01.tif",
                     plot2file = plot2file)


### Cloud EOT analysis vs ENSO ################################################
# Prepare aoi record and classifiy years as La Nina (L), El Nino (E) or 
# normal (N); weak ENSO cycles are classified as normal
enso <- aoi.list$DMI

enidENSO <- aoi.list$ONI
enidDMI <- aoi.list$DMI
st <- enidDMI$Season[1]
nd <- enidDMI$Season[nrow(enidDMI)]
fr <- grep(st, enidENSO$Season)
lr <- grep(nd, enidENSO$Season)
enidENSO <- enidENSO[fr:lr,]
enidENSO[4:15] <- enidENSO[4:15] / max(abs(enidENSO[4:15]))
enidDMI[4:15] <- enidDMI[4:15] / max(abs(enidDMI[4:15]))
enid <- cbind(enidENSO[,1:3], enidDMI[,4:15] + enidENSO[,4:15])
colnames(enidDMI)[4:15] <- colnames(enidENSO)[4:15]
enso <- enid

cloudEOT <- cloudEOT.list$EOTssn
cloudEOT.shift06m <- cloudEOT[7:(nrow(cloudEOT)-6), ]

# Compute plot for long-term normal distribution
enso$TypeClass <- "N"
yminmax = c(-0.5, 0.5)
plot.cloudEOT.seasonal.normal.all <- 
  lapply(seq(3), function(x){
    cloudEOT.shift06m <- cloudEOT[7:(nrow(cloudEOT)-6), ]
    colnames(cloudEOT.shift06m)[x+1] <-  "P_RT_NRT"
    cloudEOT.shift06m.enso.split.median <- combineAOPI(enso, cloudEOT.shift06m)
    colors <- c("black")
    plot.cloudEOT.shift06m.enso.split.median.normal <- 
      seasonPlotByAOI(cloudEOT.shift06m.enso.split.median, colors,
                      linetype = c(2), 
                      ymin = yminmax[1], ymax = yminmax[2], ylab = "Cloud cover")
  })

# Compute seasonal distribution by major aoi situation
enso$TypeClass <- "N"
enso$TypeClass[grep("L", enso$Type)] <- "C2 La Nina"
enso$TypeClass[grep("E", enso$Type)] <- "C1 El Nino"
#enso$TypeClass[grep("W", enso$Type)] <- "X"
#enso$TypeClass <- "N"
#enso$TypeClass[grep("P", enso$IOD)] <- "C3 P"
#enso$TypeClass[grep("M", enso$IOD)] <- "c4 M"

plot.cloudEOT.shift06m.enso.split.median.all <- 
  lapply(seq(3), function(x){
    plot.cloudEOT.seasonal.normal <- plot.cloudEOT.seasonal.normal.all[[x]]
    cloudEOT.shift06m <- cloudEOT[7:(nrow(cloudEOT)-6), ]
    colnames(cloudEOT.shift06m)[x+1] <- "P_RT_NRT"
    cloudEOT.shift06m.enso.split.median <- combineAOPI(enso, cloudEOT.shift06m)
    if(length(cloudEOT.shift06m.enso.split.median) == 5){
      colors <- c("blue", "lightblue", "red", "bisque", "black")
    } else if(length(cloudEOT.shift06m.enso.split.median) == 1){
      colors <- c("black")
    } else {
      colors <- c("blue", "red", "black")  
    }
    plot.cloudEOT.shift06m.enso.split.median.all <- 
      seasonPlotByAOI(cloudEOT.shift06m.enso.split.median, colors,
                      normal = plot.cloudEOT.seasonal.normal,
                      ymin = yminmax[1], ymax = yminmax[2], ylab = "Cloud cover")
    if(printToFile == TRUE){
      tiff(filename = paste0(graphicsPath, 
                             paste0("plot.cloudEOT-", as.character(x), 
                                    ".shift06m.enso.split.median.all.tif")),
           width = 2480, height = 1748 , res = 300, pointsize =  12)
      plot(plot.cloudEOT.shift06m.enso.split.median.all)
      dev.off()
    } else {
      plot(plot.cloudEOT.shift06m.enso.split.median.all)
    }
  })

# Alternative boxplot visualization
## bwplots
plot.cloudEOT.shift06m.enso.split.median.all.bw <- 
  lapply(seq(3), function(x){
    act.cloudEOT.shift06m <- cloudEOT.shift06m[,c(1,x+1)]
    colnames(act.cloudEOT.shift06m) <- c("ts", "P_RT_NRT")
    cloudEOT.shift06m.enso <- combineAOPI(enso, act.cloudEOT.shift06m, rt = "org")
    cloudEOT.shift06m.enso$month <- substr(cloudEOT.shift06m.enso$ts, 6, 7)
    cloudEOT.shift06m.enso.melt <- melt(cloudEOT.shift06m.enso, 
                                        id.vars = c("month", "TypeClass"),
                                        measure.vars = "P_RT_NRT")
    bw.plot <- ggplot(data = cloudEOT.shift06m.enso.melt, 
                      aes(y = value, x = month, fill = TypeClass, 
                          dodge = TypeClass))
    plot.cloudEOT.shift06m.enso.split.median.all.bw <- bw.plot + geom_boxplot() +
      xlab("Month") + ylab("cloudEOTtiation")
    
    if(printToFile == TRUE){
      tiff(filename = paste0(graphicsPath, 
                             paste0("plot.cloudEOT-", as.character(x), 
                                    ".shift06m.enso.split.median.all.bw.tif")),
           width = 2480, height = 1748 , res = 300, pointsize =  12)
      plot(plot.cloudEOT.shift06m.enso.split.median.all.bw)
      dev.off()
      
    } else {
      plot(plot.cloudEOT.shift06m.enso.split.median.all.bw)
    }
  })

# Create publication quality figures of correlations between enso and 
# cloudEOTitation
plot.cloudEOT.shift06m.enso.ssn_kz03k01 <- lapply(seq(3), function(x){
  act.cloudEOT.shift06m <- cloudEOT.shift06m[,c(1,x+1)]
  colnames(act.cloudEOT.shift06m) <- c("ts", "EOT")
  cloudEOT.shift06m.enso <- combineAOPI(enso, act.cloudEOT.shift06m, 
                                        parameter = "EOT", rt = "org")
  test <- cloudEOT.shift06m.enso # [cloudEOT.shift06m.enso$TypeClass == "La Nina", ]
  visCorPlotTimeSeries(df = test, 
                       x.prm = "aoi",
                       y.prm = "EOT",
                       t.prm = "ts",
                       lable.nbrs = c(c(7:12),c(1:6)),
                       p.thv = 0.05,
                       plot.filepath = paste0("plot.cloudEOT-", as.character(x), 
                                              ".shift06m.enso.ssn_kz03k01.tif"),
                       plot2file = plot2file)
})