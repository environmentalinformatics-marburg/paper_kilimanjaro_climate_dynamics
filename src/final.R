#### Set working environment ###################################################
switch(Sys.info()[["sysname"]], 
       "Linux" = setwd("/media/permanent/active/kilimanjaro"), 
       "Windows" = setwd("D:/active/kilimanjaro"))
#setwd("C:/Users/Dogbert/Desktop/kilimanjaro")) 

sourcePath <- "scripts/paper_kilimanjaro_climate_dynamics/src/"
dataPath <- "data/"
graphicsPath <- "graphics/"

printToFile = FALSE

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
source(paste0(sourcePath, "corPlotByAOI.R"))


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
precip$kz03k01 <- kz(precip$P_RT_NRT, m = 3, k = 1)

# Compute deseasoned precipitation time series and corresponding 3 month running
# mean using a Kolmogorov-Zurbenko filter with two iterations to close gaps
precip$ssn <- precip$P_RT_NRT - rep(sapply(1:12, function(i) {
  mean(precip$P_RT_NRT[seq(i, nrow(precip), 12)], na.rm = TRUE)
}), nrow(precip) / 12)
precip$ssn_kz03k01 <- kz(precip$ssn, m = 3, k = 1)
precip$ssn_kz03k02 <- kz(precip$ssn, m = 3, k = 2)

longTermDynamicsPlot(parameter = "precipitation", printToFile = printToFile)


#### Seasonal precipitation analysis ###########################################
# Split the original precipitation data set by the three main wet/dry phases
# identified in the long-term trend figure (see above)and create publication
# quality figure; to smooth the cycle while not reducing the rainfall amounts 
# significantly, use a spline prediction.
precip.seasonalwetdry <- seasonalMean(
  st = c(1975,1976,1992,2001), nd = c(2013,1992,2000,2013))

# Split the running mean data set by year and create publication quality figure;
# to smooth the cycle while not reducing the rainfall amounts significantly, 
# use a spline prediction.
precip.seasonalwetdry.split <- 
  split(precip.seasonalwetdry, precip.seasonalwetdry$season)

precip.seasonal.normal <- 
  list(lapply(precip.seasonalwetdry.split, function(x){x$p_dyn})$"1975-2013")

colors <- c("black", "blue3", "red", "cornflowerblue")

plot.precip.seasonalwetdry.all <- seasonPlotByAOI(
  lapply(precip.seasonalwetdry.split, function(x){x$p_dyn}), colors,
  linetype = c(2,1,1,1), ymin = 0, ymax = 250)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, "plot.precip.seasonalwetdry.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.seasonalwetdry.all)
  dev.off()
} else {
  plot(plot.precip.seasonalwetdry.all)
}


#### Precipitation analysis vs ENSO ############################################
# Prepare aoi record and classifiy years as La Nina (L), El Nino (E) or 
# normal (N); weak ENSO cycles are classified as normal
enso <- aoi.list$ONI
enso$TypeClass <- "N"
enso$TypeClass[grep("L", enso$Type)] <- "La Nina"
enso$TypeClass[grep("E", enso$Type)] <- "El Nino"

# Compute plot for long-term normal distribution
precip.shift06m <- precip[7:(nrow(precip)-6), ]
precip.shift06m.enso.split.median <- combineAOPI(enso, precip.shift06m)
yminmax = c(0, 250)
colors <- c("black")
plot.precip.shift06m.enso.split.median.normal <- 
  seasonPlotByAOI(precip.seasonal.normal, colors,
                  linetype = c(2), ymin = yminmax[1], ymax = yminmax[2])

# Compute seasonal distribution by major aoi situation
if(length(precip.shift06m.enso.split.median) == 5){
  colors <- c("blue", "lightblue", "red", "bisque", "black")
} else if(length(precip.shift06m.enso.split.median) == 1){
  colors <- c("black")
} else {
  colors <- c("blue", "red", "black")  
}
plot.precip.shift06m.enso.split.median.all <- 
  seasonPlotByAOI(precip.shift06m.enso.split.median, colors,
                  normal = plot.precip.shift06m.enso.split.median.normal,
                  ymin = yminmax[1], ymax = yminmax[2])
if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.enso.split.median.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.enso.split.median.all)
  dev.off()
} else {
  plot(plot.precip.shift06m.enso.split.median.all)
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
test <- precip.shift06m.enso[precip.shift06m.enso$TypeClass == "El Nino",]
corPlotByAOI(precip.shift06m.aoi = test, 
             plotFilePath = "plot.precip.shift06m.enso.ssn_kz03k01.tif",
             printToFile = printToFile)


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
corPlotByAOI(precip.shift06m.aoi = test, 
             ylable = "DMI",
             plotFilePath = "plot.precip.shift06m.iod.ssn_kz03k01.tif",
             printToFile = printToFile)


#### Precipitation analysis vs enid #############################################
# Prepare aoi record and classifiy years as enid plus (P), enid minus (M) or
# normal (N)
enidENSO <- aoi.list$ONI
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
corPlotByAOI(precip.shift06m.aoi = test, 
             ylable = "ENID",
             plotFilePath = "plot.precip.shift06m.enid.ssn_kz03k01.tif",
             printToFile = printToFile)


### Cloud EOT analysis vs ENSO ################################################
# Prepare aoi record and classifiy years as La Nina (L), El Nino (E) or 
# normal (N); weak ENSO cycles are classified as normal
enso <- aoi.list$DMIHAD
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
enso$TypeClass[grep("L", enso$Type)] <- "La Nina"
enso$TypeClass[grep("E", enso$Type)] <- "El Nino"

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
  test <- cloudEOT.shift06m.enso
  corPlotByAOI(precip.shift06m.aoi = test,
               parameter = "EOT",
               xlable = paste0("EOT", as.character(x)),
               ylable = "AOI",
               plotFilePath = paste0("plot.cloudEOT-", as.character(x), 
                                     ".shift06m.enso.ssn_kz03k01.tif"),
               printToFile = printToFile)
})




















# precip.seasonal11yr <- 
#   foreach(st = 1975:2003, nd = 1985:2013, .combine = "rbind") %do% {
#     st.process <- grep(st, precip[, 1])[7]
#     nd.process <- grep(nd, precip[, 1])[length(grep(nd, precip[, 1]))-6]
#     precip.process <- precip[st.process:nd.process, ]
#     data.frame(season = paste(st, nd, sep = "-"), 
#                p_dyn = vectorHarmonics(precip.process$kz03k01, 
#                                        st = c(st, 1), 
#                                        nd = c(nd, 12), 
#                                        m = 3))
#   }


# Compute seasonal precipitation dynamics for each year of the time series
# centered arround Dezember and create publication quality figure. 
# Compute a eleven year running mean of the seasonal precipitation dynamics
# precip.seasonal11yr <- seasonalMean(st = 1975:2003, nd = 1985:2013)
#
# Split the running mean data set by year and create publication quality figure;
# to smooth the cycle while not reducing the rainfall amounts significantly, 
# use a spline prediction.
# precip.seasonal11yr.split <- 
#   split(precip.seasonal11yr, precip.seasonal11yr$season)
# 
# colors <- colorRampPalette(
#   brewer.pal(11, "RdYlGn"))(length(precip.seasonal11yr.split))
# 
# plot.precip.seasonal11yr.all <- seasonPlotByAOI(
#   lapply(precip.seasonal11yr.split, function(x){x$p_dyn}), colors)
# 
# if(printToFile == TRUE){
#   tiff(filename = paste0(graphicsPath, "plot.precip.seasonal11yr.all.tif"),
#        width = 2480, height = 1748 , res = 300, pointsize =  12)
#   plot(plot.precip.seasonal11yr.all)
#   dev.off()
# } else {
#   plot(plot.precip.seasonal11yr.all)
# }

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

# precip.seasonalwetdry <- 
#   foreach(st = c(1975,1976,1992,2001), nd = c(2013,1992,2000,2013), .combine = "rbind") %do% {
#     st.process <- grep(st, precip[, 1])[7]
#     nd.process <- grep(nd, precip[, 1])[length(grep(nd, precip[, 1]))-6]
#     precip.process <- precip[st.process:nd.process, ]
#     data.frame(season = paste(st, nd, sep = "-"), 
#                p_dyn = vectorHarmonics(precip.process$P_RT_NRT, 
#                                        st = c(st, 1), 
#                                        nd = c(nd, 12), 
#                                        m = 3))
#   }



# Non-spline version
# plot.precip.shift06m.enso.split.median <- 
#   lapply(seq(precip.shift06m.enso.split.median), function(x){
#     xyplot(precip.shift06m.enso.split.median[[x]] ~ 
#              factor(c(7:12, 1:6), levels = c(7:12, 1:6)), type = "l", 
#            ylim = c(0, 200), col = colors[x])
#   })
# 
# plot.precip.shift06m.enso.split.median.all <- 
#   Reduce("outLayer", plot.precip.shift06m.enso.split.median)
# plot.precip.shift06m.enso.split.median.all <- 
#   Reduce("outLayer", plot.precip.shift06m.enso.split.median)

# Just for background information: the individual seasons which form the above
# combined El Nino, La Nina and normal situation spline.
# precip.shift06m.enso.split <- combineAOPI(enso, precip.shift06m, rt = "split")
# 
# precip.shift06m.enso.split.season <- lapply(precip.shift06m.enso.split, function(x){
#   split(x, x$Season)
# })
# 
# plot.precip.shift06m.enso.split.season.all <- 
#   seasonPlotByAOI(precip.shift06m.enso.split.season, colors, individual = TRUE)
# 
# if(printToFile == TRUE){
#   tiff(filename = paste0(graphicsPath, 
#                          "plot.precip.shift06m.enso.split.season.all.tif"),
#        width = 2480, height = 1748 , res = 300, pointsize =  12)
#   plot(plot.precip.shift06m.enso.split.season.all)
#   dev.off()
#   
# } else {
#   plot(plot.precip.shift06m.enso.split.season.all)
# }


# Non-spline version
# plot.precip.shift06m.iod.split.median <- 
#   lapply(seq(precip.shift06m.iod.split.median), function(x){
#     xyplot(precip.shift06m.iod.split.median[[x]] ~ 
#              factor(c(7:12, 1:6), levels = c(7:12, 1:6)), type = "l", 
#            ylim = c(0, 200), col = colors[x])
#   })
# 
# plot.precip.shift06m.iod.split.median.all <- 
#   Reduce("outLayer", plot.precip.shift06m.iod.split.median)
# plot.precip.shift06m.iod.split.median.all <- 
#   Reduce("outLayer", plot.precip.shift06m.iod.split.median)

# Just for background information: the individual seasons which form the above
# combined El Nino, La Nina and normal situation spline.
# precip.shift06m.iod.split <- combineAOPI(iod, precip.shift06m, rt = "split")
# 
# precip.shift06m.iod.split.season <- lapply(precip.shift06m.iod.split, function(x){
#   split(x, x$Season)
# })
# 
# plot.precip.shift06m.iod.split.season.all <- 
#   seasonPlotByAOI(precip.shift06m.iod.split.season, colors, individual = TRUE)
# 
# if(printToFile == TRUE){
#   tiff(filename = paste0(graphicsPath, 
#                          "plot.precip.shift06m.iod.split.season.all.tif"),
#        width = 2480, height = 1748 , res = 300, pointsize =  12)
#   plot(plot.precip.shift06m.iod.split.season.all)
#   dev.off()
#   
# } else {
#   plot(plot.precip.shift06m.iod.split.season.all)
# }



# Non-spline version
# plot.cloudEOT.shift06m.enso.split.median <- 
#   lapply(seq(cloudEOT.shift06m.enso.split.median), function(x){
#     xyplot(cloudEOT.shift06m.enso.split.median[[x]] ~ 
#              factor(c(7:12, 1:6), levels = c(7:12, 1:6)), type = "l", 
#            ylim = c(0, 200), col = colors[x])
#   })
# 
# plot.cloudEOT.shift06m.enso.split.median.all <- 
#   Reduce("outLayer", plot.cloudEOT.shift06m.enso.split.median)
# plot.cloudEOT.shift06m.enso.split.median.all <- 
#   Reduce("outLayer", plot.cloudEOT.shift06m.enso.split.median)

# Just for background information: the individual seasons which form the above
# combined El Nino, La Nina and normal situation spline.
# cloudEOT.shift06m.enso.split <- combineAOPI(enso, cloudEOT.shift06m, rt = "split")
# 
# cloudEOT.shift06m.enso.split.season <- lapply(cloudEOT.shift06m.enso.split, function(x){
#   split(x, x$Season)
# })
# 
# plot.cloudEOT.shift06m.enso.split.season.all <- 
#   seasonPlotByAOI(cloudEOT.shift06m.enso.split.season, colors, individual = TRUE,
#                   ymin = yminmax[1], ymax = yminmax[2], ylab = "Cloud cover")
# 
# if(printToFile == TRUE){
#   tiff(filename = paste0(graphicsPath, 
#                          "plot.cloudEOT.shift06m.enso.split.season.all.tif"),
#        width = 2480, height = 1748 , res = 300, pointsize =  12)
#   plot(plot.cloudEOT.shift06m.enso.split.season.all)
#   dev.off()
#   
# } else {
#   plot(plot.cloudEOT.shift06m.enso.split.season.all)
# }
# 
