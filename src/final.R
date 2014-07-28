#### Set working environment ###################################################
switch(Sys.info()[["sysname"]], 
       "Linux" = setwd("/media/permanent/active/kilimanjaro"), 
       "Windows" = setwd("D:/active/kilimanjaro"))
#setwd("C:/Users/Dogbert/Desktop/kilimanjaro")) 

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

# Read EOT
cloudEOT <- read.csv(paste0(dataPath, "cloud_cover_monthly_myd06/eot/eot_series.csv"), 
                   header = TRUE, sep = ",")
cloudEOT$Years <- rep(2003:2013, each=12)
cloudEOT$Month <- rep(1:12)
cloudEOT$ts <- as.Date(paste0(cloudEOT$Year,"-",cloudEOT$Month,"-1"))
cloudEOT <- cloudEOT[,c(6,1:3)]

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
dmi <- read.csv(paste0(dataPath, "aoi/dmi_in.txt"), header = FALSE, sep = ":")
colnames(dmi) <- c("Date", "DMI")
dmi$YearMonth <- substr(dmi$Date,1,7)
dmi <- aggregate(dmi$DMI, by=list(dmi$YearMonth), FUN = "mean")
colnames(dmi) <- c("YearMonth", "DMI")
dmi <- dmi[c((grep("1981",dmi$YearMonth)[-1]+1):(grep("2014-07",
                                                      dmi$YearMonth)[1]-1)),]
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
enso <- aoi.list$ONI
enso$TypeClass <- "N"

precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.enso.split.median <- combineAOPI(enso, precip.shift06m)

colors <- c("orange")

plot.precip.shift06m.enso.split.median.normal <- 
  seasonPlotByAOI(precip.shift06m.enso.split.median, colors,
                  ymin = 0, ymax = 300)

enso$TypeClass <- "N"
enso$TypeClass[grep("L", enso$Type)] <- "LS"
enso$TypeClass[grep("WL", enso$Type)] <- "LW"
enso$TypeClass[grep("E", enso$Type)] <- "ES"
enso$TypeClass[grep("WE", enso$Type)] <- "EW"
# enso$TypeClass[grep("W", enso$Type)] <- "N"
# enso$TypeClass[grep("M", enso$Type)] <- "X"

precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.enso.split.median <- combineAOPI(enso, precip.shift06m)

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
                  ymin = 0, ymax = 300)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.enso.split.median.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.enso.split.median.all)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.enso.split.median.all)
}

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

precip.shift06m.enso.split <- combineAOPI(enso, precip.shift06m, rt = "split")

precip.shift06m.enso.split.season <- lapply(precip.shift06m.enso.split, function(x){
  split(x, x$Season)
})

plot.precip.shift06m.enso.split.season.all <- 
  seasonPlotByAOI(precip.shift06m.enso.split.season, colors, individual = TRUE)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.enso.split.season.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.enso.split.season.all)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.enso.split.season.all)
}



# Alternative boxplot visualization
## bwplots
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
years <- unique(substr(precip$ts, 1, 4))
precip.shift06m.mat.ssn_kz03k01 <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.enso[grep(i, precip.shift06m.enso$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$ssn_kz03k01[j]
  }
}

enso.reshape.shift06m.mat <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.enso[grep(i, precip.shift06m.enso$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$aoi[j]
  }
}

plot.precip.shift06m.enso.ssn_kz03k01 <- 
  corPlotByAOI(enso.reshape.shift06m.mat, precip.shift06m.mat.ssn_kz03k01)


if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.enso.ssn_kz03k01.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.enso.ssn_kz03k01)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.enso.ssn_kz03k01)
}



#### Precipitation analysis vs IOD #############################################
# Prepare aoi record and classifiy years as IOD plus (P), IOD minus (M) or
# normal (N)
iod <- aoi.list$DMI

iod$TypeClass <- "N"

precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.iod.split.median <- combineAOPI(iod, precip.shift06m)

colors <- c("orange")

plot.precip.shift06m.iod.split.median.normal <- 
  seasonPlotByAOI(precip.shift06m.iod.split.median, colors,
                  ymin = 0, ymax = 300)

iod$TypeClass <- "N"
iod$TypeClass[grep("P", iod$IOD)] <- "I1P"
iod$TypeClass[grep("M", iod$IOD)] <- "I2M"
#iod$TypeClass[grep("E", iod$Type)] <- "X"
#iod$TypeClass[grep("L", iod$Type)] <- "X"

precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.iod.split.median <- combineAOPI(iod, precip.shift06m)

if(length(precip.shift06m.iod.split.median) == 5){
  colors <- c("blue", "lightblue", "red", "bisque", "black")
} else if (length(precip.shift06m.iod.split.median) == 4){
  colors <- c("blue", "red", "black", "bisque")
} else {
  colors <- c("blue", "red", "black")  
}

plot.precip.shift06m.iod.split.median <- 
  seasonPlotByAOI(precip.shift06m.iod.split.median, colors,
                  normal = plot.precip.shift06m.iod.split.median.normal,
                  ymin = 0, ymax = 300)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.iod.split.median.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.iod.split.median)
  dev.off()
} else {
  plot(plot.precip.shift06m.iod.split.median)
}

# Just for background information: the individual seasons which form the above
# combined El Nino, La Nina and normal situation spline.
precip.shift06m.iod.split <- combineAOPI(iod, precip.shift06m, rt = "split")

precip.shift06m.iod.split.season <- lapply(precip.shift06m.iod.split, function(x){
  split(x, x$Season)
})

plot.precip.shift06m.iod.split.season.all <- 
  seasonPlotByAOI(precip.shift06m.iod.split.season, colors, individual = TRUE)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.iod.split.season.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.iod.split.season.all)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.iod.split.season.all)
}

# Alternative boxplot visualization
## bwplots
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
years <- unique(substr(precip.shift06m.iod$ts, 1, 4))
precip.shift06m.mat.ssn_kz03k01 <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.iod[grep(i, precip.shift06m.iod$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$ssn_kz03k01[j]
  }
}

iod.reshape.shift06m.mat <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.iod[grep(i, precip.shift06m.iod$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$aoi[j]
  }
}

plot.precip.shift06m.iod.ssn_kz03k01 <- 
  corPlotByAOI(iod.reshape.shift06m.mat, precip.shift06m.mat.ssn_kz03k01)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.iod.ssn_kz03k01.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.iod.ssn_kz03k01)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.iod.ssn_kz03k01)
}




#### Precipitation analysis vs enid #############################################
# Prepare aoi record and classifiy years as enid plus (P), enid minus (M) or
# normal (N)
enidDMI <- aoi.list$DMI
enidENSO <- aoi.list$ONI
fr <- grep(enidDMI$Season[1], enidENSO$Season)
lr <- grep(enidDMI$Season[nrow(enidDMI)], enidENSO$Season)
enid <- cbind(enidENSO[fr:lr,1:3],enidDMI[,4:15] + enidENSO[fr:lr,4:15])
colnames(dmi)[4:15] <- colnames(oni)[4:15]

enid$TypeClass <- "N"

precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.enid.split.median <- combineAOPI(enid, precip.shift06m)

colors <- c("orange")

plot.precip.shift06m.enid.split.median.normal <- 
  seasonPlotByAOI(precip.shift06m.enid.split.median, colors,
                  ymin = 0, ymax = 300)

enid$TypeClass <- "N"
enid$TypeClass[grep("P", enid$enid)] <- "I1P"
enid$TypeClass[grep("M", enid$enid)] <- "I2M"
enid$TypeClass[grep("L", enid$Type)] <- "LS"
enid$TypeClass[grep("WL", enid$Type)] <- "LW"
enid$TypeClass[grep("E", enid$Type)] <- "ES"
enid$TypeClass[grep("WE", enid$Type)] <- "EW"
#enid$TypeClass[grep("E", enid$Type)] <- "X"
#enid$TypeClass[grep("L", enid$Type)] <- "X"

precip.shift06m <- precip[7:(nrow(precip)-6), ]

precip.shift06m.enid.split.median <- combineAOPI(enid, precip.shift06m)

if(length(precip.shift06m.enid.split.median) == 5){
  colors <- c("blue", "lightblue", "red", "bisque", "black")
} else if (length(precip.shift06m.enid.split.median) == 4){
  colors <- c("blue", "red", "black", "bisque")
} else {
  colors <- c("blue", "red", "black")  
}

plot.precip.shift06m.enid.split.median <- 
  seasonPlotByAOI(precip.shift06m.enid.split.median, colors,
                  normal = plot.precip.shift06m.enid.split.median.normal,
                  ymin = 0, ymax = 300)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.enid.split.median.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.enid.split.median)
  dev.off()
} else {
  plot(plot.precip.shift06m.enid.split.median)
}

# Just for background information: the individual seasons which form the above
# combined El Nino, La Nina and normal situation spline.
precip.shift06m.enid.split <- combineAOPI(enid, precip.shift06m, rt = "split")

precip.shift06m.enid.split.season <- lapply(precip.shift06m.enid.split, function(x){
  split(x, x$Season)
})

plot.precip.shift06m.enid.split.season.all <- 
  seasonPlotByAOI(precip.shift06m.enid.split.season, colors, individual = TRUE)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.enid.split.season.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.enid.split.season.all)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.enid.split.season.all)
}

# Alternative boxplot visualization
## bwplots
precip.shift06m.enid <- combineAOPI(enid, precip.shift06m, rt = "org")
precip.shift06m.enid$month <- substr(precip.shift06m.enid$ts, 6, 7)
precip.shift06m.enid.melt <- melt(precip.shift06m.enid, 
                                  id.vars = c("month", "TypeClass"),
                                  measure.vars = "P_RT_NRT")
bw.plot <- ggplot(data = precip.shift06m.enid.melt, 
                  aes(y = value, x = month, fill = TypeClass, 
                      dodge = TypeClass))
plot.precip.shift06m.enid.split.median.all.bw <- bw.plot + geom_boxplot() +
  xlab("Month") + ylab("Preciptiation")

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.enid.split.median.all.bw.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.enid.split.median.all.bw)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.enid.split.median.all.bw)
}


# Create publication quality figures of correlations between enid and 
# precipitation
years <- unique(substr(precip.shift06m.enid$ts, 1, 4))
precip.shift06m.mat.ssn_kz03k01 <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.enid[grep(i, precip.shift06m.enid$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$ssn_kz03k01[j]
  }
}

enid.reshape.shift06m.mat <- foreach(i = years, .combine = "rbind") %do% {
  dat <- precip.shift06m.enid[grep(i, precip.shift06m.enid$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$aoi[j]
  }
}

plot.precip.shift06m.enid.ssn_kz03k01 <- 
  corPlotByAOI(enid.reshape.shift06m.mat, precip.shift06m.mat.ssn_kz03k01)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.shift06m.enid.ssn_kz03k01.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.shift06m.enid.ssn_kz03k01)
  dev.off()
  
} else {
  plot(plot.precip.shift06m.enid.ssn_kz03k01)
}



### Cloud EOT analysis vs ENSO ################################################
# Prepare aoi record and classifiy years as La Nina (L), El Nino (E) or 
# normal (N); weak ENSO cycles are classified as normal
enso <- aoi.list$ONI
enso$TypeClass <- "N"

cloudEOT.shift06m <- cloudEOT[7:(nrow(cloudEOT)-6), ]
colnames(cloudEOT.shift06m)[2] <-  "P_RT_NRT"

cloudEOT.shift06m.enso.split.median <- combineAOPI(enso, cloudEOT.shift06m)

colors <- c("orange")  

plot.cloudEOT.shift06m.enso.split.median.normal <- 
  seasonPlotByAOI(cloudEOT.shift06m.enso.split.median, colors,
                  ymin = -0.2, ymax = 0.2)

enso$TypeClass <- "N"
enso$TypeClass[grep("L", enso$Type)] <- "LS"
# enso$TypeClass[grep("WL", enso$Type)] <- "LW"
enso$TypeClass[grep("E", enso$Type)] <- "ES"
# enso$TypeClass[grep("WE", enso$Type)] <- "EW"
# enso$TypeClass[grep("W", enso$Type)] <- "W"
# enso$TypeClass[grep("M", enso$Type)] <- "X"

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
                  normal = plot.cloudEOT.shift06m.enso.split.median.normal,
                  ymin = -0.2, ymax = 0.2)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.cloudEOT.shift06m.enso.split.median.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.cloudEOT.shift06m.enso.split.median.all)
  dev.off()
  
} else {
  plot(plot.cloudEOT.shift06m.enso.split.median.all)
}

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

cloudEOT.shift06m.enso.split <- combineAOPI(enso, cloudEOT.shift06m, rt = "split")

cloudEOT.shift06m.enso.split.season <- lapply(cloudEOT.shift06m.enso.split, function(x){
  split(x, x$Season)
})

plot.cloudEOT.shift06m.enso.split.season.all <- 
  seasonPlotByAOI(cloudEOT.shift06m.enso.split.season, colors, 
                  individual = TRUE,
                  ymin = -0.2, ymax = 0.2)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.cloudEOT.shift06m.enso.split.season.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.cloudEOT.shift06m.enso.split.season.all)
  dev.off()
  
} else {
  plot(plot.cloudEOT.shift06m.enso.split.season.all)
}



# Alternative boxplot visualization
## bwplots
cloudEOT.shift06m.enso <- combineAOPI(enso, cloudEOT.shift06m, rt = "org")
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
                         "plot.cloudEOT.shift06m.enso.split.median.all.bw.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.cloudEOT.shift06m.enso.split.median.all.bw)
  dev.off()
} else {
  plot(plot.cloudEOT.shift06m.enso.split.median.all.bw)
}

# Create publication quality figures of correlations between enso and 
# cloudEOT
years <- unique(substr(cloudEOT$ts, 1, 4))
cloudEOT.shift06m.mat.ssn_kz03k01 <- foreach(i = years, .combine = "rbind") %do% {
  dat <- cloudEOT.shift06m.enso[grep(i, cloudEOT.shift06m.enso$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$ssn_kz03k01[j]
  }
}

enso.reshape.shift06m.mat <- foreach(i = years, .combine = "rbind") %do% {
  dat <- cloudEOT.shift06m.enso[grep(i, cloudEOT.shift06m.enso$ts), ]
  foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
    dat$aoi[j]
  }
}

plot.cloudEOT.shift06m.enso.ssn_kz03k01 <- 
  corPlotByAOI(enso.reshape.shift06m.mat, cloudEOT.shift06m.mat.ssn_kz03k01)


if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.cloudEOT.shift06m.enso.ssn_kz03k01.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.cloudEOT.shift06m.enso.ssn_kz03k01)
  dev.off()
  
} else {
  plot(plot.cloudEOT.shift06m.enso.ssn_kz03k01)
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

