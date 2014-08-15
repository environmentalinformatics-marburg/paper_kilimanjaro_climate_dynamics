#### Set working environment ###################################################
switch(Sys.info()[["sysname"]], 
       "Linux" = setwd("/media/permanent/active/kilimanjaro"), 
       "Windows" = setwd("D:/active/kilimanjaro"))
#setwd("C:/Users/Dogbert/Desktop/kilimanjaro")) 

sourcePath <- "scripts/paper_kilimanjaro_climate_dynamics/src/"
dataPath <- "data/"
graphicsPath <- "graphics/"

printToFile <- FALSE
plot2file <- printToFile

library(kza)
library(latticeExtra)
library(Kendall)
library(foreach)
library(ade4)
library(reshape2)
library(ggplot2)
source(paste0(sourcePath, "createDataSet.R"))
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

grid31 <- function(x, y, ...) {
  update(c(x, y,
           layout = c(3, 1)),
         between = list(y = 0.3, x = 0.3))
} 

grid13 <- function(x, y, ...) {
  update(c(x, y,
           layout = c(1, 3)),
         between = list(y = 0.3, x = 0.3))
} 


#### Create data sets ##########################################################
data.set <- createDataSet()

#### Long-term temperature analysis ############################################
# Compute long-term anomalies and trends for KIA temperature and create
# publication quality figure.
# Compute 3 month running mean of original temperature values using a
# Kolmogorov-Zurbenko filter with one iteration
ta <- data.set$ta.list$KIA
ta.org <- data.set$ta.list$KIAOrg
longTermDynamicsPlot(parameter = "temperature", printToFile = printToFile)


#### Long-term precipitation analysis ##########################################
# Compute long-term anomalies and trends for KIA precipitation and create
# publication quality figure.
precip <- data.set$precip.list$KIA
longTermDynamicsPlot(parameter = "precipitation", printToFile = printToFile,
                     p.prm = "ssn_kz03k02")


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
  st = c(1975,1975,1992,2001), nd = c(2013,1991,2000,2013),
  st.shift = 7, nd.shift = 0, timespan = 12, fun = "median", prm = "P_RT_NRT")

precip.seasonalwetdry.split <- 
  split(precip.seasonalwetdry, precip.seasonalwetdry$season)
precip.seasonal.normal <- 
  list(lapply(precip.seasonalwetdry.split, function(x){x$p_dyn})$"1975-2013")

plot.precip.seasonalwetdry.all <- 
  visSeasonPlotByAOI(
  lapply(precip.seasonalwetdry.split, function(x){x$p_dyn}), colors,
  linetype = c(2,1,1,1), ymin = yminmax[1], ymax = yminmax[2], timespan = 12,
  x.text = c(450,850,500,350), 
  y.text = c(40,200,11,65),
  labels.text = c("normal", "1975-1992", "1992-2000","2001-2013"),
  colors.text = colors)

# 18 month season
# precip.18m.seasonalwetdry <- seasonalMean(
#   st = c(1975,1975,1992,2001), nd = c(2013,1992,2000,2013), 
#   st.shift = 7, nd.shift = 0, timespan = 18, fun = "median", prm = "P_RT_NRT")
# 
# precip.18m.seasonalwetdry.split <- 
#   split(precip.18m.seasonalwetdry, precip.18m.seasonalwetdry$season)
# precip.18m.seasonal.normal <- 
#   list(lapply(precip.18m.seasonalwetdry.split, function(x){x$p_dyn})$"1975-2013")
# 
# plot.precip18m.seasonalwetdry.all <- 
#   visSeasonPlotByAOI(
#     lapply(precip.18m.seasonalwetdry.split, function(x){x$p_dyn}), colors,
#     linetype = c(2,1,1,1), ymin = yminmax[1], ymax = yminmax[2], timespan = 18)
# 
# # 24 month season
# precip.24m.seasonalwetdry <- seasonalMean(
#   st = c(1975,1975,1992,2001), nd = c(2013,1992,2000,2013), 
#   st.shift = 7, nd.shift = 0, timespan = 24, fun = "median", prm = "P_RT_NRT")
# 
# precip.24m.seasonalwetdry.split <- 
#   split(precip.24m.seasonalwetdry, precip.24m.seasonalwetdry$season)
# precip.24m.seasonal.normal <- 
#   list(lapply(precip.24m.seasonalwetdry.split, function(x){x$p_dyn})$"1975-2013")

# 20 month season
precip.20m.seasonalwetdry <- seasonalMean(
  st = c(1975,1975,1992,2001), nd = c(2013,1992,2000,2013), 
  st.shift = 7, nd.shift = 0, timespan = 20, fun = "median", prm = "P_RT_NRT")

precip.20m.seasonalwetdry.split <- 
  split(precip.20m.seasonalwetdry, precip.20m.seasonalwetdry$season)
precip.20m.seasonal.normal <- 
  list(lapply(precip.20m.seasonalwetdry.split, function(x){x$p_dyn})$"1975-2013")

# plot.precip24m.seasonalwetdry.all <- 
#   visSeasonPlotByAOI(
#   lapply(precip.24m.seasonalwetdry.split, function(x){x$p_dyn}), colors,
#   linetype = c(2,1,1,1), ymin = yminmax[1], ymax = yminmax[2], timespan = 24,
#   x.text = c(1650,1100,610,400), 
#   y.text = c(25,200,11,85),
#   labels.text = c("normal", "1975-1992", "1992-2000","2001-2013"),
#   colors.text = colors)

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, "plot.precip.seasonalwetdry.all.tif"),
       width = 2480, height = 1748 , res = 300, pointsize =  12)
  plot(plot.precip.seasonalwetdry.all)
  dev.off()
#   tiff(filename = paste0(graphicsPath, "plot.precip18m.seasonalwetdry.all.tif"),
#        width = 2480, height = 1748 , res = 300, pointsize =  12)
#   plot(plot.precip18m.seasonalwetdry.all)
#   dev.off()
#   tiff(filename = paste0(graphicsPath, "plot.precip24m.seasonalwetdry.all.tif"),
#        width = 2480, height = 1748 , res = 300, pointsize =  12)
#   plot(plot.precip24m.seasonalwetdry.all)
#   dev.off()
} else {
 plot(plot.precip.seasonalwetdry.all)
#   plot(plot.precip18m.seasonalwetdry.all)
#   plot(plot.precip24m.seasonalwetdry.all)
}


#### Precipitation analysis vs ENSO ############################################
# Compute plot for long-term normal distribution
# yminmax = c(0, 250)
# #yminmax = c(-200,200)
# colors <- c("black")
# 
# # plot.precip.shift06m.normal <- 
# #   seasonPlotByAOI(precip.seasonal.normal, colors,
# #                   linetype = c(2), ymin = yminmax[1], ymax = yminmax[2])
# # 
# # plot.precip.18m.normal <- 
# #   visSeasonPlotByAOI(precip.18m.seasonal.normal, colors,
# #                      linetype = c(2), ymin = yminmax[1], ymax = yminmax[2])
# # 
# # plot.precip.24m.normal <- 
# #   visSeasonPlotByAOI(precip.24m.seasonal.normal, colors,
# #                      linetype = c(2), ymin = yminmax[1], ymax = yminmax[2],
# #                      timespan = 24)
# 
# plot.precip.20m.normal <- 
#   visSeasonPlotByAOI(precip.20m.seasonal.normal, colors,
#                      linetype = c(2), ymin = yminmax[1], ymax = yminmax[2],
#                      timespan = 20)
# 
# 
# # Prepare aoi record and classifiy years as La Nina (L), El Nino (E) or 
# # normal (N); weak ENSO cycles are classified ENSO
# # Compute seasonal distribution by major aoi situation
# red <- brewer.pal(4, "Reds")
# blue <- brewer.pal(4, "Blues")
# 
# enso <- data.set$aoi.list$ONI
# enso$TypeClass <- "C9 Normal"
# enso$TypeClass[grep("E", enso$Type)] <- "C1 El Nino"
# enso$TypeClass[grep("L", enso$Type)] <- "C2 La Nina"
# 
# colors <- c(blue[4], red[4], "black")
# linetype <- c(1,1,1)
# 
# # precip.shift06m <- precip[7:(nrow(precip)-6), ]
# # precip.shift06m.enso.split.median <- combineAOPI(enso, precip.shift06m)
# # 
# # precip.18m <- precip[7:(nrow(precip)-0), ]
# # precip.18m.enso.split.median <- mergeAOIwTS(enso, precip.18m, 
# #                                             timespan = 18,
# #                                             ts.prm = "P_RT_NRT",
# #                                             rt = "median")
# 
# precip.20m <- precip[7:(nrow(precip)-0), ]
# precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
#                                             timespan = 20,
#                                             ts.prm = "P_RT_NRT",
#                                             rt = "median")
# 
# precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20, 
#                               ts.prm = "P_RT_NRT", rt = "org")
# unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 El Nino"])
# unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C2 La Nina"])
# unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C9 Normal"])
# 
# # plot.precip.shift06m.enso.split.median.all <- 
# #   seasonPlotByAOI(precip.shift06m.enso.split.median, colors,
# #                   linetype = linetype,
# #                   normal = plot.precip.shift06m.normal,
# #                   ymin = yminmax[1], ymax = yminmax[2])
# # 
# # plot.precip.18m.enso.split.median.all <- 
# #   visSeasonPlotByAOI(precip.18m.enso.split.median, colors,
# #                   linetype = linetype,
# #                   normal = plot.precip.18m.normal,
# #                   ymin = yminmax[1], ymax = yminmax[2],
# #                   vline.pos = 501)
# # 
# # plot.precip.24m.enso.split.median.all <- 
# #   visSeasonPlotByAOI(precip.24m.enso.split.median, colors,
# #                      linetype = linetype,
# #                      normal = plot.precip.24m.normal,
# #                      ymin = yminmax[1], ymax = yminmax[2],
# #                      timespan = 24,
# #                      vline.pos = 501,
# #                      x.text = c(200,1050,950), 
# #                      y.text = c(60,180,11),
# #                      labels.text = c("El Ninos", "La Ninas", "normal"),
# #                      colors.text = colors)
# 
# plot.precip.20m.enso.all <- 
#   visSeasonPlotByAOI(precip.20m.enso.split.median, colors,
#                      linetype = linetype,
#                      normal = plot.precip.20m.normal,
#                      ymin = yminmax[1], ymax = yminmax[2],
#                      timespan = 20,
#                      vline.pos = 501,
#                      x.text = c(200,1050,950), 
#                      y.text = c(60,180,11),
#                      labels.text = c("El Ninos", "La Ninas", "normal"),
#                      colors.text = colors)
# 
# if(printToFile == TRUE){
# #   tiff(filename = paste0(graphicsPath, 
# #                          "plot.precip.shift06m.enso.split.median.all.tif"),
# #        width = 2480, height = 1748 , res = 300, pointsize =  12)
# #   plot(plot.precip.shift06m.enso.split.median.all)
# #   dev.off()
# #   tiff(filename = paste0(graphicsPath, 
# #                          "plot.precip.18m.enso.split.median.all.tif"),
# #        width = 2480, height = 1748 , res = 300, pointsize =  12)
# #   plot(plot.precip.18m.enso.split.median.all)
# #   dev.off()
#   tiff(filename = paste0(graphicsPath, 
#                          "plot.precip.20m.enso.all.tif"),
#        width = 2480, height = 1748 , res = 300, pointsize =  12)
#   plot(plot.precip.20m.enso.all)
#   dev.off()
# } else {
# #   plot(plot.precip.shift06m.enso.split.median.all)
# #   plot(plot.precip.18m.enso.split.median.all)
# #   plot(plot.precip.24m.enso.split.median.all)
#   plot(plot.precip.20m.enso.all)
# }
# 
# # Prepare aoi record and classifiy years as La Nina (L), El Nino (E) or 
# # normal (N); weak ENSO cycles are classified ENSO
# # Compute seasonal distribution by major aoi situation
# enso$TypeClass <- "C9 Normal"
# enso$TypeClass[grep("E", enso$Type)] <- "C1 El Nino"
# enso$TypeClass[grep("L", enso$Type)] <- "C2 La Nina"
# enso$TypeClass[grep("W", enso$Type)] <- "C9 Normal"
# 
# colors <- c(blue[4], red[4], "black")
# linetype <- c(1,1,1)
# 
# precip.20m <- precip[7:(nrow(precip)-0), ]
# precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
#                                             timespan = 20,
#                                             ts.prm = "P_RT_NRT",
#                                             rt = "median")
# 
# unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 El Nino"])
# unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C2 La Nina"])
# unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C9 Normal"])
# 
# plot.precip.20m.enso.ms <- 
#   visSeasonPlotByAOI(precip.20m.enso.split.median, colors,
#                      linetype = linetype,
#                      normal = plot.precip.20m.normal,
#                      ymin = yminmax[1], ymax = yminmax[2],
#                      timespan = 20,
#                      vline.pos = 501,
#                      x.text = c(200,1050,950), 
#                      y.text = c(60,180,11),
#                      labels.text = c("El Ninos", "La Ninas", "normal"),
#                      colors.text = colors)
# 
# 
# 
# if(printToFile == TRUE){
#   tiff(filename = paste0(graphicsPath, 
#                          "plot.precip.20m.enso.ms.tif"),
#        width = 2480, height = 1748 , res = 300, pointsize =  12)
#   plot(plot.precip.20m.enso.ms)
#   dev.off()
# } else {
#   plot(plot.precip.20m.enso.ms)
# }

# Prepare plot for ENSO
yminmax = c(0, 250)
plot.colors <- c("grey")

plot.precip.20m.normal <- 
  visSeasonPlotByAOI(precip.20m.seasonal.normal, plot.colors,
                     linetype = c(3), ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20)


plot.colors <- c("black", blue[3], blue[4], red[3], red[4], "darkgreen", "black")
linetype <- c(1)

# El Nino - all
enso <- data.set$aoi.list$ONI
enso$TypeClass <- "C9 Normal"
enso$TypeClass[enso$Type == "WE" | 
                 enso$Type == "ME" | 
                 enso$Type == "SE"] <- "C1 all El Nino"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 all El Nino"])
color <- plot.colors[1]
plot.precip.20m.enso.elnino.all <- 
  visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.precip.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     x.text = c(150), 
                     y.text = c(70),
                     labels.text = c("all El Ninos"),
                     colors.text = color)
plot.precip.20m.enso.elnino.all

# Pure El Nino
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WE" | enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD != "P" & enso$IOD != "M"] <- "C1 pure El Ninos"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 pure El Ninos"])
color <- plot.colors[2]
plot.precip.20m.enso.elnino.pure <- 
  visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.precip.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     x.text = c(750), 
                     y.text = c(200),
                     labels.text = c("pure El Ninos"),
                     colors.text = color)
plot.precip.20m.enso.elnino.pure

# Pure medium and stron El Nino
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD != "P" & enso$IOD != "M"] <- "C1 pure m/s El Ninos"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 pure m/s El Ninos"])
color <- plot.colors[3]
plot.precip.20m.enso.elnino.pure.ms <- 
  visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.precip.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     x.text = c(350), 
                     y.text = c(125),
                     labels.text = c("pure m/s El Ninos"),
                     colors.text = color)
plot.precip.20m.enso.elnino.pure.ms

# El Nino with IOD+
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WE" | enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD == "P" ] <- "C1 El Nino w IOD+"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 El Nino w IOD+"])
color <- plot.colors[4]
plot.precip.20m.enso.elnino.wIOD.all <- 
  visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.precip.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     x.text = c(1550), 
                     y.text = c(60),
                     labels.text = c("El Ninos w IOD+"),
                     colors.text = color)
plot.precip.20m.enso.elnino.wIOD.all

# Medium and strong El Nino with IOD+
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD == "P" ] <- "C1 m/s El Nino w IOD+"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 m/s El Nino w IOD+"])
color <- plot.colors[5]
plot.precip.20m.enso.elnino.wIOD.ms <- 
  visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.precip.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     x.text = c(900), 
                     y.text = c(5),
                     labels.text = c("m/s El Ninos w IOD+"),
                     colors.text = color)
plot.precip.20m.enso.elnino.wIOD.ms

# El Nino with IOD-
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WE" | enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD == "M" ] <- "C1 El Nino w IOD-"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 El Nino w IOD-"])

# IOD+ - all
enso$TypeClass <- "C9 Normal"
enso$TypeClass[enso$Type != "ME" & enso$Type != "SE" &
                 enso$Type != "ML" & enso$Type != "SL" &
                 enso$IOD == "P"] <- "C1 purest IOD+"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 purest IOD+"])
color <- plot.colors[6]
plot.precip.20m.iod.p.all <- 
  visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
                     linetype = 2,
                     normal = plot.precip.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     x.text = c(700), 
                     y.text = c(170),
                     labels.text = c("purest IOD+"),
                     colors.text = color)
plot.precip.20m.iod.p.all

plot.precip.20m.elnino <- 
  Reduce("outLayer", c(list(plot.precip.20m.enso.elnino.all,
                            plot.precip.20m.enso.elnino.pure,
                            plot.precip.20m.enso.elnino.pure.ms,
                            plot.precip.20m.enso.elnino.wIOD.all,
                            plot.precip.20m.enso.elnino.wIOD.ms,
                            plot.precip.20m.iod.p.all)))

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.20m.elnino.tif"),
       width = 30, height = 15, units = "cm", res = 600, pointsize =  5)
  plot(plot.precip.20m.elnino)
  dev.off()
  
  pdf(file = paste0(graphicsPath, "plot.precip.20m.elnino.pdf"),
      width = 12, height = 6, paper = "a4r")
  plot(plot.precip.20m.elnino)
  dev.off()
} else {
  plot(plot.precip.20m.elnino)
}




# Prepare plot for La Nina
yminmax = c(0, 250)
plot.colors <- c("grey")

plot.precip.20m.normal <- 
  visSeasonPlotByAOI(precip.20m.seasonal.normal, plot.colors,
                     linetype = c(2), ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20)


plot.colors <- c("black", blue[3], blue[4], red[3], red[4], "darkgreen", "black")
linetype <- c(1)

# La Nina - all
enso$TypeClass <- "C9 Normal"
enso$TypeClass[enso$Type == "WL" | 
                 enso$Type == "ML" | 
                 enso$Type == "SL"] <- "C1 all La Nina"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 all La Nina"])
color <- plot.colors[1]
plot.precip.20m.enso.lanina.all <- 
  visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.precip.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     x.text = c(200), 
                     y.text = c(25),
                     labels.text = c("all La Ninas"),
                     colors.text = color)
plot.precip.20m.enso.lanina.all

# Pure La Nina
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WL" | enso$Type == "ML" | enso$Type == "SL") &
                 enso$IOD != "P" & enso$IOD != "M"] <- "C1 pure La Ninas"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 pure La Ninas"])
color <- plot.colors[2]
plot.precip.20m.enso.lanina.pure <- 
  visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.precip.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     x.text = c(275), 
                     y.text = c(60),
                     labels.text = c("pure La Ninas"),
                     colors.text = color)
plot.precip.20m.enso.lanina.pure

# Pure medium and stron La Nina
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "ML" | enso$Type == "SL") &
                 enso$IOD != "P" & enso$IOD != "M"] <- "C1 pure m/s La Ninas"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 pure m/s La Ninas"])
color <- plot.colors[3]
plot.precip.20m.enso.lanina.pure.ms <- 
  visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.precip.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     x.text = c(1150), 
                     y.text = c(110),
                     labels.text = c("pure m/s La Ninas"),
                     colors.text = color)
plot.precip.20m.enso.lanina.pure.ms

# La Nina with IOD+
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WL" | enso$Type == "ML" | enso$Type == "SL") &
                 enso$IOD == "P" ] <- "C1 La Nina w IOD+"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 La Nina w IOD+"])
color <- plot.colors[4]
plot.precip.20m.enso.lanina.wIOD.all <- 
  visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.precip.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     x.text = c(1050), 
                     y.text = c(230),
                     labels.text = c("La Ninas w IOD+"),
                     colors.text = color)
plot.precip.20m.enso.lanina.wIOD.all

# Medium and strong La Nina with IOD+
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "ML" | enso$Type == "SL") &
                 enso$IOD == "P" ] <- "C1 m/s La Nina w IOD+"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 m/s La Nina w IOD+"])
# color <- plot.colors[5]
# plot.precip.20m.enso.lanina.wIOD.ms <- 
#   visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
#                      linetype = linetype,
#                      normal = plot.precip.20m.normal,
#                      ymin = yminmax[1], ymax = yminmax[2],
#                      timespan = 20,
#                      vline.pos = 501,
#                      x.text = c(900), 
#                      y.text = c(5),
#                      labels.text = c("m/s La Ninas w IOD+"),
#                      colors.text = color)

# La Nina with IOD-
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WL" | enso$Type == "ML" | enso$Type == "SL") &
                 enso$IOD == "M" ] <- "C1 La Nina w IOD-"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 La Nina w IOD-"])
# plot.precip.20m.enso.lanina.wIODM.all <- 
#   visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
#                      linetype = linetype,
#                      normal = plot.precip.20m.normal,
#                      ymin = yminmax[1], ymax = yminmax[2],
#                      timespan = 20,
#                      vline.pos = 501,
#                      x.text = c(900), 
#                      y.text = c(5),
#                      labels.text = c("La Ninas w IOD-"),
#                      colors.text = color)

plot.precip.20m.lanina <- 
  Reduce("outLayer", c(list(plot.precip.20m.enso.lanina.all,
                            plot.precip.20m.enso.lanina.pure,
                            plot.precip.20m.enso.lanina.pure.ms,
                            plot.precip.20m.enso.lanina.wIOD.all)))

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.precip.20m.lanina.tif"),
       width = 30, height = 15, units = "cm", res = 600, pointsize =  5)
  plot(plot.precip.20m.lanina)
  dev.off()
  
  pdf(file = paste0(graphicsPath, "plot.precip.20m.lanina.pdf"),
      width = 12, height = 6, paper = "a4r")
  plot(plot.precip.20m.lanina)
  dev.off()
} else {
  plot(plot.precip.20m.lanina)
}




# IOD
# IOD- - all
enso$TypeClass <- "C9 Normal"
enso$TypeClass[enso$Type != "ME" & enso$Type != "SE" &
                 enso$Type != "ML" & enso$Type != "SL" &
                 enso$IOD == "M"] <- "C1 purest IOD-"
precip.20m <- precip[7:(nrow(precip)-0), ]
precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                            timespan = 20,
                                            ts.prm = "P_RT_NRT",
                                            rt = "median")
precip.20m.info <- mergeAOIwTS(enso, precip.20m, timespan = 20,
                               ts.prm = "P_RT_NRT", rt = "org")
unique(precip.20m.info$Season[precip.20m.info$TypeClass == "C1 purest IOD-"])
# color <- plot.colors[5]
# plot.precip.20m.iod.m.all <- 
#   visSeasonPlotByAOI(precip.20m.enso.split.median[1], color,
#                      linetype = linetype,
#                      normal = plot.precip.20m.normal,
#                      ymin = yminmax[1], ymax = yminmax[2],
#                      timespan = 20,
#                      vline.pos = 501,
#                      x.text = c(200), 
#                      y.text = c(25),
#                      labels.text = c("purest IOD-"),
#                      colors.text = color)


# Alternative boxplot visualization
# precip.shift06m.enso <- combineAOPI(enso, precip.shift06m, rt = "org")
# precip.shift06m.enso$month <- substr(precip.shift06m.enso$ts, 6, 7)
# precip.shift06m.enso.melt <- melt(precip.shift06m.enso, 
#                                   id.vars = c("month", "TypeClass"),
#                                   measure.vars = "P_RT_NRT")
# bw.plot <- ggplot(data = precip.shift06m.enso.melt, 
#                   aes(y = value, x = month, fill = TypeClass, 
#                       dodge = TypeClass))
# plot.precip.shift06m.enso.split.median.all.bw <- bw.plot + geom_boxplot() +
#   xlab("Month") + ylab("Preciptiation")
# 
# if(printToFile == TRUE){
#   tiff(filename = paste0(graphicsPath, 
#                          "plot.precip.shift06m.enso.split.median.all.bw.tif"),
#        width = 2480, height = 1748 , res = 300, pointsize =  12)
#   plot(plot.precip.shift06m.enso.split.median.all.bw)
#   dev.off()
#   
# } else {
#   plot(plot.precip.shift06m.enso.split.median.all.bw)
# }

# Create publication quality figures of correlations between enso and 
# precipitation
#test <- precip.shift06m.enso[precip.shift06m.enso$TypeClass == "N",]
# precip.shift06m.enso <- combineAOPI(enso, precip.shift06m, rt = "org")
# precip.shift06m.enso$StartSeason <- substr(precip.shift06m.enso$Season,1,4)
# test <- precip.shift06m.enso
# m12 <- visCorPlotTimeSeries(df = test, 
#                      x.prm = "aoi",
#                      y.prm = "P_RT_NRT",
#                      t.prm = "StartSeason",
#                      lable.nbrs = c(c(7:12),c(1:6)),
#                      p.thv = 0.08,
#                      plot.filepath = "plot.precip.shift06m.enso.cor.tif",
#                      plot2file = plot2file,
#                      rt = TRUE)
# 
# precip.18m.enso <- mergeAOIwTS(enso, precip.18m, 
#                                timespan = 18, rt = "org")
# test <- precip.18m.enso # [precip.18m.enso$TypeClass == "El Nino", ]
# m18 <- visCorPlotTimeSeries(df = test, 
#                      x.prm = "aoi",
#                      y.prm = "P_RT_NRT",
#                      t.prm = "StartSeason",
#                      lable.nbrs = c(c(7:12),c(1:12)),
#                      p.thv = 0.08,
#                      plot.filepath = "plot.precip18m.enso.cor.tif",
#                      plot2file = plot2file,
#                      rt = TRUE)
# 
# precip.24m.enso <- mergeAOIwTS(enso, precip.24m, 
#                                timespan = 24, rt = "org")
# test <- precip.24m.enso # [precip.24m.enso$TypeClass == "El Nino",]
# m24 <- visCorPlotTimeSeries(df = test, 
#                      x.prm = "aoi",
#                      y.prm = "P_RT_NRT",
#                      t.prm = "StartSeason",
#                      x.lable = "ONI",
#                      y.lable = "P",
#                      lable.nbrs = c(c(7:12),c(1:12),c(1:6)),
#                      p.thv = 0.08,
#                      plot.filepath = "plot.precip24m.enso.cor.tif",
#                      plot2file = plot2file,
#                      rt = TRUE)
# 
# testm12m18 <- m12 - m18[1:12,1:12]
# testm18m24 <- m18 - m24[1:18,1:18]
# range(testm12m18, na.rm = TRUE)
# range(testm18m24, na.rm = TRUE)
enso <- data.set$aoi.list$ONI
enso$TypeClass <- "C9 Normal"
precip.20m.enso <- mergeAOIwTS(enso, precip.20m, 
                               timespan = 20, rt = "org")
plot.cor.enso <- visCorPlotTimeSeries(df = precip.20m.enso, 
                            x.prm = "aoi",
                            y.prm = "P_RT_NRT",
                            t.prm = "StartSeason",
                            x.lable = "OI",
                            y.lable = "P",
                            lable.nbrs = c(c(7:12),c(1:12),c(1:2)),
                            p.thv = 0.08,
                            x.text = 1,
                            y.text = 19,
                            labels.text = "ENSO",
                            colors.text = "black",
                            plot.filepath = "plot.precip20m.enso.cor.tif",
                            plot2file = plot2file,
                            rp = TRUE)

#### Precipitation analysis vs IOD #############################################
# Prepare aoi record and classifiy years as IOD plus (P), IOD minus (M) or
# normal (N)
iod <- data.set$aoi.list$DMIHAD
iod$TypeClass <- "C9 Normal"

# Create publication quality figures of correlations between iod and 
# precipitation
precip.20m.iod <- mergeAOIwTS(iod, precip.20m, 
                              timespan = 20, rt = "org")
plot.cor.iod <- visCorPlotTimeSeries(df = precip.20m.iod, 
                            x.prm = "aoi",
                            y.prm = "P_RT_NRT",
                            t.prm = "StartSeason",
                            y.lable = "P",
                            lable.nbrs = c(c(7:12),c(1:12),c(1:2)),
                            p.thv = 0.08,
                            x.lable = "OI",
                            x.text = 1,
                            y.text = 19,
                            labels.text = "IOD",
                            colors.text = "black",
                            plot.filepath = "plot.precip24m.iod.cor.tif",
                            plot2file = plot2file,
                            rp = TRUE)


#### Precipitation analysis vs enid #############################################
# Prepare aoi record and classifiy years as enid plus (P), enid minus (M) or
# normal (N)
enidENSO <- data.set$aoi.list$ONI
enidDMI <- data.set$aoi.list$DMIHAD
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

# Create publication quality figures of correlations between enid and 
# precipitation
precip.20m.enid <- mergeAOIwTS(enid, precip.20m, 
                              timespan = 20, rt = "org")
plot.cor.enid <- visCorPlotTimeSeries(df = test, 
                            x.prm = "aoi",
                            y.prm = "P_RT_NRT",
                            t.prm = "StartSeason",
                            y.lable = "P",
                            lable.nbrs = c(c(7:12),c(1:12),c(1:2)),
                            p.thv = 0.08,
                            x.lable = "OI",
                            x.text = 1,
                            y.text = 19,
                            labels.text = "ENSO and IOD",
                            colors.text = "black",
                            plot.filepath = "plot.precip24m.enid.tif",
                            plot2file = plot2file,
                            rp = TRUE)

# Combine cor plots
plot.cor.precip.all <- 
  Reduce("grid31", list(plot.cor.enso, plot.cor.iod, plot.cor.enid))


if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.cor.precip.all.tif"),
       width = 30, height = 15, units = "cm", res = 600, pointsize =  5)
  plot(plot.cor.precip.all)
  dev.off()
  
  pdf(file = paste0(graphicsPath, "plot.cor.precip.all.pdf"),
      width = 12, height = 6, paper = "a4r")
  plot(plot.cor.precip.all)
  dev.off()
} else {
  plot(plot.cor.precip.all)
}



### Cloud EOT analysis vs ENSO ################################################
# Prepare aoi record and classifiy years as La Nina (L), El Nino (E) or 
# normal (N); weak ENSO cycles are classified as normal
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

cloud <- data.set$cloud.list$Clust

cloud.precip <- merge(cloud, precip, by = "ts")
cor(cloud.precip$clust1,cloud.precip$P_RT_NRT, use = "pairwise.complete.obs",
    method = "kendall")
cor(cloud.precip$clust2,cloud.precip$P_RT_NRT, use = "pairwise.complete.obs",
    method = "kendall")
cor(cloud.precip$clust3,cloud.precip$P_RT_NRT, use = "pairwise.complete.obs",
    method = "kendall")

# Cluster 1 - El Nino
act.prm <- "clust1"
cloud.20m.seasonalwetdry <- seasonalMean(data = cloud,
  st = c(2003), nd = c(2013), 
  st.shift = 7, nd.shift = 0, timespan = 20, fun = "median", prm = act.prm)
cloud.20m.seasonalwetdry.split <- 
  split(cloud.20m.seasonalwetdry, cloud.20m.seasonalwetdry$season)
cloud.20m.seasonal.normal <- 
  list(lapply(cloud.20m.seasonalwetdry.split, function(x){x$p_dyn})$"2003-2013")

yminmax = c(0, 1)
plot.colors <- c("grey")

plot.cloud.20m.normal <- 
  visSeasonPlotByAOI(cloud.20m.seasonal.normal, plot.colors,
                     linetype = c(3), ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     ylab = "Cloud cover")


plot.colors <- c("black", blue[3], blue[4], red[3], red[4], "darkgreen", "black")
linetype <- c(1)


# El Nino - all
enso$TypeClass <- "C9 Normal"
enso$TypeClass[enso$Type == "WE" | 
                 enso$Type == "ME" | 
                 enso$Type == "SE"] <- "C1 all El Nino"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm =act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 all El Nino"])
color <- plot.colors[1]
plot.cloud.20m.enso.elnino.all <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(100), 
                     y.text = c(0.25),
                     labels.text = c("all El Ninos"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.all

# Pure El Nino
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WE" | enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD != "P" & enso$IOD != "M"] <- "C1 pure El Ninos"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 pure El Ninos"])
color <- plot.colors[2]
plot.cloud.20m.enso.elnino.pure <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(1300), 
                     y.text = c(0.50),
                     labels.text = c("pure El Ninos"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.pure

# Pure medium and stron El Nino
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD != "P" & enso$IOD != "M"] <- "C1 pure m/s El Ninos"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 pure m/s El Ninos"])
color <- plot.colors[3]
plot.cloud.20m.enso.elnino.pure.ms <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(200), 
                     y.text = c(0.6),
                     labels.text = c("pure m/s El Ninos"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.pure.ms

# El Nino with IOD+
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WE" | enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD == "P" ] <- "C1 El Nino w IOD+"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 El Nino w IOD+"])
color <- plot.colors[4]
plot.cloud.20m.enso.elnino.wIOD.all <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(610), 
                     y.text = c(0.6),
                     labels.text = c("El Ninos w IOD+"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.wIOD.all

# Medium and strong El Nino with IOD+
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD == "P" ] <- "C1 m/s El Nino w IOD+"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 m/s El Nino w IOD+"])
color <- plot.colors[5]
plot.cloud.20m.enso.elnino.wIOD.ms <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(250), 
                     y.text = c(0.20),
                     labels.text = c("m/s El Ninos w IOD+"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.wIOD.ms

# El Nino with IOD-
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WE" | enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD == "M" ] <- "C1 El Nino w IOD-"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 El Nino w IOD-"])

# IOD+ - all
enso$TypeClass <- "C9 Normal"
enso$TypeClass[enso$Type != "ME" & enso$Type != "SE" &
                 enso$Type != "ML" & enso$Type != "SL" &
                 enso$IOD == "P"] <- "C1 purest IOD+"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 purest IOD+"])
color <- plot.colors[6]
plot.cloud.20m.iod.p.all <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = 2,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(600), 
                     y.text = c(0.55),
                     labels.text = c("purest IOD+"),
                     colors.text = color)
plot.cloud.20m.iod.p.all

plot.cloud.20m.elnino.clust1 <- 
  Reduce("outLayer", c(list(plot.cloud.20m.enso.elnino.all,
                            plot.cloud.20m.enso.elnino.pure,
                            plot.cloud.20m.enso.elnino.pure.ms,
                            plot.cloud.20m.enso.elnino.wIOD.all,
                            plot.cloud.20m.enso.elnino.wIOD.ms,
                            plot.cloud.20m.iod.p.all)))

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.cloud.20m.elnino.clust1.tif"),
       width = 30, height = 15, units = "cm", res = 600, pointsize =  5)
  plot(plot.cloud.20m.elnino.clust1)
  dev.off()
  pdf(file = paste0(graphicsPath, "plot.cloud.20m.elnino.clust1.pdf"),
      width = 12, height = 6, paper = "a4r")
  plot(plot.cloud.20m.elnino.clust1)
  dev.off()
} else {
  plot(plot.cloud.20m.elnino.clust1)
}


# Cluster 2 - El Nino
act.prm <- "clust2"
cloud.20m.seasonalwetdry <- seasonalMean(data = cloud,
                                         st = c(2003), nd = c(2013), 
                                         st.shift = 7, nd.shift = 0, timespan = 20, fun = "median", prm = act.prm)
cloud.20m.seasonalwetdry.split <- 
  split(cloud.20m.seasonalwetdry, cloud.20m.seasonalwetdry$season)
cloud.20m.seasonal.normal <- 
  list(lapply(cloud.20m.seasonalwetdry.split, function(x){x$p_dyn})$"2003-2013")

yminmax = c(0, 1)
plot.colors <- c("grey")

plot.cloud.20m.normal <- 
  visSeasonPlotByAOI(cloud.20m.seasonal.normal, plot.colors,
                     linetype = c(3), ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20, ylab = "Cloud cover")


plot.colors <- c("black", blue[3], blue[4], red[3], red[4], "darkgreen", "black")
linetype <- c(1)


# El Nino - all
enso$TypeClass <- "C9 Normal"
enso$TypeClass[enso$Type == "WE" | 
                 enso$Type == "ME" | 
                 enso$Type == "SE"] <- "C1 all El Nino"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm =act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 all El Nino"])
color <- plot.colors[1]
plot.cloud.20m.enso.elnino.all <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(300), 
                     y.text = c(0.20),
                     labels.text = c("all El Ninos"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.all

# Pure El Nino
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WE" | enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD != "P" & enso$IOD != "M"] <- "C1 pure El Ninos"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 pure El Ninos"])
color <- plot.colors[2]
plot.cloud.20m.enso.elnino.pure <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(900), 
                     y.text = c(0.50),
                     labels.text = c("pure El Ninos"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.pure

# Pure medium and stron El Nino
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD != "P" & enso$IOD != "M"] <- "C1 pure m/s El Ninos"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 pure m/s El Ninos"])
color <- plot.colors[3]
plot.cloud.20m.enso.elnino.pure.ms <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(950), 
                     y.text = c(0.6),
                     labels.text = c("pure m/s El Ninos"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.pure.ms

# El Nino with IOD+
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WE" | enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD == "P" ] <- "C1 El Nino w IOD+"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 El Nino w IOD+"])
color <- plot.colors[4]
plot.cloud.20m.enso.elnino.wIOD.all <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(500), 
                     y.text = c(0.57),
                     labels.text = c("El Ninos w IOD+"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.wIOD.all

# Medium and strong El Nino with IOD+
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD == "P" ] <- "C1 m/s El Nino w IOD+"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 m/s El Nino w IOD+"])
color <- plot.colors[5]
plot.cloud.20m.enso.elnino.wIOD.ms <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(250), 
                     y.text = c(0.50),
                     labels.text = c("m/s El Ninos w IOD+"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.wIOD.ms

# El Nino with IOD-
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WE" | enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD == "M" ] <- "C1 El Nino w IOD-"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 El Nino w IOD-"])

# IOD+ - all
enso$TypeClass <- "C9 Normal"
enso$TypeClass[enso$Type != "ME" & enso$Type != "SE" &
                 enso$Type != "ML" & enso$Type != "SL" &
                 enso$IOD == "P"] <- "C1 purest IOD+"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 purest IOD+"])
color <- plot.colors[6]
plot.cloud.20m.iod.p.all <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = 2,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(600), 
                     y.text = c(0.50),
                     labels.text = c("purest IOD+"),
                     colors.text = color)
plot.cloud.20m.iod.p.all

plot.cloud.20m.elnino.clust2 <- 
  Reduce("outLayer", c(list(plot.cloud.20m.enso.elnino.all,
                            plot.cloud.20m.enso.elnino.pure,
                            plot.cloud.20m.enso.elnino.pure.ms,
                            plot.cloud.20m.enso.elnino.wIOD.all,
                            plot.cloud.20m.enso.elnino.wIOD.ms,
                            plot.cloud.20m.iod.p.all)))

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.cloud.20m.elnino.clust2.tif"),
       width = 30, height = 15, units = "cm", res = 600, pointsize =  5)
  plot(plot.cloud.20m.elnino.clust2)
  dev.off()
  pdf(file = paste0(graphicsPath, "plot.cloud.20m.elnino.clust2.pdf"),
      width = 12, height = 6, paper = "a4r")
  plot(plot.cloud.20m.elnino.clust2)
  dev.off()
} else {
  plot(plot.cloud.20m.elnino.clust2)
}



# Cluster 3 - El Nino
act.prm <- "clust3"
cloud.20m.seasonalwetdry <- seasonalMean(data = cloud,
                                         st = c(2003), nd = c(2013), 
                                         st.shift = 7, nd.shift = 0, timespan = 20, fun = "median", prm = act.prm)
cloud.20m.seasonalwetdry.split <- 
  split(cloud.20m.seasonalwetdry, cloud.20m.seasonalwetdry$season)
cloud.20m.seasonal.normal <- 
  list(lapply(cloud.20m.seasonalwetdry.split, function(x){x$p_dyn})$"2003-2013")

yminmax = c(0, 1)
plot.colors <- c("grey")

plot.cloud.20m.normal <- 
  visSeasonPlotByAOI(cloud.20m.seasonal.normal, plot.colors,
                     linetype = c(3), ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20, ylab = "Cloud cover")


plot.colors <- c("black", blue[3], blue[4], red[3], red[4], "darkgreen", "black")
linetype <- c(1)


# El Nino - all
enso$TypeClass <- "C9 Normal"
enso$TypeClass[enso$Type == "WE" | 
                 enso$Type == "ME" | 
                 enso$Type == "SE"] <- "C1 all El Nino"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm =act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 all El Nino"])
color <- plot.colors[1]
plot.cloud.20m.enso.elnino.all <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(100), 
                     y.text = c(0.20),
                     labels.text = c("all El Ninos"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.all

# Pure El Nino
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WE" | enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD != "P" & enso$IOD != "M"] <- "C1 pure El Ninos"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 pure El Ninos"])
color <- plot.colors[2]
plot.cloud.20m.enso.elnino.pure <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(100), 
                     y.text = c(0.25),
                     labels.text = c("pure El Ninos"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.pure

# Pure medium and stron El Nino
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD != "P" & enso$IOD != "M"] <- "C1 pure m/s El Ninos"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 pure m/s El Ninos"])
color <- plot.colors[3]
plot.cloud.20m.enso.elnino.pure.ms <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(320), 
                     y.text = c(0.35),
                     labels.text = c("pure m/s El Ninos"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.pure.ms

# El Nino with IOD+
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WE" | enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD == "P" ] <- "C1 El Nino w IOD+"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 El Nino w IOD+"])
color <- plot.colors[4]
plot.cloud.20m.enso.elnino.wIOD.all <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(300), 
                     y.text = c(0.8),
                     labels.text = c("El Ninos w IOD+"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.wIOD.all

# Medium and strong El Nino with IOD+
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD == "P" ] <- "C1 m/s El Nino w IOD+"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 m/s El Nino w IOD+"])
color <- plot.colors[5]
plot.cloud.20m.enso.elnino.wIOD.ms <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = linetype,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(375), 
                     y.text = c(0.40),
                     labels.text = c("m/s El Ninos w IOD+"),
                     colors.text = color)
plot.cloud.20m.enso.elnino.wIOD.ms

# El Nino with IOD-
enso$TypeClass <- "C9 Normal"
enso$TypeClass[(enso$Type == "WE" | enso$Type == "ME" | enso$Type == "SE") &
                 enso$IOD == "M" ] <- "C1 El Nino w IOD-"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 El Nino w IOD-"])

# IOD+ - all
enso$TypeClass <- "C9 Normal"
enso$TypeClass[enso$Type != "ME" & enso$Type != "SE" &
                 enso$Type != "ML" & enso$Type != "SL" &
                 enso$IOD == "P"] <- "C1 purest IOD+"
cloud.20m <- cloud[7:(nrow(cloud)-0), ]
cloud.20m.enso.split.median <- mergeAOIwTS(enso, cloud.20m, 
                                           timespan = 20,
                                           ts.prm = act.prm,
                                           rt = "median")
cloud.20m.info <- mergeAOIwTS(enso, cloud.20m, timespan = 20,
                              ts.prm = act.prm, rt = "org")
unique(cloud.20m.info$Season[cloud.20m.info$TypeClass == "C1 purest IOD+"])
color <- plot.colors[6]
plot.cloud.20m.iod.p.all <- 
  visSeasonPlotByAOI(cloud.20m.enso.split.median[1], color,
                     linetype = 2,
                     normal = plot.cloud.20m.normal,
                     ymin = yminmax[1], ymax = yminmax[2],
                     timespan = 20,
                     vline.pos = 501,
                     ylab = "Cloud cover",
                     x.text = c(100), 
                     y.text = c(0.70),
                     labels.text = c("purest IOD+"),
                     colors.text = color)
plot.cloud.20m.iod.p.all

plot.cloud.20m.elnino.clust3 <- 
  Reduce("outLayer", c(list(plot.cloud.20m.enso.elnino.all,
                            plot.cloud.20m.enso.elnino.pure,
                            plot.cloud.20m.enso.elnino.pure.ms,
                            plot.cloud.20m.enso.elnino.wIOD.all,
                            plot.cloud.20m.enso.elnino.wIOD.ms,
                            plot.cloud.20m.iod.p.all)))

if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.cloud.20m.elnino.clust3.tif"),
       width = 30, height = 15, units = "cm", res = 600, pointsize =  5)
  plot(plot.cloud.20m.elnino.clust3)
  dev.off()
  pdf(file = paste0(graphicsPath, "plot.cloud.20m.elnino.clust3.pdf"),
      width = 12, height = 6, paper = "a4r")
  plot(plot.cloud.20m.elnino.clust3)
  dev.off()
} else {
  plot(plot.cloud.20m.elnino.clust3)
}



plot.cloud.20m.elnino  <- 
  Reduce("grid13", list(plot.cloud.20m.elnino.clust1, 
                        plot.cloud.20m.elnino.clust2,
                        plot.cloud.20m.elnino.clust3))


if(printToFile == TRUE){
  tiff(filename = paste0(graphicsPath, 
                         "plot.cloud.20m.elnino.tif"),
       width = 30, height = 15, units = "cm", res = 600, pointsize =  5)
  plot(plot.cloud.20m.elnino)
  dev.off()
  
  pdf(file = paste0(graphicsPath, "plot.cloud.20m.elnino.pdf"),
      width = 12, height = 6, paper = "a4r")
  plot(plot.cloud.20m.elnino)
  dev.off()
} else {
  plot(plot.cloud.20m.elnino)
}











































cloudClust.20m.enso.split.median <- mergeAOIwTS(enso, cloudClust.20m, 
                                            timespan = 20,
                                            ts.prm = "clust1",
                                            rt = "median")



precip.20m.seasonal.normal <- 
  list(lapply(precip.20m.seasonalwetdry.split, function(x){x$p_dyn})$"1975-2013")





# Compute plot for long-term normal distribution
enso$TypeClass <- "N"
yminmax = c(0.0, 1.0)
plot.cloudClust.seasonal.normal.all <- 
  lapply(seq(3), function(x){
    cloudClust.shift06m <- cloudClust[7:(nrow(cloudClust)-6), ]
    colnames(cloudClust.shift06m)[x+1] <-  "P_RT_NRT"
    cloudClust.shift06m.enso.split.median <- combineAOPI(enso, cloudClust.shift06m)
    colors <- c("black")
    plot.cloudClust.shift06m.enso.split.median.normal <- 
      seasonPlotByAOI(cloudClust.shift06m.enso.split.median, colors,
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