visPLaNina <- function(yminmax = c(0,250),
                       x.text.pos = c(200, 275, 1150, 1050),
                       y.text.pos = c(25, 60, 110, 230),
                       ...){

  #### Seasonal precipitation analysis ###########################################
  # Split the original precipitation data set by the three main wet/dry phases
  # identified in the long-term trend figure (see above)and create publication
  # quality figure; to smooth the cycle while not reducing the rainfall amounts 
  # significantly, use a spline prediction.
  # Split the running mean data set by year and create publication quality figure;
  # to smooth the cycle while not reducing the rainfall amounts significantly, 
  # use a spline prediction.
  colors <- c("black", "blue3", "red", "cornflowerblue")
  #yminmax<- c(0, 250)
  #yminmax<- c(-150, 150)
  
  # 12 month season
  precip.seasonalwetdry <- seasonalMean(
    st = c(1975,1975,1992,2001), nd = c(2013,1991,2000,2013),
    st.shift = 7, nd.shift = 0, timespan = 12, fun = avrg, prm = "P_RT_NRT")
  
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
  #   st.shift = 7, nd.shift = 0, timespan = 18, fun = avrg, prm = "P_RT_NRT")
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
  #   st.shift = 7, nd.shift = 0, timespan = 24, fun = avrg, prm = "P_RT_NRT")
  # 
  # precip.24m.seasonalwetdry.split <- 
  #   split(precip.24m.seasonalwetdry, precip.24m.seasonalwetdry$season)
  # precip.24m.seasonal.normal <- 
  #   list(lapply(precip.24m.seasonalwetdry.split, function(x){x$p_dyn})$"1975-2013")
  
  # 20 month season
  precip.20m.seasonalwetdry <- NULL
  precip.20m.seasonalwetdry <- seasonalMean(
    st = c(1975,1975,1992,2001), nd = c(2013,1992,2000,2013), 
    st.shift = 7, nd.shift = 0, timespan = 20, fun = avrg, prm = "P_RT_NRT")
  
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
  red <- brewer.pal(4, "Reds")
  blue <- brewer.pal(4, "Blues")
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
  # #                                             rt = avrg)
  # 
  # precip.20m <- precip[7:(nrow(precip)-0), ]
  # precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
  #                                             timespan = 20,
  #                                             ts.prm = "P_RT_NRT",
  #                                             rt = avrg)
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
  #                                             rt = avrg)
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
  #yminmax = c(0, 500)

  plot.colors <- c("grey")
  
  plot.precip.20m.normal <- 
    visSeasonPlotByAOI(precip.20m.seasonal.normal, plot.colors,
                       linetype = c(2), ymin = yminmax[1], ymax = yminmax[2],
                       timespan = 20)
  
  
  plot.colors <- c("black", blue[3], blue[4], red[3], red[4], "darkgreen", "black")
  linetype <- c(1)
  
  # La Nina - all
  enso <- data.set$aoi.list$ONI
  enso$TypeClass <- "C9 Normal"
  enso$TypeClass[enso$Type == "WL" | 
                   enso$Type == "ML" | 
                   enso$Type == "SL"] <- "C1 all La Nina"
  precip.20m <- precip[7:(nrow(precip)-0), ]
  precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                              timespan = 20,
                                              ts.prm = "P_RT_NRT",
                                              rt = avrg)
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
                       x.text = x.text.pos[1], 
                       y.text = y.text.pos[1],
                       labels.text = c("all La Ninas"),
                       colors.text = color)
  plot.precip.20m.enso.lanina.all
  print(plot.precip.20m.enso.lanina.all)
  
  # Pure La Nina
  enso$TypeClass <- "C9 Normal"
  enso$TypeClass[(enso$Type == "WL" | enso$Type == "ML" | enso$Type == "SL") &
                   enso$IOD != "P" & enso$IOD != "M"] <- "C1 pure La Ninas"
  precip.20m <- precip[7:(nrow(precip)-0), ]
  precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                              timespan = 20,
                                              ts.prm = "P_RT_NRT",
                                              rt = avrg)
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
                       x.text = x.text.pos[2], 
                       y.text = y.text.pos[2],
                       labels.text = c("pure La Ninas"),
                       colors.text = color)
  plot.precip.20m.enso.lanina.pure
  print(plot.precip.20m.enso.lanina.pure)
  
  # Pure medium and stron La Nina
  enso$TypeClass <- "C9 Normal"
  enso$TypeClass[(enso$Type == "ML" | enso$Type == "SL") &
                   enso$IOD != "P" & enso$IOD != "M"] <- "C1 pure m/s La Ninas"
  precip.20m <- precip[7:(nrow(precip)-0), ]
  precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                              timespan = 20,
                                              ts.prm = "P_RT_NRT",
                                              rt = avrg)
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
                       x.text = x.text.pos[3], 
                       y.text = y.text.pos[3],
                       labels.text = c("pure m/s La Ninas"),
                       colors.text = color)
  plot.precip.20m.enso.lanina.pure.ms
  print(plot.precip.20m.enso.lanina.pure.ms)
  
  # La Nina with IOD+
  enso$TypeClass <- "C9 Normal"
  enso$TypeClass[(enso$Type == "WL" | enso$Type == "ML" | enso$Type == "SL") &
                   enso$IOD == "P" ] <- "C1 La Nina w IOD+"
  precip.20m <- precip[7:(nrow(precip)-0), ]
  precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                              timespan = 20,
                                              ts.prm = "P_RT_NRT",
                                              rt = avrg)
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
                       x.text = x.text.pos[4], 
                       y.text = y.text.pos[4],
                       labels.text = c("La Ninas w IOD+"),
                       colors.text = color)
  plot.precip.20m.enso.lanina.wIOD.all
  print(plot.precip.20m.enso.lanina.wIOD.all)
  
  # Medium and strong La Nina with IOD+
  enso$TypeClass <- "C9 Normal"
  enso$TypeClass[(enso$Type == "ML" | enso$Type == "SL") &
                   enso$IOD == "P" ] <- "C1 m/s La Nina w IOD+"
  precip.20m <- precip[7:(nrow(precip)-0), ]
  precip.20m.enso.split.median <- mergeAOIwTS(enso, precip.20m, 
                                              timespan = 20,
                                              ts.prm = "P_RT_NRT",
                                              rt = avrg)
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
                                              rt = avrg)
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
  
  
  plot.precip.20m.lanina <- 
    Reduce("outLayer", c(list(plot.precip.20m.enso.lanina.all,
                              plot.precip.20m.enso.lanina.pure,
                              plot.precip.20m.enso.lanina.pure.ms,
                              plot.precip.20m.enso.lanina.wIOD.all)))
  
#   if(printToFile == TRUE){
#     tiff(filename = paste0(graphicsPath, 
#                            "plot.precip.20m.lanina.tif"),
#          width = 30, height = 15, units = "cm", res = 600, pointsize =  5)
#     plot(plot.precip.20m.lanina)
#     dev.off()
#     
#     pdf(file = paste0(graphicsPath, "plot.precip.20m.lanina.pdf"),
#         width = 12, height = 6, paper = "a4r")
#     plot(plot.precip.20m.lanina)
#     dev.off()
#   } else {
#     plot(plot.precip.20m.lanina)
#   }
  
  return(plot.precip.20m.lanina)
  
}

