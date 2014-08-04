corPlotByAOI <- function(precip.shift06m.aoi,
                         parameter = "ssn_kz03k01",
                         xlable = "PR",
                         ylable = "ONI",
                         plotFilePath,
                         printToFile = FALSE){

  years <- unique(substr(precip.shift06m.aoi$ts, 1, 4))
  precip.shift06m.mat.ssn_kz03k01 <- foreach(i = years, .combine = "rbind") %do% {
    dat <- precip.shift06m.aoi[grep(i, precip.shift06m.aoi$ts), ]
    foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
      dat[, grep(parameter, colnames(dat))][j]
    }
  }
  colnames(precip.shift06m.mat.ssn_kz03k01) <- paste0(xlable, "-",
                                                      sprintf("%02d", c(7:12, 1:6)))
  
  aoi.reshape.shift06m.mat <- foreach(i = years, .combine = "rbind") %do% {
    dat <- precip.shift06m.aoi[grep(i, precip.shift06m.aoi$ts), ]
    foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
      dat$aoi[j]
    }
  }
  colnames(aoi.reshape.shift06m.mat) <- paste0(ylable, "-",
                                                sprintf("%02d", c(7:12, 1:6)))
  
#   plot.precip.shift06m.aoi.ssn_kz03k01 <- 
#     corPlotByAOI(aoi.reshape.shift06m.mat, precip.shift06m.mat.ssn_kz03k01)
  
  precip.shift06m.ssn_kz03k01.aoi.reshape.shift06m.cor <- 
    cor(aoi.reshape.shift06m.mat, precip.shift06m.mat.ssn_kz03k01, use = "complete.obs", 
        method = "kendall")
  
  print(max(precip.shift06m.ssn_kz03k01.aoi.reshape.shift06m.cor))
  print(min(precip.shift06m.ssn_kz03k01.aoi.reshape.shift06m.cor))
  for (r in 1:12) {
    for (c in 1:12) {
      if (r > c) {
        precip.shift06m.ssn_kz03k01.aoi.reshape.shift06m.cor[r,c] <- NA
      }
    }
  }
  
  colors <- colorRampPalette(brewer.pal(9, "RdBu"))
  plot.precip.shift06m.aoi.ssn_kz03k01 <- 
    levelplot(t(precip.shift06m.ssn_kz03k01.aoi.reshape.shift06m.cor), 
              col.regions = colors(100), 
              at = seq(-.99, .99, .05), xlab = "", ylab = "")
  
  if(printToFile == TRUE){
    tiff(filename = paste0(graphicsPath, plotFilePath),
         width = 2480, height = 1748 , res = 300, pointsize =  12)
    plot(plot.precip.shift06m.aoi.ssn_kz03k01)
    dev.off()
  } else {
    plot(plot.precip.shift06m.aoi.ssn_kz03k01)
  }
  
  
  
  
  
  
}
