corPlotTimeSeries <- function(df,
                              xPrm,
                              yPrm,
                              tPrm,
                              xlable = substr(xPrm, 1, 3),
                              ylable = substr(yPrm, 1, 3),
                              lableNameCounter = c(1:12),
                              method = "kendall",
                              pthv = 0.05,
                              plotFilePath = "corPlot.tif",
                              printToFile = FALSE,
                              ...){
  
  years <- unique(substr(df[, grep(tPrm, colnames(df))], 1, 4))
    
  xData <- foreach(i = years, .combine = "rbind") %do% {
    dat <- df[grep(i, df[, grep(tPrm, colnames(df))]), ]
    foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
      dat[, grep(xPrm, colnames(dat))][j]
    }
  }
  colnames(xData) <- paste0(xlable, "-", sprintf("%02d", lableNameCounter))

  yData <- foreach(i = years, .combine = "rbind") %do% {
    dat <- df[grep(i, df[, grep(tPrm, colnames(df))]), ]
    foreach(j = seq_len(nrow(dat)), .combine = "c") %do% {
      dat[, grep(yPrm, colnames(dat))][j]
    }
  }
  colnames(yData) <- paste0(ylable, "-",sprintf("%02d", lableNameCounter))
  
  # Depricated
  corEstimate <- cor(xData, yData, use = "pairwise.complete.obs", 
                     method = method)
  corEstimate <- t(corEstimate)

  cxy <- lapply(seq(ncol(xData)), function(x) {
    vx <- xData[,x]
    cy <- lapply(seq(ncol(yData)), function(y) {
      vy <- yData[,y]
      c <- cor.test(vx, vy, alternative = "two.sided", method = method)
      return(list(estimate = c$estimate,
                  pvalue = c$p.value))
    })
  })
  corEstimate <- sapply(cxy, function(x) {
                    sapply(x, function(y){y$estimate})})
  corPValue <- sapply(cxy, function(x) {
                  sapply(x, function(y){y$pvalue})})
  colnames(corEstimate) <- colnames(xData)
  rownames(corEstimate) <- colnames(yData)
  colnames(corPValue) <- colnames(xData)
  rownames(corPValue) <- colnames(yData)
  
  
  for (r in seq(nrow(corEstimate))) {
    for (c in seq(ncol(corEstimate))) {
      if (r < c) {
        corEstimate[r,c] <- NA
      }
    }
  }
  
  corEstimateSig <- corEstimate
  corEstimateSig[corPValue > pthv] <- NA
  print(max(corEstimateSig, na.rm = TRUE))
  print(min(corEstimateSig, na.rm = TRUE))

  colors <- colorRampPalette(brewer.pal(9, "RdBu"))
  plot.corEstimate <- levelplot((corEstimate), col.regions = colors(100), 
                                at = seq(-.99, .99, .05), xlab = "", ylab = "",
                                colorkey = FALSE)
  plot.corEstimateSig <- levelplot((corEstimateSig), col.regions = colors(100), 
                           at = seq(-.99, .99, .05), xlab = "", ylab = "", 
                           border = "black")
  plot.cor <- doubleYScale(plot.corEstimate, plot.corEstimateSig, 
                           add.axis = FALSE)
  
  if(printToFile == TRUE){
    tiff(filename = plotFilePath, 
         width = 2480, height = 1748 , res = 300, pointsize =  12)
    plot(plot.cor)
    dev.off()
  } else {
    plot(plot.cor)
  }
}