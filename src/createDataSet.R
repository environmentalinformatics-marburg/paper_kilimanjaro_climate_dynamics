createDataSet <- function(...){
  #### Set working environment ###################################################
  switch(Sys.info()[["sysname"]], 
         "Linux" = setwd("/media/permanent/active/kilimanjaro"), 
         "Windows" = setwd("D:/active/kilimanjaro"))
  #setwd("C:/Users/Dogbert/Desktop/kilimanjaro")) 
  
  sourcePath <- "scripts/paper_kilimanjaro_climate_dynamics/src/"
  dataPath <- "data/"
  graphicsPath <- "graphics/"
  
  source(paste0(sourcePath, "readData.R"))
  
  #### Functions #################################################################
  outLayer <- function(x, y) {
    x + as.layer(y)
  }
  
  grid31 <- function(x, y, ...) {
    update(c(x, y,
             layout = c(3, 1)),
           between = list(y = 0.3, x = 0.3))
  } 
  
  #### Read data sets ############################################################
  # Read temperature data and aggregate it to monthly values
  ta.list <- readData("temperature")
  
  # Read precipitation data and create a continous time series.
  precip.list <- readData("precipitation")
  
  # Read cloud data
  cloud.list <- readData("cloudEOT")
  
  # Read AOI
  aoi.list <- readData("aoi")
  
  #### Pre-process precipitation data ##########################################
  # Perform an outlier check on precipitation data, compute deseasoned values
  # and compute 3 month running mean of original precipitation values using a
  # Kolmogorov-Zurbenko filter with one iteration
  compute_precip <- function(precip){
    # precip.checked <- do.call(rbind, lapply(seq(1:12), function(i) {
    #   pcl <- precip[seq(i, nrow(precip), 12), ]
    #   plot(log(pcl$P_RT_NRT[with(pcl, order(P_RT_NRT))]))
    #   thv <- quantile(pcl$P_RT_NRT,  probs = c(0.999), na.rm = TRUE)
    #   pcl$P_RT_NRT[pcl$P_RT_NRT >= thv] <- NA
    #   return(pcl)
    # }))
    # precip.checked <-  precip.checked[with(precip.checked, order(ts)), ]
    # precip <- precip.checked 
    
    # plot(log(precip$P_RT_NRT[with(precip, order(P_RT_NRT))]))
    thv <- tail(sort(precip$P_RT_NRT), 3)[1]
    precip$P_RT_NRT[precip$P_RT_NRT >= thv] <- NA
    
    precip[is.na(precip$P_RT_NRT),]
    
    # precip$P_RT_NRT[is.na(precip$P_RT_NRT)] <- 
    #   abs(c(precip$P_RT_NRT[which(is.na(precip$P_RT_NRT))-1]-
    #          precip$P_RT_NRT[which(is.na(precip$P_RT_NRT))+1])/2)
    
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
    return(precip)
  }
  
  precip.list$KIA <- compute_precip(precip.list$KIA)
  precip.list$MOKIA <- compute_precip(precip.list$MOKIA)
  precip.list$MOMOSHI <- compute_precip(precip.list$MOMOSHI)
  precip.list$MOMEAN <- compute_precip(precip.list$MOMEAN)
  
  return(list(ta.list = ta.list, precip.list = precip.list, 
              cloud.list = cloud.list, aoi.list = aoi.list))
}

