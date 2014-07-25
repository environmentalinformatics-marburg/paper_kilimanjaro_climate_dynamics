seasonPlotByAOI <- function(oni, precip.shift06m){
  oni.reshape <- do.call("rbind", lapply(seq_len(nrow(oni)), function(i) {
    st <- as.Date(paste0(substr(oni[i, 3], 1, 4), "-07-01"))
    nd <- as.Date(paste0(substr(oni[i, 3], 6, 9), "-06-01"))
    data.frame(oni[i, c(1:3, ncol(oni))], 
               oni = as.numeric(oni[i, 4:(ncol(oni)-1)]), 
               month = seq(st, nd, "month"))
  }))
  
  # Adjust precipitation data records to ENSO cycle logic (i.e. start in July and
  # end in June) and combine precipitation and ONI data set
#   precip.shift06m <- precip[7:(nrow(precip)-6), ]
  
  precip.shift06m.oni <- 
    merge(precip.shift06m, oni.reshape, all.x = TRUE, by.x = "ts", by.y = "month")
  
  precip.shift06m.oni.split <- 
    split(precip.shift06m.oni, precip.shift06m.oni$TypeClass)
  
  precip.shift06m.oni.split.median <- 
    foreach(i = precip.shift06m.oni.split) %do% {
      sapply(1:12, function(j) {
        median(i$P_RT_NRT[seq(j, nrow(i), 12)], na.rm = TRUE)
      })
    }
  
  # Produce publication quality graph for mean seasonal cycles during El Nino,
  # La Nina and normal situations. To smooth the cycle while not reducing the
  # rainfall amounts significantly, use a spline prediction.
  if(length(precip.shift06m.oni.split) == 5){
    colors <- c("blue", "lightblue", "red", "bisque", "black")
  } else {
    colors <- c("blue", "red", "black")  
  }
  xhres <- seq(1, 12, 0.01)
  at <- seq(1, 1101, 100)
  labels <- c(seq(7, 12, 1), seq(1, 6, 1))
  plot.precip.shift06m.oni.split.median <- 
    lapply(seq(precip.shift06m.oni.split.median), function(x){
      sp <- predict(smooth.spline(
        precip.shift06m.oni.split.median[[x]], spar=0.01), xhres)
      xnew <- as.factor(sp[[1]])
      factor(xnew, levels = c(seq(7,12,0.01), seq(1,6,0.01)))
      xyplot(sp[[2]] ~ xnew, type = "l", lwd = "2", 
             ylim = c(0, 200), col = colors[x],
             scale=list(x=list(at = at, labels = labels)),
             xlab = "Month", ylab = "Precipitation")
    })
#   plot.precip.shift06m.oni.split.median.all <- 
#     Reduce("outLayer", plot.precip.shift06m.oni.split.median)
  
}
