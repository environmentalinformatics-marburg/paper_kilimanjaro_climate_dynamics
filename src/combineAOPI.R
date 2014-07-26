combineAOPI <- function(aoi, precip.shift06m, rt = "median"){
  aoi.reshape <- do.call("rbind", lapply(seq_len(nrow(aoi)), function(i) {
    st <- as.Date(paste0(substr(aoi[i, 3], 1, 4), "-07-01"))
    nd <- as.Date(paste0(substr(aoi[i, 3], 6, 9), "-06-01"))
    data.frame(aoi[i, c(1:3, ncol(aoi))], 
               oni = as.numeric(aoi[i, 4:(ncol(aoi)-1)]), 
               month = seq(st, nd, "month"))
  }))
  
  # Adjust precipitation data records to ENSO cycle logic (i.e. start in July and
  # end in June) and combine precipitation and aoi data set
  #   precip.shift06m <- precip[7:(nrow(precip)-6), ]
  
  precip.shift06m.aoi <- 
    merge(precip.shift06m, aoi.reshape, all.x = TRUE, by.x = "ts", by.y = "month")
  
  precip.shift06m.aoi.split <- 
    split(precip.shift06m.aoi, precip.shift06m.aoi$TypeClass)
  
  precip.shift06m.aoi.split.median <- 
    foreach(i = precip.shift06m.aoi.split) %do% {
      sapply(1:12, function(j) {
        median(i$P_RT_NRT[seq(j, nrow(i), 12)], na.rm = TRUE)
      })
    }
  if(rt == "median"){
    return(precip.shift06m.aoi.split.median)
  } else if(rt == "split") {
    return(precip.shift06m.aoi.split)
  } else {
    return(precip.shift06m.aoi)
  }
}
