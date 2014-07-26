combineAOPI <- function(oni, precip.shift06m){
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
}
