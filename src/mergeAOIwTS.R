mergeAOIwTS <- function(aoi, 
                        ts,
                        aoi.prm = "aoi",
                        ts.prm = "P_RT_NRT", 
                        rt = "median") {
  
  aoi.reshape <- do.call("rbind", lapply(seq_len(nrow(aoi)), function(i) {
    st <- as.Date(paste0(substr(aoi[i, 3], 1, 4), "-07-01"))
    nd <- as.Date(paste0(substr(aoi[i, 3], 6, 9), "-12-01"))
    data.frame(aoi[i, c(1:3, ncol(aoi))], 
               aoi = as.numeric(c(aoi[i, 4:(ncol(aoi)-1)],
                                  aoi[i+1, 4:(ncol(aoi)-7)])), 
               ts = as.character(seq(st, nd, "month")),
               stringsAsFactors=FALSE)
  }))
  aoi.reshape$month <- rep(seq(1:(nrow(aoi.reshape)/18)), each = 18)
  
  ts$ts <- as.character(substr(ts$ts,1,10))
  ts.aoi <- 
    merge(ts, aoi.reshape, all = TRUE, by.x = "ts", by.y = "ts", 
          sort = FALSE)
  ts.aoi <- ts.aoi[with(ts.aoi, order(month)), ]
  ts.aoi.split <- 
    split(ts.aoi, ts.aoi$TypeClass)
  
  ts.aoi.split.median <- 
    foreach(i = ts.aoi.split) %do% {
      sapply(1:18, function(j) {
        median(i[, grep(ts.prm, colnames(i))][seq(j, nrow(i), 18)], na.rm = TRUE)
      })
    }
  
  if(rt == "median"){
    return(ts.aoi.split.median)
  } else if(rt == "split") {
    return(ts.aoi.split)
  } else {
    return(ts.aoi)
  }
}
