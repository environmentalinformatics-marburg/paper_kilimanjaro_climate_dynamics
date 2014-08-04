combineAOPI <- function(aoi, precipfun, parameter = "P_RT_NRT", rt = "median"){
  
  aoi.reshape <- do.call("rbind", lapply(seq_len(nrow(aoi)), function(i) {
    st <- as.Date(paste0(substr(aoi[i, 3], 1, 4), "-07-01"))
    nd <- as.Date(paste0(substr(aoi[i, 3], 6, 9), "-06-01"))
    data.frame(aoi[i, c(1:3, ncol(aoi))], 
               aoi = as.numeric(aoi[i, 4:(ncol(aoi)-1)]), 
               month = seq(st, nd, "month"))
  }))
  
  # Adjust precipfunitation data records to ENSO cycle logic (i.e. start in July and
  # end in June) and combine precipfunitation and aoi data set
  #   precipfun <- precipfun[7:(nrow(precipfun)-6), ]
  
  precipfun.aoi <- 
    merge(precipfun, aoi.reshape, all.x = TRUE, by.x = "ts", by.y = "month")
  
#   precipfun.aoi <- 
#     precipfun.aoi[complete.cases(precipfun.aoi),]
#   
  precipfun.aoi.split <- 
    split(precipfun.aoi, precipfun.aoi$TypeClass)
  
  precipfun.aoi.split.median <- 
    foreach(i = precipfun.aoi.split) %do% {
      sapply(1:12, function(j) {
        median(i[, grep(parameter, colnames(i))][seq(j, nrow(i), 12)], na.rm = TRUE)
      })
    }
  
  precipfun.aoi.split.harmonics <- 
    foreach(i = precipfun.aoi.split, .combine = "rbind") %do% {
      data.frame(TypeClass = unique(i$TypeClass),
                 P_RT_NRT = vectorHarmonics(i[, grep(parameter, colnames(i))], 
                        st = c(1, 1), 
                        nd = c(nrow(i)/12, 12), 
                        m = 3))
       
   }
  
  precipfun.aoi.split.harmonics <- 
    split(precipfun.aoi.split.harmonics$P_RT_NRT, precipfun.aoi.split.harmonics$TypeClass)
  
  if(rt == "median"){
    return(precipfun.aoi.split.median)
  } else if(rt == "split") {
    return(precipfun.aoi.split)
  } else if (rt == "harmonics") {
    return(precipfun.aoi.split.harmonics)
  } else {
    return(precipfun.aoi)
  }
}
